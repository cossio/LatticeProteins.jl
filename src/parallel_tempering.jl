# Cached contact energies travel with the sequence when replicas are exchanged.
mutable struct _TemperingReplica
    seq::Vector{Int}
    energies::Vector{Float64}
    trial::Vector{Float64}
    loss::Float64
    id::Int
end

# The bundled MJ entries permit exact contact sums in units of 0.01. This enables
# Boltzmann lookup instead of an exp evaluation for every structure and proposal.
# The integer cache is a speed optimization; no significant roundoff problem was
# observed in the floating cache with its existing 100-round refreshes.
mutable struct _IntegerTemperingReplica
    seq::Vector{Int}
    energies::Vector{Int32} # units of 0.01
    trial::Vector{Int32}
    loss::Float64
    id::Int
    interactions::Matrix{Int32}
    boltzmann::Vector{Float64}
end

function _integer_energy_grid(ncontacts, interactions=MJ)
    all(isfinite, interactions) || return nothing
    maximum(abs, interactions) <= typemax(Int32) / 100 || return nothing
    grid = round.(Int32, 100 .* interactions)
    # Exact eligibility check: do not quantize a modified interaction matrix.
    all(interactions .== grid ./ 100) || return nothing
    ncontacts * maximum(abs, Int64.(grid)) <= typemax(Int32) || return nothing
    span = ncontacts * (Int64(maximum(grid)) - Int64(minimum(grid)))
    # Bound table memory; large/custom models use the general floating kernel.
    span <= 100_000 || return nothing
    return grid, exp.(-(0:span) ./ 100)
end

function _integer_energies!(energies, seq, structures, interactions)
    @inbounds for s in axes(structures, 2)
        value = Int32(0)
        for k in axes(structures, 1)
            value += interactions[seq[structures[k,s,1]], seq[structures[k,s,2]]]
        end
        energies[s] = value
    end
    return energies
end

function _integer_folding_loss(energies, cm, boltzmann)
    emin = minimum(energies)
    z = 0.0
    @inbounds @simd for e in energies
        z += boltzmann[Int(e - emin) + 1]
    end
    return (energies[cm] - emin) / 100 + log(z)
end

function _IntegerTemperingReplica(seq, cm, id, structures, interactions, boltzmann)
    E = _integer_energies!(Vector{Int32}(undef, size(structures,2)), seq, structures, interactions)
    return _IntegerTemperingReplica(copy(seq), E, similar(E),
        _integer_folding_loss(E, cm, boltzmann), id, interactions, boltzmann)
end

function _integer_mutate_energies!(replica, i, aa, sites)
    copyto!(replica.trial, replica.energies)
    old = replica.seq[i]
    delta = ntuple(b -> replica.interactions[aa,b] - replica.interactions[old,b], Val(20))
    @inbounds for t in sites.offsets[i]:(sites.offsets[i+1]-1)
        replica.trial[sites.structure[t]] += delta[replica.seq[sites.partner[t]]]
    end
    return replica.trial
end

function _TemperingReplica(seq, cm, id, structures=STRUCTURES)
    E = energies!(Vector{Float64}(undef, size(structures, 2)), seq, structures)
    return _TemperingReplica(copy(seq), E, similar(E), _folding_loss(E, cm), id)
end

function _tempering_local!(rng, replica, cm, beta, nsteps, sites)
    accepted = 0
    for _ in 1:nsteps
        i = rand(rng, 1:L)
        # Uniform among the other 19 residues: symmetric proposals, no no-ops.
        aa = rand(rng, 1:19)
        aa += aa >= replica.seq[i]
        mutate_energies!(replica.trial, replica.energies, replica.seq, i, aa, sites)
        loss = _folding_loss(replica.trial, cm)
        if beta * (loss - replica.loss) < randexp(rng)
            replica.seq[i] = aa
            replica.energies, replica.trial = replica.trial, replica.energies
            replica.loss = loss
            accepted += 1
        end
    end
    return accepted
end

function _tempering_local!(rng, replica::_IntegerTemperingReplica, cm, beta, nsteps, sites)
    accepted = 0
    for _ in 1:nsteps
        i = rand(rng, 1:L)
        aa = rand(rng, 1:19)
        aa += aa >= replica.seq[i]
        _integer_mutate_energies!(replica, i, aa, sites)
        loss = _integer_folding_loss(replica.trial, cm, replica.boltzmann)
        if beta * (loss - replica.loss) < randexp(rng)
            replica.seq[i] = aa
            replica.energies, replica.trial = replica.trial, replica.energies
            replica.loss = loss
            accepted += 1
        end
    end
    return accepted
end

function _refresh_tempering!(replica, cm, structures)
    energies!(replica.energies, replica.seq, structures)
    replica.loss = _folding_loss(replica.energies, cm)
end

function _refresh_tempering!(replica::_IntegerTemperingReplica, cm, structures)
    _integer_energies!(replica.energies, replica.seq, structures, replica.interactions)
    replica.loss = _integer_folding_loss(replica.energies, cm, replica.boltzmann)
end

# Threaded slots exclusively own their RNG, replica, and counters during this phase.
# The shared-RNG protocol is restricted to the serial path.
# The threaded loop joins before any serial swaps or endpoint bookkeeping.
function _tempering_sweep!(rngs, replicas, cm, betas, nsteps, sites, structures,
                           counts, refresh, threaded)
    function update_slot!(i)
        replica = replicas[i]
        counts.local_attempts[i] += nsteps
        counts.local_accepts[i] += _tempering_local!(rngs[i], replica, cm, betas[i], nsteps, sites)
        if refresh
            _refresh_tempering!(replica, cm, structures)
        end
    end
    if threaded && Threads.nthreads() > 1
        # Unlike :static, :dynamic also works inside another threaded caller.
        Threads.@threads :dynamic for i in eachindex(replicas)
            update_slot!(i)
        end
    else
        for i in eachindex(replicas)
            update_slot!(i)
        end
    end
    return nothing
end

# log[pi(swapped) / pi(current)] for pi_beta(s) ∝ exp(-beta * loss(s)).
_swap_logratio(beta_i, beta_j, loss_i, loss_j) = (beta_i - beta_j) * (loss_i - loss_j)

function _tempering_swap!(rng, replicas, betas, i)
    j = i + 1
    logratio = _swap_logratio(betas[i], betas[j], replicas[i].loss, replicas[j].loss)
    if -logratio < randexp(rng)
        replicas[i], replicas[j] = replicas[j], replicas[i]
        return true
    end
    return false
end

_tempering_counts(n) = (local_attempts=zeros(Int, n), local_accepts=zeros(Int, n),
    swap_attempts=zeros(Int, n - 1), swap_accepts=zeros(Int, n - 1),
    round_trips=zeros(Int, n))

function _record_endpoints!(stages, trips, replicas)
    hot, cold = first(replicas).id, last(replicas).id
    if stages[hot] == 2
        trips[hot] += 1
    end
    stages[hot] = 1
    if stages[cold] == 1
        stages[cold] = 2
    end
    return nothing
end

"""
    parallel_tempering(cm::Int; betas, nseqs=100, burnin=1000, thin=10,
                       local_steps=10, initial_seqs=nothing, rng=default_rng(),
                       structures=STRUCTURES, energy_cache=:auto, threaded=true,
                       rng_mode=:per_temperature)

Sample sequences with replica exchange, targeting `pnat(cm, seq; structures)^betas[end]`.
`betas` must contain at least two finite, nonnegative, strictly increasing inverse
selection temperatures (hot to cold). The physical folding temperature stays at 1;
all 103,346 compact structures are included by default, just as in [`log_pnat`](@ref).
Pass `structures=SHAKHNOVICH_STRUCTURES` for the 10,000-structure subset used in the
paper. `cm` indexes the selected set; use [`structure_index`](@ref) to translate it.

One **round** performs `local_steps` symmetric single-residue proposals per
replica, then proposes swaps on alternating disjoint adjacent edges. For two
replicas the only edge is attempted every round. A swap is accepted with
`min(1, exp((betas[i]-betas[j]) * (E[i]-E[j])))`, where `E = -log_pnat`.
Cached structural energies are refreshed every 100 rounds to limit roundoff drift.

`energy_cache=:auto` uses exact integer contact sums and tabulated Boltzmann factors
when every interaction is exactly representable in units of 0.01 (as in the bundled
MJ matrix). Otherwise it uses the general floating-point kernel. `:float` forces
that original kernel; `:integer` requires eligibility. No approximate rounding of
interaction parameters is performed. The target distribution is unchanged.

Discard `burnin` rounds, then save the coldest sequence after every `thin` rounds
until `nseqs` samples have been saved. Thus each replica receives
`(burnin + nseqs * thin) * local_steps` mutation proposals. `initial_seqs`, when
provided, is a `(27, length(betas))` integer matrix, ordered hot to cold; it is
copied. Otherwise each replica starts from an independent uniform sequence.
Use an explicit RNG for reproducible runs. By default, each temperature slot uses
its own `Xoshiro` stream seeded from `rng`; swaps use `rng` serially. For fixed
arguments and Julia version, results do not depend on thread count or scheduling.

With `threaded=true`, local updates and cache refreshes run on Julia's default
thread pool (start Julia with e.g. `--threads=4`). Set `threaded=false` when
parallelizing independent ladders in the caller. Both modes use the same
per-temperature streams and produce identical results.

`rng_mode=:shared` with `threaded=false` reproduces the single-stream RNG protocol
used in v2.1.0 and the saved Tubiana PT alignment. It uses `rng` for both local
updates and swaps without creating temperature streams. This also permits a
continuation call to preserve that stream by passing the same RNG and final
sequences; end chunks on complete 100-round cycles to preserve refresh and swap
phases. The default `:per_temperature` mode creates fresh local streams on each
call and is not a checkpoint interface.

Returns a named tuple with:
- `samples`: coldest-replica MSA, size `(27, nseqs)`. Samples can be correlated.
- `betas`: the fixed ladder used in the run.
- `energy_cache`: the actual cache used (`:integer` or `:float`).
- `rng_mode`: the RNG protocol used (`:per_temperature` or `:shared`).
- `log_pnat`, `replica_ids`: `(length(betas), nseqs)` traces at saved rounds.
  Rows follow temperatures; IDs follow sequences through swaps.
- `final_seqs`, `final_log_pnat`: final states at every temperature.
- `burnin_diagnostics`, `diagnostics`: separate burn-in and production counters.
  `local_attempts`/`local_accepts` are indexed by temperature;
  `swap_attempts`/`swap_accepts` by adjacent edge. `round_trips` is indexed by
  replica ID and counts completed hot→cold→hot trips entirely within that phase.

Inspect swap acceptance, round trips, and agreement between independent runs.
High folding probabilities or accepted swaps alone do not establish equilibrium.
Tune the ladder in pilot runs, then hold it fixed for production sampling.

See also [`metropolis!`](@ref), [`random_msa`](@ref).
"""
function parallel_tempering(cm::Int; betas::AbstractVector{<:Real}, nseqs::Int=100,
    burnin::Int=1000, thin::Int=10, local_steps::Int=10,
    initial_seqs::Union{Nothing,AbstractMatrix{Int}}=nothing,
    rng::AbstractRNG=default_rng(), structures::AbstractArray{Int,3}=STRUCTURES,
    energy_cache::Symbol=:auto, threaded::Bool=true,
    rng_mode::Symbol=:per_temperature)

    check_structures(structures)
    1 <= cm <= size(structures, 2) || throw(ArgumentError("invalid target structure"))
    nseqs >= 0 || throw(ArgumentError("nseqs must be nonnegative"))
    burnin >= 0 || throw(ArgumentError("burnin must be nonnegative"))
    thin > 0 || throw(ArgumentError("thin must be positive"))
    local_steps > 0 || throw(ArgumentError("local_steps must be positive"))
    energy_cache in (:auto, :float, :integer) || throw(ArgumentError("invalid energy_cache"))
    rng_mode in (:per_temperature, :shared) || throw(ArgumentError("invalid rng_mode"))
    rng_mode == :shared && threaded &&
        throw(ArgumentError("rng_mode=:shared requires threaded=false"))
    betas = Float64.(betas)
    n = length(betas)
    n >= 2 || throw(ArgumentError("at least two temperatures are required"))
    all(b -> isfinite(b) && b >= 0, betas) && all(diff(betas) .> 0) ||
        throw(ArgumentError("betas must be finite, nonnegative, and strictly increasing"))
    if initial_seqs === nothing
        initial_seqs = rand(rng, 1:20, L, n)
    else
        size(initial_seqs) == (L, n) || throw(ArgumentError("initial_seqs must have size (27, length(betas))"))
        all(a -> 1 <= a <= 20, initial_seqs) || throw(ArgumentError("residues must be in 1:20"))
    end
    grid = energy_cache == :float ? nothing : _integer_energy_grid(size(structures,1))
    energy_cache == :integer && grid === nothing &&
        throw(ArgumentError("interaction matrix/structure size is not eligible for integer caching"))
    replicas = if grid === nothing
        [_TemperingReplica(view(initial_seqs, :, i), cm, i, structures) for i in 1:n]
    else
        [_IntegerTemperingReplica(view(initial_seqs, :, i), cm, i, structures, grid...) for i in 1:n]
    end
    used_cache = grid === nothing ? :float : :integer
    rngs = rng_mode == :shared ? fill(rng, n) : [Xoshiro(rand(rng, UInt64)) for _ in 1:n]
    sites = site_contacts(structures)
    samples = Matrix{Int}(undef, L, nseqs)
    trace = Matrix{Float64}(undef, n, nseqs)
    ids = Matrix{Int}(undef, n, nseqs)
    burnin_diagnostics, diagnostics = _tempering_counts(n), _tempering_counts(n)
    round_index = 0
    for (rounds, counts, production) in ((burnin, burnin_diagnostics, false),
                                        (nseqs * thin, diagnostics, true))
        stages = zeros(Int, n)
        _record_endpoints!(stages, counts.round_trips, replicas)
        for step in 1:rounds
            round_index += 1
            _tempering_sweep!(rngs, replicas, cm, betas, local_steps, sites, structures,
                              counts, round_index % 100 == 0, threaded)
            first_edge = n == 2 || isodd(round_index) ? 1 : 2
            for i in first_edge:2:n-1
                counts.swap_attempts[i] += 1
                counts.swap_accepts[i] += _tempering_swap!(rng, replicas, betas, i)
            end
            _record_endpoints!(stages, counts.round_trips, replicas)
            if production && step % thin == 0
                sample = step ÷ thin
                samples[:, sample] .= last(replicas).seq
                for i in 1:n
                    trace[i, sample] = -replicas[i].loss
                    ids[i, sample] = replicas[i].id
                end
            end
        end
    end
    return (; samples, betas, energy_cache=used_cache, rng_mode, log_pnat=trace, replica_ids=ids,
        final_seqs=stack(r.seq for r in replicas),
        final_log_pnat=[-r.loss for r in replicas], burnin_diagnostics, diagnostics)
end
