# Cached contact energies travel with the sequence when replicas are exchanged.
mutable struct _TemperingReplica
    seq::Vector{Int}
    energies::Vector{Float64}
    trial::Vector{Float64}
    loss::Float64
    id::Int
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

# Each slot exclusively owns its RNG, replica, and counters during this phase.
# The threaded loop joins before any serial swaps or endpoint bookkeeping.
function _tempering_sweep!(rngs, replicas, cm, betas, nsteps, sites, structures,
                           counts, refresh, threaded)
    function update_slot!(i)
        replica = replicas[i]
        counts.local_attempts[i] += nsteps
        counts.local_accepts[i] += _tempering_local!(rngs[i], replica, cm, betas[i], nsteps, sites)
        if refresh
            energies!(replica.energies, replica.seq, structures)
            replica.loss = _folding_loss(replica.energies, cm)
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
                       structures=STRUCTURES, threaded=true)

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

Discard `burnin` rounds, then save the coldest sequence after every `thin` rounds
until `nseqs` samples have been saved. Thus each replica receives
`(burnin + nseqs * thin) * local_steps` mutation proposals. `initial_seqs`, when
provided, is a `(27, length(betas))` integer matrix, ordered hot to cold; it is
copied. Otherwise each replica starts from an independent uniform sequence.
Use an explicit RNG for reproducible runs. Each temperature slot uses its own
`Xoshiro` stream seeded from `rng`; swaps use `rng` serially. For fixed arguments
and a fixed Julia version, results do not depend on thread count or scheduling.
These streams differ from the sampler's original single RNG stream in v2.1.0.

With `threaded=true`, local updates and cache refreshes run on Julia's default
thread pool (start Julia with e.g. `--threads=4`). Set `threaded=false` to run
serially, for example when parallelizing independent ladders in the caller.
Both modes use the same streams and produce identical results.

Returns a named tuple with:
- `samples`: coldest-replica MSA, size `(27, nseqs)`. Samples can be correlated.
- `betas`: the fixed ladder used in the run.
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
    threaded::Bool=true)

    check_structures(structures)
    1 <= cm <= size(structures, 2) || throw(ArgumentError("invalid target structure"))
    nseqs >= 0 || throw(ArgumentError("nseqs must be nonnegative"))
    burnin >= 0 || throw(ArgumentError("burnin must be nonnegative"))
    thin > 0 || throw(ArgumentError("thin must be positive"))
    local_steps > 0 || throw(ArgumentError("local_steps must be positive"))
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
    replicas = [_TemperingReplica(view(initial_seqs, :, i), cm, i, structures) for i in 1:n]
    rngs = [Xoshiro(rand(rng, UInt64)) for _ in 1:n]
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
    return (; samples, betas, log_pnat=trace, replica_ids=ids,
        final_seqs=stack(r.seq for r in replicas),
        final_log_pnat=[-r.loss for r in replicas], burnin_diagnostics, diagnostics)
end
