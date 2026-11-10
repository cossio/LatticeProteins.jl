# Local updates of all replicas. Threaded slots exclusively own their RNG, chain, and
# counters during this phase. The shared-RNG protocol is restricted to the serial path.
# The threaded loop joins before any serial swaps or endpoint bookkeeping.
function tempering_sweep!(rngs, chains, s, betas, nsteps, sites, boltzmann, counts, threaded)
    function update_slot!(i)
        counts.local_attempts[i] += nsteps
        return counts.local_accepts[i] += local_steps!(rngs[i], chains[i], s, betas[i], nsteps, sites, boltzmann)
    end
    if threaded && Threads.nthreads() > 1
        # Unlike :static, :dynamic also works inside another threaded caller.
        Threads.@threads :dynamic for i in eachindex(chains)
            update_slot!(i)
        end
    else
        for i in eachindex(chains)
            update_slot!(i)
        end
    end
    return nothing
end

# log[pi(swapped) / pi(current)] for pi_beta(s) ∝ exp(-beta * loss(s)).
swap_logratio(beta_i, beta_j, loss_i, loss_j) = (beta_i - beta_j) * (loss_i - loss_j)

# Proposes to exchange the chains at temperatures `i` and `i + 1`. The chains (with their
# cached energies) and their identities `ids` move together.
function tempering_swap!(rng, chains, ids, betas, i)
    j = i + 1
    logratio = swap_logratio(betas[i], betas[j], chains[i].loss, chains[j].loss)
    if -logratio < randexp(rng)
        chains[i], chains[j] = chains[j], chains[i]
        ids[i], ids[j] = ids[j], ids[i]
        return true
    end
    return false
end

tempering_counts(n) = (
    local_attempts = zeros(Int, n), local_accepts = zeros(Int, n),
    swap_attempts = zeros(Int, n - 1), swap_accepts = zeros(Int, n - 1),
    round_trips = zeros(Int, n),
)

# Round-trip accounting: a replica completes a trip when it returns to the hottest
# temperature after having visited the coldest one.
function record_endpoints!(stages, trips, ids)
    hot, cold = first(ids), last(ids)
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
    ParallelTempering(; betas, burnin = 1000, thin = 10, local_steps = 10,
                      threaded = true, rng_mode = :per_temperature)

Sampling method of [`sample_msa`](@ref): replica exchange, targeting
`pnat(s, seq)^betas[end]`. `betas` must contain at least two finite, nonnegative, strictly
increasing inverse selection temperatures (hot to cold). The physical folding temperature
stays at 1, and all 103,346 compact structures compete, just as in [`logpnat`](@ref).

One **round** performs `local_steps` Metropolis steps per replica (the moves of
[`Metropolis`](@ref)), then proposes swaps on alternating disjoint adjacent edges. For
two replicas the only edge is attempted every round. A swap is accepted with
`min(1, exp((betas[i]-betas[j]) * (E[i]-E[j])))`, where `E = -logpnat`.

`burnin` rounds are discarded, then the coldest sequence is saved after every `thin`
rounds until `nseqs` samples have been saved. Thus each replica receives
`(burnin + nseqs * thin) * local_steps` mutation proposals.

    sample_msa(s, ParallelTempering(; ...); nseqs = 100, rng = default_rng(), initial_seqs = nothing)

`initial_seqs`, when provided, is a `(27, length(betas))` integer matrix of starting
sequences, ordered hot to cold; it is copied. Otherwise each replica starts from an
independent uniform sequence.

Use an explicit RNG for reproducible runs. By default, each temperature slot uses
its own `Xoshiro` stream seeded from `rng`; swaps use `rng` serially. For fixed
arguments and Julia version, results do not depend on thread count or scheduling.

With `threaded=true`, local updates run on Julia's default thread pool (start Julia
with e.g. `--threads=4`). Set `threaded=false` when parallelizing independent ladders
in the caller. Both modes use the same per-temperature streams and produce identical
results.

`rng_mode=:shared` with `threaded=false` uses `rng` for both local updates and swaps,
without creating temperature streams. This makes a run resumable: a call with
`burnin=0` that passes the same RNG and the `final_seqs` of the previous call as
`initial_seqs` continues its trajectory exactly, provided the previous call ended on an
even number of rounds (which preserves the alternation of the swap edges). It is also
the RNG protocol of v2.1.0 and of the saved Tubiana PT alignment. The default
`:per_temperature` mode creates fresh local streams on each call and is not a checkpoint
interface.

Returns a named tuple with:
- `samples`: coldest-replica MSA, size `(27, nseqs)`. Samples can be correlated.
- `betas`: the fixed ladder used in the run.
- `rng_mode`: the RNG protocol used (`:per_temperature` or `:shared`).
- `logpnat`, `replica_ids`: `(length(betas), nseqs)` traces at saved rounds.
  Rows follow temperatures; IDs follow sequences through swaps.
- `final_seqs`, `final_logpnat`: final states at every temperature.
- `burnin_diagnostics`, `diagnostics`: separate burn-in and production counters.
  `local_attempts`/`local_accepts` are indexed by temperature;
  `swap_attempts`/`swap_accepts` by adjacent edge. `round_trips` is indexed by
  replica ID and counts completed hot→cold→hot trips entirely within that phase.

Inspect swap acceptance, round trips, and agreement between independent runs.
High folding probabilities or accepted swaps alone do not establish equilibrium.
Tune the ladder in pilot runs, then hold it fixed for production sampling.

```julia
betas = [0.0; exp.(range(log(0.1), log(1000.0); length = 25))]
result = sample_msa(STRUCTURE_A, ParallelTempering(; betas); nseqs = 500)
result.samples                  # the sequences at betas[end]
result.diagnostics.round_trips  # hot → cold → hot, per replica
```
"""
struct ParallelTempering
    betas::Vector{Float64}
    burnin::Int
    thin::Int
    local_steps::Int
    threaded::Bool
    rng_mode::Symbol
    function ParallelTempering(;
            betas::AbstractVector{<:Real}, burnin::Integer = 1000, thin::Integer = 10,
            local_steps::Integer = 10, threaded::Bool = true,
            rng_mode::Symbol = :per_temperature
        )
        burnin >= 0 || throw(ArgumentError("burnin must be nonnegative"))
        thin > 0 || throw(ArgumentError("thin must be positive"))
        local_steps > 0 || throw(ArgumentError("local_steps must be positive"))
        rng_mode in (:per_temperature, :shared) || throw(ArgumentError("invalid rng_mode"))
        rng_mode == :shared && threaded &&
            throw(ArgumentError("rng_mode=:shared requires threaded=false"))
        length(betas) >= 2 || throw(ArgumentError("at least two temperatures are required"))
        all(b -> isfinite(b) && b >= 0, betas) && all(diff(betas) .> 0) ||
            throw(ArgumentError("betas must be finite, nonnegative, and strictly increasing"))
        return new(Float64.(betas), burnin, thin, local_steps, threaded, rng_mode)
    end
end

function sample_msa(
        s::Integer, method::ParallelTempering; nseqs::Integer = 100,
        rng::AbstractRNG = default_rng(),
        initial_seqs::Union{Nothing, AbstractMatrix{<:Integer}} = nothing,
        structures::StructureSet = STRUCTURES
    )
    (; betas, burnin, thin, local_steps, threaded, rng_mode) = method
    check_structure(structures, s)
    nseqs >= 0 || throw(ArgumentError("nseqs must be nonnegative"))
    n = length(betas)
    if initial_seqs === nothing
        initial_seqs = rand(rng, 1:20, L, n)
    else
        size(initial_seqs) == (L, n) || throw(ArgumentError("initial_seqs must have size (27, length(betas))"))
        all(a -> 1 <= a <= 20, initial_seqs) || throw(ArgumentError("residues must be in 1:20"))
    end
    chains = [Chain(view(initial_seqs, :, i), s, structures) for i in 1:n]
    ids = collect(1:n) # ids[i] is the replica (numbered by its initial slot) now at temperature i
    rngs = rng_mode == :shared ? fill(rng, n) : [Xoshiro(rand(rng, UInt64)) for _ in 1:n]
    sites = site_contacts(structures)
    samples = Matrix{Int}(undef, L, nseqs)
    trace = Matrix{Float64}(undef, n, nseqs)
    replica_ids = Matrix{Int}(undef, n, nseqs)
    burnin_diagnostics, diagnostics = tempering_counts(n), tempering_counts(n)
    round_index = 0
    for (rounds, counts, production) in (
            (burnin, burnin_diagnostics, false),
            (nseqs * thin, diagnostics, true),
        )
        stages = zeros(Int, n)
        record_endpoints!(stages, counts.round_trips, ids)
        for step in 1:rounds
            round_index += 1
            tempering_sweep!(
                rngs, chains, s, betas, local_steps, sites, structures.boltzmann,
                counts, threaded
            )
            first_edge = n == 2 || isodd(round_index) ? 1 : 2
            for i in first_edge:2:(n - 1)
                counts.swap_attempts[i] += 1
                counts.swap_accepts[i] += tempering_swap!(rng, chains, ids, betas, i)
            end
            record_endpoints!(stages, counts.round_trips, ids)
            if production && step % thin == 0
                sample = step ÷ thin
                samples[:, sample] .= last(chains).seq
                for i in 1:n
                    trace[i, sample] = -chains[i].loss
                    replica_ids[i, sample] = ids[i]
                end
            end
        end
    end
    return (;
        samples, betas = copy(betas), rng_mode, logpnat = trace, replica_ids,
        final_seqs = stack(chain.seq for chain in chains),
        final_logpnat = [-chain.loss for chain in chains], burnin_diagnostics, diagnostics,
    )
end
