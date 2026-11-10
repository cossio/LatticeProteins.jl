using Test: @test, @testset, @test_throws
using Random: Xoshiro
using Statistics: mean, std
using LatticeProteins: LatticeProteins, parallel_tempering, log_pnat, energies,
    STRUCTURE_A, L, potts, STRUCTURES, StructureSet, shakhnovich_structures,
    structure_index

const LP = LatticeProteins

@testset "parallel tempering" begin
    @testset "local updates keep the cached energies consistent" begin
        rng = Xoshiro(12)
        sites = LP.site_contacts(STRUCTURES)
        boltzmann = STRUCTURES.boltzmann
        centi(seq) = round.(Int32, 100 .* energies(seq))
        seq = rand(rng, 1:20, L)
        chain = LP.Chain(seq, STRUCTURE_A, STRUCTURES)
        @test chain.seq == seq && chain.seq !== seq # the sequence is copied
        @test chain.energies == centi(seq)
        @test chain.loss ≈ -log_pnat(STRUCTURE_A, seq) atol = 1.0e-13
        # All amino acids at endpoint/interior positions, including no-ops.
        for i in (1, 14, 27), aa in 1:20
            LP.mutate_energies!(chain.trial, chain.energies, chain.seq, i, aa, sites)
            mutated = copy(seq)
            mutated[i] = aa
            @test chain.trial == centi(mutated)
            @test LP.folding_loss(chain.trial, STRUCTURE_A, boltzmann) ≈ -log_pnat(STRUCTURE_A, mutated) atol = 1.0e-13
        end
        # Accepted updates and rejected proposals must leave the cache consistent:
        # integer energies never drift, however many incremental updates are applied.
        for beta in (0.0, 1.0, 3000.0)
            accepted = LP.local_steps!(rng, chain, STRUCTURE_A, beta, 200, sites, boltzmann)
            @test 0 <= accepted <= 200
            beta == 0 && @test accepted == 200
            @test chain.energies == centi(chain.seq)
            @test chain.loss == -log_pnat(STRUCTURE_A, chain.seq)
        end
    end

    @testset "serial, threaded, and nested calls agree" begin
        initial = rand(Xoshiro(70), 1:20, L, 4)
        # Include swaps and production in this comparison.
        options = (;
            betas = [0.0, 0.5, 5.0, 100.0], burnin = 99, nseqs = 4,
            thin = 2, local_steps = 2, initial_seqs = initial,
        )
        serial = parallel_tempering(STRUCTURE_A; options..., rng = Xoshiro(71), threaded = false)
        threaded = parallel_tempering(STRUCTURE_A; options..., rng = Xoshiro(71), threaded = true)
        @test serial.rng_mode == :per_temperature
        @test serial == threaded
        nested = Vector{typeof(serial)}(undef, 2)
        Threads.@threads for job in eachindex(nested)
            nested[job] = parallel_tempering(STRUCTURE_A; options..., rng = Xoshiro(71), threaded = true)
        end
        @test all(==(serial), nested)
        @test threaded.final_log_pnat == [log_pnat(STRUCTURE_A, s) for s in eachcol(threaded.final_seqs)]
    end

    @testset "sampling an exactly enumerable ensemble" begin
        # Two competing structures with one contact each, (1,4) and (1,6).
        # Only three residues affect the probability, so enumerate all 20^3
        # possibilities; the other 24 residues contribute a constant multiplicity.
        structures = StructureSet([[1 4], [1 6]])
        states = collect(Iterators.product(1:20, 1:20, 1:20))
        lp = vec([-log1p(exp(LP.MJ[a, b] - LP.MJ[a, c])) for (a, b, c) in states])
        cysteine = vec([a == 1 for (a, b, c) in states])
        betas = [0.0, 1.0, 5.0]
        expected = [sum(lp .* exp.(beta .* lp)) / sum(exp.(beta .* lp)) for beta in betas]
        weights = exp.(last(betas) .* lp)
        expected_cysteine = sum(cysteine .* weights) / sum(weights)

        r = parallel_tempering(
            1; structures, betas, rng = Xoshiro(72), threaded = false,
            burnin = 2000, nseqs = 20_000, thin = 2, local_steps = 5
        )
        # Batch means account for correlation. Each batch spans 5,000 local
        # proposals per replica; use a conservative six-standard-error bound.
        function test_mean(values, expected)
            batches = vec(mean(reshape(values, 500, :); dims = 1))
            se = std(batches) / sqrt(length(batches))
            @test mean(values) ≈ expected atol = 6se
        end
        for i in eachindex(betas)
            test_mean(r.log_pnat[i, :], expected[i])
        end
        test_mean(r.samples[1, :] .== 1, expected_cysteine)
    end

    @testset "selecting the structure ensemble" begin
        S10 = shakhnovich_structures()
        initial = repeat(potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE"); outer = (1, 2))
        full = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 3000.0],
            burnin = 0, nseqs = 0, initial_seqs = initial
        )
        subset_s = structure_index(S10, STRUCTURE_A)
        subset = parallel_tempering(
            subset_s; betas = [0.0, 3000.0],
            burnin = 0, nseqs = 0, initial_seqs = initial, structures = S10
        )
        @test full.final_log_pnat == fill(log_pnat(STRUCTURE_A, initial[:, 1]), 2)
        @test subset.final_log_pnat == fill(log_pnat(subset_s, initial[:, 1]; structures = S10), 2)
        @test all(full.final_log_pnat .< subset.final_log_pnat)
        # Exercise subset mutations, not just initialization.
        sampled = parallel_tempering(
            subset_s; betas = [0.0, 3000.0], burnin = 99,
            nseqs = 2, thin = 1, local_steps = 1, initial_seqs = initial,
            structures = S10, rng = Xoshiro(8)
        )
        @test sampled.final_log_pnat == [
            log_pnat(subset_s, s; structures = S10)
                for s in eachcol(sampled.final_seqs)
        ]
        @test sampled.log_pnat[end, :] == [
            log_pnat(subset_s, s; structures = S10)
                for s in eachcol(sampled.samples)
        ]
        @test_throws ArgumentError parallel_tempering(10_001; betas = [0.0, 1.0], structures = S10)
    end

    @testset "swap detailed balance and state identity" begin
        # Compare the acceptance ratio to the independently computed product target.
        for (bi, bj, ei, ej) in (
                (0.0, 2.0, 0.1, 0.7), (1.0, 3.0, 0.7, 0.1),
                (5.0, 3000.0, 0.02, 0.001),
            )
            logr = LP.swap_logratio(bi, bj, ei, ej)
            before = -bi * ei - bj * ej
            after = -bi * ej - bj * ei
            @test logr ≈ after - before
            @test min(0, logr) - min(0, -logr) ≈ after - before
        end
        # Exact two-state exchange problem: the stationary probability of putting
        # the lower-loss state at beta=2 is logistic(2). Detects a sign error.
        rng = Xoshiro(4)
        chains = [LP.Chain(fill(i, L), STRUCTURE_A, STRUCTURES) for i in 1:2]
        chains[1].loss, chains[2].loss = 0.0, 1.0
        originals = copy(chains)
        ids = [1, 2]
        cold_low = 0
        intact = true
        for _ in 1:20_000
            LP.tempering_swap!(rng, chains, ids, [0.0, 2.0], 1)
            cold_low += ids[2] == 1
            # chains and their identities move together
            intact &= chains[1] === originals[ids[1]] && chains[2] === originals[ids[2]]
        end
        @test intact
        @test abs(cold_low / 20_000 - 1 / (1 + exp(-2))) < 0.02
        # Round trips require a visit to both ends, and repeated endpoint residence
        # must never count as additional trips.
        ids = [1, 2]
        stages, trips = zeros(Int, 2), zeros(Int, 2)
        LP.record_endpoints!(stages, trips, ids)
        reverse!(ids)
        LP.record_endpoints!(stages, trips, ids)
        @test trips == [0, 0]
        reverse!(ids)
        LP.record_endpoints!(stages, trips, ids)
        LP.record_endpoints!(stages, trips, ids)
        @test trips == [1, 0]
    end

    @testset "stationary distribution against exact enumeration" begin
        # A deliberately small model: fold 1 contacts sites 1--4, fold 2
        # contacts sites 1--6. Only three residues affect the energy, so the
        # other 24 sites contribute the same factor 20^24 to every state.
        # Enumerate all 20^3 states independently of the sampler/loss kernels.
        structures = StructureSet([[1 4], [1 6]])
        beta = 3.0
        exact = zeros(20, 3)
        normalization = 0.0
        exact_pnat = 0.0
        for a in 1:20, b in 1:20, c in 1:20
            p = 1 / (1 + exp(LP.MJ[a, b] - LP.MJ[a, c]))
            weight = p^beta
            normalization += weight
            exact_pnat += weight * p
            exact[a, 1] += weight
            exact[b, 2] += weight
            exact[c, 3] += weight
        end
        exact ./= normalization
        exact_pnat /= normalization
        r = parallel_tempering(
            1; structures, betas = [0.0, 0.5, 1.0, 2.0, beta],
            burnin = 2000, nseqs = 40_000, thin = 2, local_steps = 5, rng = Xoshiro(20261030),
            threaded = false
        )
        observed = zeros(20, 3)
        for sample in eachcol(r.samples), (j, site) in enumerate((1, 4, 6))
            observed[sample[site], j] += 1 / size(r.samples, 2)
        end
        @test maximum(abs, observed .- exact) < 0.015
        @test abs(sum(exp, r.log_pnat[end, :]) / size(r.samples, 2) - exact_pnat) < 0.01
        @test sum(r.diagnostics.round_trips) > 100
    end

    @testset "sampling, accounting, and reproducibility" begin
        initial = rand(Xoshiro(2), 1:20, L, 4)
        initial_copy = copy(initial)
        options = (;
            betas = [0.0, 1.0, 20.0, 3000.0], burnin = 3, nseqs = 4,
            thin = 2, local_steps = 2, initial_seqs = initial,
        )
        a = parallel_tempering(STRUCTURE_A; options..., rng = Xoshiro(42))
        b = parallel_tempering(STRUCTURE_A; options..., rng = Xoshiro(42))
        @test a == b
        @test initial == initial_copy
        @test size(a.samples) == (L, 4)
        @test all(x -> 1 <= x <= 20, a.samples)
        @test a.burnin_diagnostics.local_attempts == fill(6, 4)
        @test a.diagnostics.local_attempts == fill(16, 4)
        @test a.burnin_diagnostics.swap_attempts == [2, 1, 2]
        @test a.diagnostics.swap_attempts == [4, 4, 4]
        @test all(a.diagnostics.local_accepts .<= a.diagnostics.local_attempts)
        @test a.diagnostics.local_accepts[1] == 16 # beta=0 accepts every proposal
        @test all(a.diagnostics.swap_accepts .<= a.diagnostics.swap_attempts)
        @test all(sort(ids) == 1:4 for ids in eachcol(a.replica_ids))
        @test a.log_pnat[end, :] == [log_pnat(STRUCTURE_A, s) for s in eachcol(a.samples)]
        @test a.final_log_pnat == [log_pnat(STRUCTURE_A, s) for s in eachcol(a.final_seqs)]
        @test a.samples[:, end] == a.final_seqs[:, end]

        # The two-replica swap schedule attempts its only edge every round.
        c = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 0.01], burnin = 99,
            nseqs = 2, thin = 1, local_steps = 1, rng = Xoshiro(1)
        )
        @test c.burnin_diagnostics.swap_attempts == [99]
        @test c.diagnostics.swap_attempts == [2]
        @test c.final_log_pnat == [log_pnat(STRUCTURE_A, s) for s in eachcol(c.final_seqs)]
        empty = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 1.0], nseqs = 0,
            burnin = 0, initial_seqs = initial[:, 1:2]
        )
        @test size(empty.samples) == (L, 0)
        @test empty.final_seqs == initial[:, 1:2]
        @test all(iszero, empty.diagnostics.swap_attempts)

        # Checkpoint boundaries after an even number of rounds preserve the
        # swap schedule and RNG stream of a continuous run.
        # An even ladder with thin=2 also permits reconstructing every endpoint
        # change from saved IDs, as used by the long-running figure driver.
        whole = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 0.01, 0.02, 0.03],
            burnin = 0, nseqs = 100, thin = 2, local_steps = 1, rng = Xoshiro(91),
            threaded = false, rng_mode = :shared
        )
        rng = Xoshiro(91)
        first = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 0.01, 0.02, 0.03],
            burnin = 0, nseqs = 27, thin = 2, local_steps = 1, rng,
            threaded = false, rng_mode = :shared
        )
        second = parallel_tempering(
            STRUCTURE_A; betas = [0.0, 0.01, 0.02, 0.03],
            burnin = 0, nseqs = 73, thin = 2, local_steps = 1, rng,
            initial_seqs = first.final_seqs, threaded = false, rng_mode = :shared
        )
        @test whole.rng_mode == :shared
        @test whole.samples == hcat(first.samples, second.samples)
        @test whole.log_pnat == hcat(first.log_pnat, second.log_pnat)
        @test whole.final_seqs == second.final_seqs
        ids = hcat(first.replica_ids, first.replica_ids[:, end][second.replica_ids])
        @test whole.replica_ids == ids
        stages, trips = [1, 0, 0, 0], zeros(Int, 4)
        for column in eachcol(ids)
            hot, cold = column[1], column[end]
            stages[hot] == 2 && (trips[hot] += 1)
            stages[hot] = 1
            stages[cold] == 1 && (stages[cold] = 2)
        end
        @test sum(trips) > 0
        @test trips == whole.diagnostics.round_trips
    end

    @testset "invalid arguments" begin
        @test_throws ArgumentError parallel_tempering(
            STRUCTURE_A; betas = [0.0, 1.0],
            rng_mode = :invalid
        )
        @test_throws ArgumentError parallel_tempering(
            STRUCTURE_A; betas = [0.0, 1.0],
            rng_mode = :shared, threaded = true
        )
        for betas in (
                Float64[], [1.0], [1.0, 0.0], [0.0, 0.0], [-1.0, 2.0],
                [0.0, Inf], [0.0, NaN],
            )
            @test_throws ArgumentError parallel_tempering(STRUCTURE_A; betas)
        end
        for args in (
                (; nseqs = -1), (; burnin = -1), (; thin = 0), (; local_steps = 0),
                (; initial_seqs = zeros(Int, L, 2)),
                (; initial_seqs = ones(Int, L, 3)),
            )
            @test_throws ArgumentError parallel_tempering(STRUCTURE_A; betas = [0.0, 1.0], args...)
        end
        @test_throws ArgumentError parallel_tempering(0; betas = [0.0, 1.0])
    end
end
