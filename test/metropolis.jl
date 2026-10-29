using Test: @test, @testset
using LatticeProteins: random_msa, metropolis!, pnat, potts, log_pnat, energies,
    SiteContacts, mutate_energies!, CONTACT_MAP_A, CONTACT_MAPS, L
using Random: seed!, randexp

# Reference implementation that recomputes log_pnat from scratch at every step, with the
# same sequence of random draws as metropolis!.
function naive_metropolis!(seq, cm; β, nsteps, structures)
    E = -log_pnat(cm, seq; structures)
    for _ = 1:nsteps
        i = rand(1:length(seq))
        old_aa = seq[i]
        seq[i] = rand(1:20)
        E1 = -log_pnat(cm, seq; structures)
        if β * (E1 - E) < randexp()
            E = E1
        else
            seq[i] = old_aa
        end
    end
    return seq
end

@testset "metropolis" begin
    @testset "metropolis!" begin
        seed!(42)
        seq = rand(1:20, L)

        # metropolis! returns the modified sequence
        result = metropolis!(seq, CONTACT_MAP_A; β=1, nsteps=1)
        @test result === seq  # returns the same array (in-place)
        @test length(seq) == L
        @test all(1 .≤ seq .≤ 20)

        # With nsteps=0, sequence should not change
        seq2 = rand(1:20, L)
        seq2_copy = copy(seq2)
        metropolis!(seq2, CONTACT_MAP_A; β=1, nsteps=0)
        @test seq2 == seq2_copy

        # With high β and many steps, the result should remain a valid sequence
        # and produce a valid pnat value, but the final state is not guaranteed
        # to have higher pnat than the initial state under Metropolis sampling.
        seed!(123)
        seq3 = rand(1:20, L)
        pnat_before = pnat(CONTACT_MAP_A, seq3)
        metropolis!(seq3, CONTACT_MAP_A; β=100, nsteps=200)
        pnat_after = pnat(CONTACT_MAP_A, seq3)
        @test 0 ≤ pnat_before ≤ 1
        @test 0 ≤ pnat_after ≤ 1
        @test length(seq3) == L
        @test all(1 .≤ seq3 .≤ 20)
    end

    @testset "incremental energies match recomputation" begin
        seed!(6)
        sites = SiteContacts(CONTACT_MAPS)
        seq = rand(1:20, L)
        E = energies(seq)
        E1 = similar(E)
        for _ in 1:50
            i = rand(1:L)
            new_aa = rand(1:20)
            mutate_energies!(E1, E, seq, i, new_aa, sites)
            seq[i] = new_aa
            @test E1 ≈ energies(seq) atol=1e-10
            E, E1 = E1, E
        end

        seed!(7)
        seq0 = rand(1:20, L)
        for β in (1, 10, 100)
            seed!(8)
            seq1 = metropolis!(copy(seq0), CONTACT_MAP_A; β, nsteps=300)
            seed!(8)
            seq2 = naive_metropolis!(copy(seq0), CONTACT_MAP_A; β, nsteps=300, structures=CONTACT_MAPS)
            @test seq1 == seq2
            @test seq1 != seq0
        end

        seed!(5)
        msa = random_msa(CONTACT_MAP_A; nseqs=3, β=10, nsteps=50)
        seed!(5)
        seqs = rand(1:20, L, 3)
        for n in 1:3
            naive_metropolis!(view(seqs, :, n), CONTACT_MAP_A; β=10, nsteps=50, structures=CONTACT_MAPS)
        end
        @test msa == seqs
    end

    @testset "random_msa" begin
        seed!(42)
        msa = random_msa(CONTACT_MAP_A; nseqs=5, β=1, nsteps=10)
        @test size(msa) == (L, 5)
        @test eltype(msa) == Int
        @test all(1 .≤ msa .≤ 20)

        # Single sequence
        msa1 = random_msa(CONTACT_MAP_A; nseqs=1, β=1, nsteps=10)
        @test size(msa1) == (L, 1)
    end
end

