using Test: @test, @testset
using LatticeProteins: random_msa, metropolis!, pnat, potts, log_pnat,
    CONTACT_MAP_A, L
using Random: seed!

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

        # With high β and many steps, pnat should increase (designed sequence)
        seed!(123)
        seq3 = rand(1:20, L)
        pnat_before = pnat(CONTACT_MAP_A, seq3)
        metropolis!(seq3, CONTACT_MAP_A; β=100, nsteps=200)
        pnat_after = pnat(CONTACT_MAP_A, seq3)
        @test pnat_after ≥ pnat_before
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

