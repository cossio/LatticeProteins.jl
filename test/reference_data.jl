using Test: @test, @testset, @test_throws
using LatticeProteins: pnat, potts, onehot, L, structure_index, shakhnovich_structures,
    STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D,
    Hugo_MSA, Hugo_MSA_pnat, Hugo_MSA_path

@testset "Hugo_MSA_path" begin
    for which in (:A, :B, :C, :D)
        path = Hugo_MSA_path(which)
        @test path isa String
        @test isfile(path)
    end
    @test_throws ArgumentError Hugo_MSA_path(:Z)
    @test_throws ArgumentError Hugo_MSA(:Z)
    @test_throws ArgumentError Hugo_MSA_pnat(:Z)
end

@testset "Hugo_MSA" begin
    S10 = shakhnovich_structures()
    for (which, s) in ((:A, STRUCTURE_A), (:B, STRUCTURE_B), (:C, STRUCTURE_C), (:D, STRUCTURE_D))
        msa = Hugo_MSA(which)
        pnats = Hugo_MSA_pnat(which)

        @test msa isa Vector{String}
        @test length(msa) > 0
        @test all(length.(msa) .== L)

        @test pnats isa Vector{Float64}
        @test length(pnats) == length(msa)
        @test all(0 .< pnats .≤ 1)

        # the stored folding probabilities are those of the sequences, against the subset
        # of 10,000 structures
        s10 = structure_index(S10, s)
        @test all(pnat(s10, potts(msa[i]); structures = S10) ≈ pnats[i] for i in 1:50)

        X = onehot(msa)
        @test size(X) == (21, L, length(msa))
        @test all(sum(X; dims = 1) .== 1)
    end
end
