using Test: @test, @testset, @test_throws
using LatticeProteins: pnat, potts, onehot, L, structure_index, JACQUIN_STRUCTURES,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    Hugo_MSA, Hugo_MSA_pnat, Hugo_MSA_path
using Statistics: mean

# The pnat values in the Hugo MSA files were computed against the 10,000 reference
# structures of Jacquin et al. (2016).
jacquin_pnat(cm, seq) = pnat(structure_index(JACQUIN_STRUCTURES, cm), seq; structures=JACQUIN_STRUCTURES)

@testset "Hugo_MSAs" begin
    @testset "Hugo_MSA_path" begin
        for sym in (:A, :B, :C, :D)
            path = Hugo_MSA_path(sym)
            @test isa(path, String)
            @test isfile(path)
        end
        @test_throws ArgumentError Hugo_MSA_path(:Z)
    end

    @testset "Hugo_MSA loading" begin
        for (sym, cm) in [(:A, CONTACT_MAP_A), (:B, CONTACT_MAP_B),
                          (:C, CONTACT_MAP_C), (:D, CONTACT_MAP_D)]
            msa = Hugo_MSA(sym)
            pnats = Hugo_MSA_pnat(sym)

            @test isa(msa, Vector{String})
            @test length(msa) > 0
            @test all(length.(msa) .== L)

            @test isa(pnats, Vector{Float64})
            @test length(pnats) == length(msa)
            @test all(0 .< pnats .≤ 1)

            # Verify pnat consistency for first entry
            @test jacquin_pnat(cm, potts(msa[1])) ≈ pnats[1]
        end
    end

    @testset "Hugo_MSA pnat consistency - A" begin
        msa_A = Hugo_MSA(:A)
        pnat_A = Hugo_MSA_pnat(:A)
        a10 = structure_index(JACQUIN_STRUCTURES, CONTACT_MAP_A)
        for i = 1:100
            @test pnat(a10, potts(msa_A[i]); structures=JACQUIN_STRUCTURES) ≈ pnat_A[i]
        end
    end

    @testset "Hugo_MSA pnat consistency - B" begin
        msa_B = Hugo_MSA(:B)
        pnat_B = Hugo_MSA_pnat(:B)
        b10 = structure_index(JACQUIN_STRUCTURES, CONTACT_MAP_B)
        for i = 1:100
            @test pnat(b10, potts(msa_B[i]); structures=JACQUIN_STRUCTURES) ≈ pnat_B[i]
        end
    end

    @testset "Hugo_MSA onehot" begin
        msa_A = Hugo_MSA(:A)
        X = onehot(msa_A)
        @test size(X, 1) == 21
        @test size(X, 2) == L
        @test size(X, 3) == length(msa_A)

        # Mean one-hot should have values in [0, 1]
        m = mean(X; dims=3)
        @test all(0 .≤ m .≤ 1)
    end
end
