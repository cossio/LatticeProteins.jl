using Test: @test, @testset, @test_throws
using LatticeProteins: pnat, potts, onehot, L,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    Hugo_MSA, Hugo_MSA_pnat, Hugo_MSA_path
using Statistics: mean

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
        end
    end

    @testset "Hugo_MSA pnat consistency - A" begin
        msa_A = Hugo_MSA(:A)
        pnat_A = Hugo_MSA_pnat(:A)
        for i = 1:100
            @test pnat(CONTACT_MAP_A, potts(msa_A[i])) ≈ pnat_A[i]
        end
    end

    @testset "Hugo_MSA pnat consistency - B" begin
        msa_B = Hugo_MSA(:B)
        pnat_B = Hugo_MSA_pnat(:B)
        for i = 1:100
            @test pnat(CONTACT_MAP_B, potts(msa_B[i])) ≈ pnat_B[i]
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
