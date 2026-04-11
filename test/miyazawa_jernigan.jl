using Test: @test, @testset
using LatticeProteins: MJ, load_miyazawa_jernigan_matrix, load_miyazawa_jernigan_aminoacids
using LinearAlgebra: issymmetric

@testset "miyazawa_jernigan" begin
    @testset "load_miyazawa_jernigan_matrix" begin
        M = load_miyazawa_jernigan_matrix()
        @test size(M) == (20, 20)
        @test eltype(M) == Float64

        # Matrix should be symmetric (interaction between i,j == j,i)
        @test issymmetric(M)

        # All values should be negative (attractive interactions)
        @test all(M .< 0)

        # Loading with a different type should work
        M32 = load_miyazawa_jernigan_matrix(Float32)
        @test size(M32) == (20, 20)
        @test eltype(M32) == Float32
        @test M32 ≈ Float32.(M)
    end

    @testset "load_miyazawa_jernigan_aminoacids" begin
        aas = load_miyazawa_jernigan_aminoacids()
        @test isa(aas, String)
        @test length(aas) == 20
        @test aas == "CMFILVWYAGTSNQDEHRKP"
    end

    @testset "MJ constant" begin
        @test size(MJ) == (20, 20)
        @test MJ == load_miyazawa_jernigan_matrix()
    end
end
