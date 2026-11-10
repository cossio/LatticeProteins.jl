using Test: @test, @testset
using LatticeProteins: LatticeProteins, MJ
using LinearAlgebra: issymmetric

@testset "miyazawa_jernigan" begin
    @test MJ isa Matrix{Float64}
    @test size(MJ) == (20, 20)
    # the interaction between a, b is the same as between b, a
    @test issymmetric(MJ)
    # all interactions are attractive
    @test all(MJ .< 0)
    @test MJ[1, 1] == -5.44 # cysteine-cysteine
    # the energies are multiples of 0.01, which the package computes with
    @test LatticeProteins.MJ_INT isa Matrix{Int32}
    @test LatticeProteins.MJ_INT ./ 100 == MJ
end
