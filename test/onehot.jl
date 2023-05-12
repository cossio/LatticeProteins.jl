using Test: @test, @testset, @inferred
using LinearAlgebra: I
using LatticeProteins: onehot, potts, AMINO_ACIDS

@testset "onehot" begin
    @test onehot(AMINO_ACIDS) == I
    @test potts(AMINO_ACIDS) == 1:21
    @test onehot(1:21) == onehot(AMINO_ACIDS)
end
