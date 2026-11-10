import Aqua
import LatticeProteins
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(LatticeProteins)
end
