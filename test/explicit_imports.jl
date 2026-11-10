import ExplicitImports
import LatticeProteins
using Test: @testset

@testset "ExplicitImports" begin
    ExplicitImports.test_explicit_imports(LatticeProteins)
end
