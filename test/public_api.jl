using Test: @test, @testset
import LatticeProteins

@testset "public API" begin
    api = Set(
        [
            :AMINO_ACIDS, :onehot, :potts, :aaseq,
            :L, :N_CONTACTS, :N_STRUCTURES,
            :StructureSet, :STRUCTURES, :shakhnovich_structures, :enumerate_structures,
            :contact_map, :structure_index,
            :STRUCTURE_A, :STRUCTURE_B, :STRUCTURE_C, :STRUCTURE_D,
            :MJ, :energy, :energies, :pnat, :log_pnat,
            :metropolis!, :random_msa, :parallel_tempering,
            :Hugo_MSA, :Hugo_MSA_pnat, :Hugo_MSA_path,
        ]
    )
    @test Set(names(LatticeProteins)) == union(api, [:LatticeProteins])
    # the package exports nothing, and documents everything that is public
    @test all(name -> !Base.isexported(LatticeProteins, name), api)
    @test all(name -> Docs.hasdoc(LatticeProteins, name), api)
end
