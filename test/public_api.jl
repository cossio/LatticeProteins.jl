using Test: @test, @testset
import LatticeProteins

@testset "public API" begin
    # The public API is the model with all compact structures.
    api = Set(
        [
            :L, :potts, :aaseq, :onehot,
            :N_STRUCTURES, :contact_map, :STRUCTURE_A, :STRUCTURE_B, :STRUCTURE_C, :STRUCTURE_D,
            :MJ, :energy, :pnat, :log_pnat,
            :metropolis!, :random_msa, :parallel_tempering,
        ]
    )
    @test Set(names(LatticeProteins)) == union(api, [:LatticeProteins])
    # the package exports nothing, and documents everything that is public
    @test all(name -> !Base.isexported(LatticeProteins, name), api)
    @test all(name -> Docs.hasdoc(LatticeProteins, name), api)
    # Other sets of competing structures, and the reference data sampled with one, are
    # internal.
    internal = [
        :StructureSet, :STRUCTURES, :shakhnovich_structures, :structure_index, :energies,
        :N_CONTACTS, :AMINO_ACIDS, :Hugo_MSA, :Hugo_MSA_pnat, :Hugo_MSA_path,
    ]
    @test all(name -> isdefined(LatticeProteins, name) && !Base.ispublic(LatticeProteins, name), internal)
    @test !isdefined(LatticeProteins, :enumerate_structures)
end
