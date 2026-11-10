using Test: @test, @testset
import LatticeProteins

@testset "public API" begin
    # The public API is the model with all compact structures.
    api = Set(
        [
            :L, :potts, :aaseq, :onehot,
            :enumerate_structures, :contact_map, :STRUCTURE_A, :STRUCTURE_B, :STRUCTURE_C, :STRUCTURE_D,
            :MJ, :energy, :pnat, :log_pnat,
            :metropolis!, :sample_msa, :parallel_tempering,
        ]
    )
    @test Set(names(LatticeProteins)) == union(api, [:LatticeProteins])
    # The package exports and documents everything that is public, except `L`, which is
    # public but too short a name to export.
    @test all(name -> Base.isexported(LatticeProteins, name), setdiff(api, [:L]))
    @test Base.ispublic(LatticeProteins, :L) && !Base.isexported(LatticeProteins, :L)
    @test all(name -> Docs.hasdoc(LatticeProteins, name), api)
    # Other sets of competing structures, and the reference data sampled with one, are
    # internal.
    internal = [
        :StructureSet, :STRUCTURES, :N_STRUCTURES, :shakhnovich_structures, :structure_index,
        :energies, :N_CONTACTS, :AMINO_ACIDS, :Hugo_MSA, :Hugo_MSA_pnat, :Hugo_MSA_path,
    ]
    @test all(name -> isdefined(LatticeProteins, name) && !Base.ispublic(LatticeProteins, name), internal)
end
