using Test: @test, @testset
import LatticeProteins

@testset "public API" begin
    # The public API is the model with all compact structures.
    api = Set(
        [
            :potts, :aaseq, :onehot,
            :enumerate_structures, :contact_map, :STRUCTURE_A, :STRUCTURE_B, :STRUCTURE_C, :STRUCTURE_D,
            :MJ, :energy, :pnat, :log_pnat,
            :sample_msa, :Metropolis, :ParallelTempering, :metropolis!,
        ]
    )
    @test Set(names(LatticeProteins)) == union(api, [:LatticeProteins])
    # the package exports and documents everything that is public
    @test all(name -> Base.isexported(LatticeProteins, name), api)
    @test all(name -> Docs.hasdoc(LatticeProteins, name), api)
    # Other sets of competing structures, and the reference data sampled with one, are
    # internal.
    internal = [
        :StructureSet, :STRUCTURES, :N_STRUCTURES, :shakhnovich_structures, :structure_index,
        :energies, :N_CONTACTS, :AMINO_ACIDS, :L, :Hugo_MSA, :Hugo_MSA_pnat, :Hugo_MSA_path,
    ]
    @test all(name -> isdefined(LatticeProteins, name) && !Base.ispublic(LatticeProteins, name), internal)
end
