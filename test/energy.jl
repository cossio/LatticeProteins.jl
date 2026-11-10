using Test: @test, @testset, @test_throws
using LatticeProteins: LatticeProteins, pnat, log_pnat, energy, energies, potts, contact_map,
    structure_index, StructureSet, STRUCTURES, shakhnovich_structures, MJ, N_STRUCTURES,
    STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D

# the subset of 10,000 structures of Jacquin et al. (2016), and the indices of the reference
# structures in it
const S10 = shakhnovich_structures()
const A10 = structure_index(S10, STRUCTURE_A)
const B10 = structure_index(S10, STRUCTURE_B)
const C10 = structure_index(S10, STRUCTURE_C)
const D10 = structure_index(S10, STRUCTURE_D)
jacquin_pnat(s, seq) = pnat(s, potts(seq); structures = S10)

@testset "stable folding loss" begin
    # An independent high-precision sum checks the loss over the whole range of energy
    # differences of the Boltzmann table, including near-native sequences. The tolerance
    # allows for the roundoff of summing a thousand Boltzmann factors in any order.
    boltzmann = STRUCTURES.boltzmann
    span = length(boltzmann) - 1
    @test boltzmann == exp.(-(0:span) ./ 100)
    for E in (
            Int32[-20000, -20000 + span], Int32[-500, -500], Int32[1000, -1000, 0],
            Int32[-9000; -8999; fill(-8000, 1000)],
        )
        for s in eachindex(E)
            exact = Float64(log(sum(exp((big(E[s]) - big(e)) / 100) for e in E)))
            @test LatticeProteins.folding_loss(E, s, boltzmann) ≈ exact rtol = 1.0e-12 atol = 1.0e-15
        end
    end
end

@testset "energy" begin
    seq = potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")
    @test energy(STRUCTURE_D, seq) isa Float64

    # the energy is the sum of the MJ interactions over the contacts of the structure
    for s in (STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D, 1, N_STRUCTURES)
        cm = contact_map(STRUCTURES, s)
        @test energy(s, seq) ≈ sum(MJ[seq[cm[k, 1]], seq[cm[k, 2]]] for k in axes(cm, 1)) rtol = 1.0e-14
    end
    @test energy(STRUCTURE_A, seq) != energy(STRUCTURE_D, seq)
    # the energy does not depend on the set the structure is taken from
    @test energy(D10, seq; structures = S10) == energy(STRUCTURE_D, seq)
    # sequences of other integer types
    @test energy(STRUCTURE_A, Int8.(seq)) == energy(STRUCTURE_A, seq)

    # sequence length must be 27
    @test_throws ArgumentError energy(1, [1, 2, 3])
    # amino acid indices must be in 1:20
    @test_throws ArgumentError energy(1, fill(21, 27))
    @test_throws ArgumentError energy(1, fill(0, 27))
    # the structure must be in the set
    @test_throws ArgumentError energy(0, seq)
    @test_throws ArgumentError energy(10_001, seq; structures = S10)
end

@testset "energies" begin
    seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
    E = energies(seq)
    @test E isa Vector{Float64}
    @test E == [energy(s, seq) for s in 1:N_STRUCTURES]
    @test energies(seq; structures = S10) == [energy(s, seq; structures = S10) for s in 1:10_000]
    @test_throws ArgumentError energies([1, 2, 3])
end

@testset "log_pnat and pnat" begin
    seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
    lp = log_pnat(STRUCTURE_A, seq)
    @test lp isa Float64
    @test lp < 0
    @test pnat(STRUCTURE_A, seq) == exp(lp)

    # pnat is a probability distribution over the structures of the set
    for structures in (STRUCTURES, S10), seq in (seq, potts("CMFILVWYAGTSNQDEHRKPCMFILVW"))
        E = energies(seq; structures)
        m = maximum(-E)
        p = exp.(-E .- (m + log(sum(exp.(-E .- m)))))
        @test sum(p) ≈ 1
        for s in (1, 2, length(structures) ÷ 2, length(structures))
            @test pnat(s, seq; structures) ≈ p[s]
        end
    end

    # the subset has fewer competitors, so a structure is more probable in it
    @test pnat(STRUCTURE_A, seq) < pnat(A10, seq; structures = S10)
    # the default is the full set
    @test log_pnat(STRUCTURE_A, seq; structures = STRUCTURES) == lp

    # a set of two structures with one contact each
    S = StructureSet([[1 4], [1 6]])
    @test pnat(1, seq; structures = S) ≈ 1 / (1 + exp(MJ[seq[1], seq[4]] - MJ[seq[1], seq[6]]))
    @test pnat(1, seq; structures = S) + pnat(2, seq; structures = S) ≈ 1

    @test_throws ArgumentError log_pnat(1, [1, 2, 3])
    @test_throws ArgumentError log_pnat(0, seq)
    @test_throws ArgumentError pnat(N_STRUCTURES + 1, seq)
end

@testset "pnat - known values" begin
    # Reference values computed against the 10,000 structures of Jacquin et al. (2016)
    # From Lorenzo
    @test jacquin_pnat(6685, "WRFQDPGRCEEEDERMCGPEMCRCRQK") ≈ 0.995663 rtol = 1.0e-6
    @test jacquin_pnat(6685, "LEITKLEAVDDDRARFCQENLNGCKEQ") ≈ 0.995402 rtol = 1.0e-6
    @test jacquin_pnat(6685, "FELHKLKGFDGSTASLIKESFNNWKED") ≈ 0.99407 rtol = 1.0e-6

    @test jacquin_pnat(848, "RNISIMHVVVDLLTELTIMVAEYPQFP") ≈ 1.45048e-7 rtol = 1.0e-5
    @test jacquin_pnat(848, "RMMFFVPEMSKMGDFEHFKPEEYVAKW") ≈ 0.988368 rtol = 1.0e-6
    @test jacquin_pnat(848, "TGWLEWDDFRRKETVDDHRECPIKECM") ≈ 0.993085 rtol = 1.0e-6

    @test jacquin_pnat(625, "HCEVAVLKVRIEDQYCDCHDRVHFCTE") ≈ 0.989353 rtol = 1.0e-6

    # From Hugo
    @test jacquin_pnat(A10, "EKAMPAMDPDMAHEHKKKIRAWMFEGE") ≈ 0.99463333079199023
    @test jacquin_pnat(A10, "EKALNPVDTECWKEYKKKMRPVMLDAR") ≈ 0.99281733038613029

    @test jacquin_pnat(B10, "PDRFIVQCLAQFEHFERDGDKDMRAEC") ≈ 0.9867267580946395
    @test jacquin_pnat(B10, "KEHVECDRIKPKEEVANDEDKDLKCDF") ≈ 0.99837746142507167

    @test jacquin_pnat(C10, "ICEMDADLHCPKFPSFCRDWKAEKYEF") ≈ 0.99411990031853326
    @test jacquin_pnat(C10, "LPELDFKCKWAPFHCLHREMQMCQIKI") ≈ 0.99546655562647968

    @test jacquin_pnat(D10, "FKLTPQKACEGKDKELFQAWFCSCKRD") ≈ 0.98983908221438699
    @test jacquin_pnat(D10, "FEFCKIKKLSNAKRLMCGPGFESVNAE") ≈ 0.99746046068117966
end
