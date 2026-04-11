using Test: @test, @testset, @test_throws
using LatticeProteins: pnat, potts, energy, log_pnat,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    MJ, CONTACT_MAPS, N_CONTACTS

@testset "pnat" begin
    @testset "energy" begin
        # energy should return a Float64
        seq = potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")
        E = energy(CONTACT_MAP_D, seq)
        @test isa(E, Float64)

        # energy for same sequence on different structures should differ
        E_A = energy(CONTACT_MAP_A, seq)
        E_D = energy(CONTACT_MAP_D, seq)
        @test E_A != E_D

        # energy depends on the contact map: swapping contact maps changes energy
        seq2 = potts("HCEVAVLKVRIEDQYCDCHDRVHFCTE")
        @test energy(CONTACT_MAP_A, seq2) != energy(CONTACT_MAP_B, seq2)

        # energy should be deterministic
        @test energy(CONTACT_MAP_A, seq) == energy(CONTACT_MAP_A, seq)

        # sequence length must be 27
        @test_throws AssertionError energy(1, [1, 2, 3])
    end

    @testset "log_pnat" begin
        seq = potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")
        lp = log_pnat(CONTACT_MAP_D, seq)
        @test isa(lp, Float64)
        @test lp <= 0  # log probability must be non-positive
        @test lp ≈ log(pnat(CONTACT_MAP_D, seq))

        # log_pnat for a well-designed sequence should be close to 0
        lp_good = log_pnat(CONTACT_MAP_D, seq)
        @test lp_good > -1  # high pnat means log_pnat close to 0

        # sequence length must be 27
        @test_throws AssertionError log_pnat(1, [1, 2, 3])
    end

    @testset "pnat - known values" begin
        # From Lorenzo
        @test pnat(6685, potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")) ≈ 0.995663 rtol=1e-6
        @test pnat(6685, potts("LEITKLEAVDDDRARFCQENLNGCKEQ")) ≈ 0.995402 rtol=1e-6
        @test pnat(6685, potts("FELHKLKGFDGSTASLIKESFNNWKED")) ≈ 0.994070 rtol=1e-6

        @test pnat(848, potts("RNISIMHVVVDLLTELTIMVAEYPQFP")) ≈ 1.45048e-07 rtol=1e-5
        @test pnat(848, potts("RMMFFVPEMSKMGDFEHFKPEEYVAKW")) ≈ 0.988368 rtol=1e-6
        @test pnat(848, potts("TGWLEWDDFRRKETVDDHRECPIKECM")) ≈ 0.993085 rtol=1e-6

        @test pnat(625, potts("HCEVAVLKVRIEDQYCDCHDRVHFCTE")) ≈ 0.989353 rtol=1e-6

        # From Hugo
        @test pnat(CONTACT_MAP_A, potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")) ≈ 0.99463333079199023
        @test pnat(CONTACT_MAP_A, potts("EKALNPVDTECWKEYKKKMRPVMLDAR")) ≈ 0.99281733038613029

        @test pnat(CONTACT_MAP_B, potts("PDRFIVQCLAQFEHFERDGDKDMRAEC")) ≈ 0.98672675809463950
        @test pnat(CONTACT_MAP_B, potts("KEHVECDRIKPKEEVANDEDKDLKCDF")) ≈ 0.99837746142507167

        @test pnat(CONTACT_MAP_C, potts("ICEMDADLHCPKFPSFCRDWKAEKYEF")) ≈ 0.99411990031853326
        @test pnat(CONTACT_MAP_C, potts("LPELDFKCKWAPFHCLHREMQMCQIKI")) ≈ 0.99546655562647968

        @test pnat(CONTACT_MAP_D, potts("FKLTPQKACEGKDKELFQAWFCSCKRD")) ≈ 0.98983908221438699
        @test pnat(CONTACT_MAP_D, potts("FEFCKIKKLSNAKRLMCGPGFESVNAE")) ≈ 0.99746046068117966
    end

    @testset "pnat - properties" begin
        seq = potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")
        p = pnat(CONTACT_MAP_D, seq)
        @test 0 < p ≤ 1  # probability must be in (0, 1]

        # pnat should be consistent with exp(log_pnat)
        @test p ≈ exp(log_pnat(CONTACT_MAP_D, seq))
    end
end
