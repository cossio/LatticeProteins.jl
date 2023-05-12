using Test: @test
using LatticeProteins: pnat, potts

@test pnat(6685, potts("WRFQDPGRCEEEDERMCGPEMCRCRQK")) ≈ 0.995663 rtol=1e-6
@test pnat(6685, potts("LEITKLEAVDDDRARFCQENLNGCKEQ")) ≈ 0.995402 rtol=1e-6
@test pnat(6685, potts("FELHKLKGFDGSTASLIKESFNNWKED")) ≈ 0.994070 rtol=1e-6

@test pnat(848, potts("RNISIMHVVVDLLTELTIMVAEYPQFP")) ≈ 1.45048e-07 rtol=1e-5
@test pnat(848, potts("RMMFFVPEMSKMGDFEHFKPEEYVAKW")) ≈ 0.988368 rtol=1e-6
@test pnat(848, potts("TGWLEWDDFRRKETVDDHRECPIKECM")) ≈ 0.993085 rtol=1e-6

@test pnat(625, potts("HCEVAVLKVRIEDQYCDCHDRVHFCTE")) ≈ 0.989353 rtol=1e-6
