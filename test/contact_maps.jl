using Test: @test, @testset
using LatticeProteins: CONTACT_MAPS, N_STRUCTURES, N_CONTACTS, L,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    load_contact_maps

@testset "contact_maps" begin
    @testset "constants" begin
        @test L == 27
        @test N_CONTACTS == 28
        @test N_STRUCTURES == 10_000
        @test CONTACT_MAP_A == 625
        @test CONTACT_MAP_B == 848
        @test CONTACT_MAP_C == 6683
        @test CONTACT_MAP_D == 6685
    end

    @testset "CONTACT_MAPS shape and values" begin
        @test size(CONTACT_MAPS) == (N_CONTACTS, N_STRUCTURES, 2)
        @test eltype(CONTACT_MAPS) == Int

        # Contact indices must be in valid range [1, L]
        @test all(1 .≤ CONTACT_MAPS .≤ L)

        # Each contact should be a pair of distinct positions
        for s in [CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D]
            for k in 1:N_CONTACTS
                i, j = CONTACT_MAPS[k, s, 1], CONTACT_MAPS[k, s, 2]
                @test i != j  # contacts are between different positions
            end
        end
    end

    @testset "load_contact_maps" begin
        cm = load_contact_maps()
        @test cm == CONTACT_MAPS
        @test size(cm) == (N_CONTACTS, N_STRUCTURES, 2)
    end
end
