using Test: @test, @testset
using LatticeProteins: STRUCTURES, SHAKHNOVICH_STRUCTURES, N_STRUCTURES, N_CONTACTS, L,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    load_contact_maps, structure_index

@testset "contact_maps" begin
    @testset "constants" begin
        @test L == 27
        @test N_CONTACTS == 28
        @test N_STRUCTURES == 103_346
        @test size(STRUCTURES) == (N_CONTACTS, N_STRUCTURES, 2)
        @test size(SHAKHNOVICH_STRUCTURES) == (N_CONTACTS, 10_000, 2)
        # the reference structures of Jacquin et al. (2016) are structures 625, 848, 6683
        # and 6685 of the reference set
        @test structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A) == 625
        @test structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_B) == 848
        @test structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_C) == 6683
        @test structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_D) == 6685
        @test allunique((CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D))
        @test all(1 .≤ (CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D) .≤ N_STRUCTURES)
    end

    @testset "shape and values" begin
        for S in (STRUCTURES, SHAKHNOVICH_STRUCTURES)
            @test eltype(S) == Int
            # Contact indices must be in valid range [1, L]
            @test all(1 .≤ S .≤ L)
            # Each contact should be a pair of distinct positions
            @test all(S[:, :, 1] .!= S[:, :, 2])
        end
        for s in (CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D)
            for k in 1:N_CONTACTS
                i, j = STRUCTURES[k, s, 1], STRUCTURES[k, s, 2]
                @test i != j
            end
        end
    end

    @testset "load_contact_maps" begin
        cm = load_contact_maps()
        @test cm == SHAKHNOVICH_STRUCTURES
        @test size(cm) == (N_CONTACTS, 10_000, 2)
    end
end
