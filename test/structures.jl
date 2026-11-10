using Test: @test, @testset, @test_throws
using LatticeProteins: LatticeProteins, StructureSet, STRUCTURES, shakhnovich_structures,
    contact_map, structure_index, N_STRUCTURES, N_CONTACTS, L,
    STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D

const S10 = shakhnovich_structures()

contact_key(cm) = sort([minmax(cm[k, 1], cm[k, 2]) for k in axes(cm, 1)])
# chain reversal maps contact (i, j) to (L + 1 - j, L + 1 - i)
reverse_structure(cm) = [L + 1 - cm[k, 3 - c] for k in axes(cm, 1), c in 1:2]

@testset "constants" begin
    @test L == 27
    @test N_CONTACTS == 28
    @test N_STRUCTURES == length(STRUCTURES) == 103_346
    @test length(S10) == 10_000
    @test shakhnovich_structures() === S10 # built once
    @test STRUCTURES isa StructureSet && S10 isa StructureSet
    @test sprint(show, STRUCTURES) == "StructureSet with 103346 structures of 28 contacts"
end

@testset "enumeration of the compact structures" begin
    # Shakhnovich & Gutin (1990): 103,346 compact conformations unrelated by symmetry
    @test size(STRUCTURES.contacts) == (2, N_CONTACTS, 103_346)
    @test eltype(STRUCTURES.contacts) == Int8
    @test all(1 .≤ STRUCTURES.contacts .≤ L)

    I = STRUCTURES.contacts[1, :, :]
    J = STRUCTURES.contacts[2, :, :]
    @test all(I .< J)
    # the lattice is bipartite, so contacting residues have opposite parity, and
    # residues adjacent along the chain are bonded, not in contact
    @test all(isodd.(J .- I))
    @test all(J .- I .≥ 3)
    # contacts are sorted within each structure
    @test all(issorted(collect(zip(I[:, s], J[:, s]))) for s in 1:N_STRUCTURES)

    keys = Set(contact_key(contact_map(s)) for s in 1:N_STRUCTURES)
    @test length(keys) == N_STRUCTURES # no duplicates
    # the set is closed under chain reversal; a few structures are their own reversal
    reversed = [contact_key(reverse_structure(contact_map(s))) for s in 1:N_STRUCTURES]
    @test all(in(keys), reversed)
    @test count(reversed[s] == contact_key(contact_map(s)) for s in 1:N_STRUCTURES) == 62
    # the 10,000 reference structures are among them
    @test all(contact_key(contact_map(s; structures = S10)) in keys for s in 1:10_000)
    # the enumeration is deterministic
    @test LatticeProteins.compact_contact_maps() == STRUCTURES.contacts
end

@testset "contact_map" begin
    cm = contact_map(STRUCTURE_A)
    @test cm isa Matrix{Int}
    @test size(cm) == (N_CONTACTS, 2)
    @test cm == [STRUCTURES.contacts[c, k, STRUCTURE_A] for k in 1:N_CONTACTS, c in 1:2]
    @test contact_map(STRUCTURE_A; structures = STRUCTURES) == cm
    @test_throws ArgumentError contact_map(0)
    @test_throws ArgumentError contact_map(N_STRUCTURES + 1)
    @test_throws ArgumentError contact_map(10_001; structures = S10)
end

@testset "reference structures" begin
    # structures A, B, C, D of Jacquin et al. (2016) are structures 625, 848, 6683 and 6685
    # of the subset of 10,000 structures
    for (s, s10) in zip((STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D), (625, 848, 6683, 6685))
        @test structure_index(S10, s) == s10
        @test structure_index(S10, contact_map(s)) == s10
        @test structure_index(STRUCTURES, s10; from = S10) == s
        @test structure_index(STRUCTURES, s) == s
        @test contact_key(contact_map(s)) == contact_key(contact_map(s10; structures = S10))
    end
end

@testset "structure_index" begin
    @test structure_index(S10, contact_map(7; structures = S10)) == 7
    @test structure_index(S10, 7; from = S10) == 7
    # most compact structures are not in the reference set
    @test any(structure_index(S10, s) === nothing for s in 1:100)
    # the order of the contacts and the order within each pair are ignored
    cm = contact_map(STRUCTURE_A)
    @test structure_index(STRUCTURES, reverse(cm; dims = 1)) == STRUCTURE_A
    @test structure_index(STRUCTURES, cm[:, [2, 1]]) == STRUCTURE_A
    # a contact map with a different number of contacts is never found
    @test structure_index(STRUCTURES, cm[1:27, :]) === nothing
    # the chain reversal of structure A is a different structure, absent from the reference set
    @test structure_index(S10, reverse_structure(cm)) === nothing
    @test structure_index(STRUCTURES, reverse_structure(cm)) ∉ (nothing, STRUCTURE_A)
    @test_throws ArgumentError structure_index(STRUCTURES, zeros(Int, 28, 3))
    @test_throws ArgumentError structure_index(S10, N_STRUCTURES + 1)
end

@testset "StructureSet" begin
    # a set built from contact maps holds the same structures, in the same order
    maps = [contact_map(s; structures = S10) for s in 1:50]
    S = StructureSet(maps)
    @test length(S) == 50
    @test all(contact_map(s; structures = S) == maps[s] for s in 1:50)
    @test S.contacts == S10.contacts[:, :, 1:50]
    @test structure_index(S, 625; from = S10) === nothing
    @test structure_index(S, 20; from = S10) == 20

    # any number of contacts, of any integer type
    T = StructureSet([Int8[1 4], Int8[1 6]])
    @test length(T) == 2
    @test contact_map(2; structures = T) == [1 6]
    @test sprint(show, T) == "StructureSet with 2 structures of 1 contacts"

    @test_throws ArgumentError StructureSet(Matrix{Int}[])
    @test_throws ArgumentError StructureSet([zeros(Int, 0, 2)])
    @test_throws ArgumentError StructureSet([[1 4], [1 6; 2 5]]) # different numbers of contacts
    @test_throws ArgumentError StructureSet([[1 4 6]])           # not pairs
    @test_throws ArgumentError StructureSet([[0 4]])             # residues are in 1:27
    @test_throws ArgumentError StructureSet([[1 28]])
    @test_throws ArgumentError StructureSet([[1 1000]])
    @test_throws ArgumentError StructureSet([[3 3]])             # self-contact

    # the per-site index lists every contact from both of its residues
    sites = LatticeProteins.site_contacts(T)
    @test LatticeProteins.site_contacts(T) === sites # built once
    @test sites.offsets == [1; 3; 3; 3; 4; 4; fill(5, 22)]
    @test sites.structure == [1, 2, 1, 2]
    @test sites.partner == [4, 6, 1, 1]
end
