using Test: @test, @testset
using LatticeProteins: enumerate_structures, structure_index, energy, energies, pnat, log_pnat,
    metropolis!, random_msa, potts, CONTACT_MAPS, N_STRUCTURES, N_CONTACTS, L,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D
using Random: seed!

const ALL = enumerate_structures()

contact_key(cm) = sort([minmax(cm[k, 1], cm[k, 2]) for k in axes(cm, 1)])
# chain reversal maps contact (i, j) to (L + 1 - j, L + 1 - i)
reverse_structure(cm) = [L + 1 - cm[k, 3 - c] for k in axes(cm, 1), c in 1:2]

@testset "enumerate_structures" begin
    # Shakhnovich & Gutin (1990): 103,346 compact conformations unrelated by symmetry
    @test size(ALL) == (N_CONTACTS, 103_346, 2)
    @test eltype(ALL) == Int
    @test all(1 .≤ ALL .≤ L)

    I = ALL[:, :, 1]
    J = ALL[:, :, 2]
    @test all(I .< J)
    # the lattice is bipartite, so contacting residues have opposite parity, and
    # residues adjacent along the chain are bonded, not in contact
    @test all(isodd.(J .- I))
    @test all(J .- I .≥ 3)
    # contacts are sorted within each structure
    @test all(issorted(collect(zip(I[:, s], J[:, s]))) for s in axes(ALL, 2))

    keys = Set(contact_key(ALL[:, s, :]) for s in axes(ALL, 2))
    @test length(keys) == size(ALL, 2) # no duplicates
    # the set is closed under chain reversal
    @test all(contact_key(reverse_structure(ALL[:, s, :])) in keys for s in axes(ALL, 2))
    # the 10,000 reference structures are among them
    @test all(contact_key(CONTACT_MAPS[:, s, :]) in keys for s in 1:N_STRUCTURES)
    # deterministic order
    @test enumerate_structures() == ALL
end

@testset "structure_index" begin
    for cm in (CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D)
        s = structure_index(ALL, CONTACT_MAPS[:, cm, :])
        @test s isa Int
        @test contact_key(ALL[:, s, :]) == contact_key(CONTACT_MAPS[:, cm, :])
        @test structure_index(ALL, cm) == s
        @test structure_index(CONTACT_MAPS, ALL[:, s, :]) == cm
    end
    @test structure_index(CONTACT_MAPS, 7) == 7
    # most compact structures are not in the reference set
    @test any(structure_index(CONTACT_MAPS, ALL[:, s, :]) === nothing for s in 1:100)
    # the order of the contacts and the order within each pair are ignored
    cm = CONTACT_MAPS[:, CONTACT_MAP_A, :]
    @test structure_index(CONTACT_MAPS, reverse(cm; dims=1)) == CONTACT_MAP_A
    @test structure_index(CONTACT_MAPS, cm[:, [2, 1]]) == CONTACT_MAP_A
    # the chain reversal of structure A is a different structure, absent from the reference set
    @test structure_index(CONTACT_MAPS, reverse_structure(cm)) === nothing
    @test structure_index(ALL, reverse_structure(cm)) != structure_index(ALL, cm)
end

@testset "pnat over all structures" begin
    seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
    a = structure_index(ALL, CONTACT_MAP_A)

    @test energy(a, seq; structures=ALL) == energy(CONTACT_MAP_A, seq)
    E = energies(seq; structures=ALL)
    @test E isa Vector{Float64}
    @test length(E) == size(ALL, 2)
    @test E[a] == energy(CONTACT_MAP_A, seq)
    @test E == [energy(s, seq; structures=ALL) for s in axes(ALL, 2)]
    @test energies(seq) == [energy(s, seq) for s in 1:N_STRUCTURES]

    # the reference set is a subset, so pnat can only decrease with the full set
    @test pnat(a, seq; structures=ALL) < pnat(CONTACT_MAP_A, seq)
    @test log_pnat(a, seq; structures=ALL) ≈ log(pnat(a, seq; structures=ALL))

    # normalization: pnat over the full set is a probability distribution over structures
    m = maximum(-E)
    p = exp.(-E .- (m + log(sum(exp.(-E .- m)))))
    @test sum(p) ≈ 1
    @test pnat(a, seq; structures=ALL) ≈ p[a]

    # the default is the reference set
    @test pnat(CONTACT_MAP_A, seq; structures=CONTACT_MAPS) == pnat(CONTACT_MAP_A, seq)
    @test log_pnat(CONTACT_MAP_A, seq; structures=CONTACT_MAPS) == log_pnat(CONTACT_MAP_A, seq)
end

@testset "metropolis over all structures" begin
    seed!(1)
    a = structure_index(ALL, CONTACT_MAP_A)
    seq = rand(1:20, L)
    pnat_before = pnat(a, seq; structures=ALL)
    metropolis!(seq, a; β=100, nsteps=500, structures=ALL)
    pnat_after = pnat(a, seq; structures=ALL)
    @test all(1 .≤ seq .≤ 20)
    @test pnat_after > pnat_before

    msa = random_msa(a; nseqs=2, β=100, nsteps=50, structures=ALL)
    @test size(msa) == (L, 2)
    @test all(1 .≤ msa .≤ 20)
end
