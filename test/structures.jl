using Test: @test, @testset
using LatticeProteins: enumerate_structures, structure_index, energy, energies, pnat, log_pnat,
    metropolis!, random_msa, potts, STRUCTURES, SHAKHNOVICH_STRUCTURES, N_STRUCTURES, N_CONTACTS, L,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D
using Random: seed!

contact_key(cm) = sort([minmax(cm[k, 1], cm[k, 2]) for k in axes(cm, 1)])
# chain reversal maps contact (i, j) to (L + 1 - j, L + 1 - i)
reverse_structure(cm) = [L + 1 - cm[k, 3 - c] for k in axes(cm, 1), c in 1:2]

@testset "enumerate_structures" begin
    # Shakhnovich & Gutin (1990): 103,346 compact conformations unrelated by symmetry
    @test size(STRUCTURES) == (N_CONTACTS, 103_346, 2)
    @test eltype(STRUCTURES) == Int
    @test all(1 .≤ STRUCTURES .≤ L)

    I = STRUCTURES[:, :, 1]
    J = STRUCTURES[:, :, 2]
    @test all(I .< J)
    # the lattice is bipartite, so contacting residues have opposite parity, and
    # residues adjacent along the chain are bonded, not in contact
    @test all(isodd.(J .- I))
    @test all(J .- I .≥ 3)
    # contacts are sorted within each structure
    @test all(issorted(collect(zip(I[:, s], J[:, s]))) for s in 1:N_STRUCTURES)

    keys = Set(contact_key(STRUCTURES[:, s, :]) for s in 1:N_STRUCTURES)
    @test length(keys) == N_STRUCTURES # no duplicates
    # the set is closed under chain reversal; a few structures are their own reversal
    @test all(contact_key(reverse_structure(STRUCTURES[:, s, :])) in keys for s in 1:N_STRUCTURES)
    @test count(contact_key(reverse_structure(STRUCTURES[:, s, :])) == contact_key(STRUCTURES[:, s, :]) for s in 1:N_STRUCTURES) == 62
    # the 10,000 reference structures are among them
    @test all(contact_key(SHAKHNOVICH_STRUCTURES[:, s, :]) in keys for s in 1:10_000)
    # deterministic: STRUCTURES is the output of enumerate_structures
    @test enumerate_structures() == STRUCTURES
end

@testset "structure_index" begin
    for (cm, s10) in zip((CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D), (625, 848, 6683, 6685))
        @test structure_index(SHAKHNOVICH_STRUCTURES, STRUCTURES[:, cm, :]) == s10
        @test structure_index(SHAKHNOVICH_STRUCTURES, cm) == s10
        @test structure_index(STRUCTURES, SHAKHNOVICH_STRUCTURES[:, s10, :]) == cm
        @test structure_index(STRUCTURES, cm) == cm
        @test contact_key(STRUCTURES[:, cm, :]) == contact_key(SHAKHNOVICH_STRUCTURES[:, s10, :])
    end
    @test structure_index(SHAKHNOVICH_STRUCTURES, SHAKHNOVICH_STRUCTURES[:, 7, :]) == 7
    # most compact structures are not in the reference set
    @test any(structure_index(SHAKHNOVICH_STRUCTURES, STRUCTURES[:, s, :]) === nothing for s in 1:100)
    # the order of the contacts and the order within each pair are ignored
    cm = STRUCTURES[:, CONTACT_MAP_A, :]
    @test structure_index(STRUCTURES, reverse(cm; dims=1)) == CONTACT_MAP_A
    @test structure_index(STRUCTURES, cm[:, [2, 1]]) == CONTACT_MAP_A
    # a contact map with a different number of contacts is never found
    @test structure_index(STRUCTURES, cm[1:27, :]) === nothing
    # the chain reversal of structure A is a different structure, absent from the reference set
    @test structure_index(SHAKHNOVICH_STRUCTURES, reverse_structure(cm)) === nothing
    @test structure_index(STRUCTURES, reverse_structure(cm)) != CONTACT_MAP_A
end

@testset "pnat against the reference set" begin
    seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
    a10 = structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)

    @test energy(a10, seq; structures=SHAKHNOVICH_STRUCTURES) == energy(CONTACT_MAP_A, seq)
    E = energies(seq)
    @test E isa Vector{Float64}
    @test length(E) == N_STRUCTURES
    @test E[CONTACT_MAP_A] == energy(CONTACT_MAP_A, seq)
    @test E == [energy(s, seq) for s in 1:N_STRUCTURES]
    E10 = energies(seq; structures=SHAKHNOVICH_STRUCTURES)
    @test E10 == [energy(s, seq; structures=SHAKHNOVICH_STRUCTURES) for s in 1:10_000]

    # the reference set is a subset of the full set, so pnat can only be lower with the full set
    @test pnat(CONTACT_MAP_A, seq) < pnat(a10, seq; structures=SHAKHNOVICH_STRUCTURES)
    @test log_pnat(CONTACT_MAP_A, seq) ≈ log(pnat(CONTACT_MAP_A, seq))

    # normalization: pnat is a probability distribution over structures
    m = maximum(-E)
    p = exp.(-E .- (m + log(sum(exp.(-E .- m)))))
    @test sum(p) ≈ 1
    @test pnat(CONTACT_MAP_A, seq) ≈ p[CONTACT_MAP_A]

    # the default is the full set
    @test pnat(CONTACT_MAP_A, seq; structures=STRUCTURES) == pnat(CONTACT_MAP_A, seq)
    @test log_pnat(CONTACT_MAP_A, seq; structures=STRUCTURES) == log_pnat(CONTACT_MAP_A, seq)
end

@testset "metropolis against the reference set" begin
    seed!(1)
    a10 = structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)
    seq = rand(1:20, L)
    pnat_before = pnat(a10, seq; structures=SHAKHNOVICH_STRUCTURES)
    metropolis!(seq, a10; β=100, nsteps=500, structures=SHAKHNOVICH_STRUCTURES)
    pnat_after = pnat(a10, seq; structures=SHAKHNOVICH_STRUCTURES)
    @test all(1 .≤ seq .≤ 20)
    @test pnat_after > pnat_before

    msa = random_msa(a10; nseqs=2, β=100, nsteps=50, structures=SHAKHNOVICH_STRUCTURES)
    @test size(msa) == (L, 2)
    @test all(1 .≤ msa .≤ 20)
end
