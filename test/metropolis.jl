using Test: @test, @testset, @test_throws
using LatticeProteins: LatticeProteins, random_msa, metropolis!, pnat, log_pnat, energies,
    structure_index, StructureSet, STRUCTURE_A, STRUCTURES, shakhnovich_structures, MJ, L
using Random: Xoshiro, seed!, randexp

const S10 = shakhnovich_structures()
const A10 = structure_index(S10, STRUCTURE_A)

# Reference implementation that recomputes log_pnat from scratch at every step, with the
# same sequence of random draws as metropolis!.
function naive_metropolis!(rng, seq, s; β, nsteps, structures)
    loss = -log_pnat(s, seq; structures)
    for _ in 1:nsteps
        i = rand(rng, 1:L)
        aa = rand(rng, 1:19)
        aa += aa >= seq[i] # one of the other 19 amino acids
        old = seq[i]
        seq[i] = aa
        loss1 = -log_pnat(s, seq; structures)
        if β * (loss1 - loss) < randexp(rng)
            loss = loss1
        else
            seq[i] = old
        end
    end
    return seq
end

@testset "metropolis!" begin
    seq = rand(Xoshiro(42), 1:20, L)
    seq0 = copy(seq)
    @test metropolis!(seq, STRUCTURE_A; β = 1, nsteps = 10, rng = Xoshiro(1)) === seq # in place
    @test all(1 .≤ seq .≤ 20)
    @test seq != seq0

    # with nsteps = 0 the sequence does not change
    @test metropolis!(copy(seq), STRUCTURE_A; nsteps = 0) == seq
    # at β = 0 every proposal is accepted, and proposals always change the sequence
    @test metropolis!(copy(seq), STRUCTURE_A; β = 0, nsteps = 1) != seq

    # the default RNG is the task-local one
    seed!(3)
    a = metropolis!(copy(seq0), STRUCTURE_A; β = 10, nsteps = 50)
    seed!(3)
    @test metropolis!(copy(seq0), STRUCTURE_A; β = 10, nsteps = 50) == a
    @test metropolis!(copy(seq0), STRUCTURE_A; β = 10, nsteps = 50) != a

    # views and other integer types are sampled in place too
    seqs = Int8.(rand(Xoshiro(4), 1:20, L, 2))
    column = copy(seqs[:, 2])
    metropolis!(view(seqs, :, 2), A10; β = 10, nsteps = 50, structures = S10, rng = Xoshiro(5))
    @test seqs[:, 2] == metropolis!(Int.(column), A10; β = 10, nsteps = 50, structures = S10, rng = Xoshiro(5))

    # at high β the chain finds sequences that fold into the target, on both sets
    seq = rand(Xoshiro(6), 1:20, L)
    @test pnat(STRUCTURE_A, seq) < 0.01
    metropolis!(seq, STRUCTURE_A; β = 1000, nsteps = 5000, rng = Xoshiro(7))
    @test pnat(STRUCTURE_A, seq) > 0.95
    seq = rand(Xoshiro(6), 1:20, L)
    metropolis!(seq, A10; β = 1000, nsteps = 5000, structures = S10, rng = Xoshiro(7))
    @test pnat(A10, seq; structures = S10) > 0.95

    @test_throws ArgumentError metropolis!([1, 2, 3], STRUCTURE_A)
    @test_throws ArgumentError metropolis!(fill(21, L), STRUCTURE_A)
    @test_throws ArgumentError metropolis!(copy(seq0), 0)
    @test_throws ArgumentError metropolis!(copy(seq0), 10_001; structures = S10)
    @test_throws ArgumentError metropolis!(copy(seq0), STRUCTURE_A; nsteps = -1)
end

@testset "incremental energies match recomputation" begin
    rng = Xoshiro(6)
    sites = LatticeProteins.site_contacts(STRUCTURES)
    seq = rand(rng, 1:20, L)
    E = round.(Int32, 100 .* energies(seq))
    E1 = similar(E)
    # every position, including endpoints and no-op mutations
    for i in 1:L, aa in (seq[i], rand(rng, 1:20))
        LatticeProteins.mutate_energies!(E1, E, seq, i, aa, sites)
        seq[i] = aa
        @test E1 == round.(Int32, 100 .* energies(seq))
        E, E1 = E1, E
    end

    # same trajectories as recomputing log_pnat at every step, on both structure sets
    seq0 = rand(Xoshiro(7), 1:20, L)
    for β in (1, 10, 100)
        seq1 = metropolis!(copy(seq0), A10; β, nsteps = 300, structures = S10, rng = Xoshiro(8))
        seq2 = naive_metropolis!(Xoshiro(8), copy(seq0), A10; β, nsteps = 300, structures = S10)
        @test seq1 == seq2
        @test seq1 != seq0
    end
    seq1 = metropolis!(copy(seq0), STRUCTURE_A; β = 10, nsteps = 100, rng = Xoshiro(8))
    seq2 = naive_metropolis!(Xoshiro(8), copy(seq0), STRUCTURE_A; β = 10, nsteps = 100, structures = STRUCTURES)
    @test seq1 == seq2
    @test seq1 != seq0
end

@testset "stationary distribution against exact enumeration" begin
    # Fold 1 contacts sites 1--4, fold 2 contacts sites 1--6. Only three residues affect
    # the energy, so all 20^3 states can be enumerated.
    structures = StructureSet([[1 4], [1 6]])
    β = 3.0
    exact = zeros(20, 3)
    for a in 1:20, b in 1:20, c in 1:20
        weight = (1 / (1 + exp(MJ[a, b] - MJ[a, c])))^β
        exact[a, 1] += weight
        exact[b, 2] += weight
        exact[c, 3] += weight
    end
    exact ./= sum(exact; dims = 1)
    rng = Xoshiro(2026)
    seq = rand(rng, 1:20, L)
    metropolis!(seq, 1; β, nsteps = 10_000, structures, rng)
    observed = zeros(20, 3)
    nsamples = 40_000
    for _ in 1:nsamples
        metropolis!(seq, 1; β, nsteps = 30, structures, rng)
        for (j, site) in enumerate((1, 4, 6))
            observed[seq[site], j] += 1 / nsamples
        end
    end
    @test maximum(abs, observed .- exact) < 0.015
end

@testset "random_msa" begin
    msa = random_msa(STRUCTURE_A; nseqs = 5, β = 1, nsteps = 10, rng = Xoshiro(42))
    @test msa isa Matrix{Int}
    @test size(msa) == (L, 5)
    @test all(1 .≤ msa .≤ 20)
    @test msa == random_msa(STRUCTURE_A; nseqs = 5, β = 1, nsteps = 10, rng = Xoshiro(42))
    @test size(random_msa(STRUCTURE_A; nseqs = 0)) == (L, 0)

    # independent chains started from random sequences
    rng = Xoshiro(5)
    seqs = rand(rng, 1:20, L, 3)
    for n in 1:3
        naive_metropolis!(rng, view(seqs, :, n), A10; β = 10, nsteps = 50, structures = S10)
    end
    @test random_msa(A10; nseqs = 3, β = 10, nsteps = 50, structures = S10, rng = Xoshiro(5)) == seqs

    @test_throws ArgumentError random_msa(STRUCTURE_A; nseqs = -1)
    @test_throws ArgumentError random_msa(0)
end
