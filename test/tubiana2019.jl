using Test: @test, @testset
using LatticeProteins: LatticeProteins, Hugo_MSA, Hugo_MSA_pnat, potts, onehot, log_pnat,
    metropolis!, structure_index, SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A, AMINO_ACIDS, L, N_CONTACTS
using LinearAlgebra: eigen, Symmetric, dot, norm
using Statistics: mean, cor
using Random: seed!

# Reproduces the statistics of the lattice-protein alignment of structure A shown in Fig. 4 of
#   Tubiana, Cocco & Monasson. "Learning compositional representations of interacting systems
#   with restricted Boltzmann machines: comparative study of lattice proteins."
#   Neural Computation 31.8 (2019): 1671-1717. https://doi.org/10.1162/neco_a_01210
# The alignment of that figure is generated with the protocol of Jacquin et al. (2016):
# Metropolis sampling at β = 1000 against the 10,000 reference structures, keeping one
# sequence at regular intervals. We sample a fresh alignment with this package and compare
# its site frequencies, conservation profile, first principal component and the average
# covariance between amino acids on contacting sites (panels d, e and f of the figure) with
# the reference alignment of Jacquin et al. bundled with the package.

const S = SHAKHNOVICH_STRUCTURES
const A10 = structure_index(S, CONTACT_MAP_A)

# Site frequencies (20 × 27), conservation profile, first principal component (20 × 27) and
# average covariance between amino acids on sites in contact in structure A (20 × 20).
function fig4_statistics(seqs::AbstractMatrix{Int})
    M = size(seqs, 2)
    X = Float64.(reshape(onehot(seqs)[1:20, :, :], 20 * L, M))
    f = reshape(mean(X; dims=2), 20, L)
    cons = [log(20) + sum(p > 0 ? p * log(p) : 0.0 for p in f[:, i]) for i in 1:L]
    Xc = X .- mean(X; dims=2)
    w = eigen(Symmetric(Xc * Xc' / M)).vectors[:, end]
    cm = S[:, A10, :]
    C = zeros(20, 20)
    for k in 1:N_CONTACTS
        i, j = cm[k, 1], cm[k, 2]
        fij = zeros(20, 20)
        for m in 1:M
            fij[seqs[i, m], seqs[j, m]] += 1
        end
        cov = fij ./ M .- f[:, i] * f[:, j]'
        C .+= cov .+ cov'
    end
    return f, cons, reshape(w, 20, L), C ./ N_CONTACTS
end

@testset "tubiana2019 Fig. 4" begin
    # reference alignment: the sequences of Jacquin et al. folding in A with pnat > 0.995
    msa = Hugo_MSA(:A)
    pn = Hugo_MSA_pnat(:A)
    ref = potts(msa[pn .> 0.995])
    @test size(ref, 2) > 20_000

    # Fresh alignment sampled with this package at β = 1000 against the reference structures,
    # keeping one sequence every 1,000 steps as in Jacquin et al. (2016). At this β a chain
    # rarely crosses between the two subfamilies of the alignment (the two modes of panel g
    # of the figure), so several independent chains are needed to sample both.
    seed!(2019)
    nchains, nsamples = 20, 25
    fresh = zeros(Int, L, nchains * nsamples)
    pnats = zeros(nchains * nsamples)
    for c in 1:nchains
        seq = rand(1:20, L)
        metropolis!(seq, A10; β=1000, nsteps=10_000, structures=S) # burn-in
        for t in 1:nsamples
            metropolis!(seq, A10; β=1000, nsteps=1_000, structures=S)
            k = (c - 1) * nsamples + t
            fresh[:, k] = seq
            pnats[k] = exp(log_pnat(A10, seq; structures=S))
        end
    end
    # at β = 1000 the chains sample sequences folding into A with probability ≈ 0.995
    @test mean(pnats) ≈ mean(pn) atol=0.003
    @test all(pnats .> 0.97)
    fresh = fresh[:, pnats .> 0.995]
    @test size(fresh, 2) > 150

    f_ref, cons_ref, w_ref, C_ref = fig4_statistics(ref)
    f, cons, w, C = fig4_statistics(fresh)

    # panel (d): site frequencies and conservation profile
    @test cor(vec(f), vec(f_ref)) > 0.9
    @test cor(cons, cons_ref) > 0.9
    @test maximum(abs, f .- f_ref) < 0.15
    # panel (e): first principal component (up to sign)
    @test abs(dot(vec(w), vec(w_ref))) / (norm(w) * norm(w_ref)) > 0.85
    # panel (g): both subfamilies (the two modes of the projection on PC1) are sampled
    center = vec(mean(Float64.(reshape(onehot(ref)[1:20, :, :], 20 * L, :)); dims=2))
    proj = (Float64.(reshape(onehot(fresh)[1:20, :, :], 20 * L, :)) .- center)' * vec(w_ref)
    @test 0.2 < mean(proj .> 0) < 0.85
    # panel (f): average covariance between amino acids on contacting sites
    @test cor(vec(C), vec(C_ref)) > 0.95
    # the largest covariance is cysteine-cysteine, opposite charges covary positively and
    # like charges negatively, as in the paper
    aa(c) = findfirst(==(c), AMINO_ACIDS)
    @test argmax(C) == CartesianIndex(aa('C'), aa('C'))
    @test C[aa('E'), aa('K')] > 0 && C[aa('E'), aa('R')] > 0 && C[aa('D'), aa('K')] > 0
    @test C[aa('E'), aa('E')] < 0 && C[aa('K'), aa('K')] < 0
end
