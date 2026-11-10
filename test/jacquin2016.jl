using Test: @test, @testset
using LatticeProteins: Hugo_MSA, Hugo_MSA_pnat, onehot, L, AMINO_ACIDS
using Statistics: mean

# Properties of the reference alignments of
#   Jacquin, Hugo, et al. "Benchmarking inverse statistical approaches for protein
#   structure and design with exactly solvable models."
#   PLoS Computational Biology 12.5 (2016): e1004889.
# That the stored folding probabilities match the sequences is tested in reference_data.jl,
# and sampling such alignments with this package in tubiana2019.jl.

@testset "jacquin2016" begin
    for which in (:A, :B, :C, :D)
        msa = Hugo_MSA(which)
        pnats = Hugo_MSA_pnat(which)

        # The alignments were sampled at β = 1000 (first column of the files), where
        # sequences fold with high probability: typically pnat ≈ 0.995.
        @test all(pnats .> 0.97)
        @test mean(pnats) > 0.99

        # sequences use the 20 amino acids, without gaps
        @test all(all(c in AMINO_ACIDS[1:20] for c in s) for s in msa)

        # Amino acid frequencies at each site are not uniform (1/20 = 0.05), reflecting
        # the selection for the structure, and some sites are highly conserved.
        freq = dropdims(mean(onehot(msa); dims = 3); dims = 3)
        max_freqs = [maximum(freq[1:20, i]) for i in 1:L]
        @test all(max_freqs .> 0.05)
        @test maximum(max_freqs) > 0.2
    end
end
