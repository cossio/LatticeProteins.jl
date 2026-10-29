"""
    LatticeProteins

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

A lattice protein is a simplified model where a protein sequence of 27 amino acids
folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is
determined by pairwise contacts using the Miyazawa-Jernigan interaction matrix. The
probability that a sequence folds into a structure is a Boltzmann distribution over the
competing compact structures: by default all 103,346 compact structures of the lattice
([`STRUCTURES`](@ref)), or optionally the representative subset of 10,000 structures of
Heo, Maslov & Shakhnovich (2011) used by Jacquin et al. (2016)
([`SHAKHNOVICH_STRUCTURES`](@ref)).

# References

- Jacquin, Hugo, et al. "Benchmarking inverse statistical approaches for protein
  structure and design with exactly solvable models." PLoS computational biology
  12.5 (2016): e1004889.
"""
module LatticeProteins

using LazyArtifacts: LazyArtifacts, @artifact_str
using DelimitedFiles: readdlm
using Random: AbstractRNG, default_rng, randexp

"""
    N_CONTACTS

Number of pairwise contacts in the 3×3×3 cubic lattice (28).
"""
const N_CONTACTS = 28 # number of pairwise contacts in 3x3x3 lattice

include("miyazawa_jernigan.jl")
include("contact_maps.jl")
include("structures.jl")
include("pnat.jl")
include("onehot.jl")
include("metropolis.jl")
include("parallel_tempering.jl")
include("artifacts.jl")

end # module
