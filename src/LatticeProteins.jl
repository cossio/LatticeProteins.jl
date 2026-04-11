"""
    LatticeProteins

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

A lattice protein is a simplified model where a protein sequence of 27 amino acids
folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is
determined by pairwise contacts using the Miyazawa-Jernigan interaction matrix.

# References

- Jacquin, Hugo, et al. "Benchmarking inverse statistical approaches for protein
  structure and design with exactly solvable models." PLoS computational biology
  12.5 (2016): e1004889.
"""
module LatticeProteins

using LazyArtifacts: LazyArtifacts, @artifact_str
using DelimitedFiles: readdlm
using Random: randexp
using LogExpFunctions: logsumexp

"""
    N_CONTACTS

Number of pairwise contacts in the 3×3×3 cubic lattice (28).
"""
const N_CONTACTS = 28 # number of pairwise contacts in 3x3x3 lattice

include("miyazawa_jernigan.jl")
include("contact_maps.jl")
include("pnat.jl")
include("onehot.jl")
include("metropolis.jl")
include("artifacts.jl")

end # module
