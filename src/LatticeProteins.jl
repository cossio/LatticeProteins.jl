"""
    LatticeProteins

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

A lattice protein is a simplified model where a protein sequence of 27 amino acids
folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is
determined by pairwise contacts using the Miyazawa-Jernigan interaction matrix. The
probability that a sequence folds into a structure is a Boltzmann distribution over all
103,346 compact structures of the lattice, which compete to fold the sequence.

The public API is the list of exported names (see `names(LatticeProteins)`). Everything
else is internal and can change in any release.

# References

- Jacquin, Hugo, et al. "Benchmarking inverse statistical approaches for protein
  structure and design with exactly solvable models." PLoS computational biology
  12.5 (2016): e1004889.
"""
module LatticeProteins

using LazyArtifacts: LazyArtifacts, @artifact_str
using DelimitedFiles: readdlm
using Random: AbstractRNG, Xoshiro, default_rng, randexp

# Sequences
export potts, aaseq, onehot
# Structures
export enumerate_structures, contact_map, STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D
# Energy and folding probability
export MJ, energy, pnat, log_pnat
# Sampling
export sample_msa, Metropolis, ParallelTempering, metropolis!

"""
    L

Sequence length for lattice proteins on the 3×3×3 cubic lattice (27). Internal.
"""
const L = 27

"""
    N_CONTACTS

Number of pairwise contacts of a compact structure on the 3×3×3 cubic lattice (28).
"""
const N_CONTACTS = 28

include("sequences.jl")
include("miyazawa_jernigan.jl")
include("lattice.jl")
include("structures.jl")
include("energy.jl")
include("sample_msa.jl")
include("metropolis.jl")
include("parallel_tempering.jl")
include("reference_data.jl")

end # module
