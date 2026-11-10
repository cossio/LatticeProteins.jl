"""
    LatticeProteins

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

A lattice protein is a simplified model where a protein sequence of 27 amino acids
folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is
determined by pairwise contacts using the Miyazawa-Jernigan interaction matrix. The
probability that a sequence folds into a structure is a Boltzmann distribution over a
set of competing compact structures ([`StructureSet`](@ref)): by default all 103,346
compact structures of the lattice ([`STRUCTURES`](@ref)), or optionally the representative
subset of 10,000 structures of Heo, Maslov & Shakhnovich (2011) used by Jacquin et al.
(2016) ([`shakhnovich_structures`](@ref)).

The package exports no names. Its public API is the list of names declared `public`
(see `names(LatticeProteins)`); everything else is internal and can change in any release.

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
public AMINO_ACIDS, onehot, potts, aaseq
# Lattice and structures
public L, N_CONTACTS, N_STRUCTURES
public StructureSet, STRUCTURES, shakhnovich_structures, enumerate_structures
public contact_map, structure_index
public STRUCTURE_A, STRUCTURE_B, STRUCTURE_C, STRUCTURE_D
# Energy and folding probability
public MJ, energy, energies, pnat, log_pnat
# Sampling
public metropolis!, random_msa, parallel_tempering
# Reference data
public Hugo_MSA, Hugo_MSA_pnat, Hugo_MSA_path

"""
    L

Sequence length for lattice proteins on the 3×3×3 cubic lattice (27).
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
include("metropolis.jl")
include("parallel_tempering.jl")
include("reference_data.jl")

end # module
