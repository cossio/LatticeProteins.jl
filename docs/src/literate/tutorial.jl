# # Tutorial
#
# This tutorial walks through the main features of `LatticeProteins.jl`: encoding
# amino acid sequences, computing contact energies and folding probabilities, and
# sampling sequences that fold into a target structure.

# ```@meta
# CurrentModule = LatticeProteins
# ```

# ## Setup
#
# Since the package does not export any symbols, we import the names we need
# explicitly. (Alternatively, everything is reachable through the `LatticeProteins.`
# prefix.)

using LatticeProteins: pnat, log_pnat, energy, onehot, aaseq, potts, AMINO_ACIDS,
    metropolis!, random_msa, CONTACT_MAP_A, CONTACT_MAP_B,
    STRUCTURES, SHAKHNOVICH_STRUCTURES, N_STRUCTURES, MJ, Hugo_MSA, Hugo_MSA_pnat,
    structure_index

# A lattice protein is a chain of 27 amino acids folded onto a 3×3×3 cubic lattice.
# There are 103,346 maximally compact structures (counted up to rotations and reflections
# of the cube), which the package enumerates when it is loaded. Their contact maps are
# stored in `STRUCTURES`, an array of size `(28, N_STRUCTURES, 2)`:

N_STRUCTURES

# Four of these structures are highlighted in Jacquin et al. (2016) and exposed as
# constants `CONTACT_MAP_A`, `CONTACT_MAP_B`, `CONTACT_MAP_C` and `CONTACT_MAP_D`.
# Each is simply an index into `STRUCTURES`:

CONTACT_MAP_A

# The 28 contacts of structure A, as pairs of residue positions:

STRUCTURES[:, CONTACT_MAP_A, :]

# ## Encoding sequences
#
# Sequences are written as strings over the [`AMINO_ACIDS`](@ref) alphabet (20 amino
# acids plus a gap character):

AMINO_ACIDS

# Most numerical routines work with *integer-coded* (Potts) sequences, where each
# amino acid is replaced by its index in the alphabet. Use [`potts`](@ref) to convert:

seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")

# We can convert back to a string with [`aaseq`](@ref):

aaseq(seq)

# A sequence can also be one-hot encoded into a `21 × 27` `BitMatrix` with
# [`onehot`](@ref), which is convenient for machine-learning workflows:

X = onehot("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
size(X)

# ## Energy and folding probability
#
# The contact [`energy`](@ref) of a sequence in a given structure is the sum of
# Miyazawa-Jernigan interaction energies over all contacting residue pairs:

energy(CONTACT_MAP_A, seq)

# The same sequence has a different energy in a different structure:

energy(CONTACT_MAP_B, seq)

# The probability that a sequence folds into a particular "native" structure is given
# by a Boltzmann distribution over all competing compact structures. This is
# [`pnat`](@ref) (and its logarithm [`log_pnat`](@ref)):

pnat(CONTACT_MAP_A, seq)

#-

log_pnat(CONTACT_MAP_A, seq)

# ## Sampling sequences that fold
#
# We can search for sequences that fold into a target structure using Metropolis
# Monte Carlo. [`random_msa`](@ref) starts from random sequences and evolves each one
# at inverse temperature `β`; a large `β` strongly favours sequences with high
# folding probability for the target structure.

msa = random_msa(CONTACT_MAP_A; nseqs = 20, β = 100, nsteps = 1000)
size(msa)

# Each column of the returned matrix is an integer-coded sequence. The sampled
# sequences fold into structure A with high probability:

using Statistics: mean
mean(pnat(CONTACT_MAP_A, col) for col in eachcol(msa))

# [`metropolis!`](@ref) exposes the same sampler for a single sequence, mutating it
# in place:

s = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
metropolis!(s, CONTACT_MAP_A; β = 100, nsteps = 200)
pnat(CONTACT_MAP_A, s)

# ## The subset of 10,000 structures of Jacquin et al. (2016)
#
# Jacquin et al. (2016) restricted the competing structures to a representative subset of
# 10,000 selected by the Shakhnovich group (Heo, Maslov & Shakhnovich 2011), which the
# package ships as `SHAKHNOVICH_STRUCTURES` (same layout as `STRUCTURES`):

size(SHAKHNOVICH_STRUCTURES)

# All functions accept the set of competing structures through the `structures` keyword
# argument. Structure indices refer to the set in use, and [`structure_index`](@ref)
# translates between sets. Structure A is number 625 of the reference set:

a = structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)

# The reference set has fewer competitors, so folding probabilities computed against it
# are higher than against all compact structures:

pnat(a, seq; structures = SHAKHNOVICH_STRUCTURES)

#-

pnat(CONTACT_MAP_A, seq)

# Sampling works the same way:

msa_ref = random_msa(a; nseqs = 5, β = 100, nsteps = 200, structures = SHAKHNOVICH_STRUCTURES)
mean(pnat(a, col; structures = SHAKHNOVICH_STRUCTURES) for col in eachcol(msa_ref))

# ## Reference datasets
#
# The package also bundles the reference multiple sequence alignments (MSAs) and
# their folding probabilities from Jacquin et al. (2016), accessible with
# [`Hugo_MSA`](@ref) and [`Hugo_MSA_pnat`](@ref). The folding probabilities were computed
# against the subset of 10,000 structures:

sequences = Hugo_MSA(:A)
length(sequences)

#-

pnat_values = Hugo_MSA_pnat(:A)
first(pnat_values, 5)

#-

pnat(a, potts(sequences[1]); structures = SHAKHNOVICH_STRUCTURES)

# See the [Reference](@ref) for the complete list of available functions.
