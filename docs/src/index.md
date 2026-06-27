# LatticeProteins.jl

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

In this model, a protein is represented as a chain of 27 amino acids that folds on
a 3×3×3 cubic lattice. The energy of a sequence in a given structure is determined by
pairwise contacts between amino acids using the Miyazawa-Jernigan interaction matrix.
The probability that a sequence folds into a given "native" structure is computed via a
Boltzmann distribution over 10,000 pre-computed compact structures.

This package does not export any symbols; everything is accessed either with the
`LatticeProteins.` prefix or by importing the names explicitly, e.g.
`using LatticeProteins: pnat, energy`.

## Installation

This package is registered. Install it with:

```julia
import Pkg
Pkg.add("LatticeProteins")
```

## Quick start

```julia
using LatticeProteins: pnat, energy, onehot, aaseq, potts, random_msa,
    CONTACT_MAP_A, Hugo_MSA, Hugo_MSA_pnat

# Compute the probability that a sequence folds into structure A
seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")  # convert to integer-coded sequence
p = pnat(CONTACT_MAP_A, seq)

# Compute the contact energy
E = energy(CONTACT_MAP_A, seq)

# Generate sequences that fold into structure A via Metropolis sampling
msa = random_msa(CONTACT_MAP_A; nseqs=100, β=100, nsteps=1000)
```

See the [Tutorial](@ref) for a worked walkthrough and the [Reference](@ref) for the
full list of functions.

## Contents

```@contents
Pages = ["literate/tutorial.md", "reference.md"]
Depth = 2
```

## Citation

If you use this package, please cite:

Fernandez-de-Cossio-Diaz, Jorge, Clément Roussel, Simona Cocco, and Remi Monasson.
["Accelerated Sampling with Stacked Restricted Boltzmann Machines."](https://openreview.net/forum?id=kXNJ48Hvw1)
*The Twelfth International Conference on Learning Representations* (2024).

## References

- Jacquin, Hugo, et al. ["Benchmarking inverse statistical approaches for protein
  structure and design with exactly solvable models."](https://doi.org/10.1371/journal.pcbi.1004889)
  *PLoS computational biology* 12.5 (2016): e1004889.
