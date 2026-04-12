# LatticeProteins.jl

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

In this model, a protein is represented as a chain of 27 amino acids that folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is determined by pairwise contacts between amino acids using the Miyazawa-Jernigan interaction matrix. The probability that a sequence folds into a given "native" structure is computed via a Boltzmann distribution over 10,000 pre-computed compact structures.

## Installation

This package is registered. Install it using the Julia package manager:

```julia
import Pkg
Pkg.add("LatticeProteins")
```

## Quick start

```julia
using LatticeProteins

# Compute the probability that a sequence folds into structure A
seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")  # convert to integer-coded sequence
p = pnat(CONTACT_MAP_A, seq)

# Compute the contact energy
E = energy(CONTACT_MAP_A, seq)

# One-hot encode a sequence
X = onehot("EKAMPAMDPDMAHEHKKKIRAWMFEGE")

# Convert back to amino acid string
s = aaseq(X)

# Generate sequences that fold into structure A via Metropolis sampling
msa = random_msa(CONTACT_MAP_A; nseqs=100, β=100, nsteps=1000)

# Load reference MSA datasets from Jacquin et al. (2016)
sequences = Hugo_MSA(:A)
pnat_values = Hugo_MSA_pnat(:A)
```

## API overview

### Structures and contact maps

| Symbol | Description |
|--------|-------------|
| `CONTACT_MAPS` | Pre-loaded contact maps for 10,000 compact structures |
| `N_STRUCTURES` | Total number of structures (10,000) |
| `CONTACT_MAP_A`, `CONTACT_MAP_B`, `CONTACT_MAP_C`, `CONTACT_MAP_D` | Indices of four reference structures from Jacquin et al. (2016) |

### Energy and folding probability

| Function | Description |
|----------|-------------|
| `energy(cm, seq)` | Contact energy of sequence `seq` in structure `cm` |
| `pnat(cm, seq)` | Probability that `seq` folds into structure `cm` |
| `log_pnat(cm, seq)` | Log-probability of folding |

### Sequence encoding

| Function | Description |
|----------|-------------|
| `onehot(seq)` | One-hot encode amino acid sequence(s) |
| `potts(seq)` | Convert to integer-coded (Potts) representation |
| `aaseq(X)` | Convert from one-hot or integer coding back to amino acid string(s) |
| `AMINO_ACIDS` | Amino acid alphabet: `"CMFILVWYAGTSNQDEHRKP-"` |

### Sampling

| Function | Description |
|----------|-------------|
| `metropolis!(seq, cm; β, nsteps)` | In-place Metropolis sampling of a sequence |
| `random_msa(cm; nseqs, β, nsteps)` | Generate a random MSA via Metropolis sampling |

### Reference data

| Function | Description |
|----------|-------------|
| `Hugo_MSA(which)` | Load reference MSA (`:A`, `:B`, `:C`, or `:D`) |
| `Hugo_MSA_pnat(which)` | Load reference pnat values |
| `MJ` | 20×20 Miyazawa-Jernigan interaction matrix |

## References

- Jacquin, Hugo, et al. ["Benchmarking inverse statistical approaches for protein structure and design with exactly solvable models."](https://doi.org/10.1371/journal.pcbi.1004889) *PLoS computational biology* 12.5 (2016): e1004889.