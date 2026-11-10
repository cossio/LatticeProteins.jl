# LatticeProteins.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cossio.github.io/LatticeProteins.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cossio.github.io/LatticeProteins.jl/dev/)
[![Coverage](https://codecov.io/gh/cossio/LatticeProteins.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cossio/LatticeProteins.jl)

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

In this model, a protein is represented as a chain of 27 amino acids that folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is determined by pairwise contacts between amino acids using the Miyazawa-Jernigan interaction matrix. The probability that a sequence folds into a given "native" structure is computed via a Boltzmann distribution over all 103,346 compact structures of the lattice (or, optionally, over the representative subset of 10,000 structures selected by the Shakhnovich group and used by Jacquin et al. (2016)).

## Installation

This package is registered. Install it using the Julia package manager:

```julia
import Pkg
Pkg.add("LatticeProteins")
```

## Quick start

The package does not export any names: import the ones you need, or use the
`LatticeProteins.` prefix. The names in the tables below are the public API (they are
declared `public`); anything else is internal and may change in any release.

```julia
using LatticeProteins: pnat, energy, onehot, aaseq, potts, random_msa, structure_index,
    STRUCTURE_A, shakhnovich_structures, Hugo_MSA, Hugo_MSA_pnat

# Compute the probability that a sequence folds into structure A
seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")  # convert to integer-coded sequence
p = pnat(STRUCTURE_A, seq)

# Compute the contact energy
E = energy(STRUCTURE_A, seq)

# One-hot encode a sequence
X = onehot("EKAMPAMDPDMAHEHKKKIRAWMFEGE")

# Convert back to amino acid string
s = aaseq(X)

# Generate sequences that fold into structure A via Metropolis sampling
msa = random_msa(STRUCTURE_A; nseqs=100, β=100, nsteps=1000)

# Compete only against the subset of 10,000 structures used by Jacquin et al. (2016),
# as in the reference MSAs (structure indices refer to the set in use)
S10 = shakhnovich_structures()
a = structure_index(S10, STRUCTURE_A)   # 625
p_ref = pnat(a, seq; structures=S10)
msa_ref = random_msa(a; nseqs=100, β=100, nsteps=1000, structures=S10)

# Load reference MSA datasets from Jacquin et al. (2016)
sequences = Hugo_MSA(:A)
pnat_values = Hugo_MSA_pnat(:A)
```

## API overview

### Structures

A `StructureSet` is a set of structures competing to fold a sequence. Structures are
referred to by their index in a set.

| Symbol | Description |
|--------|-------------|
| `StructureSet(contact_maps)` | A set of structures, from their contact maps (`ncontacts × 2` matrices of residue pairs) |
| `STRUCTURES` | The set of all 103,346 compact structures (the default) |
| `shakhnovich_structures()` | The representative subset of 10,000 structures (Heo, Maslov & Shakhnovich 2011) used by Jacquin et al. (2016); downloaded on first use |
| `STRUCTURE_A`, `STRUCTURE_B`, `STRUCTURE_C`, `STRUCTURE_D` | Indices in `STRUCTURES` of four reference structures from Jacquin et al. (2016) |
| `length(S)`, `contact_map(S, s)` | Number of structures of a set, and contact map of its structure `s` |
| `structure_index(S, cm)`, `structure_index(S, s; from)` | Index in `S` of a contact map, or of structure `s` of another set (by default `STRUCTURES`) |
| `enumerate_structures()` | Regenerate the set of all compact structures |
| `L`, `N_CONTACTS`, `N_STRUCTURES` | Sequence length (27), contacts per compact structure (28), number of compact structures (103,346) |

### Energy and folding probability

| Function | Description |
|----------|-------------|
| `energy(s, seq)` | Contact energy of sequence `seq` in structure `s` |
| `energies(seq)` | Contact energies of `seq` in every structure |
| `pnat(s, seq)` | Probability that `seq` folds into structure `s` |
| `log_pnat(s, seq)` | Log-probability of folding |
| `MJ` | 20×20 Miyazawa-Jernigan interaction matrix |

All of these functions take a `structures` keyword argument, the `StructureSet` of competing structures (default: `STRUCTURES`).

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
| `metropolis!(seq, s; β, nsteps, structures, rng)` | In-place Metropolis sampling of a sequence |
| `random_msa(s; nseqs, β, nsteps, structures, rng)` | Generate a random MSA via Metropolis sampling |
| `parallel_tempering(s; betas, nseqs, burnin, thin, local_steps, structures, rng, threaded)` | Replica-exchange sampling with swap and round-trip diagnostics |

For strong selection, run replicas over a fixed ladder from weak to strong selection:

```julia
using LatticeProteins: parallel_tempering, STRUCTURE_A
using Random: Xoshiro

# Starting ladder for a pilot; inspect diagnostics before a production run.
betas = [0.0; exp.(range(log(0.1), log(3000.0); length=31))]
betas[end] = 3000.0
result = parallel_tempering(STRUCTURE_A; betas, burnin=5000,
    nseqs=500, thin=10, local_steps=10, rng=Xoshiro(42))
msa = result.samples                         # only the beta=3000 samples
swap_rates = result.diagnostics.swap_accepts ./ result.diagnostics.swap_attempts
result.diagnostics.round_trips               # hot → cold → hot, per replica ID
```

Burn-in and thinning are measured in exchange rounds, each containing
`local_steps` mutation proposals per replica. The example uses 50,000 burn-in
proposals per replica. The target distribution is proportional to `pnat^beta`;
the default folding model includes all 103,346 compact structures. To use the
10,000-structure subset of Jacquin et al. (2016), pass
`structures=shakhnovich_structures()` and translate the target with
`structure_index(shakhnovich_structures(), STRUCTURE_A)`. Inspect temperature
traversal and compare independent runs: a high `pnat` alone does not demonstrate
equilibration, and saved samples are not necessarily independent.

Both samplers keep the energies of the current sequence in every structure and update
them after each mutation. The MJ interactions are multiples of 0.01, so these energies
are exact integers in that unit: they never drift, and Boltzmann factors are read from a
table instead of being exponentiated for every structure and proposal.

The [Tubiana figure scripts](docs/figures/tubiana2019/README.md) run four independent
ladders, resume from checkpoints, and continue until the monitored convergence
criteria pass, including sequence frequencies and principal components.
The full-structure figure uses the resulting 316,000 saved PT draws, with 172,717
retained at `P_nat > 0.99`; the alignment and convergence report are included in
`docs/src/assets/` so the figure can be redrawn directly from these samples.

Start Julia with `--threads=4` (or another thread count) to run local replica updates
in parallel. Temperature-specific RNG streams make results identical to
`threaded=false` for the same seed and Julia version. Use `threaded=false` when
parallelizing multiple independent ladders externally. Seeded trajectories differ
from the single-stream implementation in version 2.1.0.
To reproduce that original protocol, as used by the saved Tubiana alignment, select
`rng_mode=:shared` together with `threaded=false`. The figure drivers do this explicitly;
this mode also lets a run be continued exactly from a checkpoint.

For a reproducible comparison against independent Metropolis chains with the same
total mutation budget, run:

```bash
julia --project=. --threads=4 examples/parallel_tempering.jl /tmp/pt-results
```

The script saves sequences, temperature traces, exchange rates, and run parameters;
see its header for shorter pilot settings. The output directory must not already
exist. It defaults to the full model; an optional final `shakhnovich` argument
selects the 10,000-structure subset and records that choice in the results.

### Reference data

| Function | Description |
|----------|-------------|
| `Hugo_MSA(which)` | Load reference MSA (`:A`, `:B`, `:C`, or `:D`) of Jacquin et al. (2016) |
| `Hugo_MSA_pnat(which)` | Load reference pnat values |
| `Hugo_MSA_path(which)` | Path of the file with both |

## Citation

If you use this package, please cite:

Fernandez-de-Cossio-Diaz, Jorge, Clément Roussel, Simona Cocco, and Remi Monasson. ["Accelerated Sampling with Stacked Restricted Boltzmann Machines."](https://openreview.net/forum?id=kXNJ48Hvw1) *The Twelfth International Conference on Learning Representations* (2024).

Citation metadata is available in [CITATION.cff](CITATION.cff).

## References

- Jacquin, Hugo, et al. ["Benchmarking inverse statistical approaches for protein structure and design with exactly solvable models."](https://doi.org/10.1371/journal.pcbi.1004889) *PLoS computational biology* 12.5 (2016): e1004889.
