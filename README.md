# LatticeProteins.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cossio.github.io/LatticeProteins.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cossio.github.io/LatticeProteins.jl/dev/)

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

In this model, a protein is represented as a chain of 27 amino acids that folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is determined by pairwise contacts between amino acids using the Miyazawa-Jernigan interaction matrix. The probability that a sequence folds into a given "native" structure is computed via a Boltzmann distribution over all 103,346 compact structures of the lattice (or, optionally, over the representative subset of 10,000 structures selected by the Shakhnovich group and used by Jacquin et al. (2016)).

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

# Compete only against the subset of 10,000 structures used by Jacquin et al. (2016),
# as in the reference MSAs (structure indices refer to the set in use)
a = structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)   # 625
p_ref = pnat(a, seq; structures=SHAKHNOVICH_STRUCTURES)
msa_ref = random_msa(a; nseqs=100, β=100, nsteps=1000, structures=SHAKHNOVICH_STRUCTURES)

# Load reference MSA datasets from Jacquin et al. (2016)
sequences = Hugo_MSA(:A)
pnat_values = Hugo_MSA_pnat(:A)
```

## API overview

### Structures and contact maps

| Symbol | Description |
|--------|-------------|
| `STRUCTURES` | Contact maps of all 103,346 compact structures (the default set) |
| `N_STRUCTURES` | Number of compact structures (103,346) |
| `SHAKHNOVICH_STRUCTURES` | Contact maps of the representative subset of 10,000 structures (Heo, Maslov & Shakhnovich 2011) used by Jacquin et al. (2016) |
| `CONTACT_MAP_A`, `CONTACT_MAP_B`, `CONTACT_MAP_C`, `CONTACT_MAP_D` | Indices in `STRUCTURES` of four reference structures from Jacquin et al. (2016) |
| `structure_index(structures, cm)` | Index of a contact map within a set of structures |
| `enumerate_structures()` | Regenerate the contact maps of all compact structures |

### Energy and folding probability

| Function | Description |
|----------|-------------|
| `energy(cm, seq)` | Contact energy of sequence `seq` in structure `cm` |
| `energies(seq)` | Contact energies of `seq` in every structure |
| `pnat(cm, seq)` | Probability that `seq` folds into structure `cm` |
| `log_pnat(cm, seq)` | Log-probability of folding |

All of these take a `structures` keyword argument selecting the set of competing structures (default: `STRUCTURES`).

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
| `metropolis!(seq, cm; β, nsteps, structures)` | In-place Metropolis sampling of a sequence |
| `random_msa(cm; nseqs, β, nsteps, structures)` | Generate a random MSA via Metropolis sampling |
| `parallel_tempering(cm; betas, nseqs, burnin, thin, local_steps, structures)` | Replica-exchange sampling with swap and round-trip diagnostics |

For strong selection, run replicas over a fixed ladder from weak to strong selection:

```julia
using LatticeProteins: parallel_tempering, CONTACT_MAP_A
using Random: Xoshiro

# Starting ladder for a pilot; inspect diagnostics before a production run.
betas = [0.0; exp.(range(log(0.1), log(3000.0); length=31))]
betas[end] = 3000.0
result = parallel_tempering(CONTACT_MAP_A; betas, burnin=5000,
    nseqs=500, thin=10, local_steps=10, rng=Xoshiro(42))
msa = result.samples                         # only the beta=3000 samples
swap_rates = result.diagnostics.swap_accepts ./ result.diagnostics.swap_attempts
result.diagnostics.round_trips               # hot → cold → hot, per replica ID
```

Burn-in and thinning are measured in exchange rounds, each containing
`local_steps` mutation proposals per replica. The example uses 50,000 burn-in
proposals per replica. The target distribution is proportional to `pnat^beta`;
the default folding model includes all 103,346 compact structures. To use the
paper's 10,000-structure subset, pass `structures=SHAKHNOVICH_STRUCTURES` and translate
the target with `structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)`. Inspect temperature
traversal and compare independent runs: a high `pnat` alone does not demonstrate
equilibration, and saved samples are not necessarily independent.

The bundled MJ interactions allow exact integer contact-energy caching and a
Boltzmann lookup table, used automatically by `parallel_tempering`. This speeds up
full-structure sampling without changing the target distribution. Interactions
that do not lie on the supported grid use the floating-point kernel; set
`energy_cache=:float` to select it explicitly. The result reports the cache used.

The [Tubiana figure scripts](docs/figures/tubiana2019/README.md) run four independent
ladders, resume from checkpoints, and continue until the monitored convergence
criteria pass, including sequence frequencies and principal components.
The full-structure figure uses the resulting 316,000 saved PT draws, with 172,717
retained at `P_nat > 0.99`; the alignment and convergence report are included in
`docs/src/assets/` so the figure can be redrawn directly from these samples.

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
| `Hugo_MSA(which)` | Load reference MSA (`:A`, `:B`, `:C`, or `:D`) |
| `Hugo_MSA_pnat(which)` | Load reference pnat values |
| `MJ` | 20×20 Miyazawa-Jernigan interaction matrix |

## Citation

If you use this package, please cite:

Fernandez-de-Cossio-Diaz, Jorge, Clément Roussel, Simona Cocco, and Remi Monasson. ["Accelerated Sampling with Stacked Restricted Boltzmann Machines."](https://openreview.net/forum?id=kXNJ48Hvw1) *The Twelfth International Conference on Learning Representations* (2024).

Citation metadata is available in [CITATION.cff](CITATION.cff).

## References

- Jacquin, Hugo, et al. ["Benchmarking inverse statistical approaches for protein structure and design with exactly solvable models."](https://doi.org/10.1371/journal.pcbi.1004889) *PLoS computational biology* 12.5 (2016): e1004889.
