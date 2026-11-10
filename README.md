# LatticeProteins.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cossio.github.io/LatticeProteins.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cossio.github.io/LatticeProteins.jl/dev/)
[![Coverage](https://codecov.io/gh/cossio/LatticeProteins.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cossio/LatticeProteins.jl)

Simulate the lattice proteins model on a 3×3×3 cubic lattice in Julia.

In this model, a protein is represented as a chain of 27 amino acids that folds on a 3×3×3 cubic lattice. The energy of a sequence in a given structure is determined by pairwise contacts between amino acids using the Miyazawa-Jernigan interaction matrix. The probability that a sequence folds into a given "native" structure is computed via a Boltzmann distribution over all 103,346 compact structures of the lattice.

## Installation

This package is registered. Install it using the Julia package manager:

```julia
import Pkg
Pkg.add("LatticeProteins")
```

## Quick start

The names in the tables below are the public API, and are exported. Anything else is
internal and may change in any release.

```julia
using LatticeProteins

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
msa = sample_msa(STRUCTURE_A, Metropolis(β=100, nsteps=1000); nseqs=100).samples

# The 28 contacts of structure A, as pairs of residue positions
cm = contact_map(STRUCTURE_A)
```

## API overview

### Structures

A structure is referred to by its index, from 1 to 103,346.

| Symbol | Description |
|--------|-------------|
| `enumerate_structures()` | The set of all compact structures; `length(enumerate_structures())` is 103,346 |
| `contact_map(s)` | Contact map of structure `s`: a `28 × 2` matrix of pairs of residue positions |
| `STRUCTURE_A`, `STRUCTURE_B`, `STRUCTURE_C`, `STRUCTURE_D` | Indices of four reference structures from Jacquin et al. (2016) |

### Energy and folding probability

| Function | Description |
|----------|-------------|
| `energy(s, seq)` | Contact energy of sequence `seq` in structure `s` |
| `pnat(s, seq)` | Probability that `seq` folds into structure `s` |
| `logpnat(s, seq)` | Log-probability of folding |

### Sequence encoding

| Function | Description |
|----------|-------------|
| `onehot(seq)` | One-hot encode amino acid sequence(s) |
| `potts(seq)` | Convert to integer-coded (Potts) representation: indices in the alphabet `"CMFILVWYAGTSNQDEHRKP-"` |
| `aaseq(X)` | Convert from one-hot or integer coding back to amino acid string(s) |

### Sampling

| Function | Description |
|----------|-------------|
| `sample_msa(s, method; nseqs, rng)` | Sample an MSA of sequences folding into structure `s`, with the given method; the result's `samples` field is the alignment |
| `Metropolis(; β, nsteps)` | Method: independent Metropolis chains, one per sequence |
| `ParallelTempering(; betas, burnin, thin, local_steps, threaded)` | Method: replica exchange, with swap and round-trip diagnostics |

For strong selection, run replicas over a fixed ladder from weak to strong selection:

```julia
using LatticeProteins
using Random: Xoshiro

# Starting ladder for a pilot; inspect diagnostics before a production run.
betas = [0.0; exp.(range(log(0.1), log(3000.0); length=31))]
betas[end] = 3000.0
method = ParallelTempering(; betas, burnin=5000, thin=10, local_steps=10)
result = sample_msa(STRUCTURE_A, method; nseqs=500, rng=Xoshiro(42))
msa = result.samples                         # only the beta=3000 samples
swap_rates = result.diagnostics.swap_accepts ./ result.diagnostics.swap_attempts
result.diagnostics.round_trips               # hot → cold → hot, per replica ID
```

Burn-in and thinning are measured in exchange rounds, each containing
`local_steps` mutation proposals per replica. The example uses 50,000 burn-in
proposals per replica. The target distribution is proportional to `pnat^beta`.
Inspect temperature
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

### Internals

Internally, the functions above take a keyword argument `structures` that replaces the
set of competing structures. The documentation and the tests use it with the subset of
10,000 structures of Jacquin et al. (2016) (`LatticeProteins.shakhnovich_structures()`),
to compare with published results and with the reference alignments of that paper
(`LatticeProteins.Hugo_MSA`). None of this is part of the public API; see the
[reference](https://cossio.github.io/LatticeProteins.jl/dev/reference/#Internals).

## Citation

If you use this package, please cite:

Fernandez-de-Cossio-Diaz, Jorge, Clément Roussel, Simona Cocco, and Remi Monasson. ["Accelerated Sampling with Stacked Restricted Boltzmann Machines."](https://openreview.net/forum?id=kXNJ48Hvw1) *The Twelfth International Conference on Learning Representations* (2024).

Citation metadata is available in [CITATION.cff](CITATION.cff).

## References

- Jacquin, Hugo, et al. ["Benchmarking inverse statistical approaches for protein structure and design with exactly solvable models."](https://doi.org/10.1371/journal.pcbi.1004889) *PLoS computational biology* 12.5 (2016): e1004889.
