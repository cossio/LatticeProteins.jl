# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.0.0

The package now has a small declared public API, for the model in which all 103,346
compact structures compete, and a single energy kernel shared by all samplers.

### Breaking

- **Public API.** The package now exports its public names: `potts`, `aaseq`, `onehot`;
  `enumerate_structures`, `contact_map`, `STRUCTURE_A`, ..., `STRUCTURE_D`; `MJ`, `energy`,
  `pnat`, `log_pnat`; `sample_msa`, `Metropolis`, `ParallelTempering`, `metropolis!`.
  Everything else is internal and can change in any release. Since 2.x exported nothing,
  `using LatticeProteins` now brings these names into scope.
- **Structures.** `enumerate_structures()` returns the set of all compact structures,
  computed once, instead of recomputing an array at each call.
  `length(enumerate_structures())` is their number and `contact_map(s)` gives the
  contacts of structure `s`; these replace `N_STRUCTURES` and indexing `STRUCTURES`.
- **Reference structures.** `CONTACT_MAP_A`, ..., `CONTACT_MAP_D` are renamed
  `STRUCTURE_A`, ..., `STRUCTURE_D`: they are indices of structures, not contact maps.
- **Sampling.** `sample_msa(s, method; nseqs, rng)` replaces both `random_msa` and
  `parallel_tempering`. The method is `Metropolis(; β, nsteps)` or
  `ParallelTempering(; betas, burnin, thin, local_steps, threaded, rng_mode)`, and holds
  the parameters of the sampler. The result is a named tuple whose field `samples` is
  the alignment; with `ParallelTempering` it has the fields that `parallel_tempering`
  returned. `nseqs` defaults to 100 for both methods.
- **The subset of 10,000 structures** of Jacquin et al. (2016) is not part of the public
  API, nor is anything that only serves it: the `structures` keyword argument,
  `structure_index`, and the reference alignments `Hugo_MSA`, `Hugo_MSA_pnat`,
  `Hugo_MSA_path`, which were sampled with that subset. The documentation and the tests
  still use them. The subset is `LatticeProteins.shakhnovich_structures()` instead of
  the constant `SHAKHNOVICH_STRUCTURES`, a `LatticeProteins.StructureSet` instead of a
  bare `(ncontacts, nstructures, 2)` array, and is downloaded the first time it is
  requested instead of when the package is installed.
- **Metropolis moves** are those of parallel tempering, in `metropolis!` and in
  `sample_msa` with either method: a step
  proposes one of the 19 *other* amino acids at a random position, where it used to
  propose any of the 20 (leaving the sequence unchanged one time in 20). The stationary
  distribution is the same; seeded trajectories differ from 2.x.
- **Errors.** Invalid arguments throw an `ArgumentError`: sequences of the wrong length
  or with amino acids outside `1:20` (previously an `AssertionError`), and structure
  indices that are out of range (previously a `BoundsError`).
- **`potts`** throws an `ArgumentError` for characters outside the alphabet, which it
  used to code silently as cysteine (index 1), and for sequences of unequal lengths.
- **Removed** `load_contact_maps`, `load_miyazawa_jernigan_matrix` and
  `load_miyazawa_jernigan_aminoacids`: use `MJ`.

| 2.x | 3.0 |
|:--|:--|
| `CONTACT_MAP_A` | `STRUCTURE_A` |
| `random_msa(s; nseqs, β, nsteps)` | `sample_msa(s, Metropolis(; β, nsteps); nseqs).samples` |
| `parallel_tempering(s; betas, nseqs, ...)` | `sample_msa(s, ParallelTempering(; betas, ...); nseqs)` |
| `STRUCTURES[:, s, :]` | `contact_map(s)` |
| `N_STRUCTURES`, `size(STRUCTURES, 2)` | `length(enumerate_structures())` |
| `SHAKHNOVICH_STRUCTURES` | `LatticeProteins.shakhnovich_structures()` (internal) |
| `structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)` | `LatticeProteins.structure_index(LatticeProteins.shakhnovich_structures(), STRUCTURE_A)` (internal) |

### Added

- `contact_map`.
- `metropolis!` and `sample_msa` take a random number generator `rng`.
- `ParallelTempering` has an explicit serial `rng_mode=:shared`, which uses a single
  RNG stream: a run can then be continued exactly from the final states of a previous
  call. It is also the RNG protocol of v2.1.0 and of the saved Tubiana alignment.
  Per-temperature RNGs and threaded updates remain the defaults.
- Checkpointed continuation and convergence diagnostics for the full-model Tubiana
  figure, including independent-ladder R-hat, effective sample sizes,
  conditional-frequency uncertainties, and downloadable trajectories. The
  full-structure figure uses the PT sample that passes all monitored convergence
  criteria: 316,000 draws, with 172,717 at `P_nat > 0.99`.

### Changed

- **One energy kernel.** Contact energies are computed as exact integers in units of
  0.01, the resolution of the MJ matrix, and Boltzmann factors are read from a table.
  Energies updated incrementally during sampling therefore never drift, and the
  periodic recomputation of parallel tempering is gone. The floating-point kernel is
  removed. `energy`, `energies`, `log_pnat` and `pnat` still return `Float64`, equal to
  their previous values up to roundoff.
- Metropolis and parallel tempering are about three times faster on the full
  structure set. A seeded parallel tempering run makes the same
  proposals as in 2.1.1, and can differ only if an acceptance test is decided by
  roundoff; the saved Tubiana alignment is reproduced exactly (checked on its first
  44,000 draws).
- Structures are stored as `Int8`: the precompiled package is 8 MB instead of 52 MB.
- `potts` of strings and `onehot` of a vector of strings are roughly 100 times faster
  on an alignment of 46,000 sequences.

### Fixed

- Parallel tempering no longer reads out of bounds when given (internally) structures
  with residues outside `1:27`.

## 2.1.1

- Parallelize local replica updates with temperature-specific RNG streams, producing
  identical results across thread counts and with `threaded=false`. Seeded trajectories
  differ from v2.1.0. Share the stable folding-loss calculation across samplers and
  probability evaluation, and add an exact-distribution sampling check.

## 2.1.0

- Add full-model parallel tempering with reproducible RNGs, cached mutation
  energies, temperature traces, and separate burn-in/production exchange and
  round-trip diagnostics. Include an equal-budget Metropolis comparison script.

## 2.0.0

- **Breaking:** folding probabilities are computed against all 103,346 compact
  structures of the lattice by default (`STRUCTURES`, generated by
  `enumerate_structures`), instead of the subset of 10,000 structures used by Jacquin
  et al. (2016). `energy`, `pnat`, `log_pnat`, `metropolis!` and `random_msa` take a
  `structures` keyword argument; pass `structures = SHAKHNOVICH_STRUCTURES` (the array
  previously named `CONTACT_MAPS`) to compete against the subset as before.
  `CONTACT_MAP_A`, ..., `CONTACT_MAP_D` and `N_STRUCTURES` now refer to the full set;
  `structure_index` translates structure indices between sets.
- Add `energies`, the contact energies of a sequence in every structure. `metropolis!`
  updates them incrementally instead of recomputing them at every step.
- Add the documentation, with a tutorial and a reproduction of Fig. 4 of Tubiana et al.
  (2019).
- Require Julia 1.12.

## 1.0.1

- Add docstrings and tests. `random_msa` accepts a non-integer `β`, and `metropolis!`
  checks that amino acid indices are in `1:20`.

## 1.0.0 - 2023-08-03

- No changes.

## 0.1.0

- Register at https://github.com/cossio/CossioJuliaRegistry.
