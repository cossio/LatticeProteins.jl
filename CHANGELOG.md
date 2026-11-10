# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.0.0

The package now has a small declared public API, for the model in which all 103,346
compact structures compete, and a single energy kernel shared by all samplers.

### Breaking

- **Public API.** The public names are declared `public`: `AMINO_ACIDS`, `L`, `potts`,
  `aaseq`, `onehot`; `N_STRUCTURES`, `contact_map`, `STRUCTURE_A`, ..., `STRUCTURE_D`;
  `MJ`, `energy`, `pnat`, `log_pnat`; `metropolis!`, `random_msa`, `parallel_tempering`.
  Everything else is internal and can change in any release.
- **Structures** are referred to by their index in `1:N_STRUCTURES`, and
  `contact_map(s)` gives the contacts of structure `s`. The array `STRUCTURES` and the
  function `enumerate_structures` are no longer part of the API.
- **Reference structures.** `CONTACT_MAP_A`, ..., `CONTACT_MAP_D` are renamed
  `STRUCTURE_A`, ..., `STRUCTURE_D`: they are indices of structures, not contact maps.
- **Other sets of competing structures are internal.** This covers the `structures`
  keyword argument, the subset of 10,000 structures of Jacquin et al. (2016), now
  `LatticeProteins.shakhnovich_structures()` instead of the constant
  `SHAKHNOVICH_STRUCTURES`, `structure_index`, and the reference alignments
  `Hugo_MSA`, `Hugo_MSA_pnat`, `Hugo_MSA_path`, which were sampled with that subset.
  They remain available, and are used by the documentation and the tests, but without
  stability guarantees. A set of structures is a `LatticeProteins.StructureSet`
  instead of a bare `(ncontacts, nstructures, 2)` array. The subset is downloaded the
  first time it is requested, no longer when the package is installed.
- **`metropolis!` and `random_msa`** make the same moves as `parallel_tempering`: a step
  proposes one of the 19 *other* amino acids at a random position, where it used to
  propose any of the 20 (leaving the sequence unchanged one time in 20). The stationary
  distribution is the same; seeded trajectories differ from 2.x.
- **Errors.** Invalid arguments throw an `ArgumentError`: sequences of the wrong length
  or with amino acids outside `1:20` (previously an `AssertionError`), and structure
  indices outside `1:N_STRUCTURES` (previously a `BoundsError`).
- **`potts`** throws an `ArgumentError` for characters outside `AMINO_ACIDS`, which it
  used to code silently as cysteine (index 1), and for sequences of unequal lengths.
- **Removed** `load_contact_maps`, `load_miyazawa_jernigan_matrix` and
  `load_miyazawa_jernigan_aminoacids`: use `MJ` and `AMINO_ACIDS`.
- `energies` and `N_CONTACTS` are internal.

| 2.x | 3.0 |
|:--|:--|
| `CONTACT_MAP_A` | `STRUCTURE_A` |
| `STRUCTURES[:, s, :]` | `contact_map(s)` |
| `size(STRUCTURES, 2)` | `N_STRUCTURES` |
| `SHAKHNOVICH_STRUCTURES` | `LatticeProteins.shakhnovich_structures()` (internal) |
| `structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)` | `LatticeProteins.structure_index(LatticeProteins.shakhnovich_structures(), STRUCTURE_A)` (internal) |

### Added

- `contact_map`.
- `metropolis!` and `random_msa` take a random number generator `rng`.
- `parallel_tempering` has an explicit serial `rng_mode=:shared`, which uses a single
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
  periodic recomputation of `parallel_tempering` is gone. The floating-point kernel is
  removed. `energy`, `energies`, `log_pnat` and `pnat` still return `Float64`, equal to
  their previous values up to roundoff.
- `metropolis!` and `parallel_tempering` are about three times faster on the full
  structure set. A seeded `parallel_tempering` run makes the same
  proposals as in 2.1.1, and can differ only if an acceptance test is decided by
  roundoff; the saved Tubiana alignment is reproduced exactly (checked on its first
  44,000 draws).
- Structures are stored as `Int8`: the precompiled package is 8 MB instead of 52 MB.
- `potts` of strings and `onehot` of a vector of strings are roughly 100 times faster
  on an alignment of 46,000 sequences.

### Fixed

- `parallel_tempering` no longer reads out of bounds when given (internally) structures
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
