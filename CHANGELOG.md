# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 3.0.0

The package now has a declared public API, built around a `StructureSet` type, and a
single energy kernel shared by all samplers.

### Breaking

- **Structure sets.** A set of competing structures is a `StructureSet`, built from a
  vector of contact maps (`ncontacts × 2` matrices), instead of a bare
  `(ncontacts, nstructures, 2)` array. `STRUCTURES` is a `StructureSet`, as is the value
  of `enumerate_structures()`, and the `structures` keyword argument of every function
  takes one. Use `length(S)` for the number of structures and `contact_map(S, s)` for
  what was `S[:, s, :]`.
- **The subset of 10,000 structures** is `shakhnovich_structures()` instead of the
  constant `SHAKHNOVICH_STRUCTURES`. It is downloaded the first time it is requested,
  no longer when the package is installed. `load_contact_maps` is removed.
- **Reference structures.** `CONTACT_MAP_A`, ..., `CONTACT_MAP_D` are renamed
  `STRUCTURE_A`, ..., `STRUCTURE_D`: they are indices of structures, not contact maps.
- **`metropolis!` and `random_msa`** make the same moves as `parallel_tempering`: a step
  proposes one of the 19 *other* amino acids at a random position, where it used to
  propose any of the 20 (leaving the sequence unchanged one time in 20). The stationary
  distribution is the same; seeded trajectories differ from 2.x.
- **Errors.** Invalid arguments throw an `ArgumentError`: sequences of the wrong length
  or with amino acids outside `1:20` (previously an `AssertionError`), structure indices
  outside the set (previously a `BoundsError`), and contact maps with residues outside
  `1:27`, which `StructureSet` rejects on construction.
- **`potts`** throws an `ArgumentError` for characters outside `AMINO_ACIDS`, which it
  used to code silently as cysteine (index 1), and for sequences of unequal lengths.
- **Removed** `load_miyazawa_jernigan_matrix` and `load_miyazawa_jernigan_aminoacids`:
  use `MJ` and `AMINO_ACIDS`.
- **Public API.** The names documented in the reference are declared `public`.
  Everything else is internal and can change in any release.

| 2.x | 3.0 |
|:--|:--|
| `SHAKHNOVICH_STRUCTURES` | `shakhnovich_structures()` |
| `CONTACT_MAP_A` | `STRUCTURE_A` |
| `STRUCTURES[:, s, :]` | `contact_map(STRUCTURES, s)` |
| `size(STRUCTURES, 2)` | `length(STRUCTURES)` |
| `structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)` | `structure_index(shakhnovich_structures(), STRUCTURE_A)` |
| `structure_index(STRUCTURES, SHAKHNOVICH_STRUCTURES[:, s, :])` | `structure_index(STRUCTURES, s; from = shakhnovich_structures())` |
| `pnat(cm, seq; structures = array)` | `pnat(s, seq; structures = StructureSet(contact_maps))` |

### Added

- `StructureSet`, `contact_map`, `shakhnovich_structures`, and the `from` keyword of
  `structure_index`.
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
- A structure set stores its per-site contact index, so sampling repeatedly on a set
  other than `STRUCTURES` no longer rebuilds it on every call.
- Structures are stored as `Int8`: the precompiled package is 8 MB instead of 52 MB.
- `potts` of strings and `onehot` of a vector of strings are roughly 100 times faster
  on an alignment of 46,000 sequences.

### Fixed

- `parallel_tempering` no longer reads out of bounds for structures with residues
  outside `1:27` (`StructureSet` rejects them).

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
