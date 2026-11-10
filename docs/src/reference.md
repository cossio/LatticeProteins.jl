# Reference

```@meta
CurrentModule = LatticeProteins
```

```@docs
LatticeProteins
```

## Public API

The public API is the model in which all 103,346 compact structures of the lattice compete
to fold a sequence. It consists of the names below, which are exported.

### Sequences

Sequences are strings of amino acids, integer-coded vectors ([`potts`](@ref)) or one-hot
encoded arrays ([`onehot`](@ref)). Energies and samplers work with integer-coded
sequences, of length 27.

```@docs
potts
aaseq
onehot
```

### Structures

A structure is referred to by its index among all compact structures.

```@docs
enumerate_structures
contact_map
STRUCTURE_A
STRUCTURE_B
STRUCTURE_C
STRUCTURE_D
```

### Energy and folding probability

```@docs
energy
pnat
logpnat
```

### Sampling

```@docs
sample_msa
Metropolis
ParallelTempering
```

## Internals

!!! warning "Not part of the public API"
    The names in this section are internal: they are documented because the documentation,
    the figure scripts and the tests use them, but they can change or disappear in any
    release.

### Other sets of competing structures

Internally, the structures that compete to fold a sequence form a [`StructureSet`](@ref).
Every function of the public API that takes a structure index also accepts a keyword
argument `structures`, by default [`enumerate_structures`](@ref)`()`; the index then refers
to that set.
This is how the documentation reproduces the model of Jacquin et al. (2016), where only a
subset of 10,000 structures compete ([`shakhnovich_structures`](@ref)).

```@docs
StructureSet
shakhnovich_structures
structure_index
energies
N_CONTACTS
```

### Sequences, interactions and moves

```@docs
AMINO_ACIDS
L
MJ
metropolis!
```

### Continuing a parallel-tempering run

The figure scripts extend a [`ParallelTempering`](@ref) run from a checkpoint, with three
internal additions to its interface:

- `ParallelTempering(; ..., threaded = false, rng_mode = :shared)` uses the `rng` given
  to `sample_msa` for both local updates and swaps, without creating a stream per
  temperature. The default is `rng_mode = :per_temperature`, which creates fresh local
  streams on each call and cannot be resumed.
- `sample_msa(s, method; ..., initial_seqs)` starts the replicas from the columns of
  `initial_seqs`, a `(27, length(betas))` integer matrix ordered hot to cold.
- The result has the fields `final_seqs` and `final_logpnat`, the final states at every
  temperature (and `betas`, `rng_mode`, which repeat the method's).

A call with `burnin = 0` that passes the same RNG and the `final_seqs` of the previous
call as `initial_seqs` continues its trajectory exactly, provided the previous call
ended on an even number of rounds (which preserves the alternation of the swap edges).
The shared stream is also the RNG protocol of v2.1.0 and of the saved Tubiana alignment.

### Reference alignments of Jacquin et al. (2016)

These alignments were sampled against the subset of 10,000 structures.

```@docs
Hugo_MSA
Hugo_MSA_pnat
Hugo_MSA_path
```
