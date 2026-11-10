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
MJ
energy
pnat
log_pnat
```

### Sampling

```@docs
metropolis!
sample_msa
parallel_tempering
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

### Sequences

```@docs
AMINO_ACIDS
L
```

### Reference alignments of Jacquin et al. (2016)

These alignments were sampled against the subset of 10,000 structures.

```@docs
Hugo_MSA
Hugo_MSA_pnat
Hugo_MSA_path
```
