# Reference

```@meta
CurrentModule = LatticeProteins
```

This page documents the public API of `LatticeProteins.jl`: the names declared `public`
in the package. None of them are exported, so they must be accessed with the
`LatticeProteins.` prefix or imported explicitly (e.g. `using LatticeProteins: pnat`).
Names that are not listed here are internal and may change in any release.

```@docs
LatticeProteins
```

## Sequences

Sequences are strings over the alphabet [`AMINO_ACIDS`](@ref), integer-coded vectors
([`potts`](@ref)) or one-hot encoded arrays ([`onehot`](@ref)). Energies and samplers
work with integer-coded sequences.

```@docs
AMINO_ACIDS
potts
onehot
aaseq
```

## Structures

A [`StructureSet`](@ref) is a set of structures competing to fold a sequence. Functions
take it as the keyword argument `structures`, and refer to a structure by its index in
the set.

```@docs
StructureSet
STRUCTURES
shakhnovich_structures
contact_map
structure_index
enumerate_structures
STRUCTURE_A
STRUCTURE_B
STRUCTURE_C
STRUCTURE_D
L
N_CONTACTS
N_STRUCTURES
```

## Energy and folding probability

```@docs
MJ
energy
energies
pnat
log_pnat
```

## Sampling

```@docs
metropolis!
random_msa
parallel_tempering
```

## Reference data

```@docs
Hugo_MSA
Hugo_MSA_pnat
Hugo_MSA_path
```
