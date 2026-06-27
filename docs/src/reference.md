# Reference

```@meta
CurrentModule = LatticeProteins
```

This page documents all symbols provided by `LatticeProteins.jl`. None of them are
exported, so they must be accessed with the `LatticeProteins.` prefix or imported
explicitly (e.g. `using LatticeProteins: pnat`).

## Structures and contact maps

```@docs
CONTACT_MAPS
N_STRUCTURES
N_CONTACTS
L
CONTACT_MAP_A
CONTACT_MAP_B
CONTACT_MAP_C
CONTACT_MAP_D
load_contact_maps
```

## Energy and folding probability

```@docs
energy
pnat
log_pnat
```

## Sequence encoding

```@docs
AMINO_ACIDS
onehot
potts
aaseq
```

## Sampling

```@docs
metropolis!
random_msa
```

## Miyazawa-Jernigan matrix

```@docs
MJ
load_miyazawa_jernigan_matrix
load_miyazawa_jernigan_aminoacids
```

## Reference data

```@docs
Hugo_MSA
Hugo_MSA_pnat
Hugo_MSA_path
```
