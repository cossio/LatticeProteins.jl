# Validates an integer-coded sequence.
function check_sequence(seq::AbstractVector{<:Integer})
    Base.require_one_based_indexing(seq)
    length(seq) == L || throw(ArgumentError("sequence must have length $L"))
    return all(a -> 1 ≤ a ≤ 20, seq) || throw(ArgumentError("amino acid indices must be in 1:20"))
end

# Contact energy of `seq` in structure `s` of `S`, in units of 0.01 (see `MJ_INT`). Assumes
# `seq` has been validated with `check_sequence` and that `s` is a valid structure index.
function structure_energy(seq::AbstractVector{<:Integer}, S::StructureSet, s::Integer)
    contacts = S.contacts
    e = Int32(0)
    @inbounds for k in axes(contacts, 2)
        e += MJ_INT[seq[contacts[1, k, s]], seq[contacts[2, k, s]]]
    end
    return e
end

# Writes the energies of `seq` in all structures of `S` into `E`, in units of 0.01. Assumes
# `seq` has been validated with `check_sequence`.
function energies!(E::Vector{Int32}, seq::AbstractVector{<:Integer}, S::StructureSet)
    length(E) == length(S) || throw(DimensionMismatch("one energy per structure is needed"))
    for s in eachindex(E)
        E[s] = structure_energy(seq, S, s)
    end
    return E
end

# -log(pnat) of structure `s`, given the energies `E` of a sequence in all structures of a
# set (in units of 0.01) and the Boltzmann factors `boltzmann` of that set. Energies are
# shifted by their minimum before summing the Boltzmann factors, and the native energy is
# centred before adding log(z), which avoids cancellation when pnat is close to one.
function folding_loss(E::Vector{Int32}, s::Integer, boltzmann::Vector{Float64})
    emin = minimum(E)
    z = 0.0
    @inbounds @simd for e in E
        z += boltzmann[Int(e - emin) + 1]
    end
    return (E[s] - emin) / 100 + log(z)
end

"""
    energy(s, seq; structures = STRUCTURES)

Contact energy of amino acid sequence `seq` in structure `s`: the sum of the
Miyazawa-Jernigan interaction energies ([`MJ`](@ref)) over the pairs of residues in
contact in that structure.

`seq` is a vector of length 27 of amino acid indices in `1:20` (see [`potts`](@ref)). `s`
is the index of the structure in the [`StructureSet`](@ref) `structures`, which defaults
to all compact structures [`STRUCTURES`](@ref).
"""
function energy(s::Integer, seq::AbstractVector{<:Integer}; structures::StructureSet = STRUCTURES)
    check_sequence(seq)
    check_structure(structures, s)
    return structure_energy(seq, structures, s) / 100
end

"""
    energies(seq; structures = STRUCTURES)

Contact energies of `seq` in every structure of `structures`, as a `Vector{Float64}` of
length `length(structures)`. Equivalent to `[energy(s, seq; structures) for s in
1:length(structures)]`, but faster.

See also [`energy`](@ref).
"""
function energies(seq::AbstractVector{<:Integer}; structures::StructureSet = STRUCTURES)
    check_sequence(seq)
    return energies!(Vector{Int32}(undef, length(structures)), seq, structures) ./ 100
end

"""
    log_pnat(s, seq; structures = STRUCTURES)

Log-probability that sequence `seq` folds into structure `s` (the "native" structure), as
opposed to any of the other structures in `structures`.

The probability is a Boltzmann distribution at temperature 1 over the structures in
`structures`:

```math
\\log P_\\text{nat}(s \\mid \\text{seq}) = -\\log \\sum_{s'} \\exp\\bigl(E(s, \\text{seq}) - E(s', \\text{seq})\\bigr)
```

By default `structures` is [`STRUCTURES`](@ref), all 103,346 compact structures of the
lattice. Pass `structures = shakhnovich_structures()` to sum only over the representative
subset of 10,000 structures used by Jacquin et al. (2016) (see
[`shakhnovich_structures`](@ref)); note that structure indices refer to the set in use.

See also [`pnat`](@ref), [`energy`](@ref).
"""
function log_pnat(s::Integer, seq::AbstractVector{<:Integer}; structures::StructureSet = STRUCTURES)
    check_sequence(seq)
    check_structure(structures, s)
    E = energies!(Vector{Int32}(undef, length(structures)), seq, structures)
    return -folding_loss(E, s, structures.boltzmann)
end

"""
    pnat(s, seq; structures = STRUCTURES)

Probability that sequence `seq` folds into structure `s`.
Equivalent to `exp(log_pnat(s, seq; structures))`.

See also [`log_pnat`](@ref).
"""
pnat(s::Integer, seq::AbstractVector{<:Integer}; structures::StructureSet = STRUCTURES) = exp(log_pnat(s, seq; structures))
