"""
    energy(cm::Int, seq::AbstractVector{Int}; structures = STRUCTURES)

Compute the contact energy of amino acid sequence `seq` in structure `cm`.

The energy is the sum of Miyazawa-Jernigan interaction energies over all contacting
residue pairs in the given structure. `seq` must be a vector of length 27 with
integer amino acid indices in `1:20`. `cm` is the index of the structure in `structures`,
which defaults to all compact structures [`STRUCTURES`](@ref) and can be set to the
subset [`SHAKHNOVICH_STRUCTURES`](@ref).
"""
function energy(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = STRUCTURES)
    check_sequence(seq)
    check_structures(structures)
    checkbounds(structures, :, cm, :)
    return structure_energy(seq, structures, cm)
end

# Contact energy of `seq` in structure `s` of `structures`. Assumes `seq` has been validated
# with `check_sequence` and that `s` is a valid structure index.
function structure_energy(seq::AbstractVector{Int}, structures::AbstractArray{Int,3}, s::Int)
    e = 0.0
    for k in axes(structures, 1)
        @inbounds i, j = structures[k, s, 1], structures[k, s, 2]
        a, b = seq[i], seq[j]
        @inbounds e += MJ[a, b] # a, b ∈ 1:20 by check_sequence
    end
    return e
end

"""
    energies(seq::AbstractVector{Int}; structures = STRUCTURES)

Contact energies of `seq` in every structure of `structures`, as a `Vector{Float64}` of
length `size(structures, 2)`. Equivalent to `[energy(s, seq; structures) for s in
1:size(structures, 2)]`, but faster.

See also [`energy`](@ref).
"""
function energies(seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = STRUCTURES)
    check_sequence(seq)
    check_structures(structures)
    return energies!(Vector{Float64}(undef, size(structures, 2)), seq, structures)
end

# Writes the energies of `seq` in all `structures` into `E`. Assumes `seq` and `structures`
# have been validated with `check_sequence` and `check_structures`.
function energies!(E::Vector{Float64}, seq::AbstractVector{Int}, structures::AbstractArray{Int,3})
    @assert length(E) == size(structures, 2)
    for s in axes(structures, 2)
        E[s] = structure_energy(seq, structures, s)
    end
    return E
end

# -log(pnat), shared by probability evaluation and both samplers. Shift before
# summing exponentials, and center the native energy before adding log(z) to
# avoid cancellation when the native probability is close to one.
function _folding_loss(E, cm)
    emin = minimum(E)
    z = 0.0
    for e in E
        z += exp(emin - e)
    end
    return (E[cm] - emin) + log(z)
end

"""
    log_pnat(cm::Int, seq::AbstractVector{Int}; structures = STRUCTURES)

Compute the log-probability that sequence `seq` folds into structure `cm` (the
"native" structure), as opposed to any of the other structures in `structures`.

The probability is computed using a Boltzmann distribution at temperature 1 over
the structures in `structures`:

```math
\\log P_\\text{nat}(\\text{cm} \\mid \\text{seq}) = -\\log \\sum_s \\exp\\bigl(E(\\text{cm}, \\text{seq}) - E(s, \\text{seq})\\bigr)
```

By default `structures` is [`STRUCTURES`](@ref), all 103,346 compact structures of the
lattice. Pass `structures = SHAKHNOVICH_STRUCTURES` to sum only over the representative
subset of 10,000 structures used by Jacquin et al. (2016) (see
[`SHAKHNOVICH_STRUCTURES`](@ref)); note that structure indices refer to the set in use.

See also [`pnat`](@ref), [`energy`](@ref).
"""
function log_pnat(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = STRUCTURES)
    E = energies(seq; structures)
    return -_folding_loss(E, cm)
end

"""
    pnat(cm::Int, seq::AbstractVector{Int}; structures = STRUCTURES)

Probability that sequence `seq` folds into structure `cm`.
Equivalent to `exp(log_pnat(cm, seq; structures))`.

See also [`log_pnat`](@ref).
"""
pnat(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = STRUCTURES) = exp(log_pnat(cm, seq; structures))
