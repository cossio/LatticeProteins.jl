"""
    energy(cm::Int, seq::AbstractVector{Int}; structures = CONTACT_MAPS)

Compute the contact energy of amino acid sequence `seq` in structure `cm`.

The energy is the sum of Miyazawa-Jernigan interaction energies over all contacting
residue pairs in the given structure. `seq` must be a vector of length 27 with
integer amino acid indices in `1:20`. `cm` is the index of the structure in `structures`,
which defaults to the 10,000 reference structures [`CONTACT_MAPS`](@ref) and can be set
to the output of [`enumerate_structures`](@ref).
"""
function energy(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = CONTACT_MAPS)
    check_sequence(seq)
    check_structures(structures)
    checkbounds(structures, :, cm, :)
    E = 0.0
    for k in axes(structures, 1)
        @inbounds i, j = structures[k, cm, 1], structures[k, cm, 2]
        a, b = seq[i], seq[j]
        @inbounds E += MJ[a, b] # a, b ∈ 1:20 by check_sequence
    end
    return E
end

"""
    energies(seq::AbstractVector{Int}; structures = CONTACT_MAPS)

Contact energies of `seq` in every structure of `structures`, as a `Vector{Float64}` of
length `size(structures, 2)`. Equivalent to `[energy(s, seq; structures) for s in
1:size(structures, 2)]`, but faster.

See also [`energy`](@ref).
"""
function energies(seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = CONTACT_MAPS)
    check_sequence(seq)
    check_structures(structures)
    return energies!(Vector{Float64}(undef, size(structures, 2)), seq, structures)
end

# Writes the energies of `seq` in all `structures` into `E`. Assumes `seq` and `structures`
# have been validated with `check_sequence` and `check_structures`.
function energies!(E::Vector{Float64}, seq::AbstractVector{Int}, structures::AbstractArray{Int,3})
    @assert length(E) == size(structures, 2)
    for s in axes(structures, 2)
        e = 0.0
        for k in axes(structures, 1)
            @inbounds i, j = structures[k, s, 1], structures[k, s, 2]
            a, b = seq[i], seq[j]
            @inbounds e += MJ[a, b] # a, b ∈ 1:20 by check_sequence
        end
        E[s] = e
    end
    return E
end

# log(sum(exp.(-E))), computed stably and without allocating.
function logsumexp_neg(E::Vector{Float64})
    m = -minimum(E)
    s = 0.0
    for e in E
        s += exp(-e - m)
    end
    return m + log(s)
end

"""
    log_pnat(cm::Int, seq::AbstractVector{Int}; structures = CONTACT_MAPS)

Compute the log-probability that sequence `seq` folds into structure `cm` (the
"native" structure), as opposed to any of the other structures in `structures`.

The probability is computed using a Boltzmann distribution at temperature 1 over
the structures in `structures`:

```math
\\log P_\\text{nat}(\\text{cm} \\mid \\text{seq}) = -\\log \\sum_s \\exp\\bigl(E(\\text{cm}, \\text{seq}) - E(s, \\text{seq})\\bigr)
```

By default `structures` is the set of 10,000 reference structures [`CONTACT_MAPS`](@ref)
used by Jacquin et al. (2016), which is a subset of the 103,346 compact structures of the
lattice. To sum over all of them, pass `structures = enumerate_structures()` (see
[`enumerate_structures`](@ref)).

See also [`pnat`](@ref), [`energy`](@ref).
"""
function log_pnat(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = CONTACT_MAPS)
    E = energies(seq; structures)
    return -E[cm] - logsumexp_neg(E)
end

"""
    pnat(cm::Int, seq::AbstractVector{Int}; structures = CONTACT_MAPS)

Probability that sequence `seq` folds into structure `cm`.
Equivalent to `exp(log_pnat(cm, seq; structures))`.

See also [`log_pnat`](@ref).
"""
pnat(cm::Int, seq::AbstractVector{Int}; structures::AbstractArray{Int,3} = CONTACT_MAPS) = exp(log_pnat(cm, seq; structures))
