"""
    energy(cm::Int, seq::AbstractVector{Int})

Compute the contact energy of amino acid sequence `seq` in structure `cm`.

The energy is the sum of Miyazawa-Jernigan interaction energies over all contacting
residue pairs in the given structure. `seq` must be a vector of length 27 with
integer amino acid indices in `1:20`. `cm` is the structure index in `1:10000`
(see [`N_STRUCTURES`](@ref)).
"""
function energy(cm::Int, seq::AbstractVector{Int})
    @assert length(seq) == 27
    E = 0.0
    for (i,j) in eachrow(CONTACT_MAPS[:,cm,:])
        E += MJ[seq[i], seq[j]]
    end
    return E
end

"""
    log_pnat(cm::Int, seq::AbstractVector{Int})

Compute the log-probability that sequence `seq` folds into structure `cm` (the
"native" structure), as opposed to any of the other pre-computed structures.

The probability is computed using a Boltzmann distribution at temperature 1 over
the 10,000 pre-computed compact structures:

```math
\\log P_\\text{nat}(\\text{cm} \\mid \\text{seq}) = -\\log \\sum_s \\exp\\bigl(E(\\text{cm}, \\text{seq}) - E(s, \\text{seq})\\bigr)
```

See also [`pnat`](@ref), [`energy`](@ref).
"""
function log_pnat(cm::Int, seq::AbstractVector{Int})
    @assert length(seq) == 27
    E = energy(cm, seq)
    # note that this considers only a pre-defined set of 10000 contact maps, instead of all possible structures
    competing_energies = [energy(s, seq) for s in axes(CONTACT_MAPS, 2)]
    return -logsumexp(E .- competing_energies)
end

"""
    pnat(cm::Int, seq::AbstractVector{Int})

Probability that sequence `seq` folds into structure `cm`.
Equivalent to `exp(log_pnat(cm, seq))`.

See also [`log_pnat`](@ref).
"""
pnat(cm::Int, seq::AbstractVector{Int}) = exp(log_pnat(cm, seq))
