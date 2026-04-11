"""
    metropolis!(seq::AbstractVector{Int}, cm::Int; β=1, nsteps=1)

Perform Metropolis Monte Carlo sampling on the amino acid sequence `seq`, targeting
structure `cm`. The sequence is mutated in place.

At each step, a random position is selected and mutated to a random amino acid (1–20).
The move is accepted or rejected using the Metropolis criterion with inverse temperature `β`,
where the energy is `-log_pnat(cm, seq)`.

# Arguments
- `seq`: integer-coded amino acid sequence of length 27 (values in `1:20`).
- `cm`: target structure index.
- `β`: inverse temperature (higher values favor sequences that fold into `cm`).
- `nsteps`: number of Monte Carlo steps to perform.

Returns the mutated `seq`.
"""
function metropolis!(seq::AbstractVector{Int}, cm::Int; β::Real=1, nsteps::Int=1)
    @assert length(seq) == 27 && all(a -> 1 ≤ a ≤ 20, seq)
    E = -log_pnat(cm, seq)
    for _ = 1:nsteps
        i = rand(1:length(seq)) # randomly select a position in the sequence
        old_aa = seq[i]
        seq[i] = rand(1:20)
        E1 = -log_pnat(cm, seq)
        if β * (E1 - E) < randexp()
            E = E1 # accept move
        else
            seq[i] = old_aa # revert change
        end
    end
    return seq
end

"""
    random_msa(cm::Int; nseqs=1, β=1, nsteps=1000)

Generate a random multiple sequence alignment (MSA) of `nseqs` sequences that fold into
structure `cm`, using Metropolis Monte Carlo sampling.

Each sequence is initialized randomly and then evolved for `nsteps` Metropolis steps
at inverse temperature `β`. Returns a `Matrix{Int}` of size `(27, nseqs)`, where each
column is an integer-coded amino acid sequence.

See also [`metropolis!`](@ref).
"""
function random_msa(cm::Int; nseqs::Int = 1, β::Real=1, nsteps::Int=1000)
    seqs = rand(1:20, L, nseqs)
    for n = 1:nseqs
        metropolis!(view(seqs,:,n), cm; β, nsteps)
    end
    return seqs
end
