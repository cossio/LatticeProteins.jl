function metropolis!(seq::AbstractVector{Int}, cm::Int; β::Real=1, nsteps::Int=1)
    @assert length(seq) == 27
    i = rand(1:length(seq))   # randomly select a position in the sequence
    E = -log_pnat(cm, seq)
    for _ = 1:nsteps
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

function random_msa(cm::Int; nseqs::Int = 1, β::Int=1, nsteps::Int=1000)
    seqs = rand(1:20, L, nseqs)
    for n in 1:nseqs
        metropolis!(view(seqs,:,n), cm; β, nsteps)
    end
    return seqs
end
