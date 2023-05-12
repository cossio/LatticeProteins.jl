function metropolis_step!(seq::AbstractVector{Int}; β::Real=1)
    i = rand(1:length(seq))   # randomly select a position in the sequence
    j = rand(1:20)  # randomly select a new amino acid
    ΔE = MJ[seq[i], j] - MJ[seq[i], seq[i]] # compute the energy difference
    # accept or reject the new amino acid
    if rand() < exp(-β*ΔE)
        seq[i] = j
    end
end
