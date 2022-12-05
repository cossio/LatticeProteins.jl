function energy(cmap, seq)
    E = 0.0
    for (i,j) in eachrow(cmap)
        E += MJ[seq[i], seq[j]]
    end
    return E
end

function log_pnat(cmap, seq)
    # note that this considers only a pre-defined set of 10000 contact maps,
    # instead of all possible structures
    competing_energies = [energy(CONTACT_MAPS[:,s,:], seq) for s in axes(CONTACT_MAPS, 2)]
    E = energy(cmap, seq)
    return logsumexp(E .- competing_energies)
end

pnat(cmap, seq) = exp(log_pnat(cmap, seq))
