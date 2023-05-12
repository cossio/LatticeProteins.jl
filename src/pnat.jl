function energy(cm::Int, seq)
    E = 0.0
    for (i,j) in eachrow(CONTACT_MAPS[:,cm,:])
        E += MJ[seq[i], seq[j]]
    end
    return E
end

function log_pnat(cm::Int, seq)
    E = energy(cm, seq)
    # note that this considers only a pre-defined set of 10000 contact maps, instead of all possible structures
    competing_energies = [energy(s, seq) for s in axes(CONTACT_MAPS, 2)]
    return -logsumexp(E .- competing_energies)
end

pnat(cm, seq) = exp(log_pnat(cm, seq))
