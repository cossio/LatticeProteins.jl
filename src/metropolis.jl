# Contacts of a set of structures grouped by residue position: for each position `i`, the
# structures having a contact at `i` and the partner residue in that contact (stored in
# CSR layout: the entries of position `i` are `offsets[i]:(offsets[i+1]-1)`). This lets
# `metropolis!` update the energies of all structures after a single-site mutation by
# visiting only the contacts that involve the mutated position.
struct SiteContacts
    offsets::Vector{Int}
    structure::Vector{Int32}
    partner::Vector{Int32}
end

function SiteContacts(structures::AbstractArray{Int,3})
    check_structures(structures)
    counts = zeros(Int, L)
    for s in axes(structures, 2), k in axes(structures, 1)
        counts[structures[k, s, 1]] += 1
        counts[structures[k, s, 2]] += 1
    end
    offsets = Vector{Int}(undef, L + 1)
    offsets[1] = 1
    for i in 1:L
        offsets[i + 1] = offsets[i] + counts[i]
    end
    structure = Vector{Int32}(undef, offsets[end] - 1)
    partner = similar(structure)
    fill!(counts, 0)
    for s in axes(structures, 2), k in axes(structures, 1)
        i, j = structures[k, s, 1], structures[k, s, 2]
        t = offsets[i] + counts[i]
        structure[t] = s
        partner[t] = j
        counts[i] += 1
        t = offsets[j] + counts[j]
        structure[t] = s
        partner[t] = i
        counts[j] += 1
    end
    return SiteContacts(offsets, structure, partner)
end

"""
    metropolis!(seq::AbstractVector{Int}, cm::Int; β=1, nsteps=1, structures=CONTACT_MAPS)

Perform Metropolis Monte Carlo sampling on the amino acid sequence `seq`, targeting
structure `cm`. The sequence is mutated in place.

At each step, a random position is selected and mutated to a random amino acid (1–20).
The move is accepted or rejected using the Metropolis criterion with inverse temperature `β`,
where the energy is `-log_pnat(cm, seq; structures)`.

The energies of `seq` in all structures are kept and updated incrementally after each
mutation, so a step costs much less than a call to [`log_pnat`](@ref).

# Arguments
- `seq`: integer-coded amino acid sequence of length 27 (values in `1:20`).
- `cm`: target structure index (in `structures`).
- `β`: inverse temperature (higher values favor sequences that fold into `cm`).
- `nsteps`: number of Monte Carlo steps to perform.
- `structures`: set of competing structures, by default the 10,000 reference structures
  [`CONTACT_MAPS`](@ref). Pass the output of [`enumerate_structures`](@ref) to compete
  against all compact structures.

Returns the mutated `seq`.
"""
function metropolis!(seq::AbstractVector{Int}, cm::Int; β::Real=1, nsteps::Int=1, structures::AbstractArray{Int,3}=CONTACT_MAPS)
    check_sequence(seq)
    return metropolis!(seq, cm, structures, SiteContacts(structures); β, nsteps)
end

# Assumes `seq` has been validated with `check_sequence`, and `sites == SiteContacts(structures)`.
function metropolis!(seq::AbstractVector{Int}, cm::Int, structures::AbstractArray{Int,3}, sites::SiteContacts; β::Real, nsteps::Int)
    E = energies!(Vector{Float64}(undef, size(structures, 2)), seq, structures) # energies of seq in all structures
    E1 = similar(E) # energies of the proposed sequence
    F = E[cm] + logsumexp_neg(E) # -log_pnat
    for _ = 1:nsteps
        i = rand(1:length(seq)) # randomly select a position in the sequence
        old_aa = seq[i]
        new_aa = rand(1:20)
        copyto!(E1, E)
        # only contacts involving position i change energy
        @inbounds for t in sites.offsets[i]:(sites.offsets[i + 1] - 1)
            b = seq[sites.partner[t]]
            E1[sites.structure[t]] += MJ[new_aa, b] - MJ[old_aa, b]
        end
        F1 = E1[cm] + logsumexp_neg(E1)
        if β * (F1 - F) < randexp()
            seq[i] = new_aa # accept move
            F = F1
            E, E1 = E1, E
        end
    end
    return seq
end

"""
    random_msa(cm::Int; nseqs=1, β=1, nsteps=1000, structures=CONTACT_MAPS)

Generate a random multiple sequence alignment (MSA) of `nseqs` sequences that fold into
structure `cm`, using Metropolis Monte Carlo sampling.

Each sequence is initialized randomly and then evolved for `nsteps` Metropolis steps
at inverse temperature `β`. Returns a `Matrix{Int}` of size `(27, nseqs)`, where each
column is an integer-coded amino acid sequence. `structures` is the set of competing
structures (see [`metropolis!`](@ref)).

See also [`metropolis!`](@ref).
"""
function random_msa(cm::Int; nseqs::Int = 1, β::Real=1, nsteps::Int=1000, structures::AbstractArray{Int,3}=CONTACT_MAPS)
    sites = SiteContacts(structures)
    seqs = rand(1:20, L, nseqs)
    for n = 1:nseqs
        metropolis!(view(seqs, :, n), cm, structures, sites; β, nsteps)
    end
    return seqs
end
