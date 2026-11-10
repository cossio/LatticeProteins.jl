# State of a Monte Carlo chain targeting structure `s` of a structure set: the sequence,
# its energies in all the structures of the set (in units of 0.01), a buffer for the
# energies of a proposed sequence, and `loss = -log_pnat(s, seq)`. Both samplers of the
# package advance chains with `local_steps!`.
mutable struct Chain
    seq::Vector{Int}
    energies::Vector{Int32}
    trial::Vector{Int32}
    loss::Float64
end

# Assumes `seq` has been validated with `check_sequence`, and `s` with `check_structure`.
# The sequence is copied.
function Chain(seq::AbstractVector{<:Integer}, s::Integer, S::StructureSet)
    E = energies!(Vector{Int32}(undef, length(S)), seq, S)
    return Chain(collect(Int, seq), E, similar(E), folding_loss(E, s, S.boltzmann))
end

# Given the energies `E` of `seq` in all structures, writes into `E1` the energies of the
# sequence obtained by mutating position `i` of `seq` to amino acid `aa`, visiting only
# the contacts that involve position `i`. `sites` is the per-site contact index of the
# structure set of `E`.
function mutate_energies!(E1::Vector{Int32}, E::Vector{Int32}, seq::Vector{Int}, i::Int, aa::Int, sites::SiteContacts)
    copyto!(E1, E)
    old = seq[i]
    # energy change of a contact between position i and an amino acid b
    ΔE = ntuple(b -> MJ_INT[aa, b] - MJ_INT[old, b], Val(20))
    @inbounds for t in sites.offsets[i]:(sites.offsets[i + 1] - 1)
        E1[sites.structure[t]] += ΔE[seq[sites.partner[t]]]
    end
    return E1
end

# Performs `nsteps` Metropolis steps on `chain` at inverse temperature `β`, with energy
# `-log_pnat` of structure `s`, and returns the number of accepted moves. A step proposes
# to mutate a uniformly chosen position to one of the other 19 amino acids, chosen
# uniformly. `sites` and `boltzmann` are those of the structure set of the chain.
function local_steps!(
        rng::AbstractRNG, chain::Chain, s::Integer, β::Real, nsteps::Integer,
        sites::SiteContacts, boltzmann::Vector{Float64}
    )
    accepted = 0
    for _ in 1:nsteps
        i = rand(rng, 1:L)
        aa = rand(rng, 1:19)
        aa += aa >= chain.seq[i]
        mutate_energies!(chain.trial, chain.energies, chain.seq, i, aa, sites)
        loss = folding_loss(chain.trial, s, boltzmann)
        if β * (loss - chain.loss) < randexp(rng)
            chain.seq[i] = aa
            chain.energies, chain.trial = chain.trial, chain.energies
            chain.loss = loss
            accepted += 1
        end
    end
    return accepted
end

"""
    metropolis!(seq, s; β = 1, nsteps = 1, rng = default_rng())

Perform Metropolis Monte Carlo sampling on the amino acid sequence `seq`, targeting
structure `s`. The sequence is mutated in place, and returned.

Each step proposes to mutate a random position to a different amino acid, chosen
uniformly among the other 19. The move is accepted or rejected using the Metropolis
criterion with inverse temperature `β`, where the energy is
`-log_pnat(s, seq)`. The stationary distribution is thus proportional to
`pnat(s, seq)^β`.

The energies of `seq` in all structures are kept and updated incrementally after each
mutation, so a step costs much less than a call to [`log_pnat`](@ref). They are computed
from scratch at the start of each call: prefer one call with a large `nsteps` to many
short ones.

# Arguments
- `seq`: integer-coded amino acid sequence of length 27 (values in `1:20`).
- `s`: index of the target structure, in `1:N_STRUCTURES`.
- `β`: inverse temperature (higher values favor sequences that fold into `s`).
- `nsteps`: number of Monte Carlo steps to perform.
- `rng`: random number generator.

See also [`random_msa`](@ref), [`parallel_tempering`](@ref).
"""
function metropolis!(
        seq::AbstractVector{<:Integer}, s::Integer; β::Real = 1, nsteps::Integer = 1,
        structures::StructureSet = STRUCTURES, rng::AbstractRNG = default_rng()
    )
    check_sequence(seq)
    check_structure(structures, s)
    nsteps ≥ 0 || throw(ArgumentError("nsteps must be nonnegative"))
    chain = Chain(seq, s, structures)
    local_steps!(rng, chain, s, β, nsteps, site_contacts(structures), structures.boltzmann)
    copyto!(seq, chain.seq)
    return seq
end

"""
    random_msa(s; nseqs = 1, β = 1, nsteps = 1000, rng = default_rng())

Generate a random multiple sequence alignment (MSA) of `nseqs` sequences that fold into
structure `s`, using Metropolis Monte Carlo sampling.

Each sequence is initialized randomly and then evolved independently for `nsteps`
Metropolis steps at inverse temperature `β` (see [`metropolis!`](@ref)). Returns a
`Matrix{Int}` of size `(27, nseqs)`, where each column is an integer-coded amino acid
sequence.

See also [`metropolis!`](@ref), [`parallel_tempering`](@ref).
"""
function random_msa(
        s::Integer; nseqs::Integer = 1, β::Real = 1, nsteps::Integer = 1000,
        structures::StructureSet = STRUCTURES, rng::AbstractRNG = default_rng()
    )
    nseqs ≥ 0 || throw(ArgumentError("nseqs must be nonnegative"))
    seqs = rand(rng, 1:20, L, nseqs)
    for n in 1:nseqs
        metropolis!(view(seqs, :, n), s; β, nsteps, structures, rng)
    end
    return seqs
end
