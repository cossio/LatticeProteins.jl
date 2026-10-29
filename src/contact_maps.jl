"""
    L

Sequence length for lattice proteins on the 3×3×3 cubic lattice (27).
"""
const L = 27 # sequence length

"""
    load_contact_maps()

Load the pre-computed contact maps of the 10,000 reference compact structures used by
Jacquin et al. (2016). Returns an `Array{Int, 3}` of size `(28, 10000, 2)`, where
`cm[:, s, :]` gives the 28 pairs of contacting residue indices of structure `s`.
See [`JACQUIN_STRUCTURES`](@ref).
"""
function load_contact_maps()
    cm = readdlm(joinpath(artifact"ContactMaps", "contact_maps_10000.dat"), Int)[:, 2:end]
    return reshape(cm, 28, 10_000, 2) # array[:, s, :] gives contacts of structure 's'
end

"""
    JACQUIN_STRUCTURES

Contact maps of the 10,000 reference compact structures used by Jacquin et al. (2016), a
subset of [`STRUCTURES`](@ref). An `Array{Int, 3}` of size `(28, 10000, 2)`; use
`JACQUIN_STRUCTURES[:, s, :]` to get the contact pairs of structure `s`.
See [`load_contact_maps`](@ref).

Pass `structures = JACQUIN_STRUCTURES` to [`pnat`](@ref), [`metropolis!`](@ref) and related
functions to compete only against this set, as in Jacquin et al. (2016) and in the
[`Hugo_MSA_pnat`](@ref) values. Structure indices refer to the set in use; use
[`structure_index`](@ref) to translate between the two sets. For example,
`structure_index(JACQUIN_STRUCTURES, CONTACT_MAP_A)` gives the index (625) of reference
structure A in this set.
"""
const JACQUIN_STRUCTURES = load_contact_maps()
