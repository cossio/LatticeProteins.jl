"""
    L

Sequence length for lattice proteins on the 3×3×3 cubic lattice (27).
"""
const L = 27 # sequence length

"""
    load_contact_maps()

Load the pre-computed contact maps of the 10,000 reference compact structures on the 3×3×3
cubic lattice used by Jacquin et al. (2016).
Returns an `Array{Int, 3}` of size `(28, 10000, 2)`, where `CONTACT_MAPS[:, s, :]` gives the
28 pairs of contacting residue indices for structure `s`.
"""
function load_contact_maps()
    cm = readdlm(joinpath(artifact"ContactMaps", "contact_maps_10000.dat"), Int)[:, 2:end]
    return reshape(cm, 28, 10_000, 2) # array[:, s, :] gives contacts of structure 's'
end

"""
    CONTACT_MAPS

Pre-loaded contact maps of the 10,000 reference compact structures. An `Array{Int, 3}` of
size `(28, 10000, 2)`. Use `CONTACT_MAPS[:, s, :]` to get the contact pairs for structure `s`.
See [`load_contact_maps`](@ref).

These are the structures used by Jacquin et al. (2016), and the default set of competing
structures in [`pnat`](@ref) and [`metropolis!`](@ref). They are a subset of the 103,346
compact structures of the lattice, which can be generated with
[`enumerate_structures`](@ref).
"""
const CONTACT_MAPS = load_contact_maps()

"""
    N_STRUCTURES

Number of reference compact structures in [`CONTACT_MAPS`](@ref) (10,000).
"""
const N_STRUCTURES = size(CONTACT_MAPS, 2)

"""
    CONTACT_MAP_A

Index of reference structure A (625) from Jacquin et al. (2016).
"""
const CONTACT_MAP_A = 625

"""
    CONTACT_MAP_B

Index of reference structure B (848) from Jacquin et al. (2016).
"""
const CONTACT_MAP_B = 848

"""
    CONTACT_MAP_C

Index of reference structure C (6683) from Jacquin et al. (2016).
"""
const CONTACT_MAP_C = 6683

"""
    CONTACT_MAP_D

Index of reference structure D (6685) from Jacquin et al. (2016).
"""
const CONTACT_MAP_D = 6685
