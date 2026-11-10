# Contacts of a set of structures grouped by residue position: for each position `i`, the
# structures having a contact at `i` and the partner residue in that contact (stored in
# CSR layout: the entries of position `i` are `offsets[i]:(offsets[i+1]-1)`). This lets the
# samplers update the energies of all structures after a single-site mutation by visiting
# only the contacts that involve the mutated position.
struct SiteContacts
    offsets::Vector{Int}
    structure::Vector{Int32}
    partner::Vector{Int8}
end

function SiteContacts(contacts::Array{Int8, 3})
    counts = zeros(Int, L)
    for i in contacts
        counts[i] += 1
    end
    offsets = Vector{Int}(undef, L + 1)
    offsets[1] = 1
    for i in 1:L
        offsets[i + 1] = offsets[i] + counts[i]
    end
    structure = Vector{Int32}(undef, offsets[end] - 1)
    partner = Vector{Int8}(undef, offsets[end] - 1)
    fill!(counts, 0)
    for s in axes(contacts, 3), k in axes(contacts, 2)
        i, j = contacts[1, k, s], contacts[2, k, s]
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
    StructureSet(contact_maps)

A set of structures competing to fold a sequence. `contact_maps` is a vector with one
contact map per structure: a matrix of size `(ncontacts, 2)` whose rows are the pairs of
residue positions (in `1:27`) in contact. All structures must have the same number of
contacts.

The functions that take a structure index ([`contact_map`](@ref), [`energy`](@ref),
[`pnat`](@ref), [`metropolis!`](@ref), [`sample_msa`](@ref), ...) accept the set of
competing structures as the keyword argument `structures`, and the index then refers to
that set. The package has two sets:

- [`enumerate_structures`](@ref)`()`: all 103,346 compact structures of the lattice (the
  default);
- [`shakhnovich_structures`](@ref)`()`: the subset of 10,000 structures used by Jacquin et
  al. (2016).

`length(S)` is the number of structures of a set `S`, and [`structure_index`](@ref) finds
a structure in a set, which translates indices between sets.

This type and the `structures` keyword are internal: the public API is the model with all
compact structures.

# Example

```julia
# two structures with a single contact each
S = StructureSet([[1 4], [1 6]])
pnat(1, seq; structures = S)
```
"""
mutable struct StructureSet
    # contacts[:, k, s] are the two residues of contact k of structure s
    const contacts::Array{Int8, 3}
    # exp(-ΔE) for every energy difference ΔE = 0, 0.01, 0.02, ... that can occur between
    # two structures of the set
    const boltzmann::Vector{Float64}
    # built on first use by `site_contacts`
    @atomic sites::Union{Nothing, SiteContacts}

    function StructureSet(contacts::Array{Int8, 3})
        size(contacts, 1) == 2 && size(contacts, 2) > 0 && size(contacts, 3) > 0 ||
            throw(ArgumentError("a structure set needs at least one structure, with at least one contact"))
        all(i -> 1 ≤ i ≤ L, contacts) ||
            throw(ArgumentError("contacts must be between residue positions in 1:$L"))
        all(contacts[1, k, s] != contacts[2, k, s] for s in axes(contacts, 3), k in axes(contacts, 2)) ||
            throw(ArgumentError("a residue cannot be in contact with itself"))
        span = size(contacts, 2) * (Int(maximum(MJ_INT)) - Int(minimum(MJ_INT)))
        return new(contacts, exp.(-(0:span) ./ 100), nothing)
    end
end

function StructureSet(contact_maps::AbstractVector{<:AbstractMatrix{<:Integer}})
    isempty(contact_maps) && throw(ArgumentError("a structure set needs at least one structure"))
    foreach(Base.require_one_based_indexing, contact_maps)
    n = size(first(contact_maps), 1)
    all(cm -> size(cm) == (n, 2), contact_maps) ||
        throw(ArgumentError("contact maps must be matrices of size (ncontacts, 2), with the same number of contacts"))
    all(cm -> all(i -> 1 ≤ i ≤ L, cm), contact_maps) ||
        throw(ArgumentError("contacts must be between residue positions in 1:$L"))
    contacts = Array{Int8, 3}(undef, 2, n, length(contact_maps))
    for (s, cm) in enumerate(contact_maps), k in 1:n
        contacts[1, k, s] = cm[k, 1]
        contacts[2, k, s] = cm[k, 2]
    end
    return StructureSet(contacts)
end

Base.length(S::StructureSet) = size(S.contacts, 3)

ncontacts(S::StructureSet) = size(S.contacts, 2)

Base.show(io::IO, S::StructureSet) =
    print(io, "StructureSet with ", length(S), " structures of ", ncontacts(S), " contacts")

# Per-site contact index of `S`, built on first use. Threads racing to build it compute
# the same index, so whichever is stored is correct.
function site_contacts(S::StructureSet)
    sites = @atomic S.sites
    sites === nothing || return sites
    sites = SiteContacts(S.contacts)
    @atomic S.sites = sites
    return sites
end

# Validates the index of a structure of `S`.
function check_structure(S::StructureSet, s::Integer)
    return 1 ≤ s ≤ length(S) || throw(ArgumentError("structure index $s is not in 1:$(length(S))"))
end

"""
    contact_map(s)

Contact map of structure `s`, one of the 103,346 compact structures of
[`enumerate_structures`](@ref): a `Matrix{Int}` of size `(28, 2)` whose rows are the pairs
`(i, j)` of residue positions in contact, with `i < j`, sorted lexicographically.

```julia
contact_map(STRUCTURE_A)  # the 28 contacts of structure A
```
"""
function contact_map(s::Integer; structures::StructureSet = STRUCTURES)
    check_structure(structures, s)
    return Int[structures.contacts[c, k, s] for k in 1:ncontacts(structures), c in 1:2]
end

# Canonical representation of a contact map, independent of the order of the contacts and
# of the order within each pair.
contact_key(cm::AbstractMatrix{<:Integer}) = sort!([minmax(Int(cm[k, 1]), Int(cm[k, 2])) for k in axes(cm, 1)])

"""
    structure_index(S::StructureSet, cm::AbstractMatrix)
    structure_index(S::StructureSet, s::Integer; from = STRUCTURES)

Internal. Index within the set `S` of the structure with contact map `cm`, or `nothing` if
it is not there. `cm` is a matrix of contact pairs as returned by [`contact_map`](@ref); the order
of the contacts and the order within each pair are ignored.

The second method looks for structure `s` of the set `from` (by default
[`enumerate_structures`](@ref)`()`), and so translates structure indices from one set to
another.

# Example

```julia
S10 = shakhnovich_structures()
a = structure_index(S10, STRUCTURE_A)  # 625: index of structure A in the subset
pnat(a, seq; structures = S10)
```
"""
function structure_index(S::StructureSet, cm::AbstractMatrix{<:Integer})
    size(cm, 2) == 2 || throw(ArgumentError("a contact map is a matrix of size (ncontacts, 2)"))
    key = contact_key(cm)
    length(key) == ncontacts(S) || return nothing
    candidate = similar(key)
    for s in 1:length(S)
        for k in 1:ncontacts(S)
            candidate[k] = minmax(Int(S.contacts[1, k, s]), Int(S.contacts[2, k, s]))
        end
        sort!(candidate)
        candidate == key && return s
    end
    return nothing
end

structure_index(S::StructureSet, s::Integer; from::StructureSet = STRUCTURES) = structure_index(S, contact_map(s; structures = from))

# Enumerated when the package is precompiled, which takes about a second.
const STRUCTURES = StructureSet(compact_contact_maps())

"""
    enumerate_structures()

The set of all compact structures of a 27-mer on the 3×3×3 cubic lattice: all the
self-avoiding walks visiting the 27 lattice sites, modulo the 48 rotations and reflections
of the cube. There are `length(enumerate_structures())` = 103,346 of them, and a structure
is referred to by its index in this set; [`contact_map`](@ref) gives its contacts.

The chain direction is kept: reversing the chain maps contact `(i, j)` to
`(28 - j, 28 - i)`, which in general gives a different contact map, so a structure and its
chain reversal count as distinct structures (62 of them are their own reversal).

The structures are enumerated once, in a fixed order, when the package is precompiled:
calling this function costs nothing, and always returns the same object.

# References

- Shakhnovich, E. I., and A. M. Gutin. "Enumeration of all compact conformations of
  copolymers with random sequence of links." *The Journal of Chemical Physics* 93.8
  (1990): 5967-5971.
"""
enumerate_structures() = STRUCTURES

# Number of compact structures (103,346).
const N_STRUCTURES = length(STRUCTURES)

# The subset is read from an artifact that is downloaded on first use, so it is not loaded
# with the package.
const shakhnovich_set = Base.OncePerProcess{StructureSet}() do
    # one contact per row, "structure i j", the 28 contacts of a structure on consecutive rows
    data = readdlm(joinpath(artifact"ContactMaps", "contact_maps_10000.dat"), Int)
    contacts = reshape(Int8.(data[:, 2:3]), N_CONTACTS, :, 2) # [contact, structure, residue]
    return StructureSet(permutedims(contacts, (3, 1, 2)))
end

"""
    shakhnovich_structures()

The [`StructureSet`](@ref) of the representative subset of 10,000 compact structures
selected by the Shakhnovich group (Heo, Maslov & Shakhnovich, 2011) out of the 103,346 of
[`enumerate_structures`](@ref), and used as the set of competing structures by Jacquin et al. (2016)
and Tubiana et al. (2019).

Pass `structures = shakhnovich_structures()` to [`pnat`](@ref), [`metropolis!`](@ref) and
related functions to compete only against this set, as in those papers and in the
[`Hugo_MSA_pnat`](@ref) values. Structure indices refer to the set in use; use
[`structure_index`](@ref) to translate between the two sets. For example,
`structure_index(shakhnovich_structures(), STRUCTURE_A)` gives the index (625) of reference
structure A in this set.

The contact maps are downloaded the first time this function is called, and the set is
built once per session. Internal: the public API is the model with all compact structures.

# References

- Heo, M., Maslov, S., and Shakhnovich, E. "Topology of protein interaction network shapes
  protein abundances and strengths of their functional and nonspecific interactions."
  *Proceedings of the National Academy of Sciences* 108.10 (2011): 4258-4263.
- Jacquin, H., et al. "Benchmarking inverse statistical approaches for protein structure
  and design with exactly solvable models." *PLoS Computational Biology* 12.5 (2016):
  e1004889.
"""
shakhnovich_structures() = shakhnovich_set()

# The reference structures of Jacquin et al. (2016) are given there as indices in the
# subset of 10,000 structures. Their indices in STRUCTURES are written out here, so that
# loading the package does not need the subset; they depend on the order of the
# enumeration, and the tests check them against the subset.

"""
    STRUCTURE_A

Index of reference structure A of Jacquin et al. (2016) among the compact structures
([`enumerate_structures`](@ref)). It is structure 625 of the subset of 10,000 structures used in that paper.
"""
const STRUCTURE_A = 22451

"""
    STRUCTURE_B

Index of reference structure B of Jacquin et al. (2016) among the compact structures
([`enumerate_structures`](@ref)). It is structure 848 of the subset of 10,000 structures used in that paper.
"""
const STRUCTURE_B = 65036

"""
    STRUCTURE_C

Index of reference structure C of Jacquin et al. (2016) among the compact structures
([`enumerate_structures`](@ref)). It is structure 6683 of the subset of 10,000 structures used in that paper.
"""
const STRUCTURE_C = 66913

"""
    STRUCTURE_D

Index of reference structure D of Jacquin et al. (2016) among the compact structures
([`enumerate_structures`](@ref)). It is structure 6685 of the subset of 10,000 structures used in that paper.
"""
const STRUCTURE_D = 90697
