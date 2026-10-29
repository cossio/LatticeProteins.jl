# Coordinates of lattice site `s` (sites are numbered 1:27) and the inverse map.
lattice_coords(s::Int) = ((s - 1) % 3, (s - 1) ÷ 3 % 3, (s - 1) ÷ 9)
lattice_site(x::Int, y::Int, z::Int) = 1 + x + 3y + 9z

# Lattice neighbours of each site.
function lattice_neighbors()
    nbrs = [Int[] for _ in 1:L]
    for s in 1:L, d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
        x, y, z = lattice_coords(s) .+ d
        if 0 ≤ x ≤ 2 && 0 ≤ y ≤ 2 && 0 ≤ z ≤ 2
            push!(nbrs[s], lattice_site(x, y, z))
        end
    end
    return nbrs
end

# The 48 rotations and reflections of the cube, as permutations of the 27 sites.
function lattice_symmetries()
    syms = Vector{Int}[]
    for p in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)), f in 0:7
        g = Vector{Int}(undef, L)
        for s in 1:L
            c = lattice_coords(s)
            o = ntuple(a -> ((f >> (a - 1)) & 1) == 1 ? 2 - c[p[a]] : c[p[a]], 3)
            g[s] = lattice_site(o...)
        end
        push!(syms, g)
    end
    return syms
end

# Whether the first `depth` sites of `path` are lexicographically no larger than their image
# under each of the `symmetries`. The lexicographic comparison of a walk with its image is
# decided by the first site where they differ, so a prefix that is larger than its image
# cannot be extended to a canonical walk.
function is_canonical(path::Vector{Int}, depth::Int, symmetries::Vector{Vector{Int}})
    for g in symmetries
        for i in 1:depth
            a, b = path[i], g[path[i]]
            b < a && return false
            b > a && break
        end
    end
    return true
end

# Depth-first extension of `path` (which occupies `depth` sites so far, and is canonical)
# over all self-avoiding continuations, pruning the ones that are not canonical. Complete
# walks are converted to a contact map (flattened as `[i1, j1, i2, j2, ...]`) and appended
# to `out`.
function walk!(out::Vector{Vector{Int}}, path::Vector{Int}, pos::Vector{Int}, used::UInt32,
               depth::Int, nbrs::Vector{Vector{Int}}, stabilizer::Vector{Vector{Int}})
    if depth == L
        contacts = Tuple{Int,Int}[]
        for s in 1:L, u in nbrs[s]
            u > s || continue
            i, j = minmax(pos[s], pos[u])
            j - i > 1 && push!(contacts, (i, j)) # adjacent sites not bonded along the chain
        end
        @assert length(contacts) == N_CONTACTS
        sort!(contacts)
        cm = Vector{Int}(undef, 2 * N_CONTACTS)
        for (k, (i, j)) in enumerate(contacts)
            cm[2k - 1] = i
            cm[2k] = j
        end
        push!(out, cm)
        return
    end
    for u in nbrs[path[depth]]
        if (used & (UInt32(1) << u)) == 0
            path[depth + 1] = u
            pos[u] = depth + 1
            if is_canonical(path, depth + 1, stabilizer)
                walk!(out, path, pos, used | (UInt32(1) << u), depth + 1, nbrs, stabilizer)
            end
        end
    end
end

"""
    enumerate_structures()

Enumerate all maximally compact structures of a 27-mer on the 3×3×3 cubic lattice, that is,
all self-avoiding walks visiting the 27 lattice sites, modulo the 48 rotations and
reflections of the cube. There are 103,346 of them. The chain direction is kept: reversing
the chain maps contact `(i, j)` to `(28 - j, 28 - i)`, which in general gives a different
contact map, so a structure and its chain reversal count as distinct structures (62 of
them are their own reversal).

Returns an `Array{Int, 3}` of size `(28, 103346, 2)`: `structures[:, s, :]` gives the 28
contact pairs of structure `s`, as rows `(i, j)` with `i < j`, sorted lexicographically.
The enumeration order is deterministic. Takes about a second.

The result is available pre-computed as [`STRUCTURES`](@ref), which is the default set of
competing structures in [`pnat`](@ref), [`metropolis!`](@ref) and related functions; this
function is only needed to regenerate it.

# References

- Shakhnovich, E. I., and A. M. Gutin. "Enumeration of all compact conformations of
  copolymers with random sequence of links." *The Journal of Chemical Physics* 93.8
  (1990): 5967-5971.
"""
function enumerate_structures()
    nbrs = lattice_neighbors()
    syms = lattice_symmetries()
    out = Vector{Int}[]
    # A walk through the 27 sites must start on the majority sublattice (corners or face
    # centres). Every walk is the image under a cube symmetry of one starting from the
    # corner (0,0,0) or from the face centre (1,1,0), so it suffices to start from these two
    # sites and keep the walks that are minimal among their images under the symmetries
    # fixing the starting site (the stabilizer).
    for s0 in (lattice_site(0, 0, 0), lattice_site(1, 1, 0))
        stabilizer = [g for g in syms if g[s0] == s0 && g != 1:L]
        path = zeros(Int, L)
        pos = zeros(Int, L)
        path[1] = s0
        pos[s0] = 1
        walk!(out, path, pos, UInt32(1) << s0, 1, nbrs, stabilizer)
    end
    cms = reshape(reduce(hcat, out), 2, N_CONTACTS, length(out))
    return permutedims(cms, (2, 3, 1))
end

# Validates an array of structures: `structures[:, s, :]` must be the contact pairs of structure `s`.
function check_structures(structures::AbstractArray{Int,3})
    @assert size(structures, 3) == 2 "structures must be an array of size (ncontacts, nstructures, 2)"
end

# Validates an integer-coded sequence.
function check_sequence(seq::AbstractVector{Int})
    @assert length(seq) == L "sequence must have length $L"
    @assert all(a -> 1 ≤ a ≤ 20, seq) "amino acid indices must be in 1:20"
end

# Canonical representation of a contact map, independent of the order of the contacts and
# of the order within each pair.
contact_key(cm::AbstractMatrix{Int}) = sort!([minmax(cm[k, 1], cm[k, 2]) for k in axes(cm, 1)])

"""
    structure_index(structures, cm::AbstractMatrix{Int})
    structure_index(structures, cm::Int)

Index within `structures` of the structure with contact map `cm`, or `nothing` if it is not
there. `structures` is an array of structures such as [`STRUCTURES`](@ref) or
[`JACQUIN_STRUCTURES`](@ref). `cm` is a `28 × 2` matrix of contact pairs, such as
`STRUCTURES[:, CONTACT_MAP_A, :]`; the order of the contacts and the order within each
pair are ignored. If `cm` is an integer, it refers to structure `STRUCTURES[:, cm, :]` of
the default set.

# Example

```julia
a = structure_index(JACQUIN_STRUCTURES, CONTACT_MAP_A)  # 625: index of structure A in the reference set
pnat(a, seq; structures = JACQUIN_STRUCTURES)
```
"""
function structure_index(structures::AbstractArray{Int,3}, cm::AbstractMatrix{Int})
    check_structures(structures)
    key = contact_key(cm)
    length(key) == size(structures, 1) || return nothing
    candidate = similar(key)
    for s in axes(structures, 2)
        for (t, k) in enumerate(axes(structures, 1))
            candidate[t] = minmax(structures[k, s, 1], structures[k, s, 2])
        end
        sort!(candidate)
        candidate == key && return s
    end
    return nothing
end

structure_index(structures::AbstractArray{Int,3}, cm::Int) = structure_index(structures, STRUCTURES[:, cm, :])

"""
    STRUCTURES

Contact maps of all 103,346 compact structures of a 27-mer on the 3×3×3 cubic lattice, as
generated by [`enumerate_structures`](@ref). An `Array{Int, 3}` of size `(28, 103346, 2)`;
use `STRUCTURES[:, s, :]` to get the contact pairs of structure `s`.

This is the default set of competing structures in [`pnat`](@ref), [`metropolis!`](@ref)
and related functions. See [`JACQUIN_STRUCTURES`](@ref) for the subset of 10,000 reference
structures used by Jacquin et al. (2016).
"""
const STRUCTURES = enumerate_structures()

"""
    N_STRUCTURES

Number of compact structures in [`STRUCTURES`](@ref) (103,346).
"""
const N_STRUCTURES = size(STRUCTURES, 2)

# Index in STRUCTURES of structure `s` of JACQUIN_STRUCTURES.
jacquin_structure(s::Int) = something(structure_index(STRUCTURES, JACQUIN_STRUCTURES[:, s, :]))

"""
    CONTACT_MAP_A

Index in [`STRUCTURES`](@ref) of reference structure A of Jacquin et al. (2016), which is
structure 625 of [`JACQUIN_STRUCTURES`](@ref).
"""
const CONTACT_MAP_A = jacquin_structure(625)

"""
    CONTACT_MAP_B

Index in [`STRUCTURES`](@ref) of reference structure B of Jacquin et al. (2016), which is
structure 848 of [`JACQUIN_STRUCTURES`](@ref).
"""
const CONTACT_MAP_B = jacquin_structure(848)

"""
    CONTACT_MAP_C

Index in [`STRUCTURES`](@ref) of reference structure C of Jacquin et al. (2016), which is
structure 6683 of [`JACQUIN_STRUCTURES`](@ref).
"""
const CONTACT_MAP_C = jacquin_structure(6683)

"""
    CONTACT_MAP_D

Index in [`STRUCTURES`](@ref) of reference structure D of Jacquin et al. (2016), which is
structure 6685 of [`JACQUIN_STRUCTURES`](@ref).
"""
const CONTACT_MAP_D = jacquin_structure(6685)
