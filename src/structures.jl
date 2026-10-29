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

# Whether `path` is lexicographically no larger than each of its images under `symmetries`.
function is_canonical(path::Vector{Int}, symmetries::Vector{Vector{Int}})
    for g in symmetries
        for i in 1:L
            a, b = path[i], g[path[i]]
            b < a && return false
            b > a && break
        end
    end
    return true
end

# Depth-first extension of `path` (which occupies `depth` sites so far) over all
# self-avoiding continuations. Complete walks that are canonical are converted to a
# contact map (flattened as `[i1, j1, i2, j2, ...]`) and appended to `out`.
function walk!(out::Vector{Vector{Int}}, path::Vector{Int}, pos::Vector{Int}, used::UInt32,
               depth::Int, nbrs::Vector{Vector{Int}}, stabilizer::Vector{Vector{Int}})
    if depth == L
        is_canonical(path, stabilizer) || return
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
            walk!(out, path, pos, used | (UInt32(1) << u), depth + 1, nbrs, stabilizer)
        end
    end
end

"""
    enumerate_structures()

Enumerate all maximally compact structures of a 27-mer on the 3×3×3 cubic lattice, that is,
all self-avoiding walks visiting the 27 lattice sites, modulo the 48 rotations and
reflections of the cube. There are 103,346 of them. The chain direction is kept: a structure
and its chain reversal are distinct structures, since they place residue 1 at different ends
of the walk and have different contact maps.

Returns an `Array{Int, 3}` of size `(28, 103346, 2)` with the same layout as
[`CONTACT_MAPS`](@ref): `structures[:, s, :]` gives the 28 contact pairs of structure `s`,
as rows `(i, j)` with `i < j`, sorted lexicographically. The enumeration order is
deterministic. Takes a few seconds.

The 10,000 reference structures of [`CONTACT_MAPS`](@ref) are a subset of these; use
[`structure_index`](@ref) to locate one of them. Pass the result as the `structures`
keyword argument of [`energy`](@ref), [`pnat`](@ref), [`log_pnat`](@ref),
[`metropolis!`](@ref) or [`random_msa`](@ref) to compute folding probabilities over all
compact structures instead of the 10,000 reference ones:

```julia
ALL = enumerate_structures()
a = structure_index(ALL, CONTACT_MAP_A)
pnat(a, seq; structures = ALL)
```

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
    # fixing the starting site.
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

# Canonical representation of a contact map, independent of the order of the contacts and
# of the order within each pair.
contact_key(cm::AbstractMatrix{Int}) = sort!([minmax(cm[k, 1], cm[k, 2]) for k in axes(cm, 1)])

"""
    structure_index(structures, cm::AbstractMatrix{Int})
    structure_index(structures, cm::Int)

Index within `structures` of the structure with contact map `cm`, or `nothing` if it is not
there. `structures` is an array of structures such as [`CONTACT_MAPS`](@ref) or the output
of [`enumerate_structures`](@ref). `cm` is a `28 × 2` matrix of contact pairs, such as
`CONTACT_MAPS[:, CONTACT_MAP_A, :]`; the order of the contacts and the order within each
pair are ignored. If `cm` is an integer, it refers to the reference structure
`CONTACT_MAPS[:, cm, :]`.

# Example

```julia
ALL = enumerate_structures()
a = structure_index(ALL, CONTACT_MAP_A)  # index of reference structure A in ALL
pnat(a, seq; structures = ALL)
```
"""
function structure_index(structures::AbstractArray{Int,3}, cm::AbstractMatrix{Int})
    check_structures(structures)
    key = contact_key(cm)
    for s in axes(structures, 2)
        contact_key(view(structures, :, s, :)) == key && return s
    end
    return nothing
end

structure_index(structures::AbstractArray{Int,3}, cm::Int) = structure_index(structures, CONTACT_MAPS[:, cm, :])

# Validates an array of structures: `structures[:, s, :]` must be the contact pairs of structure `s`.
function check_structures(structures::AbstractArray{Int,3})
    @assert size(structures, 3) == 2 "structures must be an array of size (ncontacts, nstructures, 2)"
end

# Validates an integer-coded sequence.
function check_sequence(seq::AbstractVector{Int})
    @assert length(seq) == L "sequence must have length $L"
    @assert all(a -> 1 ≤ a ≤ 20, seq) "amino acid indices must be in 1:20"
end
