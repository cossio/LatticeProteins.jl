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
function walk!(
        out::Vector{Vector{Int8}}, path::Vector{Int}, pos::Vector{Int}, used::UInt32,
        depth::Int, nbrs::Vector{Vector{Int}}, stabilizer::Vector{Vector{Int}}
    )
    if depth == L
        contacts = Tuple{Int, Int}[]
        for s in 1:L, u in nbrs[s]
            u > s || continue
            i, j = minmax(pos[s], pos[u])
            j - i > 1 && push!(contacts, (i, j)) # adjacent sites not bonded along the chain
        end
        @assert length(contacts) == N_CONTACTS
        sort!(contacts)
        cm = Vector{Int8}(undef, 2 * N_CONTACTS)
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
    return
end

# Contact maps of all maximally compact structures of a 27-mer, enumerated in a fixed
# order (see the docstring of `contact_map` for what is enumerated). Returns an array of
# size `(2, N_CONTACTS, nstructures)`, the layout of `StructureSet`.
function compact_contact_maps()
    nbrs = lattice_neighbors()
    syms = lattice_symmetries()
    out = Vector{Int8}[]
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
    return reshape(reduce(hcat, out), 2, N_CONTACTS, length(out))
end
