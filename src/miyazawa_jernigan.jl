const MJ_DIR = joinpath(dirname(pathof(LatticeProteins)), "..", "MJ")

"""
    load_miyazawa_jernigan_matrix([T=Float64])

Load the 20×20 Miyazawa-Jernigan contact energy matrix from the bundled data file.
Element `(i, j)` gives the interaction energy between amino acid types `i` and `j`,
where the amino acid indexing follows the order in [`AMINO_ACIDS`](@ref) (first 20 characters,
excluding the gap character).

Returns a `Matrix{T}`.
"""
function load_miyazawa_jernigan_matrix(::Type{T} = Float64) where {T}
    readdlm(joinpath(MJ_DIR, "MJ.txt"), T)
end

"""
    load_miyazawa_jernigan_aminoacids()

Load the amino acid ordering used in the Miyazawa-Jernigan matrix file.
Returns a `String` of single-letter amino acid codes.
"""
function load_miyazawa_jernigan_aminoacids()
    readline(joinpath(MJ_DIR, "AAs.txt"))
end

"""
    MJ

The 20×20 Miyazawa-Jernigan contact energy matrix (preloaded as `Matrix{Float64}`).
See [`load_miyazawa_jernigan_matrix`](@ref).
"""
const MJ = load_miyazawa_jernigan_matrix()
