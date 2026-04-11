"""
    AMINO_ACIDS

Amino acid alphabet string used for encoding: `"CMFILVWYAGTSNQDEHRKP-"` (20 amino acids
plus a gap character). This ordering follows the convention used in the lattice proteins
literature.
"""
const AMINO_ACIDS = "CMFILVWYAGTSNQDEHRKP-"

"""
    onehot(seq::AbstractString)

One-hot encode an amino acid sequence string. Returns a `BitMatrix` of size `(21, L)`,
where `L` is the length of `seq`. Each column is a one-hot vector over the
[`AMINO_ACIDS`](@ref) alphabet.

    onehot(seqs::AbstractVector{<:AbstractString})

One-hot encode a vector of amino acid sequence strings. All sequences must have the same
length `L`. Returns a `BitArray{3}` of size `(21, L, N)`, where `N = length(seqs)`.

    onehot(A::Union{AbstractVector{Int}, AbstractMatrix{Int}})

One-hot encode integer-coded amino acid sequences (values in `1:21`).
"""
function onehot(seq::AbstractString)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(AMINO_ACIDS)
end

function onehot(seqs::AbstractVector{<:AbstractString})
    L = only(unique(length.(seqs)))
    return reshape(mapreduce(onehot, hcat, seqs), 21, L, length(seqs))
end

onehot(A::Union{AbstractVector{Int}, AbstractMatrix{Int}}) = reshape(A, 1, size(A)...) .== 1:21

"""
    aaseq(X::Union{BitMatrix, BitArray{3}})

Convert a one-hot encoded representation back to amino acid sequence string(s).
`X` must have 21 rows (one per amino acid in [`AMINO_ACIDS`](@ref)).

    aaseq(P::AbstractVector{Int})

Convert an integer-coded sequence to an amino acid string.

    aaseq(P::AbstractMatrix{Int})

Convert a matrix of integer-coded sequences (each column is one sequence) to a vector
of amino acid strings.
"""
function aaseq(X::Union{BitMatrix, BitArray{3}})
    @assert size(X, 1) == 21
    return aaseq(potts(X))
end

aaseq(P::AbstractVector{Int}) = join([AMINO_ACIDS[i] for i in P])
aaseq(P::AbstractMatrix{Int}) = [aaseq(view(P,:,n)) for n in axes(P, 2)]

"""
    potts(X::BitMatrix)
    potts(X::BitArray{3})

Convert a one-hot encoded representation to integer-coded (Potts) representation by
taking the `argmax` along the first dimension.

    potts(s::Union{AbstractString, AbstractVector{<:AbstractString}})

Convert amino acid string(s) to integer-coded (Potts) representation, using the
[`AMINO_ACIDS`](@ref) alphabet ordering.
"""
potts(X::BitMatrix) = vec(Int.(first.(Tuple.(argmax(X; dims=1)))))
potts(X::BitArray{3}) = reshape(Int.(first.(Tuple.(argmax(X; dims=1)))), size(X,2), size(X,3))
potts(s::Union{AbstractString, AbstractVector{<:AbstractString}}) = potts(onehot(s))
