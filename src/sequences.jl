"""
    AMINO_ACIDS

Amino acid alphabet string used for encoding: `"CMFILVWYAGTSNQDEHRKP-"` (20 amino acids
plus a gap character). This ordering follows the convention used in the lattice proteins
literature. Internal: [`potts`](@ref) and [`aaseq`](@ref) convert between amino acids and
their indices in this alphabet.
"""
const AMINO_ACIDS = "CMFILVWYAGTSNQDEHRKP-"

"""
    onehot(seq::AbstractString)

One-hot encode an amino acid sequence string. Returns a `BitMatrix` of size `(21, L)`,
where `L` is the length of `seq`. Each column is a one-hot vector over the alphabet
`"CMFILVWYAGTSNQDEHRKP-"` (20 amino acids and a gap): row `a` stands for its `a`-th
character.

    onehot(seqs::AbstractVector{<:AbstractString})

One-hot encode a vector of amino acid sequence strings. All sequences must have the same
length `L`. Returns a `BitArray{3}` of size `(21, L, N)`, where `N = length(seqs)`.

    onehot(A::Union{AbstractVector{<:Integer}, AbstractMatrix{<:Integer}})

One-hot encode integer-coded amino acid sequences (values in `1:21`).
"""
function onehot(seq::AbstractString)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(AMINO_ACIDS)
end

function onehot(seqs::AbstractVector{<:AbstractString})
    L = only(unique(length.(seqs)))
    X = falses(21, L, length(seqs))
    for (n, seq) in enumerate(seqs), (i, c) in enumerate(seq)
        a = findfirst(==(c), AMINO_ACIDS)
        a === nothing || (X[a, i, n] = true)
    end
    return X
end

onehot(A::Union{AbstractVector{<:Integer}, AbstractMatrix{<:Integer}}) = reshape(A, 1, size(A)...) .== 1:21

"""
    aaseq(X::Union{BitMatrix, BitArray{3}})

Convert a one-hot encoded representation back to amino acid sequence string(s).
`X` must have 21 rows, one per character of the alphabet `"CMFILVWYAGTSNQDEHRKP-"`.

    aaseq(P::AbstractVector{<:Integer})

Convert an integer-coded sequence to an amino acid string.

    aaseq(P::AbstractMatrix{<:Integer})

Convert a matrix of integer-coded sequences (each column is one sequence) to a vector
of amino acid strings.
"""
function aaseq(X::Union{BitMatrix, BitArray{3}})
    size(X, 1) == 21 || throw(ArgumentError("a one-hot encoding has 21 rows"))
    return aaseq(potts(X))
end

aaseq(P::AbstractVector{<:Integer}) = join(AMINO_ACIDS[i] for i in P)
aaseq(P::AbstractMatrix{<:Integer}) = [aaseq(view(P, :, n)) for n in axes(P, 2)]

"""
    potts(X::BitMatrix)
    potts(X::BitArray{3})

Convert a one-hot encoded representation to integer-coded (Potts) representation by
taking the `argmax` along the first dimension.

    potts(s::Union{AbstractString, AbstractVector{<:AbstractString}})

Convert amino acid string(s) to integer-coded (Potts) representation: each character is
replaced by its index in the alphabet `"CMFILVWYAGTSNQDEHRKP-"` (20 amino acids and a gap),
so `potts("CM-") == [1, 2, 21]`. A vector of strings (all of the same length) gives a
matrix with one sequence per column. Throws an `ArgumentError` for characters that are
not in the alphabet.
"""
potts(X::BitMatrix) = vec(Int.(first.(Tuple.(argmax(X; dims = 1)))))
potts(X::BitArray{3}) = reshape(Int.(first.(Tuple.(argmax(X; dims = 1)))), size(X, 2), size(X, 3))
potts(s::AbstractString) = [aa_index(c) for c in s]

function potts(seqs::AbstractVector{<:AbstractString})
    isempty(seqs) && throw(ArgumentError("cannot convert an empty collection of sequences"))
    L = length(first(seqs))
    all(s -> length(s) == L, seqs) || throw(ArgumentError("all sequences must have the same length"))
    P = Matrix{Int}(undef, L, length(seqs))
    for (n, seq) in enumerate(seqs), (i, c) in enumerate(seq)
        P[i, n] = aa_index(c)
    end
    return P
end

# Index of amino acid `c` in AMINO_ACIDS (which is ASCII, so byte and character indices agree).
function aa_index(c::AbstractChar)
    a = findfirst(==(c), AMINO_ACIDS)
    a === nothing && throw(ArgumentError("'$c' is not in the amino acid alphabet \"$AMINO_ACIDS\""))
    return a
end
