"""
    Hugo_MSA_path(which::Symbol)

Return the file path to one of the reference MSA datasets from Jacquin et al. (2016).
Internal, like the other `Hugo_MSA` functions: these alignments were sampled against the
subset of 10,000 structures ([`shakhnovich_structures`](@ref)), not the full model.
`which` must be one of `:A`, `:B`, `:C`, or `:D`, corresponding to the four reference
structures ([`STRUCTURE_A`](@ref), etc.). Each line of the file holds the inverse
temperature used for sampling, the folding probability and the sequence.

The files are downloaded the first time one of them is requested.
"""
function Hugo_MSA_path(which::Symbol)
    which in (:A, :B, :C, :D) || throw(ArgumentError("Unknown MSA $which; must be one of :A, :B, :C, :D"))
    return joinpath(artifact"MSAs_Hugo_PlosCb", "align_$which.msa")
end

"""
    Hugo_MSA(which::Symbol)

Load a reference MSA dataset from Jacquin et al. (2016) as a `Vector{String}` of amino acid
sequences. `which` must be one of `:A`, `:B`, `:C`, or `:D`.

See also [`Hugo_MSA_pnat`](@ref), [`Hugo_MSA_path`](@ref).
"""
function Hugo_MSA(which::Symbol)
    path = Hugo_MSA_path(which)
    return string.(readdlm(path)[:, 3])
end

"""
    Hugo_MSA_pnat(which::Symbol)

Load the native fold probabilities (pnat values) for the reference MSA dataset
from Jacquin et al. (2016). `which` must be one of `:A`, `:B`, `:C`, or `:D`.
Returns a `Vector{Float64}`. These values were computed against the subset of 10,000
structures, i.e. `pnat(s, seq; structures = shakhnovich_structures())` with `s` the index
of the structure in [`shakhnovich_structures`](@ref).

See also [`Hugo_MSA`](@ref).
"""
function Hugo_MSA_pnat(which::Symbol)
    path = Hugo_MSA_path(which)
    return Float64.(readdlm(path)[:, 2])
end
