"""
    Hugo_MSA_path(which::Symbol)

Return the file path to one of the reference MSA datasets from Jacquin et al. (2016).
`which` must be one of `:A`, `:B`, `:C`, or `:D`, corresponding to the four reference
structures ([`CONTACT_MAP_A`](@ref), etc.).
"""
function Hugo_MSA_path(which::Symbol)
    if which == :A
        return joinpath(artifact"MSAs_Hugo_PlosCb", "align_A.msa")
    elseif which == :B
        return joinpath(artifact"MSAs_Hugo_PlosCb", "align_B.msa")
    elseif which == :C
        return joinpath(artifact"MSAs_Hugo_PlosCb", "align_C.msa")
    elseif which == :D
        return joinpath(artifact"MSAs_Hugo_PlosCb", "align_D.msa")
    else
        throw(ArgumentError("Unknown MSA $which; must be one of :A, :B, :C, :D"))
    end
end

"""
    Hugo_MSA(which::Symbol)

Load a reference MSA dataset from Jacquin et al. (2016) as a `Vector{String}` of amino acid
sequences. `which` must be one of `:A`, `:B`, `:C`, or `:D`.

See also [`Hugo_MSA_pnat`](@ref), [`Hugo_MSA_path`](@ref).
"""
function Hugo_MSA(which::Symbol)
    path = Hugo_MSA_path(which)
    return string.(readdlm(path)[:,3])
end

"""
    Hugo_MSA_pnat(which::Symbol)

Load the native fold probabilities (pnat values) for the reference MSA dataset
from Jacquin et al. (2016). `which` must be one of `:A`, `:B`, `:C`, or `:D`.
Returns a `Vector{Float64}`.

See also [`Hugo_MSA`](@ref).
"""
function Hugo_MSA_pnat(which::Symbol)
    path = Hugo_MSA_path(which)
    return Float64.(readdlm(path)[:,2])
end
