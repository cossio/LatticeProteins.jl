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

function Hugo_MSA(which::Symbol)
    path = Hugo_MSA_path(which)
    return string.(readdlm(path)[:,3])
end

function Hugo_MSA_pnat(which::Symbol)
    path = Hugo_MSA_path(which)
    return Float64.(readdlm(path)[:,2])
end
