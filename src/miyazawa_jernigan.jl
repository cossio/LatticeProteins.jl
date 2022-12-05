const MJ_DIR = joinpath(dirname(pathof(LatticeProteins)), "..", "MJ")

function load_miyazawa_jernigan_matrix(::Type{T} = Float64) where {T}
    readdlm(joinpath(MJ_DIR, "MJ.txt"), T)
end

function load_miyazawa_jernigan_aminoacids()
    readline(joinpath(MJ_DIR, "AAs.txt"))
end

const MJ = load_miyazawa_jernigan_matrix()
