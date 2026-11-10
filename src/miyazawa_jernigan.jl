const MJ_PATH = joinpath(@__DIR__, "..", "MJ", "MJ.txt")
Base.include_dependency(MJ_PATH)

"""
    MJ

The 20×20 Miyazawa-Jernigan contact energy matrix (a `Matrix{Float64}`). Element `(a, b)`
is the interaction energy between amino acids `a` and `b`, indexed as in integer-coded
sequences (see [`potts`](@ref)): the amino acids are in the order `"CMFILVWYAGTSNQDEHRKP"`.
Internal: [`energy`](@ref) sums these interactions over the contacts of a structure.
"""
const MJ = readdlm(MJ_PATH, Float64)

# The package computes with energies in units of 0.01, the resolution of the MJ table.
# Sums of contact energies are then exact integers: energies updated incrementally during
# sampling never drift, and Boltzmann factors can be read from a table (see `StructureSet`)
# instead of evaluating an exponential for every structure.
const MJ_INT = round.(Int32, 100 .* MJ)
@assert MJ_INT ./ 100 == MJ "the MJ energies must be multiples of 0.01"
