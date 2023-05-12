module LatticeProteins

using LazyArtifacts: LazyArtifacts, @artifact_str
using DelimitedFiles: readdlm
using Random: randexp
using LogExpFunctions: logsumexp

const N_CONTACTS = 28 # number of pairwise contacts in 3x3x3 lattice

include("miyazawa_jernigan.jl")
include("contact_maps.jl")
include("pnat.jl")
include("onehot.jl")
include("metropolis.jl")

end # module
