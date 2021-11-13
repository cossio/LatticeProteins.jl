using Documenter, LatticeProteins

makedocs(
    modules=[LatticeProteins],
    sitename="LatticeProteins.jl"
)

deploydocs(
    repo = "github.com/cossio/LatticeProteins.jl.git",
    devbranch = "main"
)
