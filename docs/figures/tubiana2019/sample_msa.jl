# Sample an alignment of sequences folding into structure A with the protocol of
# Jacquin et al. (2016): Metropolis at β = 1000, one sequence kept every 1,000 steps after a
# burn-in of 10,000 steps, from several independent chains started at random sequences.
#
# Usage: julia -t <threads> sample_msa.jl <shakhnovich|all> <chains> <samples per chain> <output file> [β] [burn-in]
#
# β defaults to 1000 and the burn-in to 10,000 steps.
#
# `shakhnovich` competes against the subset of 10,000 structures (the model of Jacquin et al.
# 2016 and of Tubiana et al. 2019), `all` against all 103,346 compact structures. The output has one
# line per sequence: "1000.0 <pnat> <sequence>", the format of the Hugo_MSA files.
using LatticeProteins: LatticeProteins, STRUCTURES, SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A,
    structure_index, aaseq, log_pnat, L
using Random: seed!

which, nchains, per_chain, outfile = ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), ARGS[4]
which in ("all", "shakhnovich") || error("first argument must be 'shakhnovich' or 'all'")
β = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 1000.0
burnin = length(ARGS) ≥ 6 ? parse(Int, ARGS[6]) : 10_000
interval = 1000
S = which == "all" ? STRUCTURES : SHAKHNOVICH_STRUCTURES
cm = which == "all" ? CONTACT_MAP_A : structure_index(SHAKHNOVICH_STRUCTURES, CONTACT_MAP_A)
sites = LatticeProteins.site_contacts(S)

results = Vector{Vector{Tuple{Float64,String}}}(undef, nchains)
seed!(2026)
t0 = time()
Threads.@threads for c in 1:nchains
    seq = rand(1:20, L)
    LatticeProteins.metropolis!(seq, cm, S, sites; β, nsteps=burnin)
    out = Tuple{Float64,String}[]
    for t in 1:per_chain
        LatticeProteins.metropolis!(seq, cm, S, sites; β, nsteps=interval)
        push!(out, (exp(log_pnat(cm, seq; structures=S)), aaseq(seq)))
    end
    results[c] = out
end
open(outfile, "w") do io
    for out in results, (p, s) in out
        println(io, "1000.0 ", p, " ", s)
    end
end
println("wrote $(sum(length, results)) sequences from $nchains chains (β = $β, burn-in $burnin) in $(round(time() - t0)) s")
