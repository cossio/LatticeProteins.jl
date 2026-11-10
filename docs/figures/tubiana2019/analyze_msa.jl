# Statistics of an MSA for structure A, as in Fig. 4 of Tubiana, Cocco & Monasson (2019).
# Usage: julia --project=<package dir> analyze_msa.jl <msa file: "beta pnat seq" per line> <output dir> [pnat threshold]
# The competing-structure set is irrelevant here: the contact map of structure A is the same in both sets.
using LatticeProteins: potts, onehot, contact_map, structure_index, shakhnovich_structures,
    AMINO_ACIDS, STRUCTURE_A, L, N_CONTACTS
using LinearAlgebra: eigen, Symmetric
using Statistics: mean
using DelimitedFiles: writedlm
using Random: Xoshiro, randperm
using TOML

msafile, outdir = ARGS[1], ARGS[2]
thr = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.995
mkpath(outdir)
rows = [split(l) for l in eachline(msafile) if !isempty(strip(l))]
pn = [parse(Float64, r[2]) for r in rows]; seqs = [String(r[3]) for r in rows]
keep = pn .> thr
println("$(length(seqs)) sequences, $(count(keep)) with pnat > $thr")
seqs = seqs[keep]; pn = pn[keep]
M = length(seqs)
M >= 2 || error("need at least two sequences above the requested threshold $thr")
X = Float64.(reshape(onehot(seqs)[1:20, :, :], 20 * L, M))   # 540 × M one-hot
f = reshape(mean(X; dims = 2), 20, L)                           # site frequencies f_i(a)
cons = [log(20) + sum(p > 0 ? p * log(p) : 0.0 for p in f[:, i]) for i in 1:L]
# PCA: top eigenvector of the covariance of the one-hot encoding
Xc = X .- mean(X; dims = 2)
Cov = Symmetric(Xc * Xc' / M)
ev = eigen(Cov)
w = ev.vectors[:, end]
W = reshape(w, 20, L)
charged = [findfirst(==(c), AMINO_ACIDS) for c in "DEKRH"]
sum(W[charged, :]) < 0 && (w .*= -1)                         # W is a view of w
proj = Xc' * w
println("top eigenvalues: ", round.(ev.values[end:-1:(end - 3)]; sigdigits = 3))
# average covariance between pairs of sites in contact in structure A (20×20)
S10 = shakhnovich_structures()
a10 = structure_index(S10, STRUCTURE_A)
cmA = contact_map(a10; structures = S10)
P = potts(seqs)                                                # 27 × M
C = zeros(20, 20)
for k in 1:N_CONTACTS
    i, j = cmA[k, 1], cmA[k, 2]
    fij = zeros(20, 20)
    for m in 1:M
        fij[P[i, m], P[j, m]] += 1
    end
    fij ./= M
    cov = fij .- f[:, i] * f[:, j]'
    C .+= cov .+ cov'
end
C ./= N_CONTACTS
writedlm(joinpath(outdir, "freq.csv"), f, ',')
writedlm(joinpath(outdir, "conservation.csv"), cons, ',')
writedlm(joinpath(outdir, "pc1.csv"), W, ',')
writedlm(joinpath(outdir, "proj.csv"), proj, ',')
writedlm(joinpath(outdir, "contact_cov.csv"), C, ',')
writedlm(joinpath(outdir, "pnat.csv"), pn, ',')
open(joinpath(outdir, "sample_seqs.txt"), "w") do io
    # Draw panel (c) across the whole alignment, not the beginning of one chain.
    for i in randperm(Xoshiro(2019), M)[1:min(60, M)]
        println(io, seqs[i])
    end
end
open(joinpath(outdir, "analysis.toml"), "w") do io
    TOML.print(
        io, Dict(
            "pnat_threshold" => thr, "input_sequences" => length(rows),
            "retained_sequences" => M, "mean_pnat" => mean(pn),
            "eigenvalues" => ev.values[end:-1:(end - 4)], "panel_c_seed" => 2019
        )
    )
end
# 3D walks of structures A and G (6863 in the reference set) for drawing. The 27 lattice
# sites are numbered so that site `1 + x + 3y + 9z` has coordinates `(x, y, z)`.
site_coords(s) = ((s - 1) % 3, (s - 1) ÷ 3 % 3, (s - 1) ÷ 9)
function site_neighbors(s)
    nbrs = Int[]
    for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
        x, y, z = site_coords(s) .+ d
        0 ≤ x ≤ 2 && 0 ≤ y ≤ 2 && 0 ≤ z ≤ 2 && push!(nbrs, 1 + x + 3y + 9z)
    end
    return nbrs
end
# A self-avoiding walk through the 27 sites whose contacts are those of `cm`.
function find_walk(cm)
    key = Set(zip(cm[:, 1], cm[:, 2]))
    nbrs = site_neighbors.(1:L)
    path = zeros(Int, L); pos = zeros(Int, L); found = Ref{Vector{Int}}()
    function dfs(depth, used)
        isassigned(found) && return
        if depth == L
            contacts = Set{Tuple{Int, Int}}()
            for s in 1:L, u in nbrs[s]
                u > s || continue
                i, j = minmax(pos[s], pos[u]); j - i > 1 && push!(contacts, (i, j))
            end
            contacts == key && (found[] = copy(path)); return
        end
        for u in nbrs[path[depth]]
            if (used & (UInt32(1) << u)) == 0
                path[depth + 1] = u; pos[u] = depth + 1
                dfs(depth + 1, used | (UInt32(1) << u))
            end
        end
        return
    end
    for s0 in (1, 5) # a corner, (0,0,0), or a face centre, (1,1,0)
        path[1] = s0; pos[s0] = 1; dfs(1, UInt32(1) << s0)
        isassigned(found) && break
    end
    return [site_coords(s)[c] for s in found[], c in 1:3]
end
for (name, s) in (("A", a10), ("G", 6863))
    cm = contact_map(s; structures = S10)
    writedlm(joinpath(outdir, "walk_$name.csv"), find_walk(cm), ',')
    writedlm(joinpath(outdir, "contacts_$name.csv"), cm, ',')
end
println("conservation (first 10 sites): ", round.(cons[1:10]; digits = 2))
println("done")
