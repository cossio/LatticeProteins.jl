### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 984b1c6a-2d74-11f0-081d-4da16915f95e
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ fcad2697-1102-4ef7-bd4a-31e0f7314f57
using ProgressLogging: @progress, @withprogress, @logprogress

# ╔═╡ ab9b600a-4b75-4798-aa60-cc60d69e87b2
using Statistics: mean, cov

# ╔═╡ 5cbdc67d-7d58-4daa-8513-485576dfb2ee
using Random: randexp

# ╔═╡ 39d6c4a0-7cef-44a7-ba9d-0e8e1e8471e3
import LatticeProteins

# ╔═╡ 4e91c459-b6e1-4fff-a0bf-fcd478b827ff
import PlutoUI, Makie, CairoMakie

# ╔═╡ ded31a02-5c40-49fa-a33e-d6e52d7e162c
PlutoUI.TableOfContents()

# ╔═╡ 7bb78f2e-6676-4207-a4a8-391765c0ef60
md"# Functions"

# ╔═╡ ede43b37-13e6-4822-89c1-dfd4ca3d18da
tri(m::AbstractMatrix) = [m[i,j] for i = axes(m,1) for j = axes(m,2) if i ≤ j]

# ╔═╡ 4b53d605-4ca1-4b6c-95b6-f0ef3776a415
md"# Simulate"

# ╔═╡ 1c5175c8-f3b2-44af-8417-e5a6f68532df
hugo_msa_A = LatticeProteins.Hugo_MSA(:A)

# ╔═╡ a44bb919-6ea0-4410-812b-700721e3eef0
hugo_pnat_A = LatticeProteins.Hugo_MSA_pnat(:A)

# ╔═╡ 7f713ac9-8d32-4666-a640-2dd2b25ab196
mean(hugo_pnat_A)

# ╔═╡ bb5919b0-5cec-47be-ac05-89cbf5cfc435
function sample_msa(cm::Int; nseqs::Int = 1, β::Real=1, nsteps::Int=1000)
	L = 27
    seqs = rand(1:20, L, nseqs)
    @progress for n = 1:nseqs
        LatticeProteins.metropolis!(view(seqs,:,n), cm; β, nsteps)
    end
    return seqs
end

# ╔═╡ f4a32a91-b467-4368-aabf-540f1de89cd8
sampled_msa_100 = sample_msa(LatticeProteins.CONTACT_MAP_A; nseqs=1000, β=100, nsteps=1000)

# ╔═╡ fc1f387f-6bb4-4966-8d97-029b164d2c88
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=300, height=300, xlabel="<s_i>_Hugo", ylabel="<s_i>_Metropolis")
	Makie.ablines!(ax, 0, 1; color=:black)
	Makie.scatter!(ax, vec(mean(LatticeProteins.onehot(hugo_msa_A); dims=3)), vec(mean(LatticeProteins.onehot(sampled_msa_100); dims=3)); markersize=6, color=:blue, label="i after")
	Makie.axislegend(ax; position=:lt)

	ax = Makie.Axis(fig[1,2]; width=300, height=300, xlabel="<s_i s_j>_Hugo", ylabel="<s_i s_j>_Metropolis")
	Makie.ablines!(ax, 0, 1; color=:black)
	Makie.scatter!(ax, tri(cov(reshape(LatticeProteins.onehot(hugo_msa_A), 21 * 27, :); dims=2, corrected=false)), tri(cov(reshape(LatticeProteins.onehot(sampled_msa_100), 21 * 27, :); dims=2, corrected=false)); markersize=6, color=:blue, label="i after")
	Makie.axislegend(ax; position=:lt)

	ax = Makie.Axis(fig[1,3]; width=500, height=300, xlabel="Pnat", ylabel="frequency")
	Makie.hist!(ax, [LatticeProteins.pnat(LatticeProteins.CONTACT_MAP_A, seq) for seq = eachcol(sampled_msa_100)]; color=:blue, label="i after", normalization=:pdf)
	Makie.hist!(ax, hugo_pnat_A; color=:green, label="Hugo", normalization=:pdf)
	Makie.axislegend(ax; position=:ct)
	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ e002546a-da33-46c1-8ebf-5399b851d251
sampled_msa_1000 = sample_msa(LatticeProteins.CONTACT_MAP_A; nseqs=10000, β=1000, nsteps=1000)

# ╔═╡ c937b300-c274-4a8e-955b-a125724f591b
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=300, height=300, xlabel="<s_i>_Hugo", ylabel="<s_i>_Metropolis")
	Makie.ablines!(ax, 0, 1; color=:black)
	Makie.scatter!(ax, vec(mean(LatticeProteins.onehot(hugo_msa_A); dims=3)), vec(mean(LatticeProteins.onehot(sampled_msa_1000); dims=3)); markersize=6, color=:blue, label="i after")
	Makie.axislegend(ax; position=:lt)

	ax = Makie.Axis(fig[1,2]; width=300, height=300, xlabel="<s_i s_j>_Hugo", ylabel="<s_i s_j>_Metropolis")
	Makie.ablines!(ax, 0, 1; color=:black)
	Makie.scatter!(ax, tri(cov(reshape(LatticeProteins.onehot(hugo_msa_A), 21 * 27, :); dims=2, corrected=false)), tri(cov(reshape(LatticeProteins.onehot(sampled_msa_1000), 21 * 27, :); dims=2, corrected=false)); markersize=6, color=:blue, label="i after")
	Makie.axislegend(ax; position=:lt)

	ax = Makie.Axis(fig[1,3]; width=500, height=300, xlabel="Pnat", ylabel="frequency")
	Makie.hist!(ax, [LatticeProteins.pnat(LatticeProteins.CONTACT_MAP_A, seq) for seq = eachcol(sampled_msa_1000)]; color=:blue, label="i after", normalization=:pdf)
	Makie.hist!(ax, hugo_pnat_A; color=:green, label="Hugo", normalization=:pdf)
	Makie.axislegend(ax; position=:ct)
	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╠═984b1c6a-2d74-11f0-081d-4da16915f95e
# ╠═39d6c4a0-7cef-44a7-ba9d-0e8e1e8471e3
# ╠═4e91c459-b6e1-4fff-a0bf-fcd478b827ff
# ╠═fcad2697-1102-4ef7-bd4a-31e0f7314f57
# ╠═ab9b600a-4b75-4798-aa60-cc60d69e87b2
# ╠═5cbdc67d-7d58-4daa-8513-485576dfb2ee
# ╠═ded31a02-5c40-49fa-a33e-d6e52d7e162c
# ╠═7bb78f2e-6676-4207-a4a8-391765c0ef60
# ╠═ede43b37-13e6-4822-89c1-dfd4ca3d18da
# ╠═4b53d605-4ca1-4b6c-95b6-f0ef3776a415
# ╠═1c5175c8-f3b2-44af-8417-e5a6f68532df
# ╠═a44bb919-6ea0-4410-812b-700721e3eef0
# ╠═7f713ac9-8d32-4666-a640-2dd2b25ab196
# ╠═bb5919b0-5cec-47be-ac05-89cbf5cfc435
# ╠═f4a32a91-b467-4368-aabf-540f1de89cd8
# ╠═fc1f387f-6bb4-4966-8d97-029b164d2c88
# ╠═e002546a-da33-46c1-8ebf-5399b851d251
# ╠═c937b300-c274-4a8e-955b-a125724f591b
