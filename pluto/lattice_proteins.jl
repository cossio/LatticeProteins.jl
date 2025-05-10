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

# ╔═╡ 39d6c4a0-7cef-44a7-ba9d-0e8e1e8471e3
import LatticeProteins

# ╔═╡ 4e91c459-b6e1-4fff-a0bf-fcd478b827ff
import PlutoUI, Makie, CairoMakie

# ╔═╡ ded31a02-5c40-49fa-a33e-d6e52d7e162c
PlutoUI.TableOfContents()

# ╔═╡ 4b53d605-4ca1-4b6c-95b6-f0ef3776a415
md"# Simulate"

# ╔═╡ 1c5175c8-f3b2-44af-8417-e5a6f68532df
hugo_msa_A = LatticeProteins.Hugo_MSA(:A)

# ╔═╡ a44bb919-6ea0-4410-812b-700721e3eef0
hugo_pnat_A = LatticeProteins.Hugo_MSA_pnat(:A)

# ╔═╡ 7f713ac9-8d32-4666-a640-2dd2b25ab196
mean(hugo_pnat_A)

# ╔═╡ bb5919b0-5cec-47be-ac05-89cbf5cfc435


# ╔═╡ Cell order:
# ╠═984b1c6a-2d74-11f0-081d-4da16915f95e
# ╠═39d6c4a0-7cef-44a7-ba9d-0e8e1e8471e3
# ╠═4e91c459-b6e1-4fff-a0bf-fcd478b827ff
# ╠═fcad2697-1102-4ef7-bd4a-31e0f7314f57
# ╠═ab9b600a-4b75-4798-aa60-cc60d69e87b2
# ╠═ded31a02-5c40-49fa-a33e-d6e52d7e162c
# ╠═4b53d605-4ca1-4b6c-95b6-f0ef3776a415
# ╠═1c5175c8-f3b2-44af-8417-e5a6f68532df
# ╠═a44bb919-6ea0-4410-812b-700721e3eef0
# ╠═7f713ac9-8d32-4666-a640-2dd2b25ab196
# ╠═bb5919b0-5cec-47be-ac05-89cbf5cfc435
