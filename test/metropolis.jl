import Makie
import CairoMakie
using LatticeProteins: random_msa
using LatticeProteins: pnat
using LatticeProteins: CONTACT_MAP_A
using LatticeProteins: Hugo_MSA
using LatticeProteins: potts
using DelimitedFiles: readdlm

msa = random_msa(CONTACT_MAP_A; nseqs=100, β=100, nsteps=1000)
#pnat(CONTACT_MAP_A, msa[:,1])
