using Test: @test
using LatticeProteins: pnat
using LatticeProteins: CONTACT_MAP_A
using LatticeProteins: CONTACT_MAP_B
using LatticeProteins: Hugo_MSA
using LatticeProteins: Hugo_MSA_pnat
using LatticeProteins: potts
using LatticeProteins: onehot
using Statistics: mean
using Statistics: cov

msa_A = Hugo_MSA(:A)
pnat_A = Hugo_MSA_pnat(:A)

msa_B = Hugo_MSA(:B)
pnat_B = Hugo_MSA_pnat(:B)

for i = 1:100
    @test pnat(CONTACT_MAP_A, potts(msa_A[i])) ≈ pnat_A[i]
    @test pnat(CONTACT_MAP_B, potts(msa_B[i])) ≈ pnat_B[i]
end

mean(onehot(msa_A); dims=3)
