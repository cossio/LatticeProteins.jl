using Test: @test
using LatticeProteins: pnat, CONTACT_MAP_A, CONTACT_MAP_B, Hugo_MSA, Hugo_MSA_pnat, potts

msa_A = Hugo_MSA(:A)
pnat_A = Hugo_MSA_pnat(:A)

msa_B = Hugo_MSA(:B)
pnat_B = Hugo_MSA_pnat(:B)

for i in 1:100
    @test pnat(CONTACT_MAP_A, potts(msa_A[i])) ≈ pnat_A[i]
    @test pnat(CONTACT_MAP_B, potts(msa_B[i])) ≈ pnat_B[i]
end
