using LatticeProteins: random_msa, pnat, CONTACT_MAP_A

msa = random_msa(CONTACT_MAP_A; nseqs=1, β=100, nsteps=10000)
pnat(CONTACT_MAP_A, msa[:,1])
