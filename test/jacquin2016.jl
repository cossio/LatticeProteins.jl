using Test: @test, @testset
using LatticeProteins: Hugo_MSA, Hugo_MSA_pnat, pnat, potts, energy, log_pnat, onehot,
    CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D,
    MJ, CONTACT_MAPS, L, N_STRUCTURES, N_CONTACTS, AMINO_ACIDS,
    random_msa, metropolis!
using Statistics: mean
using Random: seed!

# Tests based on properties from:
#   Jacquin, Hugo, et al. "Benchmarking inverse statistical approaches for protein
#   structure and design with exactly solvable models."
#   PLoS Computational Biology 12.5 (2016): e1004889.

@testset "jacquin2016" begin
    @testset "pnat normalization" begin
        # pnat defines a probability distribution over structures for any sequence.
        # The sum of pnat over all 10,000 structures must equal 1.
        for seq_str in ["EKAMPAMDPDMAHEHKKKIRAWMFEGE",  # Hugo MSA A
                        "PDRFIVQCLAQFEHFERDGDKDMRAEC",  # Hugo MSA B
                        "WRFQDPGRCEEEDERMCGPEMCRCRQK",  # Hugo MSA D
                        "CMFILVWYAGTSNQDEHRKPCMFILVW"]  # arbitrary sequence
            seq = potts(seq_str)
            pnat_sum = sum(pnat(s, seq) for s in 1:N_STRUCTURES)
            @test pnat_sum ≈ 1.0 atol=1e-10
        end
    end

    @testset "Hugo MSA pnat distribution" begin
        # The Hugo MSAs were sampled at β = 1000 (see first column of .msa files).
        # At this inverse temperature, sequences should have very high pnat (≥ 0.97).
        # The mean pnat should be approximately 0.995.
        for (sym, cm) in [(:A, CONTACT_MAP_A), (:B, CONTACT_MAP_B),
                          (:C, CONTACT_MAP_C), (:D, CONTACT_MAP_D)]
            pnats = Hugo_MSA_pnat(sym)

            # All sequences in the Hugo MSA should have high pnat
            @test all(pnats .> 0.97)

            # Mean pnat should be above 0.99 (paper: typical pnat ≈ 0.995 at β = 1000)
            @test mean(pnats) > 0.99
        end
    end

    @testset "Hugo MSA pnat consistency" begin
        # Verify that pnat values stored in the Hugo MSA files match our computation.
        for (sym, cm) in [(:A, CONTACT_MAP_A), (:B, CONTACT_MAP_B),
                          (:C, CONTACT_MAP_C), (:D, CONTACT_MAP_D)]
            msa = Hugo_MSA(sym)
            pnats = Hugo_MSA_pnat(sym)

            # Check first 10 sequences for each structure
            for i in 1:10
                @test pnat(cm, potts(msa[i])) ≈ pnats[i]
            end
        end
    end

    @testset "Hugo MSA sequence properties" begin
        # Sequences in the Hugo MSA should be valid: length 27, using only the 20 standard AAs.
        for sym in [:A, :B, :C, :D]
            msa = Hugo_MSA(sym)
            @test all(length(s) == L for s in msa)
            @test all(all(c in AMINO_ACIDS[1:20] for c in s) for s in msa)
        end
    end

    @testset "Hugo MSA amino acid frequencies" begin
        # At β = 1000, sequences are not random: amino acid frequencies should be
        # non-uniform, reflecting the structure-specific selection.
        for sym in [:A, :B, :C, :D]
            msa = Hugo_MSA(sym)
            X = onehot(msa)
            # Site-averaged amino acid frequencies (21 × 27 → marginal over sequences)
            freq = dropdims(mean(X; dims=3); dims=3)

            # No gap characters should appear
            @test all(freq[21, :] .== 0)

            # Frequencies at each position should be non-uniform (max > 1/20 = 0.05)
            for pos in 1:L
                @test maximum(freq[1:20, pos]) > 0.05
            end

            # Some positions should be highly conserved (max frequency > 0.2)
            max_freqs = [maximum(freq[1:20, pos]) for pos in 1:L]
            @test maximum(max_freqs) > 0.2
        end
    end

    @testset "Metropolis sampling convergence" begin
        # At high β, Metropolis sampling should produce sequences with high pnat,
        # similar to the Hugo MSAs.
        seed!(42)
        seq = rand(1:20, L)
        initial_pnat = pnat(CONTACT_MAP_A, seq)

        # After sufficient steps at β = 1000, pnat should be high
        metropolis!(seq, CONTACT_MAP_A; β=1000, nsteps=5000)
        final_pnat = pnat(CONTACT_MAP_A, seq)
        @test final_pnat > 0.95
    end

    @testset "energy is sum of MJ contacts" begin
        # The energy of a sequence in a structure is the sum of MJ interaction
        # energies over the 28 contact pairs.
        seq = potts("EKAMPAMDPDMAHEHKKKIRAWMFEGE")
        for cm in [CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D]
            E = energy(cm, seq)
            E_manual = sum(
                MJ[seq[CONTACT_MAPS[k, cm, 1]], seq[CONTACT_MAPS[k, cm, 2]]]
                for k in 1:N_CONTACTS
            )
            @test E ≈ E_manual
        end
    end

    @testset "contact map structure" begin
        # Each structure has 28 contacts between distinct residues in 1:27.
        for cm in [CONTACT_MAP_A, CONTACT_MAP_B, CONTACT_MAP_C, CONTACT_MAP_D]
            contacts = CONTACT_MAPS[:, cm, :]
            @test size(contacts) == (N_CONTACTS, 2)

            for k in 1:N_CONTACTS
                i, j = contacts[k, 1], contacts[k, 2]
                # Contacts are between valid residue indices
                @test 1 ≤ i ≤ L
                @test 1 ≤ j ≤ L
                # Contacts are between distinct residues
                @test i ≠ j
                # Contacts are ordered (i < j)
                @test i < j
            end

            # All 28 contacts should be unique
            contact_pairs = Set((contacts[k,1], contacts[k,2]) for k in 1:N_CONTACTS)
            @test length(contact_pairs) == N_CONTACTS
        end
    end
end
