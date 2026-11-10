using Test: @test, @testset, @test_throws
using LinearAlgebra: I
using LatticeProteins: onehot, potts, aaseq, AMINO_ACIDS

@testset "sequences" begin
    @testset "onehot - single string" begin
        # Full alphabet produces identity matrix
        @test onehot(AMINO_ACIDS) == I
        # Single amino acid produces a one-hot column
        X = onehot("C")
        @test size(X) == (21, 1)
        @test X[1, 1] == true
        @test sum(X) == 1
        # Short sequence
        X = onehot("CM")
        @test size(X) == (21, 2)
        @test X[1, 1] == true  # C is index 1
        @test X[2, 2] == true  # M is index 2
        @test sum(X) == 2
        # Gap character
        X = onehot("-")
        @test size(X) == (21, 1)
        @test X[21, 1] == true
    end

    @testset "onehot - multiple strings" begin
        seqs = ["CMFILVWYAGTSNQDEHRKP-", "CMFILVWYAGTSNQDEHRKP-"]
        X = onehot(seqs)
        @test size(X) == (21, 21, 2)
        @test X[:, :, 1] == I
        @test X[:, :, 2] == I

        # Different sequences
        seqs2 = ["CCC", "MMM"]
        X2 = onehot(seqs2)
        @test size(X2) == (21, 3, 2)
        @test all(X2[1, :, 1])   # all C in first seq
        @test all(X2[2, :, 2])   # all M in second seq

        # same encoding as one sequence at a time
        seqs3 = ["CMFILVWYAG", "TSNQDEHRKP", "-PKRHEDQNS"]
        X3 = onehot(seqs3)
        @test X3 isa BitArray{3}
        @test all(X3[:, :, n] == onehot(seqs3[n]) for n in eachindex(seqs3))
        @test_throws ArgumentError onehot(["CMF", "IL"])
    end

    @testset "onehot - integer input" begin
        @test onehot(1:21) == onehot(AMINO_ACIDS)

        # Vector input
        v = [1, 2, 3]
        X = onehot(v)
        @test size(X) == (21, 3)
        @test X[1, 1] == true
        @test X[2, 2] == true
        @test X[3, 3] == true

        @test onehot(Int8[1, 2, 3]) == X

        # Matrix input
        M = [1 2; 3 4]
        X = onehot(M)
        @test size(X) == (21, 2, 2)
        @test X[1, 1, 1] == true  # position (1,1) has amino acid 1
        @test X[2, 1, 2] == true  # position (1,2) has amino acid 2
        @test X[3, 2, 1] == true  # position (2,1) has amino acid 3
        @test X[4, 2, 2] == true  # position (2,2) has amino acid 4
    end

    @testset "potts" begin
        @test potts(AMINO_ACIDS) == 1:21

        # potts from BitMatrix (single sequence)
        X = onehot("CMF")
        @test potts(X) == [1, 2, 3]

        # potts from BitArray{3} (multiple sequences)
        seqs = ["CMF", "ILV"]
        X = onehot(seqs)
        P = potts(X)
        @test size(P) == (3, 2)
        @test P[:, 1] == [1, 2, 3]
        @test P[:, 2] == [4, 5, 6]

        # potts from string
        @test potts("CMFILVWYAGTSNQDEHRKP-") == collect(1:21)

        # potts from vector of strings
        P = potts(["CMF", "ILV"])
        @test P isa Matrix{Int}
        @test P[:, 1] == [1, 2, 3]
        @test P[:, 2] == [4, 5, 6]
        @test potts(["CMF", "ILV"]) == potts(onehot(["CMF", "ILV"]))
        @test potts(SubString("XCMFX", 2, 4)) == [1, 2, 3]

        # characters outside the alphabet have no integer code
        @test_throws ArgumentError potts("CMX")
        @test_throws ArgumentError potts("cmf")
        @test_throws ArgumentError potts(["CMF", "ILX"])
        # sequences of a matrix must have the same length
        @test_throws ArgumentError potts(["CMF", "IL"])
        @test_throws ArgumentError potts(String[])
    end

    @testset "aaseq" begin
        # roundtrip: string -> onehot -> aaseq
        seq = "CMFILVWYAGTSNQDEHRKP-"
        @test aaseq(onehot(seq)) == seq

        # aaseq from integer vector
        @test aaseq([1, 2, 3]) == "CMF"
        @test aaseq(collect(1:21)) == AMINO_ACIDS

        # aaseq from integer matrix (multiple sequences)
        P = [1 4; 2 5; 3 6]
        result = aaseq(P)
        @test result == ["CMF", "ILV"]

        # aaseq from BitArray{3} (multiple sequences)
        seqs = ["CMF", "ILV"]
        X = onehot(seqs)
        result = aaseq(X)
        @test result == seqs

        @test aaseq(Int8[1, 2, 3]) == "CMF"
        # a one-hot encoding has 21 rows
        @test_throws ArgumentError aaseq(falses(20, 3))
    end

    @testset "roundtrip consistency" begin
        # string -> onehot -> potts -> aaseq -> string
        seq = "CMFILVWYAGTSNQDEHRKP-"
        @test aaseq(potts(onehot(seq))) == seq

        # string -> potts -> onehot -> aaseq
        @test aaseq(onehot(potts(seq))) == seq

        # Multiple sequences roundtrip
        seqs = ["CMFIL", "VWYAG"]
        X = onehot(seqs)
        @test aaseq(X) == seqs
    end
end
