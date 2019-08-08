include("../src/Bioinformatics.jl")
using Test

@testset "Tests" begin

    @testset "Alignments" begin
        @test Bioinformatics.globalAlign("AAGT", "AT") == ("AAGT", "-A-T")
        @test Bioinformatics.globalAlign("ATTATCT", "TTTCTA") == ("ATTATCT-", "-TT-TCTA")
        @test Bioinformatics.localAlign("GGTTGACTA", "TGTTACGG") == ("GTTGAC", "GTT-AC")
    end;

    @testset "Distances" begin
        @test Bioinformatics.hammingDist("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7
        @test Bioinformatics.editDist("TGCATAT", "ATCCGAT") == 4
    end;

    @testset "IO" begin
        @test Bioinformatics.countBases(Bioinformatics.readStringFromFile(joinpath(@__DIR__, "data", "rosalind_dna.txt"))) == 
                (201, 182, 223, 195)
        @test Bioinformatics.gcContent(Bioinformatics.readFASTA(joinpath(@__DIR__, "data", "rosalind_gc.txt")))["Rosalind_1093"] == 
                48.46560846560847
    end;

    @testset "Numeric" begin
        @test Bioinformatics.countBases("AGTCGACG") == (2, 2, 3, 1)
        @test_throws ErrorException Bioinformatics.countBases("AGTGCFF")
        @test round(Bioinformatics.gcContent("CCACCCTCGTGGTATGGCTAGGCATTCAG"), digits = 2) == 58.62
        @test Bioinformatics.gcContent(Dict("DNA_1" => "AGTC", "DNA_2" => "AATG")) == Dict("DNA_1" => 50.0, "DNA_2" => 25.0)
    end;

    @testset "String" begin
        @test Bioinformatics.verifyDna("ACGT") == true
        @test_throws ErrorException Bioinformatics.verifyDna("ACFY")
        @test Bioinformatics.transcribeDnaToRna("GATGGAACTT") == "GAUGGAACUU"
        @test Bioinformatics.translateRNA(Bioinformatics.transcribeDnaToMRna("AAACTTGAATAAACGT")) == "FELIC"
        @test Bioinformatics.reverseComplement("AAAACCCGGT") == "ACCGGGTTTT"
        @test Bioinformatics.findMotif("GATATATGCATATACTT", "ATAT") == [2, 4, 10]
        @test Bioinformatics.profileMatrix(Dict("Rosalind_1" => "ATCCAGCT", 
            "Rosalind_2" => "GGGCAACT", 
            "Rosalind_3" => "ATGGATCT", 
            "Rosalind_4" => "AAGCAACC", 
            "Rosalind_5" => "TTGGAACT", 
            "Rosalind_6" => "ATGCCATT", 
            "Rosalind_7" => "ATGGCACT")) ==
            [[5.0  1.0  0.0  0.0  5.0  5.0  0.0  0.0]
             [0.0  0.0  1.0  4.0  2.0  0.0  6.0  1.0]
             [1.0  1.0  6.0  3.0  0.0  1.0  0.0  0.0]
             [1.0  5.0  0.0  0.0  0.0  1.0  1.0  6.0]]
        @test round(Bioinformatics.proteinMass("SKADYEK"), digits = 2) == 821.39
        @test Bioinformatics.longestCommonSubstring(Dict("Rosalind_1" => "AATCCACT",
             "Rosalind_2" => "GGCAACT",
             "Rosalind_3" => "AATCT",
             "Rosalind_4" => "AAAACC")) == "AA"
        @test Bioinformatics.allKmers("GTAGAGCTGT", 5) == 
            Set(["GTAGA", "TAGAG", "AGAGC", "GAGCT", "AGCTG", "GCTGT"])
        @test Bioinformatics.kmerCounts("AAGTCGTCGACGT", 4) == ["GTCG" => 2]
        @test Bioinformatics.ltClumps("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4) == 
            Set(["GAAGA", "CGACA"])
    end;

end;