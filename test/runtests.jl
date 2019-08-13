include("../src/Bioinformatics.jl")
using Test

@testset "Tests" begin

    @testset "distances.jl" begin
        s1 = "kitten"
        s2 = "sitting"
        @test Bioinformatics.edit_dist(s1, s2) == 3

        s3 = "karolin"
        s4 = "kathrin"
        @test Bioinformatics.hamming_dist(s3, s4) == 3
    end

    @testset "sequence.jl" begin
        seq = Bioinformatics.Sequence("ATGACAGAT", "DNA")

        @test Bioinformatics.transcription(seq) == Bioinformatics.Sequence(
            "AUGACAGAU",
            "RNA"
        )
        @test Bioinformatics.reverse_complement(seq) == Bioinformatics.Sequence(
            "ATCTGTCAT",
            "DNA"
        )

        @test Bioinformatics.translation(seq) == Bioinformatics.Sequence(
            "MTD",
            "AA"
        )

        @test Bioinformatics.kmers(seq, 8) == ["ATGACAGA", "TGACAGAT"]

        @test collect(values(Bioinformatics.possible_proteins(seq)))[1] == Bioinformatics.Sequence("MTD", "AA")
    end

    @testset "stats.jl" begin
        seq = Bioinformatics.Sequence("atagataactcgcatag", "DNA")

        freqs = Bioinformatics.frequency(seq)
        @test freqs['A'] == 7
        @test freqs['C'] == 3
        @test freqs['G'] == 3
        @test freqs['T'] == 4

        seq = Bioinformatics.Sequence(
            "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGG",
            "DNA"
        )
        gc_content = Bioinformatics.gc_content(seq)
        @test round(gc_content, digits = 2) == 0.61
    end

end
