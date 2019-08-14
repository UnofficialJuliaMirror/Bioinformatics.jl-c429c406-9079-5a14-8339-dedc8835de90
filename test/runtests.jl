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

        @test collect(values(Bioinformatics.possible_proteins(seq)))[1] == Bioinformatics.Sequence(
            "MTD",
            "AA"
        )
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

        seq = Bioinformatics.Sequence(
            collect(values(Bioinformatics.readFASTA("./example_data/P35858.fasta")))[1],
            "AA"
        )
        seq_stats = Bioinformatics.protparam(seq)
        @test seq_stats["Number of amino acids"] == 605
        @test seq_stats["Molecular weight"] - 66035.0 < 0.01
        @test Bioinformatics.isoelectric_point(seq) - 6.35 < 0.01
        @test seq_stats["# of negatively charged residues"] == 59
        @test seq_stats["# of positively charged residues"] == 54
        @test seq_stats["Extinction coefficient"] == 61555
        @test seq_stats["Instability index"] - 46.09 < 0.01
        @test seq_stats["Aliphatic index"] - 111.652 < 0.01
        @test seq_stats["Grand average of hydropathicity (GRAVY)"] -
              0.018843 < 0.01
    end

end
