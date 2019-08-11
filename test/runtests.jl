include("../src/Bioinformatics.jl")
using Test

@testset "Tests" begin

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
