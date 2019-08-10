include("../src/Bioinformatics.jl")
using Test

@testset "Tests" begin

    @testset "stats.jl" begin
        seq = "atagataactcgcatag"
        
        freqs = Bioinformatics.frequency(seq)
        @test freqs['a'] == 7
        @test freqs['c'] == 3
        @test freqs['g'] == 3
        @test freqs['t'] == 4

        seq = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
        gc_content = Bioinformatics.gc_content(seq)
        @test round(gc_content, digits = 2) == 0.61
    end;

end;