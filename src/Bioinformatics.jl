__precompile__()

module Bioinformatics

export Sequence, frequency, gc_content

include("io.jl")
include("sequence.jl")
include("stats.jl")

end
