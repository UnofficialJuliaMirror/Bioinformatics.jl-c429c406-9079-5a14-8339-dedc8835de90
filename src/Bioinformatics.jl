__precompile__()

module Bioinformatics

export readStringFromFile,
       readFASTA,
       Sequence,
       transcription,
       reverse_complement,
       translation,
       frequency,
       gc_content

include("io.jl")
include("sequence.jl")
include("stats.jl")

end
