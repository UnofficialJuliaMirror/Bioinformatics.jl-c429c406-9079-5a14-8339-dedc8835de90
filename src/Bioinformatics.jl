__precompile__()

module Bioinformatics

export frequency,
       gc_content,
       readFASTA,
       reading_frames,
       readStringFromFile,
       reverse_complement,
       Sequence,
       transcription,
       translation


include("io.jl")
include("sequence.jl")
include("stats.jl")

end
