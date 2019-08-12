__precompile__()

module Bioinformatics

export edit_dist,
       frequency,
       gc_content,
       hamming_dist,
       possible_proteins,
       protein_mass,
       readFASTA,
       reading_frames,
       readStringFromFile,
       reverse_complement,
       Sequence,
       transcription,
       translation

include("distances.jl")
include("io.jl")
include("sequence.jl")
include("stats.jl")

end
