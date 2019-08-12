__precompile__()

module Bioinformatics

export dotmatrix,
       edit_dist,
       frequency,
       gc_content,
       hamming_dist,
       plot_dotmatrix,
       plot_gc_content,
       possible_proteins,
       protein_mass,
       readFASTA,
       reading_frames,
       readStringFromFile,
       reverse_complement,
       Sequence,
       transcription,
       translation

include("sequence.jl")
include("alignments.jl")
include("distances.jl")
include("io.jl")
include("plots.jl")
include("stats.jl")

end
