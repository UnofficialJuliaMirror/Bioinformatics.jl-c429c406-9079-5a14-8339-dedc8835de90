__precompile__()
module Bioinformatics

export aliphatic_index,
       dotmatrix,
       edit_dist,
       extinction_coeff,
       frequency,
       gravy,
       gc_content,
       hamming_dist,
       instability_index,
       plot_dotmatrix,
       plot_gc_content,
       possible_proteins,
       protein_mass,
       protparam,
       readFASTA,
       reading_frames,
       readStringFromFile,
       reverse_complement,
       Sequence,
       transcription,
       translation

include("sequence.jl")
include("data.jl")
include("alignments.jl")
include("distances.jl")
include("io.jl")
include("plots.jl")
include("stats.jl")

end
