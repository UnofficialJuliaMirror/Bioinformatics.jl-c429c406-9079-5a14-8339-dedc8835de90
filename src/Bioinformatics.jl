__precompile__()
module Bioinformatics

export dotmatrix,
       edit_dist,
       frequency,
       gravy,
       gc_content,
       hamming_dist,
       instability_index,
       minimum_skew,
       plot_dotmatrix,
       plot_gc_content,
       possible_proteins,
       protein_mass,
       protparam,
       readFASTA,
       reading_frames,
       reverse_complement,
       Sequence,
       skew,
       skew_plot,
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
