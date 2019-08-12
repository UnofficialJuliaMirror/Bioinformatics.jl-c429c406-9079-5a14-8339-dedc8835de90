include("sequence.jl")

using Plots

function plot_gc_content(seq::Sequence, window_size::Int64)
    x = 1:(length(seq)-window_size)
    y = gc_content(seq, window_size)
    plot(x, y, title = "GC-Content Distribution", labels=["GC-Content"])
end
