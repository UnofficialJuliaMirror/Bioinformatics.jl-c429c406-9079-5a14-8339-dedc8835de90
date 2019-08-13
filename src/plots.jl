using Plots

"""
    function plot_gc_content(seq::Sequence, window_size::Int64)

Plot GC-content of given sequence.
"""
function plot_gc_content(seq::Sequence, window_size::Int64)
    x = 1:(length(seq)-window_size)
    y = gc_content(seq, window_size)
    plot(x, y, title = "GC-Content Distribution", labels=["GC-Content"])
end

"""
    function plot_dotmatrix(dotmatrix::Array{Int8})

Plot dotplot for given dotmatrix.
"""
function plot_dotmatrix(dotmatrix::Array{Int8})
    return heatmap(dotmatrix, xlabel="Sequence 1", ylabel="Sequence 2")
end
