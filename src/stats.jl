"""
    frequency(seq::String)

    Calculate frequencies of symbols.
"""
function frequency(seq::String)
    freqs = Dict()
    for s in seq
        if haskey(freqs, s)
            freqs[s] += 1
        else
            freqs[s] = 1
        end
    end
    return freqs
end

"""
    gc_content(seq::String)

    Calculate GC content of a DNA sequence.
"""
function gc_content(seq::String)
    gc_count = 0
    for s in seq
        if s in "GCgc"
            gc_count += 1
        end
    end
    return gc_count / length(seq)
end