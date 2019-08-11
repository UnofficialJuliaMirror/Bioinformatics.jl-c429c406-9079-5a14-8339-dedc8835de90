function frequency(seq::Sequence)
    freqs = Dict()
    for s in seq.seq
        if haskey(freqs, s)
            freqs[s] += 1
        else
            freqs[s] = 1
        end
    end
    return freqs
end

function gc_content(seq::Sequence)
    gc_count = 0
    for s in seq.seq
        if s in "GC"
            gc_count += 1
        end
    end
    return gc_count / length(seq)
end
