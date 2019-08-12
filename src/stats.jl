const aa_mass = Dict(
    "monoisotopic" => Dict(
        'A' => 71.03711,
        'R' => 156.10111,
        'N' => 114.04293,
        'D' => 115.02694,
        'C' => 103.00919,
        'E' => 129.04259,
        'Q' => 128.05858,
        'G' => 57.02146,
        'H' => 137.05891,
        'I' => 113.08406,
        'L' => 113.08406,
        'K' => 128.09496,
        'M' => 131.04049,
        'F' => 147.06841,
        'P' => 97.05276,
        'S' => 87.03203,
        'T' => 101.04768,
        'W' => 186.07931,
        'Y' => 163.06333,
        'V' => 99.06841
    ),
    "average" => Dict(
        'A' => 71.0788,
        'R' => 156.1875,
        'N' => 114.1038,
        'D' => 115.0886,
        'C' => 103.1388,
        'E' => 129.1155,
        'Q' => 128.1307,
        'G' => 57.0519,
        'H' => 137.1411,
        'I' => 113.1594,
        'L' => 113.1594,
        'K' => 128.1741,
        'M' => 131.1926,
        'F' => 147.1766,
        'P' => 97.1167,
        'S' => 87.0782,
        'T' => 101.1051,
        'W' => 186.2132,
        'Y' => 163.1760,
        'V' => 99.1326
    )
)

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

function protein_mass(seq::Sequence, type="monoisotopic")
    if seq.type != "AA"
        error("Sequence must be a protein sequence.")
    end
    if type=="monoisotopic"
        mass = 	18.01524
    else
        mass = 18.01056
    end
    for aa in seq.seq
        mass += aa_mass[type][aa]
    end
    return mass
end
