function frequency(seq::Sequence)
    freqs = Dict([c => 0 for c in alphabets[seq.type]])
    for s in seq.seq
        freqs[s] += 1
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

function gc_content(seq::Sequence, window_size::Int64)
    gc_contents = Float64[]
    for i in 1:(length(seq)-window_size)
        window = Sequence(seq[i:(i+window_size)], seq.type)
        push!(gc_contents, gc_content(window))
    end
    return gc_contents
end

function protein_mass(seq::Sequence, type = "monoisotopic")
    if seq.type != "AA"
        error("Sequence must be a protein sequence.")
    end
    if type == "monoisotopic"
        mass = 18.01524
    else
        mass = 18.01056
    end
    for aa in seq.seq
        mass += aa_mass[type][aa]
    end
    return mass
end

function extinction_coeff(n_tyr, n_trp, n_cys::Int64)
    return 1490 * n_tyr + 5500 * n_trp + 125 * n_cys
end

function instability_index(seq::Sequence)
    ii = 10 / length(seq) *
         sum([instability_values[seq[i:(i+1)]] for i in 1:(length(seq)-1)])
    return ii
end

function aliphatic_index(x_ala, x_val, x_ile, x_leu::Float64)
    return x_ala + 2.9 * x_val + 3.9 * (x_ile + x_leu)
end

function gravy(seq::Sequence)
    return sum([hydropathicity[aa] for aa in seq.seq]) / length(seq)
end

function protparam(seq::Sequence)
    statistics = Dict()
    statistics["Number of amino acids"] = length(seq)
    statistics["Molecular weight"] = protein_mass(seq, "average")
    statistics["Amino acid composition"] = frequency(seq)
    statistics["Total number of negatively charged residues (Asp + Glu)"] = statistics["Amino acid composition"]['D'] +
                                                                            statistics["Amino acid composition"]['E']
    statistics["Total number of positively charged residues (Arg + Lys)"] = statistics["Amino acid composition"]['R'] +
                                                                            statistics["Amino acid composition"]['K']
    statistics["Extinction coefficient"] = extinction_coeff(
        statistics["Amino acid composition"]['Y'],
        statistics["Amino acid composition"]['W'],
        statistics["Amino acid composition"]['C']
    )
    statistics["Instability index"] = instability_index(seq)
    statistics["Aliphatic index"] = aliphatic_index(
        100 * statistics["Amino acid composition"]['A'] /
        statistics["Number of amino acids"],
        100 * statistics["Amino acid composition"]['V'] /
        statistics["Number of amino acids"],
        100 * statistics["Amino acid composition"]['I'] /
        statistics["Number of amino acids"],
        100 * statistics["Amino acid composition"]['L'] /
        statistics["Number of amino acids"]
    )
    statistics["Grand average of hydropathicity (GRAVY)"] = gravy(seq)
    return statistics
end
