import Base.==

const alphabets = Dict(
    "DNA" => "ACGTMRWSYKVHDBN",
    "RNA" => "ACGUMRWSYKVHDBN",
    "AA" => "ACDEFGHIKLMNPQRSTVWY-"
)

const complements = Dict(
    "DNA" => Dict(
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'W' => 'W',
        'S' => 'S',
        'M' => 'K',
        'K' => 'M',
        'R' => 'Y',
        'Y' => 'R',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N'
    ),
    "RNA" => Dict(
        'A' => 'U',
        'C' => 'G',
        'G' => 'C',
        'U' => 'A',
        'W' => 'W',
        'S' => 'S',
        'M' => 'K',
        'K' => 'M',
        'R' => 'Y',
        'Y' => 'R',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N'
    )
)

const codons = Dict(
    "GCT" => 'A',
    "GCC" => 'A',
    "GCA" => 'A',
    "GCG" => 'A',
    "TGT" => 'C',
    "TGC" => 'C',
    "GAT" => 'D',
    "GAC" => 'D',
    "GAA" => 'E',
    "GAG" => 'E',
    "TTT" => 'F',
    "TTC" => 'F',
    "GGT" => 'G',
    "GGC" => 'G',
    "GGA" => 'G',
    "GGG" => 'G',
    "CAT" => 'H',
    "CAC" => 'H',
    "ATA" => 'I',
    "ATT" => 'I',
    "ATC" => 'I',
    "AAA" => 'K',
    "AAG" => 'K',
    "TTA" => 'L',
    "TTG" => 'L',
    "CTT" => 'L',
    "CTC" => 'L',
    "CTA" => 'L',
    "CTG" => 'L',
    "ATG" => 'M',
    "AAT" => 'N',
    "AAC" => 'N',
    "CCT" => 'P',
    "CCC" => 'P',
    "CCA" => 'P',
    "CCG" => 'P',
    "CAA" => 'Q',
    "CAG" => 'Q',
    "CGT" => 'R',
    "CGC" => 'R',
    "CGA" => 'R',
    "CGG" => 'R',
    "AGA" => 'R',
    "AGG" => 'R',
    "TCT" => 'S',
    "TCC" => 'S',
    "TCA" => 'S',
    "TCG" => 'S',
    "AGT" => 'S',
    "AGC" => 'S',
    "ACT" => 'T',
    "ACC" => 'T',
    "ACA" => 'T',
    "ACG" => 'T',
    "GTT" => 'V',
    "GTC" => 'V',
    "GTA" => 'V',
    "GTG" => 'V',
    "TGG" => 'W',
    "TAT" => 'Y',
    "TAC" => 'Y',
    "TAA" => '-',
    "TAG" => '-',
    "TGA" => '-'
)

struct Sequence
    seq::String
    type::String
    Sequence(seq, type) =
        all(base -> base in alphabets[type], uppercase(seq)) ?
        new(uppercase(seq), type) : error("String is not a valid $type string!")
end

==(s1::Sequence, s2::Sequence) = (s1.seq == s2.seq) && (s1.type == s2.type)

Base.getindex(s::Sequence, i::Int64) = getindex(s.seq, i)
Base.getindex(s::Sequence, r::UnitRange{Int64}) = getindex(s.seq, r)
Base.length(s::Sequence) = length(s.seq)
Base.replace(s::Sequence, p::Pair{Char,Char}) = replace(s.seq, p)
Base.reverse(s::Sequence) = reverse(s.seq)

function transcription(dna_seq::Sequence)
    if dna_seq.type != "DNA"
        error("Only DNA sequences can be transcribed.")
    end
    return Sequence(replace(dna_seq, 'T' => 'U'), "RNA")
end

function reverse_complement(s::Sequence)
    if (s.type == "AA")
        error("Amino acid sequence doesn't have reverse complement")
    end
    len = length(s)
    rev = reverse(s)
    reverse_complement = Char['N' for i in 1:len]
    for i in 1:len
        reverse_complement[i] = complements[s.type][rev[i]]
    end
    return Sequence(string(reverse_complement...), s.type)
end

function translation(s::Sequence, start_pos::Int64 = 1)
    if (s.type == "AA")
        error("Amino acid sequence cannot be translated.")
    end
    len = length(s)
    translated_seq = Char[]
    for i in start_pos:3:(len - 2)
        cod = s[i:i+2]
        push!(translated_seq, codons[cod])
    end
    return Sequence(string(translated_seq...), "AA")
end

function reading_frames(s::Sequence)
    if (s.type == "AA")
        error("Amino acid sequence cannot be translated.")
    end
    reading_frames = Dict{String,Sequence}()
    rc = reverse_complement(s)
    for i in 1:3
        reading_frames["5'3' Frame $i"] = translation(s, i)
        reading_frames["3'5' Frame $i"] = translation(rc, i)
    end
    return reading_frames
end

function possible_proteins(s::Sequence)
    if (s.type == "AA")
        error("Amino acid sequence cannot be translated.")
    end
    regex = r"M[ACDEFGHIKLMNPQRSTVWY]+(?=-|$)"
    rfs = reading_frames(s)
    possible_proteins = Dict{String,Sequence}()
    for k in keys(rfs)
        m = match(regex, rfs[k].seq)
        if m != nothing
            possible_proteins[k] = Sequence(m.match, "AA")
        end
    end
    return possible_proteins
end
