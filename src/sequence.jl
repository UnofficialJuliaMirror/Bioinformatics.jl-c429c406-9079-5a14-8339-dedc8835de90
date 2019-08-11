const alphabets = Dict(
    "DNA" => "ACGTMRWSYKVHDBN",
    "RNA" => "ACGUMRWSYKVHDBN",
    "AA" => "ARNDCQEGHILKMFPSTWYVBZX"
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

struct Sequence
    seq::String
    type::String
    Sequence(seq, type) =
        all(base -> base in alphabets[type], uppercase(seq)) ?
        new(uppercase(seq), type) : error("String is not a valid $type string!")
end

==(s1::Sequence, s2::Sequence) = (s1.seq == s2.seq) && (s1.type == s2.type)

Base.getindex(s::Sequence, i::Int64) = getindex(s.seq, i)
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
    if (s.type != "DNA") && (s.type != "RNA")
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
