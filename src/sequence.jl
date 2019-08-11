const alphabets = Dict(
    "DNA" => "ACGTMRWSYKVHDBN",
    "RNA" => "ACGUMRWSYKVHDBN",
    "AA" => "ARNDCQEGHILKMFPSTWYVBZX"
)

struct Sequence
    seq::String
    type::String
    Sequence(seq, type) =
        all(base -> base in alphabets[type], uppercase(seq)) ?
        new(uppercase(seq), type) : error("String is not a valid $type string!")
end

==(s1::Sequence, s2::Sequence) = (s1.seq == s2.seq) && (s1.type == s2.type)

Base.length(s::Sequence) = length(s.seq)
