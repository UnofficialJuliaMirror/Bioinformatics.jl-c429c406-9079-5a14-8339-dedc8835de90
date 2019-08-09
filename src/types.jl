const dnaAlphabet = ['A', 'C', 'G', 'T']
const rnaAlphabet = ['A', 'C', 'G', 'U']

const proteinAlphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

abstract type Seq end   

struct DNASeq <: Seq
    seq::String
    DNASeq(seq) = all(base->base in dnaAlphabet, seq) ? new(seq) : error("String is not a valid DNA string!")
end

struct RNASeq <: Seq
    seq::String
    RNASeq(seq) = all(base->base in rnaAlphabet, seq) ? new(seq) : error("String is not a valid RNA string!")
end

struct ProteinSeq <: Seq
    seq::String
    ProteinSeq(seq) = all(base->base in proteinAlphabet, seq) ? new(seq) : error("String is not a valid protein string!")
end
