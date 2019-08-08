const dnaAlphabet = ['A', 'C', 'G', 'T']
const rnaAlphabet = ['A', 'C', 'G', 'U']

const dnaComplements = Dict('A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C')
const rnaComplements = Dict('A' => 'U', 'U' => 'A', 'C' => 'G', 'G' => 'C')

const proteinAlphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                         'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

const dnaCodonTable = Dict("TTT" => 'F',      "CTT" => 'L',      "ATT" => 'I',
                           "GTT" => 'V',      "TTC" => 'F',      "CTC" => 'L',
                           "ATC" => 'I',      "GTC" => 'V',      "TTA" => 'L',
                           "CTA" => 'L',      "ATA" => 'I',      "GTA" => 'V',
                           "TTG" => 'L',      "CTG" => 'L',      "ATG" => 'M',
                           "GTG" => 'V',      "TCT" => 'S',      "CCT" => 'P',
                           "ACT" => 'T',      "GCT" => 'A',      "TCC" => 'S',
                           "CCC" => 'P',      "ACC" => 'T',      "GCC" => 'A',
                           "TCA" => 'S',      "CCA" => 'P',      "ACA" => 'T',
                           "GCA" => 'A',      "TCG" => 'S',      "CCG" => 'P',
                           "ACG" => 'T',      "GCG" => 'A',      "TAT" => 'Y',
                           "CAT" => 'H',      "AAT" => 'N',      "GAT" => 'D',
                           "TAC" => 'Y',      "CAC" => 'H',      "AAC" => 'N',
                           "GAC" => 'D',      "TAA" => "Stop",   "CAA" => 'Q',
                           "AAA" => 'K',      "GAA" => 'E',      "TAG" => "Stop",
                           "CAG" => 'Q',      "AAG" => 'K',      "GAG" => 'E',
                           "TGT" => 'C',      "CGT" => 'R',      "AGT" => 'S',
                           "GGT" => 'G',      "TGC" => 'C',      "CGC" => 'R',
                           "AGC" => 'S',      "GGC" => 'G',      "TGA" => "Stop",
                           "CGA" => 'R',      "AGA" => 'R',      "GGA" => 'G',
                           "TGG" => 'W',      "CGG" => 'R',      "AGG" => 'R',
                           "GGG" => 'G')

const rnaCodonTable = Dict("UUU" => 'F',      "CUU" => 'L',      "AUU" => 'I',
                           "GUU" => 'V',      "UUC" => 'F',      "CUC" => 'L',
                           "AUC" => 'I',      "GUC" => 'V',      "UUA" => 'L',
                           "CUA" => 'L',      "AUA" => 'I',      "GUA" => 'V',
                           "UUG" => 'L',      "CUG" => 'L',      "AUG" => 'M',
                           "GUG" => 'V',      "UCU" => 'S',      "CCU" => 'P',
                           "ACU" => 'T',      "GCU" => 'A',      "UCC" => 'S',
                           "CCC" => 'P',      "ACC" => 'T',      "GCC" => 'A',
                           "UCA" => 'S',      "CCA" => 'P',      "ACA" => 'T',
                           "GCA" => 'A',      "UCG" => 'S',      "CCG" => 'P',
                           "ACG" => 'T',      "GCG" => 'A',      "UAU" => 'Y',
                           "CAU" => 'H',      "AAU" => 'N',      "GAU" => 'D',
                           "UAC" => 'Y',      "CAC" => 'H',      "AAC" => 'N',
                           "GAC" => 'D',      "UAA" => "Stop",   "CAA" => 'Q',
                           "AAA" => 'K',      "GAA" => 'E',      "UAG" => "Stop",
                           "CAG" => 'Q',      "AAG" => 'K',      "GAG" => 'E',
                           "UGU" => 'C',      "CGU" => 'R',      "AGU" => 'S',
                           "GGU" => 'G',      "UGC" => 'C',      "CGC" => 'R',
                           "AGC" => 'S',      "GGC" => 'G',      "UGA" => "Stop",
                           "CGA" => 'R',      "AGA" => 'R',      "GGA" => 'G',
                           "UGG" => 'W',      "CGG" => 'R',      "AGG" => 'R',
                           "GGG" => 'G')

"""
    verifyDna(dna::String)

    Verifies that the string contains only 'A', 'T', 'G', 'C' characters.
"""
function verifyDna(dna::String)
    if all(base->base in dnaAlphabet, dna)
        return true
    else
        error("String is not a valid DNA string!")
    end
end

"""
    reverseComplement(dna::String)

    Calculates a reverse complement of a DNA string.
"""
function reverseComplement(dna::String)
    reversed = reverse(dna)
    reverseComplement = join([dnaComplements[c] for c in reversed])
    return reverseComplement
end

"""
    reversePalindrome(dna::String)

    Finds reverse palindrome substrings of a DNA string.
"""
function reversePalindrome(dna::String)
    res = Dict()
    for i = 1:length(dna)
        for j = 4:length(dna)
            if i + j > length(dna) + 1
                continue
            else
                subStr = dna[i:i + j - 1]
                lenSubStr = length(subStr)
                revComp = reverseComplement(subStr)
                if subStr == revComp
                    res[i] = subStr
                end
            end
        end
    end
    return res
end

"""
    allKmers(dna::String, k::Int)

    Returns all the possible substrings of length k that are contained in
    a given DNA string.
"""
function allKmers(dna::String, k::Int)
    kmers = Set{String}()
    for i = 1:length(dna) - k + 1
        push!(kmers, dna[i:i + k - 1])
    end
    return kmers
end

"""
    kmerCounts(dna::String, k::Int)

    Returns counts of all k-mers which more than one.
"""
function kmerCounts(dna::String, k::Int64)
    patterns = Dict()
    for i = 1:length(dna) - k
        kmer = dna[i:i + k - 1]
        if kmer in keys(patterns)
            continue
        else
            indexes = findMotif(dna, kmer)
            patterns[kmer] = length(indexes)
        end
    end
    mostFreqPatterns = filter((k, v)->v > 1, patterns)
    return sort(collect(mostFreqPatterns), by = x->x[2], rev = true)
end

"""
    transcribeDnaToRna(dna::String)

    Transcibes a DNA string to RNA.
"""
function transcribeDnaToRna(dna::String)
    rna = ""
    for c in dna
        if c == 'T'
            rna = string(rna, 'U')
        else
            rna = string(rna, c)
        end
    end
    return rna
end

"""
    transcribeDnaToMRna(dna::String)

    Transcibes a DNA string to mRNA.
"""
function transcribeDnaToMRna(dna::String)
    rna = transcribeDnaToRna(dna)
    mRNA = ""
    for b in rna
        mRNA = string(mRNA, rnaComplements[b])
    end
    return mRNA
end

"""
    translateDNA(dna::String)

    Translates a DNA string to protein string.
"""
function translateDNA(dna::String)
    proteinString = ""
    for i = 1:3:length(dna) - 2
        aa = dnaCodonTable[dna[i:i + 2]]
        if aa == "Stop"
            continue
        else
            proteinString = string(proteinString, aa)
        end
    end
    return proteinString
end

"""
    translateRNA(rna::String)

    Translates a RNA string to protein string.
"""
function translateRNA(rna::String)
    proteinString = ""
    for i = 1:3:length(rna) - 2
        aa = rnaCodonTable[rna[i:i + 2]]
        if aa == "Stop"
            continue
        else
            proteinString = string(proteinString, aa)
        end
    end
    return proteinString
end

"""
    findMotif(s::String, t::String)

    Returns all locations of t as a substring of s.
"""
function findMotif(s::String, t::String)
    indexes = []
    for i = 1:(length(s) - length(t) + 1)
        if s[i:i + length(t) - 1] == t
            append!(indexes, i)
        end
    end
    return indexes
end

"""
    profileMatrix(dna::Dict)

    Calculates profile matrix of a collection of DNA strings.
"""
function profileMatrix(dna::Dict)
    colSize = length(collect(Iterators.take(dna, 1))[1][2])
    profile = zeros(4, colSize)
    for k in keys(dna)
        str = (dna[k])
        for i = 1:length(str)
            if str[i] == 'A'
                profile[1, i] += 1
            elseif str[i] == 'C'
                profile[2, i] += 1
            elseif str[i] == 'G'
                profile[3, i] += 1
            else
                profile[4, i] += 1
            end
        end
    end
    return profile
end

"""
    consensusString(profileMatrix::Array)

    Calculates consensus string from a profile matrix.
"""
function consensusString(profileMatrix::Array)
    dict = Dict([1] => 'A', [2] => 'C',
              [3] => 'G', [4] => 'T',
              [1,2] => 'M', [1,3] => 'R',
              [1,4] => 'W', [2,3] => 'S',
              [2,4] => 'Y', [3,4] => 'K',
              [1,2,3] => 'V', [1,2,4] => 'H',
              [1,3,4] => 'D', [2,3,4] => 'B',
              [1,2,3,4] => 'N')
    rowSize, colSize = size(profileMatrix)
    consensus = ""
    for i = 1:colSize
        ind = find(b->b == maximum(profileMatrix[:, i]), profileMatrix[:, i])
        b = dict[ind]
        consensus = string(consensus, b)
    end
    return consensus
end

"""
    longestCommonSubstring(dna::Dict)

    Calculates longest common substring of the collection.
"""
function longestCommonSubstring(dna::Dict)
    v = collect(values(dna))
    sort!(v, by = length)
    motif = ""
    shortestStr = v[1]
    otherStr = v[2:end]
    len = length(shortestStr)
    for i = 1:len
        for j = i:len
            m = shortestStr[i:j]
            found = false
            for dna in otherStr
                if occursin(m, dna)
                    found = true
                else
                    found = false
                    break
                end
            end
            if found & (length(m) > length(motif))
                motif = m
            end
        end
    end
    return motif
end

"""
    ltClumps(dna::String, k::Int64, L::Int64, t::Int64)

    Given integers L and k, a k-mer Pattern forms an (L, t)-clump inside
    a longer string if there is an interval of length L in which this k-mer
    appears at least t times.[^1]

[^1]: Compeau, P., & Pevzner, P. (2015). Bioinformatics algorithms:
      An active learning approach. La Jolla: Active Learning.
"""
function ltClumps(dna::String, k::Int64, L::Int64, t::Int64)
    result = Set()
    for i = 1:length(dna) - L
        window = dna[i:i + L - 1]
        for j = 1:length(window) - k
            kmer = window[j:j + k - 1]
            if length(findMotif(window, kmer)) >= t
                push!(result, kmer)
            end
        end
    end
    return result
end
