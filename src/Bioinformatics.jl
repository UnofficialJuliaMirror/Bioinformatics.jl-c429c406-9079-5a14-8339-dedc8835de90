__precompile__()

module Bioinformatics
export allKmers,
       consensusString,
       countBases,
       editDist,
       findMotif,
       gcContent,
       globalAlign,
       hammingDist,
       kmerCounts,
       localAlign,
       longestCommonSubstring,
       ltClumps,
       minimumSkew,
       profileMatrix,
       proteinMass,
       readFASTA,
       readStringFromFile,
       reverseComplement,
       reversePalindrome,
       transcribeDnaToMRna,
       transcribeDnaToRna,
       translateDNA,
       translateRNA,
       verifyDna

include("align.jl")
include("distances.jl")
include("io.jl")
include("numeric.jl")
include("string.jl")

end
