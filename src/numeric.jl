const proteinMassTable = Dict('A' => 71.03711, 'C' => 103.00919,
                              'D' => 115.02694, 'E' => 129.04259,
                              'F' => 147.06841, 'G' => 57.02146,
                              'H' => 137.05891, 'I' => 113.08406,
                              'K' => 128.09496, 'L' => 113.08406,
                              'M' => 131.04049, 'N' => 114.04293,
                              'P' => 97.05276, 'Q' => 128.05858,
                              'R' => 156.10111, 'S' => 87.03203,
                              'T' => 101.04768, 'V' => 99.06841,
                              'W' => 186.07931, 'Y' => 163.06333)

"""
  countBases(dna::String)

  Counts the amount of each base in a DNA chain.
"""
function countBases(dna::String)
  countA = 0
  countC = 0
  countG = 0
  countT = 0
  for c in dna
      if c == 'A'
          countA += 1
      elseif c == 'C'
          countC += 1
      elseif c == 'G'
          countG += 1
      elseif c == 'T'
          countT += 1
      else
          error("Char $c is not a valid base!")
      end
  end
  return countA, countC, countG, countT
end

"""
    gcContent(dna::String)

    Calculates G+C ratio.
"""
function gcContent(dna::String)
    n = length(dna)
    m = 0
    for b in dna
        if b == 'G' || b == 'C'
            m += 1
        end
    end
    return 100 * m / n
end

"""
    gcContent(dna::String)

    Calculates G+C ratio.
"""
function gcContent(dna::Dict)
    data = Dict()
    for (k, v) in dna
        data[k] = gcContent(v)
    end
    return data
end

"""
    proteinMass(protein::String)

    Calculates mass of given protein string.
"""
function proteinMass(protein::String)
    mass = 0.0
    for aa in protein
        mass += proteinMassTable[aa]
    end
    return mass
end
