# Functions

## Sequence

* Sequence constructor:

  ```s1 = Bioinformatics.Sequence("CGATATAGATT", "DNA")```

* ```function transcription(dna_seq::Sequence)```

  Transcribe DNA to RNA.

* ```function reverse_complement(s::Sequence)```

  Reverse complement of given DNA/RNA sequence.

* ```function translation(s::Sequence, start_pos::Int64 = 1)```

  Translate given DNA sequence to Amino Acid sequence.

* ```function reading_frames(s::Sequence)```

  Find all reading frames of given sequence.

* ```function possible_proteins(s::Sequence)```

  Find possible proteins of given sequence.

## Alignments

* ```dotmatrix(s1::Sequence, s2::Sequence)```

  Calculate the dotplot matrix for given two sequences.

## Distances

* ```edit_dist(s1::String, s2::String)```

  Calculate edit distance between two strings.

* ```hamming_dist(s1::String, s2::String)```

  Calculate Hamming distance between two strings.

## IO

* ```function readFASTA(filename::String)```

  Parse a FASTA formatted file.

## Plots

* ```function plot_gc_content(seq::Sequence, window_size::Int64)```

  Plot GC-content of given sequence.

* ```function plot_dotmatrix(dotmatrix::Array{Int8})```

  Plot dotplot for given dotmatrix.

## Stats

* ```function frequency(seq::Sequence)```

  Count bases / amino acids for given sequence.

* ```function gc_content(seq::Sequence)```

  Calculate GC-content of given sequence.

* ```function gc_content(seq::Sequence, window_size::Int64)```

  Calculate GC-content of given sequence using a window.

* ```function protein_mass(seq::Sequence, type = "monoisotopic")```

  Calculate mass of given amino acid sequence.

* ```function extinction_coeff(n_tyr, n_trp, n_cys::Int64)```

  The extinction coefficient indicates how much light a protein absorbs at a certain wavelength. It is useful to have an estimation of this coefficient for following a protein which a spectrophotometer when purifying it.

* ```function instability_index(seq::Sequence)```

  The instability index provides an estimate of the stability of your protein in a test tube.

* ```function aliphatic_index(x_ala, x_val, x_ile, x_leu::Float64)```

  The aliphatic index of a protein is defined as the relative volume occupied by aliphatic side chains (alanine, valine, isoleucine, and leucine).

* ```function gravy(seq::Sequence)```

  The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values of all the amino acids, divided by the number of residues in the sequence.

* ```function isoelectric_point(seq::Sequence)```

  The isoelectric point, is the pH at which a molecule carries no net electrical charge or is electrically neutral in the statistical mean.

* ```function protparam(seq::Sequence)```

  Computes various physico-chemical properties that can be deduced from a protein sequence.

# References

* [ExPASy](https://web.expasy.org/protparam/protparam-doc.html)
* [Isoelectric](http://isoelectric.org/algorithms.html)
