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
