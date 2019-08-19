var documenterSearchIndex = {"docs":
[{"location":"#Bioinformatics.jl-Documentation-1","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.jl Documentation","text":"","category":"section"},{"location":"#","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.jl Documentation","text":"CurrentModule = Bioinformatics\nDocTestSetup = quote\n    using Bioinformatics\nend","category":"page"},{"location":"#","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.jl Documentation","text":"Modules = [Bioinformatics]\nPrivate = false","category":"page"},{"location":"#Bioinformatics.Sequence","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.Sequence","text":"Structure for sequence data.\n\n\n\n\n\n","category":"type"},{"location":"#Bioinformatics.dotmatrix-Tuple{Sequence,Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.dotmatrix","text":"function dotmatrix(s1::Sequence, s2::Sequence)\n\nCalculate the dotplot matrix for given two sequences.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.edit_dist-Tuple{String,String}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.edit_dist","text":"function edit_dist(s1::String, s2::String)\n\nCalculate edit distance between two strings.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.frequency-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.frequency","text":"function frequency(seq::Sequence)\n\nCount bases / amino acids for given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.gc_content-Tuple{Sequence,Int64}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.gc_content","text":"function gc_content(seq::Sequence, window_size::Int64)\n\nCalculate GC-content of given sequence using a window.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.gc_content-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.gc_content","text":"function gc_content(seq::Sequence)\n\nCalculate GC-content of given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.global_alignment_linear_gap-Tuple{Sequence,Sequence,Any,Any}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.global_alignment_linear_gap","text":"function global_alignment_linear_gap(seq1, seq2, sm, d)\n\nNeedleman-Wunsch algorithm with linear gap penalty.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.gravy-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.gravy","text":"function gravy(seq::Sequence)\n\nThe GRAVY value for a peptide or protein is calculated as the sum of hydropathy values of all the amino acids, divided by the number of residues in the sequence.\n\nSee also: https://web.expasy.org/protparam/protparam-doc.html.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.hamming_dist-Tuple{String,String}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.hamming_dist","text":"function hamming_dist(s1::String, s2::String)\n\nCalculate Hamming distance between two strings.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.instability_index-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.instability_index","text":"function instability_index(seq::Sequence)\n\nThe instability index provides an estimate of the stability of your protein in a test tube.\n\nSee also: https://web.expasy.org/protparam/protparam-doc.html.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.minimum_skew-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.minimum_skew","text":"function minimum_skew(seq::Sequence)\n\nFind genome positions where the skew is minimum.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.plot_dotmatrix-Tuple{Array{Int8,N} where N}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.plot_dotmatrix","text":"function plot_dotmatrix(dotmatrix::Array{Int8})\n\nPlot dotplot for given dotmatrix.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.plot_gc_content-Tuple{Sequence,Int64}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.plot_gc_content","text":"function plot_gc_content(seq::Sequence, window_size::Int64)\n\nPlot GC-content of given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.possible_proteins-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.possible_proteins","text":"function possible_proteins(s::Sequence)\n\nFind possible proteins of given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.protein_mass","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.protein_mass","text":"function protein_mass(seq::Sequence, type = \"monoisotopic\")\n\nCalculate mass of given amino acid sequence.\n\n\n\n\n\n","category":"function"},{"location":"#Bioinformatics.protparam-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.protparam","text":"function protparam(seq::Sequence)\n\nComputes various physico-chemical properties that can be deduced from a protein sequence.\n\nSee also: https://web.expasy.org/protparam/protparam-doc.html.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.readFASTA-Tuple{String}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.readFASTA","text":"function readFASTA(filename::String)\n\nParse a FASTA formatted file.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.reading_frames-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.reading_frames","text":"function reading_frames(s::Sequence)\n\nFind all reading frames of given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.reverse_complement-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.reverse_complement","text":"function reverse_complement(s::Sequence)\n\nReverse complement of given DNA/RNA sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.skew-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.skew","text":"function skew(seq::Sequence)\n\nCalculate the skew values of sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.skew_plot-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.skew_plot","text":"function skew_plot(seq::Sequence)\n\nPlot skew diagram for given sequence.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.transcription-Tuple{Sequence}","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.transcription","text":"function transcription(dna_seq::Sequence)\n\nTranscribe DNA to RNA.\n\n\n\n\n\n","category":"method"},{"location":"#Bioinformatics.translation","page":"Bioinformatics.jl Documentation","title":"Bioinformatics.translation","text":"function translation(s::Sequence, start_pos::Int64 = 1)\n\nTranslate given DNA sequence to Amino Acid sequence.\n\n\n\n\n\n","category":"function"}]
}