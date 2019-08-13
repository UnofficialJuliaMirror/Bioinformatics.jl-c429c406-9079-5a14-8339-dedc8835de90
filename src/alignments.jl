"""
    function dotmatrix(s1::Sequence, s2::Sequence)

Calculate the dotplot matrix for given two sequences.
"""
function dotmatrix(s1::Sequence, s2::Sequence)
    mat = zeros(Int8, (length(s1), length(s2)))
    for i in 1:length(s1)
        for j in 1:length(s2)
            if s1[i] == s2[j]
                mat[i, j] = 1
            end
        end
    end
    return mat
end
