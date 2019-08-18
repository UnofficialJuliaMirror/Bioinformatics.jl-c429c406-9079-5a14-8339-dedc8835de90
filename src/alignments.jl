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

function score_pos(c1, c2, sub_mat, gap_penalty)
    if c1 == '-' || c2 == '-'
        return gap_penalty
    else
        return sub_mat[(c1, c2)]
    end
end

function score_align(seq1, seq2, sub_mat, gap_penalty)
    res = 0
    for i in 1:(length(seq1))
        res += score_pos(seq1[i], seq2[i], sub_mat, gap_penalty)
    end
    return res
end

function score_affine_gap(seq1, seq2, sub_mat, gap_open_penalty, gap_extend_penalty)
    res = 0
    ingap1 = false
    ingap2 = false
    for i in 1:(length(seq1))
        if seq1[i] == '-'
            if ingap1
                res += gap_extend_penalty
            else
                ingap1 = true
                res += gap_open_penalty
            end
        elseif seq2[i] == '-'
            if ingap2
                res += gap_extend_penalty
            else
                ingap2 = true
                res += gap_open_penalty
            end
        else
            if ingap1
                ingap1 = false
            end
            if ingap2
                ingap2 = false
            end
            res += sub_mat[(seq1[i], seq2[i])]
        end
    end
    return res
end
