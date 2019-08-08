"""
    function globalAlign(str1::String, str2::String, match::Int64=1,
                         mismatch::Int64=-1, gap::Int64=-1)

    Applies global alignment to given two strings.
    (Needleman-Wunsch algorithm)
"""
function globalAlign(str1::String, str2::String,
                     match::Int64=1, mismatch::Int64=-1, gap::Int64=-1)

    function w(x, y)
        return (x == y) ? match : mismatch
    end

    m = length(str1)
    n = length(str2)
    D = zeros(m+1, n+1)

    for i = 1:m+1
        D[i,1] = gap*(i-1)
    end

    for j = 1:n+1
        D[1,j] = gap*(j-1)
    end

    for i = 2:m+1
        for j = 2:n+1
            D[i,j] = max(D[i-1, j-1] + w(str1[i-1], str2[j-1]),
                         D[i-1, j] + gap,
                         D[i, j-1] + gap)
        end
    end

    align1 = ""
    align2 = ""

    i = m+1
    j = n+1

    while (i > 1 && j > 1)
        score_current = D[i, j]
        score_diag    = D[i-1, j-1]
        score_up      = D[i, j-1]
        score_left    = D[i-1, j]
        if score_current == score_diag + w(str1[i-1], str2[j-1])
            align1 = string(str1[i-1], align1)
            align2 = string(str2[j-1], align2)
            i -= 1
            j -= 1
        elseif score_current == (score_left + gap)
            align1 = string(str1[i-1], align1)
            align2 = string("-", align2)
            i -= 1
        else
            align1 = string("-", align1)
            align2 = string(str2[j-1], align2)
            j -= 1
        end
    end

    while i > 1
        align1 = string(str1[i-1], align1)
        align2 = string("-", align2)
        i -= 1
    end

    while j > 1
        align1 = string("-", align1)
        align2 = string(str2[j-1], align2)
        j -= 1
    end

    return align1, align2

end

"""
    function localAlign(str1::String, str2::String, match::Int64=3,
                        mismatch::Int64=-3, gap::Int64=-2)

    Applies local alignment to given two strings.
    (Smith-Waterman algorithm)
"""
function localAlign(str1::String, str2::String, match::Int64=3,
                    mismatch::Int64=-3, gap::Int64=-2)

    function w(x, y)
        return (x == y) ? match : mismatch
    end

    m = length(str1)
    n = length(str2)
    D = zeros(m+1, n+1)

    for i = 1:m+1
        D[i,1] = 0
    end

    for j = 1:n+1
        D[1,j] = 0
    end

    for i = 2:m+1
        for j = 2:n+1
            D[i,j] = max(D[i-1, j-1] + w(str1[i-1], str2[j-1]),
                         D[i-1, j] + gap,
                         D[i, j-1] + gap,
                         0)
        end
    end

    align1 = ""
    align2 = ""

    i, j = Tuple(argmax(D))

    while (D[i, j] != 0)
        score_current = D[i, j]
        score_diag    = D[i-1, j-1]
        score_up      = D[i, j-1]
        score_left    = D[i-1, j]
        if score_current == score_diag + w(str1[i-1], str2[j-1])
            align1 = string(str1[i-1], align1)
            align2 = string(str2[j-1], align2)
            i -= 1
            j -= 1
        elseif score_current == (score_left + gap)
            align1 = string(str1[i-1], align1)
            align2 = string("-", align2)
            i -= 1
        else
            align1 = string("-", align1)
            align2 = string(str2[j-1], align2)
            j -= 1
        end
    end

    return align1, align2

end

"""
    prettyPrint(str1::String, str2::String)

    Prints two aligned strings to console.
"""
function prettyPrint(str1::String, str2::String)
    alignment = ""
    for i = 1:length(str1)
        if str1[i] == str2[i]
            alignment = string(alignment, "|")
        else
            alignment = string(alignment, " ")
        end
    end
    println("Alignment: ")
    println(str1)
    println(alignment)
    println(str2)
end
