"""
    hammingDist(s::String, t::String)

    Calculates Hamming distance between two strings.
"""
function hammingDist(s::String, t::String)
    @assert length(s) == length(t)
    dist = 0
    for i = 1:length(s)
        if s[i] != t[i]
            dist += 1
        end
    end
    return dist
end

"""
    editDist(s::String, t::String)

    Edit distance is a way of quantifying how dissimilar two strings are
    to one another by counting the minimum number of operations required
    to transform one string into the other.[^1]

[^1]: https://www.wikiwand.com/en/Edit_distance
"""
function editDist(s::String, t::String)
    m = length(s)
    n = length(t)
    d = zeros(m+1, n+1)
    for i = 1:m+1
        d[i, 1] = i - 1
    end
    for j = 1:n+1
        d[1, j] = j - 1
    end
    for j = 2:n+1
        for i = 2:m+1
            if s[i-1] == t[j-1]
                d[i, j] = d[i-1, j-1]
            else
                d[i, j] = min(d[i-1, j] + 1,
                              d[i, j-1] + 1,
                              d[i-1, j-1] + 1)
            end
        end
    end
    return d[m, n]
end
