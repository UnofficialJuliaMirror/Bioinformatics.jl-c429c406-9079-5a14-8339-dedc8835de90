"""
    readStringFromFile(filename::String)

    Loads a text file that contains a DNA/RNA string.
"""
function readStringFromFile(filename::String)
    s = open(filename) do f
        read(f, String)
    end
    return String(strip(s))
end

"""
    readFASTA(filename::String)

    Reads a FASTA formatted file.
"""
function readFASTA(filename::String)
    data = Dict()
    lines = open(filename) do f
        readlines(f)
    end
    key = ""
    val = ""
    for i = 1:length(lines)
        line = lines[i]
        if startswith(line, '>')
            key = match(r"([^>])+", line).match
        else
            val = string(val, line)
            if i == length(lines)
                data[key] = val
                continue
            end
            nextLine = lines[i+1]
            if startswith(nextLine, '>')
                data[key] = val
                key = ""
                val = ""
            end
        end
    end
    return data
end
