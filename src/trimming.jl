

"""
    Read a block of fastq reads into Vector{String}[]
"""
function read_block!(io, data::Vector{Vector{String}}, blocksize=length(data))
    r = 0
    for i = 1:blocksize
        eof(io) && break
        data[i][1] = readline(io)
        data[i][2] = readline(io)
        data[i][3] = readline(io)
        data[i][4] = readline(io)
        r += 1
    end
    r
end


"""
    Write a block of fastq reads from Vector{String}[]
    Writes trimmed sequence and quality string
"""
function write_block(io, data, rn, trim_data)
    for i = 1:rn
        ind = 1:trim_data[2, i]
        write(io, data[i][1], "\n")
        write(io, data[i][2][ind], "\n")
        write(io, data[i][3], "\n")
        write(io, data[i][4][ind], "\n")
    end
end


"""
    Count the number of Ambiguous bases in a sequence
    i.e. count number of N's
"""
count_n(seq) = count(isambiguous, seq)
