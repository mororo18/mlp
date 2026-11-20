#
# Port of Data.hpp -- reads the instance from "../distance_matrix"
#
# The distance matrix is kept as a plain Julia matrix. The algorithm code
# works with 0-based node ids (just like the C++ version), so every access
# adds 1 to the indices.
#

function load_data()
    dimension = 0
    cost = zeros(Float64, 0, 0)
    rnd = Int[]

    open("../distance_matrix", "r") do file
        line = readline(file)
        dimension = parse(Int, line)
        cost = zeros(Float64, dimension, dimension)

        flag = false
        rnd_size = -1
        i = 0
        while !eof(file)
            line = readline(file)
            j = i + 1
            while occursin(" ", line)
                pos = findfirst(' ', line)
                value = parse(Float64, line[1:pos-1])
                line = line[pos+1:end]
                cost[i+1, j+1] = value
                cost[j+1, i+1] = value
                j += 1
            end

            i += 1

            if flag
                if rnd_size < 0
                    rnd_size = parse(Int, line)
                    sizehint!(rnd, rnd_size)
                else
                    push!(rnd, parse(Int, line))
                end
            end

            if line == "RND"
                flag = true
            end
        end
    end

    return dimension, cost, rnd
end
