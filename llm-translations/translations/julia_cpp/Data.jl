dimension = nothing
matrix = nothing

function get_instance_info()

    f = open("../distance_matrix")

    line = readline(f)
    dimension = parse(Int, line)
    matrix = zeros(dimension, dimension)

    rand_values = Array{Int, 1}(undef, 0)
    rnd = false
    sz = -1
    name = false
    inst = ""

    i = 1
    for line in eachline(f)
        if sz > 0 
            push!(rand_values, parse(Int, line))
            continue
        end

        if rnd == true
            sz = parse(Int, line)
            continue
        end

        if line == "RND"
            rnd = true
            continue
        end

        if name == true
            inst = line
            name = false
            continue
        end

        if line == "EOF"
            name = true
            continue
        end

        j = i+1
        while length(line) > 0
            index = findfirst(isequal(' '), line)
            d = parse(Float64, line[begin:index])
            matrix[i,j] = d
            matrix[j,i] = d

            line = line[index+1:end]

            j += 1
        end

        i += 1
    end

    close(f)

    return dimension, matrix, rand_values

end
