dimension = nothing
matrix = nothing

function get_instance_info()

    f = open("../distance_matrix")

    line = readline(f)
    dimension = parse(Int, line)
    #matrix = Matrix{Float64}(undef, dimension,dimension)
    matrix = zeros(dimension, dimension)

    """
    matrix = Array{Array{Float64, 1},1}()
    for i in 1:dimension
        l = Array{Float64, 1}()
        for j in 1:dimension
            push!(l, 0.0)
        end
        #println(l)

        push!(matrix, l)
    end
    """

    name = false
    inst = ""

    i = 1
    for line in eachline(f)
        #println(line)
        if name == true
            println(line)
            inst = line
            break
        end

        if line == "EOF"
            #break
            name = true
            continue
        end

        j = i+1
        while length(line) > 0
            index = findfirst(isequal(' '), line)
            d = parse(Float64, line[begin:index])
            #setindex!(matrix, d, i, j)
            matrix[i,j] = d
            matrix[j,i] = d

            line = line[index+1:end]

            j += 1
        end

        i += 1
    end

    """
    println("opa")
    for i in 1:dimension
        for j in 1:dimension
            print(matrix[i,j])
            print(" ")
        end
        println()
    end
    """

    close(f)

    return dimension, matrix, inst

end

#get_instance_info()
