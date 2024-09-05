function s_print(s)
    for i=1,#s do
        io.write(s[i], " ")
    end
    print()
end

function readData(c, rnd) 
    local filename = "../distance_matrix"
    local fh = assert(io.open(filename, "rb"))
    local contents = assert(fh:read("*a"))
    fh:close()

    local i = 0
    local index = 1
    local dimension = 0
    --local c = {}
    local lines = {}

    for l in string.gmatch(contents, "[^\n]+") do
        lines[#lines + 1] = l
    end

    --print(lines[1])
    dimension = tonumber(lines[index])
    --print(dimension)
    index = index + 1
    for i=1,dimension do
        c[i] = {}
        for j=1,dimension do
            c[i][j] = 0
        end
    end


    for i=1, dimension-1 do
        local j = i + 1
        for cost in string.gmatch(lines[index], "[^%s]+") do
            --io.write(cost, " ");
            if cost ~= "\n" then
                --io.write("\"", cost, "\" ")
                c[i][j] = tonumber(cost)
                c[j][i] = tonumber(cost)
                j = j + 1
            end
            --print()
        end
        index = index + 1
    end

    index = index + 3
    
    local rnd_size = tonumber(lines[index])
    --print(lines[index], #rnd)

    local final = 0;
    for i = 1, rnd_size do
        rnd[#rnd+1] = lines[index+i]
        final = i
    end
    --print(rnd[rnd_size])

    --[[
    for index, l in ipairs(lines) do
        --print(lines[index]);
        print(index);
        if i == 0 then
            dimension = tonumber(l)
            --print("dimension", dimension) 
            for i=1,dimension do
                c[i] = {}
                for j=1,dimension do
                    c[i][j] = 0
                    --io.write(c[i][j], " ")
                end
                --print()
            end
        elseif  l ~= "EOF" then
            local j = i + 1
            for cost in string.gmatch(l, "[^%s]+") do
                --io.write(cost, " ");
                if cost ~= "\n" then
                    --io.write("\"", cost, "\" ")
                    c[i][j] = tonumber(cost)
                    c[j][i] = tonumber(cost)
                    j = j + 1
                end
                --print()
            end

        end
        i = i + 1
    end
    ]]--

    --[[
    for i=1,dimension do
        for j=1,dimension do
            io.write(c[i][j], " ")
        end
        print()
    end
    --]]

    return dimension, rnd_size
end
