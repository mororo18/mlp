function readData(c) 
    local filename = "../distance_matrix"
    local fh = assert(io.open(filename, "rb"))
    local contents = assert(fh:read("*a"))
    fh:close()

    local i = 0
    --local c = {}
    local dimension = 0
    local lines = {}

    for l in string.gmatch(contents, "[^\n]+") do
        lines[#lines + 1] = l
    end

    for index, l in ipairs(lines) do
        --print(lines[index]);
        --print(l);
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

    --[[
    for i=1,dimension do
        for j=1,dimension do
            io.write(c[i][j], " ")
        end
        print()
    end
    --]]

    return dimension
end
