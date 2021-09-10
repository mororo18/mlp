dofile("Data.lua")

function main() 
    local info = {
        c = {}
    }

    info.dimension = readData(info.c)

    print(info.dimension)
    for i=1,info.dimension do
        for j=1,info.dimension do
            io.write(info.c[i][j], " ")
        end
        print()
    end
end

main()

--[[
--comments
]]--

