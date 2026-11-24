-- Data.lua
-- Translated from Data.hpp

local function matrix_allocate(sz)
    local ptr = {}
    for i = 0, sz - 1 do
        ptr[i] = {}
    end
    return ptr
end

local function loadData()
    local dimension
    local matrix
    local rnd = {}
    local rnd_n = 0

    local file = io.open("../distance_matrix", "r")

    if file then
        local line = file:read("*l")
        dimension = tonumber(line)
        matrix = matrix_allocate(dimension)

        local flag = false
        local rnd_size = -1
        local i = 0

        for text_line in file:lines() do
            local j = i + 1
            local rest = text_line
            while rest:find(" ") do
                local pos = rest:find(" ")
                local value = tonumber(rest:sub(1, pos - 1))
                rest = rest:sub(pos + 1)
                matrix[i][j] = value
                matrix[j][i] = value
                j = j + 1
            end

            i = i + 1

            if flag then
                if rnd_size < 0 then
                    rnd_size = tonumber(rest)
                else
                    rnd[rnd_n] = tonumber(rest)
                    rnd_n = rnd_n + 1
                end
            end

            if rest == "RND" then flag = true end
        end

        file:close()
    end

    return dimension, matrix, rnd, rnd_n
end

local Data = {}
Data.__index = Data

function Data.new()
    return setmetatable({ dimen = nil, cost = nil, rnd = nil, rnd_n = 0 }, Data)
end

function Data:load()
    self.dimen, self.cost, self.rnd, self.rnd_n = loadData()
end

function Data:getRndVec()
    return self.rnd
end

function Data:getDimen()
    return self.dimen
end

function Data:getCostPtr()
    return self.cost
end

return {
    matrix_allocate = matrix_allocate,
    loadData = loadData,
    Data = Data,
}
