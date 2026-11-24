-- readData.lua
-- Translated from readData.cpp (TSPLIB-format instance reader).
-- Note: as in the original C++ project, this function is not called from
-- main; the program's actual input path is Data.lua's loadData(), which
-- reads "../distance_matrix".

local function calcDistEuc(X, Y, I, J)
    return math.sqrt((X[I] - X[J]) ^ 2 + (Y[I] - Y[J]) ^ 2)
end

local function calcDistAtt(X, Y, I, J)
    local rij = math.sqrt(((X[I] - X[J]) ^ 2 + (Y[I] - Y[J]) ^ 2) / 10)
    local tij = math.floor(rij + 0.5)
    local dij
    if tij < rij then
        dij = tij + 1
    else
        dij = tij
    end
    return dij
end

local function calcLatLong(X, Y, n)
    local PI = 3.141592
    local latit, longit = {}, {}
    for i = 1, n do
        local deg = math.floor(X[i])
        local min = X[i] - deg
        latit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0
    end
    for i = 1, n do
        local deg = math.floor(Y[i])
        local min = Y[i] - deg
        longit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0
    end
    return latit, longit
end

local function calcDistGeo(latit, longit, I, J)
    local RRR = 6378.388
    local q1 = math.cos(longit[I] - longit[J])
    local q2 = math.cos(latit[I] - latit[J])
    local q3 = math.cos(latit[I] + latit[J])
    return math.floor(RRR * math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0)
end

-- Whitespace tokenizer mimicking C++ ifstream's `>>` operator.
local function makeReader(filename)
    local f = io.open(filename, "r")
    if not f then return nil end
    local content = f:read("*a")
    f:close()

    local tokens = {}
    for tok in content:gmatch("%S+") do
        tokens[#tokens + 1] = tok
    end

    local pos = 0
    local function nextToken()
        pos = pos + 1
        return tokens[pos]
    end
    local function nextNumber()
        return tonumber(nextToken())
    end
    return nextToken, nextNumber
end

local function readData(argv)
    if #argv < 1 then
        print("\nFaltando parametros")
        print(" ./exec [Instancia] ")
        os.exit(1)
    end

    if #argv > 1 then
        print("\nMuitos parametros")
        print(" ./exec [Instancia] ")
        os.exit(1)
    end

    local instancia = argv[1]
    local nextToken, nextNumber = makeReader(instancia)
    if not nextToken then
        print("arquivo nao pode ser aberto")
        os.exit(1)
    end

    local arquivo = ""
    while arquivo ~= "DIMENSION:" and arquivo ~= "DIMENSION" do
        arquivo = nextToken()
    end
    if arquivo == "DIMENSION" then arquivo = nextToken() end

    local N = nextNumber()

    while arquivo ~= "EDGE_WEIGHT_TYPE:" and arquivo ~= "EDGE_WEIGHT_TYPE" do
        arquivo = nextToken()
    end
    if arquivo == "EDGE_WEIGHT_TYPE" then arquivo = nextToken() end

    local ewt = nextToken()

    local x, y = {}, {}
    local dist = {}
    for i = 0, N do
        dist[i] = {}
    end

    if ewt == "EXPLICIT" then
        while arquivo ~= "EDGE_WEIGHT_FORMAT:" and arquivo ~= "EDGE_WEIGHT_FORMAT" do
            arquivo = nextToken()
        end
        if arquivo == "EDGE_WEIGHT_FORMAT" then arquivo = nextToken() end
        local ewf = nextToken()

        if ewf == "FUNCTION" then
            print("FUNCTION - Nao suportado!")

        elseif ewf == "FULL_MATRIX" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for i = 1, N do
                for j = 1, N do
                    dist[i][j] = nextNumber()
                end
            end

        elseif ewf == "UPPER_ROW" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for i = 1, N - 1 do
                for j = i + 1, N do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end
            for i = 1, N do dist[i][i] = 0 end

        elseif ewf == "LOWER_ROW" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for i = 2, N do
                for j = 1, i - 1 do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end
            for i = 1, N do dist[i][i] = 0 end

        elseif ewf == "UPPER_DIAG_ROW" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for i = 1, N do
                for j = i, N do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end

        elseif ewf == "LOWER_DIAG_ROW" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for i = 1, N do
                for j = 1, i do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end

        elseif ewf == "UPPER_COL" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for j = 2, N do
                for i = 1, j - 1 do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end
            for i = 1, N do dist[i][i] = 0 end

        elseif ewf == "LOWER_COL" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            -- Mirrors the original C++ loop verbatim: it increments `j`
            -- instead of `i` in the inner loop, so `i` never advances.
            -- This branch is unreachable in practice (readData() is never
            -- invoked by main.lua) so the bug is preserved rather than fixed.
            local j = 1
            while j < N do
                local i = j + 1
                while i < N + 1 do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                    j = j + 1
                end
            end
            for i = 1, N do dist[i][i] = 0 end

        elseif ewf == "UPPER_DIAG_COL" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for j = 1, N do
                for i = 1, j do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end

        elseif ewf == "LOWER_DIAG_COL" then
            while arquivo ~= "EDGE_WEIGHT_SECTION" do arquivo = nextToken() end
            for j = 1, N do
                for i = j, N do
                    dist[i][j] = nextNumber()
                    dist[j][i] = dist[i][j]
                end
            end
        end

    elseif ewt == "EUC_2D" then
        while arquivo ~= "NODE_COORD_SECTION" do arquivo = nextToken() end
        for i = 1, N do
            local _tempCity = nextNumber()
            x[i] = nextNumber()
            y[i] = nextNumber()
        end
        for i = 1, N do
            for j = 1, N do
                dist[i][j] = math.floor(calcDistEuc(x, y, i, j) + 0.5)
            end
        end

    elseif ewt == "EUD_3D" then
        print("EUC_3D - Nao suportado!")
    elseif ewt == "MAX_2D" then
        print("MAX_2D - Nao suportado!")
    elseif ewt == "MAX_3D" then
        print("MAX_3D - Nao Suportado!")
    elseif ewt == "MAN_2D" then
        print("MAN_2D - Nao suportado!")
    elseif ewt == "MAN_3D" then
        print("MAN_3D - Nao Suportado!")

    elseif ewt == "CEIL_2D" then
        while arquivo ~= "NODE_COORD_SECTION" do arquivo = nextToken() end
        for i = 1, N do
            local _tempCity = nextNumber()
            x[i] = nextNumber()
            y[i] = nextNumber()
        end
        for i = 1, N do
            for j = 1, N do
                dist[i][j] = math.ceil(calcDistEuc(x, y, i, j))
            end
        end

    elseif ewt == "GEO" then
        while arquivo ~= "NODE_COORD_SECTION" do arquivo = nextToken() end
        for i = 1, N do
            local _tempCity = nextNumber()
            x[i] = nextNumber()
            y[i] = nextNumber()
        end
        local latitude, longitude = calcLatLong(x, y, N)
        for i = 1, N do
            for j = 1, N do
                dist[i][j] = calcDistGeo(latitude, longitude, i, j)
            end
        end

    elseif ewt == "ATT" then
        while arquivo ~= "NODE_COORD_SECTION" do arquivo = nextToken() end
        for i = 1, N do
            local _tempCity = nextNumber()
            x[i] = nextNumber()
            y[i] = nextNumber()
        end
        for i = 1, N do
            for j = 1, N do
                dist[i][j] = calcDistAtt(x, y, i, j)
            end
        end

    elseif ewt == "XRAY1" then
        print("XRAY1 - Nao suportado!")
    elseif ewt == "XRAY2" then
        print("XRAY2 - Nao suportado!")
    elseif ewt == "SPECIAL" then
        print("SPECIAL - Nao suportado!")
    end

    return N, dist
end

return { readData = readData }
