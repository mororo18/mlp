-- Classic random-access comparison across different data-structure shapes
-- for the same logical data (T/C/W triple per (i,j) pair, i<=j):
--   A) flat 1D, square index (current code's actual layout)
--   B) flat 1D, packed triangular index (no gaps)
--   C) nested tables t[i][j][k]  (true 3-level sub-tables, pre-flattening style)
--   D) "array of structs": t[i][j] = {T=,C=,W=} (2-level nesting + named leaf table)
--   E) struct of arrays: three separate flat 1D packed tables (T[], C[], W[])
--
-- Random (shuffled) access order defeats hardware prefetch and exposes
-- the real cost of pointer-chasing through nested tables vs a single
-- flat index -- a much harsher, more decisive test than sequential scan.

local clock = os.clock
local DIM = 300 -- smaller than the sequential test: nested-table builds are heavier

local function row_off_square(i) return DIM * (i - 1) end

local function build_flat_square()
    local t = {}
    for i = 1, DIM do
        local ro = row_off_square(i)
        for j = i, DIM do
            local base = 3 * (ro + (j - 1))
            t[base + 1] = (i + j + 1) * 0.001
            t[base + 2] = (i + j + 2) * 0.001
            t[base + 3] = (i + j + 3) * 0.001
        end
    end
    return t
end

local function build_flat_packed()
    local t = {}
    local row_off = {}
    local acc = 0
    for i = 1, DIM do
        row_off[i] = acc
        acc = acc + (DIM - i + 1)
    end
    for i = 1, DIM do
        local ro = row_off[i]
        for j = i, DIM do
            local base = 3 * (ro + (j - i))
            t[base + 1] = (i + j + 1) * 0.001
            t[base + 2] = (i + j + 2) * 0.001
            t[base + 3] = (i + j + 3) * 0.001
        end
    end
    return t, row_off
end

local function build_nested_ijk()
    local t = {}
    for i = 1, DIM do
        t[i] = {}
        for j = i, DIM do
            t[i][j] = { (i + j + 1) * 0.001, (i + j + 2) * 0.001, (i + j + 3) * 0.001 }
        end
    end
    return t
end

local function build_array_of_structs()
    local t = {}
    for i = 1, DIM do
        t[i] = {}
        for j = i, DIM do
            t[i][j] = { T = (i + j + 1) * 0.001, C = (i + j + 2) * 0.001, W = (i + j + 3) * 0.001 }
        end
    end
    return t
end

local function build_soa()
    local Tt, Ct, Wt = {}, {}, {}
    local row_off = {}
    local acc = 0
    for i = 1, DIM do
        row_off[i] = acc
        acc = acc + (DIM - i + 1)
    end
    for i = 1, DIM do
        local ro = row_off[i]
        for j = i, DIM do
            local idx = ro + (j - i) + 1
            Tt[idx] = (i + j + 1) * 0.001
            Ct[idx] = (i + j + 2) * 0.001
            Wt[idx] = (i + j + 3) * 0.001
        end
    end
    return Tt, Ct, Wt, row_off
end

print("Building structures...")
collectgarbage("collect")
local flatSq = build_flat_square()
collectgarbage("collect")
local flatPk, pkOff = build_flat_packed()
collectgarbage("collect")
local nested = build_nested_ijk()
collectgarbage("collect")
local aos = build_array_of_structs()
collectgarbage("collect")
local soaT, soaC, soaW, soaOff = build_soa()
collectgarbage("collect")

-- Common random order of valid (i,j) pairs, shared across all structures
-- so every variant does the exact same logical accesses.
local pairs_ij = {}
for i = 1, DIM do
    for j = i, DIM do
        pairs_ij[#pairs_ij + 1] = { i, j }
    end
end
math.randomseed(1234)
for i = #pairs_ij, 2, -1 do
    local r = math.random(i)
    pairs_ij[i], pairs_ij[r] = pairs_ij[r], pairs_ij[i]
end
print(("DIM=%d, %d (i,j) pairs, random access order shared across all variants"):format(DIM, #pairs_ij))

local REPS = 8

local function time_it(label, fn)
    collectgarbage("collect")
    local t0 = clock()
    local sum = 0
    for _ = 1, REPS do
        sum = sum + fn()
    end
    local t1 = clock()
    print(("%-38s %.4f s  (sum=%.3f)"):format(label, t1 - t0, sum))
end

print()
print(("=== RANDOM access, %d reps ==="):format(REPS))

time_it("A) flat 1D square (current code)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (row_off_square(i) + (j - 1))
        sum = sum + flatSq[base + 1] + flatSq[base + 2] + flatSq[base + 3]
    end
    return sum
end)

time_it("B) flat 1D packed triangular", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (pkOff[i] + (j - i))
        sum = sum + flatPk[base + 1] + flatPk[base + 2] + flatPk[base + 3]
    end
    return sum
end)

time_it("C) nested t[i][j][k]", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local cell = nested[i][j]
        sum = sum + cell[1] + cell[2] + cell[3]
    end
    return sum
end)

time_it("D) array of structs t[i][j]={T=,C=,W=}", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local cell = aos[i][j]
        sum = sum + cell.T + cell.C + cell.W
    end
    return sum
end)

time_it("E) struct of arrays (3 flat packed)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local idx = soaOff[i] + (j - i) + 1
        sum = sum + soaT[idx] + soaC[idx] + soaW[idx]
    end
    return sum
end)

print()
print("=== SEQUENTIAL access (same pairs, original ascending order), for comparison ===")

table.sort(pairs_ij, function(a, b)
    if a[1] ~= b[1] then return a[1] < b[1] end
    return a[2] < b[2]
end)

time_it("A) flat 1D square (current code)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (row_off_square(i) + (j - 1))
        sum = sum + flatSq[base + 1] + flatSq[base + 2] + flatSq[base + 3]
    end
    return sum
end)

time_it("B) flat 1D packed triangular", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (pkOff[i] + (j - i))
        sum = sum + flatPk[base + 1] + flatPk[base + 2] + flatPk[base + 3]
    end
    return sum
end)

time_it("C) nested t[i][j][k]", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local cell = nested[i][j]
        sum = sum + cell[1] + cell[2] + cell[3]
    end
    return sum
end)

time_it("D) array of structs t[i][j]={T=,C=,W=}", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local cell = aos[i][j]
        sum = sum + cell.T + cell.C + cell.W
    end
    return sum
end)

time_it("E) struct of arrays (3 flat packed)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local idx = soaOff[i] + (j - i) + 1
        sum = sum + soaT[idx] + soaC[idx] + soaW[idx]
    end
    return sum
end)
