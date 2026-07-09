-- LuaJIT-only (needs `ffi`, not available in lua5.3). Tests whether the
-- LuaJIT Numerical Computing Performance Guide's headline claim (FFI
-- numeric arrays much faster than Lua tables) holds up against our
-- ALREADY-flat table layout (not a naive nested-table strawman), at the
-- same DIM=1500 cache-busting scale and access patterns as bench_04, so
-- results are directly comparable.
local ok, ffi = pcall(require, "ffi")
if not ok then
    print("ffi not available (expected on lua5.3, this test is luajit-only): " .. tostring(ffi))
    os.exit(0)
end

local clock = os.clock
local DIM = 1500

local function mem_kb()
    collectgarbage("collect")
    collectgarbage("collect")
    return collectgarbage("count")
end

local function build_table_square()
    local t = {}
    for i = 1, DIM do
        local ro = DIM * (i - 1)
        for j = i, DIM do
            local base = 3 * (ro + (j - 1))
            t[base + 1] = (i + j + 1) * 0.001
            t[base + 2] = (i + j + 2) * 0.001
            t[base + 3] = (i + j + 3) * 0.001
        end
    end
    return t
end

local function build_ffi_square()
    local n = 3 * DIM * DIM
    local arr = ffi.new("double[?]", n) -- zero-initialized by ffi.new
    for i = 1, DIM do
        local ro = DIM * (i - 1)
        for j = i, DIM do
            local base = 3 * (ro + (j - 1))
            arr[base + 0] = (i + j + 1) * 0.001 -- FFI arrays are 0-indexed
            arr[base + 1] = (i + j + 2) * 0.001
            arr[base + 2] = (i + j + 3) * 0.001
        end
    end
    return arr
end

print(("DIM=%d, %d slots (~%.1f MB as double[])"):format(DIM, 3 * DIM * DIM, 3 * DIM * DIM * 8 / 1e6))

print()
print("=== Build time + memory footprint ===")
collectgarbage("collect")
local base = mem_kb()
local t0 = clock()
local tbl = build_table_square()
local t1 = clock()
local m1 = mem_kb()
print(("Lua table (current layout): build=%.4f s  +%.1f KB"):format(t1 - t0, m1 - base))

collectgarbage("collect")
base = mem_kb()
t0 = clock()
local arr = build_ffi_square()
t1 = clock()
m1 = mem_kb()
print(("FFI double[] array:         build=%.4f s  +%.1f KB (Lua-heap only; cdata itself is outside the GC heap)"):format(t1 - t0, m1 - base))

-- shared random order of valid (i,j) pairs
local pairs_ij = {}
for i = 1, DIM do
    for j = i, DIM do
        pairs_ij[#pairs_ij + 1] = { i, j }
    end
end
math.randomseed(2026)
for i = #pairs_ij, 2, -1 do
    local r = math.random(i)
    pairs_ij[i], pairs_ij[r] = pairs_ij[r], pairs_ij[i]
end

local REPS = 5
local function time_it(label, fn)
    collectgarbage("collect")
    local t0 = clock()
    local sum = 0
    for _ = 1, REPS do
        sum = sum + fn()
    end
    local t1 = clock()
    print(("%-28s %.4f s  (sum=%.3f)"):format(label, t1 - t0, sum))
end

print()
print(("=== FULL RANDOM access (%d pairs), %d reps ==="):format(#pairs_ij, REPS))

time_it("Lua table (current)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (DIM * (i - 1) + (j - 1))
        sum = sum + tbl[base + 1] + tbl[base + 2] + tbl[base + 3]
    end
    return sum
end)

time_it("FFI double[]", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (DIM * (i - 1) + (j - 1))
        sum = sum + arr[base + 0] + arr[base + 1] + arr[base + 2]
    end
    return sum
end)

-- windowed/neighborhood pattern, mirrors search_reinsertion's real shape
print()
print("=== WINDOWED neighborhood access ===")
local WINDOW = 60
local OUTER_STEP = 11
local WREPS = 60

time_it("Lua table (current)", function()
    local sum = 0
    for _ = 1, WREPS do
        for i = 1, DIM - WINDOW, OUTER_STEP do
            local ro = DIM * (i - 1)
            for j = i, i + WINDOW do
                local base = 3 * (ro + (j - 1))
                sum = sum + tbl[base + 1] + tbl[base + 2] + tbl[base + 3]
            end
        end
    end
    return sum
end)

time_it("FFI double[]", function()
    local sum = 0
    for _ = 1, WREPS do
        for i = 1, DIM - WINDOW, OUTER_STEP do
            local ro = DIM * (i - 1)
            for j = i, i + WINDOW do
                local base = 3 * (ro + (j - 1))
                sum = sum + arr[base + 0] + arr[base + 1] + arr[base + 2]
            end
        end
    end
    return sum
end)
