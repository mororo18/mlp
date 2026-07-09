-- Larger-scale, factorial follow-up to table_structs_random_test.lua.
-- That test showed nested tables (t[i][j][k]) and per-entry named
-- structs (t[i][j]={T=,C=,W=}) clearly lose to any flat representation,
-- so this run drops them and focuses on the flat survivors, crossed as
-- a 2x2 factorial to separate the two independent questions:
--   layout:  interleaved (current code, *3+k)      vs  SoA (3 arrays)
--   packing: square (has triangular gaps)          vs  packed (no gaps)
--
-- DIM is picked so the data clearly exceeds L3 (6 MiB on this machine,
-- `lscpu` checked) -- at DIM=1500 the square layout spans 3*1500*1500 =
-- 6.75M slots, tens of MB, well past L1/L2/L3, so any real cache-locality
-- difference between layouts should show up in wall time instead of
-- being hidden by everything fitting in cache.

local clock = os.clock
local unpack = table.unpack or unpack
local DIM = 1500

local function mem_kb()
    collectgarbage("collect")
    collectgarbage("collect")
    return collectgarbage("count")
end

-- ---- builders ----

local function build_interleaved_square()
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

local function build_interleaved_packed()
    local t, row_off = {}, {}
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

local function build_soa_square()
    local Tt, Ct, Wt = {}, {}, {}
    for i = 1, DIM do
        local ro = DIM * (i - 1)
        for j = i, DIM do
            local idx = ro + (j - 1) + 1
            Tt[idx] = (i + j + 1) * 0.001
            Ct[idx] = (i + j + 2) * 0.001
            Wt[idx] = (i + j + 3) * 0.001
        end
    end
    return Tt, Ct, Wt
end

local function build_soa_packed()
    local Tt, Ct, Wt, row_off = {}, {}, {}, {}
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

print(("DIM=%d -> square range=%d slots (~%.1f MB @16B/slot), packed range=%d slots (~%.1f MB)")
    :format(DIM, 3 * DIM * DIM, 3 * DIM * DIM * 16 / 1e6,
        3 * (DIM * (DIM + 1) / 2), 3 * (DIM * (DIM + 1) / 2) * 16 / 1e6))

print()
print("=== Build time + memory footprint ===")

local function build_and_report(label, fn)
    collectgarbage("collect")
    local base = mem_kb()
    local t0 = clock()
    local results = { fn() }
    local t1 = clock()
    local m1 = mem_kb()
    print(("%-28s build=%.4f s  +%.1f KB"):format(label, t1 - t0, m1 - base))
    return unpack(results)
end

local intSq = build_and_report("interleaved+square (current)", build_interleaved_square)
local intPk, intPkOff = build_and_report("interleaved+packed", build_interleaved_packed)
local soaSqT, soaSqC, soaSqW = build_and_report("SoA+square", build_soa_square)
local soaPkT, soaPkC, soaPkW, soaPkOff = build_and_report("SoA+packed", build_soa_packed)

-- Shared random order of valid (i,j) pairs
local pairs_ij = {}
for i = 1, DIM do
    for j = i, DIM do
        pairs_ij[#pairs_ij + 1] = { i, j }
    end
end
math.randomseed(777)
for i = #pairs_ij, 2, -1 do
    local r = math.random(i)
    pairs_ij[i], pairs_ij[r] = pairs_ij[r], pairs_ij[i]
end
print(("\n%d (i,j) pairs total"):format(#pairs_ij))

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
print(("=== FULL RANDOM access (all %d pairs, shuffled), %d reps ==="):format(#pairs_ij, REPS))

time_it("interleaved+square (current)", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (DIM * (i - 1) + (j - 1))
        sum = sum + intSq[base + 1] + intSq[base + 2] + intSq[base + 3]
    end
    return sum
end)

time_it("interleaved+packed", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local base = 3 * (intPkOff[i] + (j - i))
        sum = sum + intPk[base + 1] + intPk[base + 2] + intPk[base + 3]
    end
    return sum
end)

time_it("SoA+square", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local idx = DIM * (i - 1) + (j - 1) + 1
        sum = sum + soaSqT[idx] + soaSqC[idx] + soaSqW[idx]
    end
    return sum
end)

time_it("SoA+packed", function()
    local sum = 0
    for p = 1, #pairs_ij do
        local i, j = pairs_ij[p][1], pairs_ij[p][2]
        local idx = soaPkOff[i] + (j - i) + 1
        sum = sum + soaPkT[idx] + soaPkC[idx] + soaPkW[idx]
    end
    return sum
end)

-- Windowed/neighborhood access: closer to what search_reinsertion's
-- profiled hot lines actually do (fixed outer i, scanning a bounded
-- range of j/k around it, repeated for many different outer positions)
-- rather than one giant uniform shuffle over the entire matrix.
print()
print("=== WINDOWED neighborhood access (mimics search_reinsertion) ===")
local WINDOW = 60
local OUTER_STEP = 11
local WREPS = 60

local function windowed_interleaved(t, offFn)
    local sum = 0
    for i = 1, DIM - WINDOW, OUTER_STEP do
        local ro = offFn(i)
        for j = i, i + WINDOW do
            local base = 3 * (ro + (j - i))
            sum = sum + t[base + 1] + t[base + 2] + t[base + 3]
        end
    end
    return sum
end

local function windowed_soa(Tt, Ct, Wt, offFn)
    local sum = 0
    for i = 1, DIM - WINDOW, OUTER_STEP do
        local ro = offFn(i)
        for j = i, i + WINDOW do
            local idx = ro + (j - i) + 1
            sum = sum + Tt[idx] + Ct[idx] + Wt[idx]
        end
    end
    return sum
end

local function windowed_square_explicit(t)
    local sum = 0
    for i = 1, DIM - WINDOW, OUTER_STEP do
        local ro = DIM * (i - 1)
        for j = i, i + WINDOW do
            local base = 3 * (ro + (j - 1))
            sum = sum + t[base + 1] + t[base + 2] + t[base + 3]
        end
    end
    return sum
end
time_it("interleaved+square (current)", function()
    local sum = 0
    for _ = 1, WREPS do sum = sum + windowed_square_explicit(intSq) end
    return sum
end)

time_it("interleaved+packed", function()
    local sum = 0
    for _ = 1, WREPS do
        sum = sum + windowed_interleaved(intPk, function(i) return intPkOff[i] end)
    end
    return sum
end)

local function windowed_soa_square_explicit(Tt, Ct, Wt)
    local sum = 0
    for i = 1, DIM - WINDOW, OUTER_STEP do
        local ro = DIM * (i - 1)
        for j = i, i + WINDOW do
            local idx = ro + (j - 1) + 1
            sum = sum + Tt[idx] + Ct[idx] + Wt[idx]
        end
    end
    return sum
end
time_it("SoA+square", function()
    local sum = 0
    for _ = 1, WREPS do sum = sum + windowed_soa_square_explicit(soaSqT, soaSqC, soaSqW) end
    return sum
end)

time_it("SoA+packed", function()
    local sum = 0
    for _ = 1, WREPS do
        sum = sum + windowed_soa(soaPkT, soaPkC, soaPkW, function(i) return soaPkOff[i] end)
    end
    return sum
end)
