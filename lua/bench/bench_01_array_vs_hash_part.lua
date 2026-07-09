-- Empirical probe of Lua table internals: array part vs hash part.
-- Goal: validate (not assume) whether a table filled the way `seq` is
-- filled in main.canc.lua (subseq_fill: triangular i<=j pattern, so the
-- flat index sequence has gaps) ends up array-backed (contiguous) or
-- falls into the hash part, and whether pre-sizing helps.

local clock = os.clock

local function mem_kb()
    collectgarbage("collect")
    collectgarbage("collect")
    return collectgarbage("count")
end

local DIM = 100 -- mirrors dimension+1 for a ~99-city instance (rat99)

-- 1) Dense sequential fill: t[1..N] = 0.0, N = same total slots as the
--    real seq array would span (3 * DIM * DIM), fully populated (no gaps).
local function build_dense()
    local t = {}
    local N = 3 * DIM * DIM
    for i = 1, N do
        t[i] = 0.0
    end
    return t, N
end

-- 2) Triangular fill: exact same loop shape as subseq_fill() in
--    main.canc.lua (i=1..DIM, j=i..DIM, k=1..3), which leaves gaps in the
--    flat index space wherever j < i would have mapped.
local function build_triangular()
    local t = {}
    local count = 0
    for i = 1, DIM do
        for j = i, DIM do
            for k = 1, 3 do
                local idx = 3 * (DIM * (i - 1) + (j - 1)) + k
                t[idx] = 0.0
                count = count + 1
            end
        end
    end
    return t, count
end

-- 3) Same triangular key SET as (2), but inserted in a shuffled order,
--    to see whether insertion order (not just final key set) matters to
--    Lua's array-part promotion heuristic.
local function build_triangular_shuffled()
    local keys = {}
    for i = 1, DIM do
        for j = i, DIM do
            for k = 1, 3 do
                keys[#keys + 1] = 3 * (DIM * (i - 1) + (j - 1)) + k
            end
        end
    end
    -- Fisher-Yates shuffle with a fixed seed for reproducibility
    math.randomseed(42)
    for i = #keys, 2, -1 do
        local j = math.random(i)
        keys[i], keys[j] = keys[j], keys[i]
    end
    local t = {}
    for _, idx in ipairs(keys) do
        t[idx] = 0.0
    end
    return t, #keys, keys
end

print("=== Memory footprint (KB, GC-settled) ===")
collectgarbage("collect")
local base = mem_kb()

local dense, denseN = build_dense()
local m1 = mem_kb()
print(string.format("dense  1..%d (%d slots, fully populated): +%.1f KB (%.4f KB/slot)",
    denseN, denseN, m1 - base, (m1 - base) / denseN))
dense = nil
collectgarbage("collect")

local base2 = mem_kb()
local tri, triN = build_triangular()
local m2 = mem_kb()
print(string.format("triangular (subseq_fill pattern), %d real keys over range 1..%d: +%.1f KB (%.4f KB/slot)",
    triN, 3 * DIM * DIM, m2 - base2, (m2 - base2) / triN))

local base3 = mem_kb()
local trish, trishN, keys = build_triangular_shuffled()
local m3 = mem_kb()
print(string.format("triangular SHUFFLED insertion order, %d keys: +%.1f KB (%.4f KB/slot)",
    trishN, m3 - base3, (m3 - base3) / trishN))

print()
print("=== Timing: fill (write) cost ===")
collectgarbage("collect")
local t0 = clock()
local d2 = build_dense()
local t1 = clock()
print(string.format("dense fill:              %.4f s", t1 - t0))
d2 = nil

collectgarbage("collect")
t0 = clock()
local tr2 = build_triangular()
t1 = clock()
print(string.format("triangular fill (real pattern): %.4f s", t1 - t0))

collectgarbage("collect")
t0 = clock()
local trs2 = build_triangular_shuffled()
t1 = clock()
print(string.format("triangular fill SHUFFLED order: %.4f s", t1 - t0))

print()
print("=== Timing: read cost (sequential order vs shuffled order) over the SAME triangular table ===")
-- Read the triangular table's real keys in ascending order (matches how
-- subseq_load/search_* actually walk it) vs in shuffled order, many
-- repetitions to get a stable signal. If this is a true contiguous
-- array part, ascending order should benefit from cache locality and be
-- faster than shuffled order touching the same key set. If it's hash
-- part, order should barely matter (hash lookups don't care about
-- "logical" adjacency, and are already cache-unfriendly regardless).
local REPS = 200
local ascKeys = {}
for i = 1, DIM do
    for j = i, DIM do
        for k = 1, 3 do
            ascKeys[#ascKeys + 1] = 3 * (DIM * (i - 1) + (j - 1)) + k
        end
    end
end

local sum = 0
t0 = clock()
for _ = 1, REPS do
    for i = 1, #ascKeys do
        sum = sum + tri[ascKeys[i]]
    end
end
t1 = clock()
print(string.format("read ascending order:  %.4f s (sum=%s)", t1 - t0, tostring(sum)))

sum = 0
t0 = clock()
for _ = 1, REPS do
    for i = 1, #keys do
        sum = sum + tri[keys[i]]
    end
end
t1 = clock()
print(string.format("read shuffled order:    %.4f s (sum=%s)", t1 - t0, tostring(sum)))

print()
print("=== table.new (LuaJIT only) pre-sizing test, if available ===")
local ok, tnew = pcall(require, "table.new")
if ok and tnew then
    collectgarbage("collect")
    local base4 = mem_kb()
    t0 = clock()
    local pre = tnew(3 * DIM * DIM, 0)
    for i = 1, DIM do
        for j = i, DIM do
            for k = 1, 3 do
                local idx = 3 * (DIM * (i - 1) + (j - 1)) + k
                pre[idx] = 0.0
            end
        end
    end
    t1 = clock()
    local m4 = mem_kb()
    print(string.format("table.new preallocated + triangular fill: %.4f s, +%.1f KB", t1 - t0, m4 - base4))
else
    print("table.new not available in this interpreter (expected on lua5.3): " .. tostring(tnew))
end
