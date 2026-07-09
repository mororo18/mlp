-- v2: isolate memory-layout/cache effect from index-arithmetic cost.
-- v1 was confounded: idx_packed recomputed a division-heavy triangular
-- offset on every single call, so the "packed" slowdown could just be
-- the index math, not cache locality. Here both schemes precompute a
-- per-row offset ONCE per outer i (row_off[i]), so the inner-loop index
-- arithmetic is the same shape/cost (one mul, one add, one sub) for
-- both -- isolating whatever's left to genuine memory layout.

local clock = os.clock
local DIM = 700

-- current scheme: square, offset per row is just DIM*(i-1), no division
local function build_square()
    local t = {}
    local row_off = {}
    for i = 1, DIM do
        row_off[i] = DIM * (i - 1)
    end
    for i = 1, DIM do
        local ro = row_off[i]
        for j = i, DIM do
            local base = 3 * (ro + (j - 1))
            t[base + 1] = (i + j + 1) * 0.001
            t[base + 2] = (i + j + 2) * 0.001
            t[base + 3] = (i + j + 3) * 0.001
        end
    end
    return t, row_off
end

-- packed scheme: real triangular offset, precomputed ONCE per row
-- (division happens DIM times total at setup, not per access)
local function build_packed()
    local t = {}
    local row_off = {}
    local acc = 0
    for i = 1, DIM do
        row_off[i] = acc
        acc = acc + (DIM - i + 1) -- number of j values for this i (j=i..DIM)
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

collectgarbage("collect")
local sq, sq_off = build_square()
collectgarbage("collect")
local pk, pk_off = build_packed()
collectgarbage("collect")

print(("DIM=%d, square range=%d slots, packed range=%d slots (%.1f%% of square)")
    :format(DIM, 3 * DIM * DIM, 3 * (DIM * (DIM + 1) / 2), 100 * (DIM + 1) / (2 * DIM)))

local REPS = 15

-- both loops now do EXACTLY the same shape of arithmetic per access:
-- one precomputed row offset lookup, one sub, one mul-by-3, one add.
local function full_scan_square(t, row_off)
    local sum = 0
    for i = 1, DIM do
        local ro = row_off[i]
        for j = i, DIM do
            local base = 3 * (ro + (j - 1))
            sum = sum + t[base + 1] + t[base + 2] + t[base + 3]
        end
    end
    return sum
end

local function full_scan_packed(t, row_off)
    local sum = 0
    for i = 1, DIM do
        local ro = row_off[i]
        for j = i, DIM do
            local base = 3 * (ro + (j - i))
            sum = sum + t[base + 1] + t[base + 2] + t[base + 3]
        end
    end
    return sum
end

print()
print(("=== Full repeated scan, %d reps (index arithmetic now equal-cost) ==="):format(REPS))

local t0 = clock()
local s1 = 0
for _ = 1, REPS do
    s1 = s1 + full_scan_square(sq, sq_off)
end
local t1 = clock()
print(("square (current, gapped):  %.4f s  (sum=%.3f)"):format(t1 - t0, s1))

t0 = clock()
local s2 = 0
for _ = 1, REPS do
    s2 = s2 + full_scan_packed(pk, pk_off)
end
t1 = clock()
print(("packed (triangular, dense): %.4f s  (sum=%.3f)"):format(t1 - t0, s2))
