-- Validates whether table.insert()/table.remove() (confirmed via
-- `luajit -jv` to trigger trace "stitch" events in the real main.canc.lua
-- code -- i.e. NOT inlined by the JIT, bails to the interpreter every
-- call) are actually slower than the hand-written equivalent, which
-- traces/compiles inline.
local clock = os.clock
local N = 2000000

collectgarbage("collect")
local t0 = clock()
local t = {}
for i = 1, N do
    table.insert(t, i)
end
local t1 = clock()
print(("table.insert x%d:      %.4f s"):format(N, t1 - t0))

collectgarbage("collect")
t0 = clock()
local u = {}
for i = 1, N do
    u[#u + 1] = i
end
t1 = clock()
print(("manual u[#u+1]=x x%d:  %.4f s"):format(N, t1 - t0))

-- table.remove from the end (matches neighbd_list's usage pattern:
-- small list, remove-by-index) vs manual shift
local SMALL = 5
local REPS = 500000

collectgarbage("collect")
t0 = clock()
for _ = 1, REPS do
    local list = { 1, 2, 3, 4, 5 }
    table.remove(list, 3)
end
t1 = clock()
print(("table.remove(list,3) x%d: %.4f s"):format(REPS, t1 - t0))

collectgarbage("collect")
t0 = clock()
for _ = 1, REPS do
    local list = { 1, 2, 3, 4, 5 }
    -- manual remove-by-index (shift down), same semantics as table.remove
    local idx = 3
    for k = idx, #list - 1 do
        list[k] = list[k + 1]
    end
    list[#list] = nil
end
t1 = clock()
print(("manual shift-remove x%d:  %.4f s"):format(REPS, t1 - t0))
