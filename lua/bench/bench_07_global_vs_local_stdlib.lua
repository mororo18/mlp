-- Validates the LuaJIT Numerical Computing Performance Guide's advice:
-- "cache non-FFI functions in upvalues: local sin = math.sin". Our real
-- code accesses math.floor/math.random/math.huge/os.clock via the global
-- `math`/`os` table on every call, including inside search_reinsertion
-- (60% of profiled time) and friends -- never localized. Checking
-- whether that actually costs anything measurable in this environment.

local clock = os.clock
local N = 5000000

collectgarbage("collect")
local t0 = clock()
local acc = 0
for i = 1, N do
    acc = acc + math.floor(i / 2) -- global table + field lookup, every call
end
local t1 = clock()
print(("global math.floor  x%d: %.4f s (acc=%d)"):format(N, t1 - t0, acc))

collectgarbage("collect")
local floor = math.floor -- localize once, outside any loop
t0 = clock()
acc = 0
for i = 1, N do
    acc = acc + floor(i / 2)
end
t1 = clock()
print(("local  floor        x%d: %.4f s (acc=%d)"):format(N, t1 - t0, acc))

-- math.huge: not a function call, just a field read (math.huge), but
-- still a global+field lookup done once per search_* call (thousands of
-- times) in the real code, e.g. `local cost_best = math.huge` at the top
-- of search_swap/search_reinsertion/search_two_opt.
collectgarbage("collect")
t0 = clock()
local s = 0
for i = 1, N do
    local cost_best = math.huge
    s = s + (cost_best > i and 1 or 0)
end
t1 = clock()
print(("global math.huge   x%d: %.4f s"):format(N, t1 - t0))

collectgarbage("collect")
local HUGE = math.huge
t0 = clock()
s = 0
for i = 1, N do
    local cost_best = HUGE
    s = s + (cost_best > i and 1 or 0)
end
t1 = clock()
print(("local  HUGE        x%d: %.4f s"):format(N, t1 - t0))
