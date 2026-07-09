-- Illustrative pattern for main.lua: FFI array on luajit, plain Lua
-- table fallback on lua5.3, same [] indexing syntax either way so the
-- rest of the code (subseq_fill, subseq_load, search_*) doesn't need to
-- know or care which one it got.

local ffi_ok, ffi = pcall(require, "ffi")

-- Factory: allocate a flat numeric array of size n (1-based indices
-- 1..n). FFI branch on luajit, plain table on lua5.3.
local function new_numeric_array(n)
    if ffi_ok then
        -- ffi.new zero-initializes automatically; size n+1 so index 0
        -- exists but goes unused, keeping every existing 1-based index
        -- formula in the codebase unchanged (no -1 rewrites needed).
        return ffi.new("double[?]", n + 1)
    else
        local t = {}
        for i = 1, n do
            t[i] = 0.0
        end
        return t
    end
end

-- Usage mirrors current subseq_fill/subseq_load exactly -- seq[idx] read
-- and write work identically whether seq is a table or a cdata array:
local function subseq_fill_example(seq, info)
    for i = 1, info.dimension + 1 do
        for j = i, info.dimension + 1 do
            for k = 1, 3 do
                seq[(3 * ((info.dimension + 1) * (i - 1) + (j - 1)) + k)] = 0.0
            end
        end
    end
end

-- At construction time (main.lua's main(), where `solut.seq = {}` is set
-- up today):
local function make_solut(info)
    local n = 3 * (info.dimension + 1) * (info.dimension + 1)
    return {
        seq = new_numeric_array(n),
        s = {}, -- tour/solution sequence: stays a plain table, needs
                -- table.insert/remove/# which FFI arrays don't support
        cost = 0,
    }
end

-- quick sanity check this file's logic actually works standalone
local ffi_ok2, ffi2 = pcall(require, "ffi")
print("ffi available:", ffi_ok2)
local info = { dimension = 4 }
local solut = make_solut(info)
subseq_fill_example(solut.seq, info)
print("type check:", type(solut.seq), "sample value at [1]:", solut.seq[1])
