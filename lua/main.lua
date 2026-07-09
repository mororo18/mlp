dofile("Data.lua")

-- Localized stdlib references: avoids a global-table lookup on every
-- call in the hot loop (search_reinsertion et al.). Real, validated win
-- in lua5.3 (~25-35%, pure interpreter pays the lookup every time);
-- neutral in luajit (the JIT already optimizes the lookup away during
-- trace compilation) -- see lua/IMPLEMENTATION.md and
-- lua/bench/bench_07_global_vs_local_stdlib.lua.
local floor = math.floor
local random = math.random
local randomseed = math.randomseed
local huge = math.huge
local min = math.min
local clock = os.clock

-- FFI array for seq (LuaJIT only; falls back to a plain Lua table on
-- lua5.3, where `ffi` doesn't exist). `[]` indexing is syntactically
-- identical for a Lua table and an FFI cdata array, so every function
-- that reads/writes seq[idx] (subseq_fill, subseq_load, search_*) needs
-- no changes -- only this factory branches. n+1 keeps index 0 unused but
-- present, so every existing 1-based index formula stays unchanged. See
-- lua/IMPLEMENTATION.md and lua/bench/bench_08_ffi_vs_table.lua,
-- lua/bench/bench_09_ffi_conditional_pattern.lua.
local ffi_ok, ffi = pcall(require, "ffi")
local function new_numeric_array(n)
    if ffi_ok then
        return ffi.new("double[?]", n + 1)
    else
        local t = {}
        for i = 1, n do
            t[i] = 0.0
        end
        return t
    end
end

function s_print(solut)
    for i=1,#solut.s do
        io.write(solut.s[i], " ")
    end
    print()
end

function table_print(tbl)
    for i=1,#tbl do
        io.write(tbl[i], " ")
    end
    print()
end

function matrix_print(info)
    for i = 1, #info.c do
        for j = 1, #info.c do
            io.write(info.c[i][j], " ")
        end
        print()
    end

end

function seq_print(solut, info)
    --[[
    for i=1,#solut.seq do
        for j=1,#solut.seq do
            for s=1,3 do
            end
        end
    end
    ]]--

    for i=1,info.dimension+1 do
        for j=i,info.dimension+1 do
            for k=1,3 do
                io.write(solut.seq[(3 * ((info.dimension+1) * (i - 1) + (j - 1)) + k)] , " ") --= {0, 0, 0}
            end
            io.write("| ")
        end
        print()
    end
end

function table.clone(org)
    if type(jit) == 'table' then
        --print(jit.version)  --LuaJIT 2.0.2
        return {unpack(org)}
    else
        return {table.unpack(org)}
    end

end

function solut_clone(src, dest)
    --local cpy = {}

    --cpy.seq = table.clone(solut.seq)
    for i=1, #src.s do
        dest.s[i] = src.s[i]
    end
    --cpy.s = table.clone(solut.s)
    dest.cost = src.cost

    return dest
end


function subseq_fill(seq, info)
    for i=1,info.dimension+1 do
        for j=i,info.dimension+1 do
            for k=1,3 do
                seq[(3 * ((info.dimension+1) * (i - 1) + (j - 1)) + k)] = 0.0
            end
        end
    end
end

function subseq_load(solut, info)
    local dimen = info.dimension
    local seq = solut.seq

    local T = info.T
    local W = info.W
    local C = info.C

    local c = info.c
    local s = solut.s
    local s_size = dimen+1



    for i=1,s_size do
        local k = 1 - i - (i ~= 1 and 0 or 1)

        seq[(3 * ((s_size) * (i - 1) + (i - 1)) + T)] = 0.0
        seq[(3 * ((s_size) * (i - 1) + (i - 1)) + C)] = 0.0
        seq[(3 * ((s_size) * (i - 1) + (i - 1)) + W)] = (i ~= 1 and 1 or 0)

        local some_T = 0.0
        local some_C = 0.0

        local s_j_prev = s[i]

        for j=i+1,s_size do
            local j_prev = j-1

            local s_j = s[j]

            -- set T
            local cost_T= c[s_j_prev][s_j] + some_T
            seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)] = cost_T

            -- set C
            local cost_C = cost_T + some_C
            seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)]  = cost_C

            -- set W
            seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)] = j + k


            s_j_prev = s_j
            some_T = cost_T
            some_C = cost_C
        end
    end

    solut.cost = seq[(3 * ((dimen+1) * (1 - 1) + (dimen+1 - 1)) + C)] - info.EPSILON
    --print(solut.cost)
end

function quicksort(arr, left, right, r, info)
    if left < right then
        local pivotIndex = partition(arr, left, right, r, info)
        quicksort(arr, left, pivotIndex - 1, r, info)
        quicksort(arr, pivotIndex + 1, right, r, info)
    end
end

function partition(arr, left, right, r, info)
    local pivotValue = arr[right]
    local i = left - 1

    for j = left, right - 1 do
        if info.c[r][arr[j]] < info.c[r][pivotValue] then
            i = i + 1
            arr[i], arr[j] = arr[j], arr[i]
        end
    end

    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1
end

function sort(arr, r, info)
    quicksort(arr, 1, #arr, r, info)
end

function construction(alpha, info) 
    local s = {1}

    local cList = {}
    for i=2,info.dimension do
        cList[#cList + 1] = i
        --io.write(cList[i-1], " ")
    end
    --print()

    local r = 1

    --print("alpha", alpha)
    while #cList > 0 do
    --while table.getn(cList) > 0 do
        --table.sort(cList, function(i, j) return info.c[r][i] <= info.c[r][j] end)
        sort(cList, r, info)

        --print(alpha)
        local range = floor(#cList * alpha) + 1
        local i = random(1, range)
        i = info.rnd[info.rnd_index] + 1
        info.rnd_index = info.rnd_index + 1

        local c = table.remove(cList, i)
        table.insert(s, c)
        --print(c, info.c[r][c])
        r = c
    end
    table.insert(s, 1)

    return table.clone(s)
end

function swap(s, i, j)
    s[j], s[i] = s[i], s[j]
end

function reverse(s, i, j)
    local l = j
    for k = i,floor((j+i)/2) do
        --print(k, l)
        swap(s, k, l)
        l = l - 1
    end
end

function reinsert(s, i, j, pos)
    if i < pos then
        for k = i,j do
            local e = s[i]
            table.insert(s, pos, e)
            table.remove(s, i)
        end
    else
        for k = i,j do
            local e = s[j]
            table.remove(s, j)
            table.insert(s, pos, e)
        end
    end
end

function search_swap(solut, info)
    local cost_best = huge
    local I = -1
    local J = -1
    local dimen = info.dimension

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_concat_4 = 0.0
    local cost_new = 0.0

    local s = solut.s
    local seq = solut.seq
    local c = info.c

    local s_size = dimen+1

    local T = info.T
    local C = info.C
    local W = info.W

    local EP = info.EPSILON

    for i = 2, dimen-1 do
        local i_prev = i - 1
        local i_next = i + 1

        local s_i = s[i]
        local s_i_next = s[i_next]

        cost_concat_1 =                 seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)] + c[s[i_prev]][s_i_next]
        cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + T)] + c[s_i][     s[i_next+1]]

        cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)]                                                    +           --       1st subseq
        seq[(3 * ((s_size) * (i - 1) + (i_next - 1)) + W)]               * (cost_concat_1) + c[s_i_next][s_i]  +           -- concat 2nd subseq
        seq[(3 * ((s_size) * (i_next+1 - 1) + (s_size - 1)) + W)]   * (cost_concat_2) + seq[(3 * ((s_size) * (i_next+1 - 1) + (s_size - 1)) + C)]   -- concat 3rd subseq

        if cost_new < cost_best then
            cost_best = cost_new - EP
            I = i
            J = i_next
        end

        for j = i_next+1,dimen do
            local j_prev = j-1
            local j_next = j+1
            local s_j = s[j]


            cost_concat_1 =                 seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)]       + c[s[i_prev]][s_j]
            cost_concat_2 = cost_concat_1                     + c[s_j][s_i_next]
            cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + T)]  + c[s[j_prev]][s_i]
            cost_concat_4 = cost_concat_3                           + c[s_i][s[j_next]]


            cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)]                                                 +      -- 1st subseq
            cost_concat_1 +                                                             -- concat 2nd subseq (single node)
            seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + W)]      * cost_concat_2 + seq[(3 * ((s_size) * (i_next - 1) + (j_prev - 1)) + C)] +      -- concat 3rd subseq
            cost_concat_3 +                                                             -- concat 4th subseq (single node)
            seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_4 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)]   -- concat 5th subseq

            if(cost_new < cost_best) then
                cost_best = cost_new - EP
                I = i;
                J = j;
            end

        end
    end

    if cost_best < solut.seq[(3 * ((s_size) * (1 - 1) + (s_size - 1)) + info.C)] - info.EPSILON then
        --print("swap")
        --print(cost_best)
        swap(solut.s, I, J)
        subseq_load(solut, info)
        --print(solut.seq[1][info.dimension+1][info.C])
        return true
    end

    return false
end

function search_two_opt(solut, info)
    local cost_best = huge
    local I = -1
    local J = -1
    local dimen = info.dimension

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_new = 0.0

    local s = solut.s
    local seq = solut.seq
    local c = info.c

    local T = info.T
    local C = info.C
    local W = info.W

    local EP = info.EPSILON

    local s_size = dimen+1

    for i = 2,dimen-1 do
        local i_prev = i - 1
        local rev_seq_cost = seq[(3 * ((s_size) * (i - 1) + (i+1 - 1)) + T)]

        local s_j_prev = s[i+1]

        for j = i+2,dimen do
            local j_next = j+1
            local s_j = s[j]

            rev_seq_cost = rev_seq_cost + c[s_j_prev][s_j] * (seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)]-1.0)

            cost_concat_1 =                 seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)]   + c[s_j][     s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)]        + c[s[j_next]][s[i]]

            cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)]                                                        +   --   1st subseq
                    seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)]                * cost_concat_1 + rev_seq_cost                  +   -- concat 2nd subseq (reversed seq)
                    seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_2 + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)]       -- concat 3rd subseq

            if cost_new < cost_best then
                cost_best = cost_new - EP
                I = i
                J = j
            end

            s_j_prev = s_j
        end
    end

    if cost_best < solut.seq[(3 * ((s_size) * (1 - 1) + (s_size - 1)) + C)] - EP then
        reverse(solut.s, I, J)
        subseq_load(solut, info)
        return true
    end

    return false
end

function search_reinsertion(solut, info, opt)
    local cost_best = huge
    local I = -1
    local J = -1
    local POS = -1
    local dimen = info.dimension

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_new = 0.0

    local s = solut.s
    local seq = solut.seq
    local c = info.c
    local s_size = dimen+1

    local T = info.T
    local C = info.C
    local W = info.W

    local EP = info.EPSILON

    for i = 2, dimen-opt+1 do
        local j = opt+i-1
        local i_prev = i-1
        local j_next = j+1

        local s_i = s[i]
        local s_j = s[j]

        local s_i_prev = s[i_prev]
        local s_j_next = s[j_next]

        for k = 1, i_prev-1 do
            local k_next = k+1

            cost_concat_1 =                 seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + T)]            + c[s[k]][s_i];
            cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)]            + c[s_j][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + T)]  + c[s_i_prev][s_j_next];

            cost_new = seq[(3 * ((s_size) * (1 - 1) + (k - 1)) + C)]                                                             +   --       1st subseq
            seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)]                * cost_concat_1 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)]                  +   -- concat 2nd subseq (reinserted seq)
            seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + W)]      * cost_concat_2 + seq[(3 * ((s_size) * (k_next - 1) + (i_prev - 1)) + C)]        +   -- concat 3rd subseq
            seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + W)] * cost_concat_3     + seq[(3 * ((s_size) * (j_next - 1) + (s_size - 1)) + C)];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - EP
                I = i
                J = j
                POS = k
                --test = table.clone(test_a)
            end
        end

        for k = i+opt,dimen do
            local k_next = k+1

            cost_concat_1 =                 seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + T)]   + c[s_i_prev][s_j_next];
            cost_concat_2 = cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + T)]   + c[s[k]][     s_i];
            cost_concat_3 = cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + T)]        + c[s_j][     s[k_next]];

            cost_new = seq[(3 * ((s_size) * (1 - 1) + (i_prev - 1)) + C)]                                                        +   --       1st subseq
                    seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + W)]           * cost_concat_1 + seq[(3 * ((s_size) * (j_next - 1) + (k - 1)) + C)]             +   -- concat 2nd subseq
                    seq[(3 * ((s_size) * (i - 1) + (j - 1)) + W)]                * cost_concat_2 + seq[(3 * ((s_size) * (i - 1) + (j - 1)) + C)]                  +   -- concat 3rd subseq (reinserted seq)
                    seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + W)] * cost_concat_3     + seq[(3 * ((s_size) * (k_next - 1) + (s_size - 1)) + C)];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - EP
                I = i
                J = j
                POS = k
                --test = table.clone(test_b)
            end

        end
    end

    if cost_best < solut.cost then
        --print("reinsert", I, J, POS+1)
        --print(cost_best)
        --s_print(solut)
        reinsert(solut.s, I, J, POS+1)
        --s_print(solut)
        subseq_load(solut, info)
        --print(solut.cost)

        if cost_best ~= solut.cost then
            print("ERROR")
            os.exit(1)
        end
        return true
    end

    return false
end

function RVND(solut, info)
    local SWAP        = 0  
    local REINSERTION = 1
    local OR_OPT_2    = 2
    local OR_OPT_3    = 3
    local TWO_OPT     = 4

    local neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}

    --s_print(solut)
    --reinsert(solut.s, 2, 8, 13)
    --s_print(solut)

    --os.exit()

    while #neighbd_list > 0 do
        local index = random(1, #neighbd_list)

        index = info.rnd[info.rnd_index] + 1
        info.rnd_index = info.rnd_index + 1

        local neighbd = neighbd_list[index]

        local improve = false

        if neighbd == SWAP then
            improve = search_swap(solut, info)
            --print("swap")
        elseif neighbd == REINSERTION then
            improve = search_reinsertion(solut, info, REINSERTION)
            --print("reinsertion")
        elseif neighbd == OR_OPT_2 then
            improve = search_reinsertion(solut, info, OR_OPT_2)
            --print("or-opt2")
        elseif neighbd == OR_OPT_3 then
            improve = search_reinsertion(solut, info, OR_OPT_3)
            --print("or-opt3")
        elseif neighbd == TWO_OPT then
            improve = search_two_opt(solut, info)
            --print("two-opt")
        end


        if improve == true then
            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
            --print(solut.cost, "atual", info.rnd[info.rnd_index], "proximo", info.rnd[info.rnd_index+1])
        else 
            table.remove(neighbd_list, index)
        end
        --table_print(neighbd_list)
        
    end
end

function perturb(solut, info)
    local s = table.clone(solut.s)

    local A_start = 1
    local A_end = 1
    local B_start = 1
    local B_end = 1

    local size_max = floor(#s/10)
    size_max = size_max >= 2 and size_max or 2
    local size_min = 2

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do
        A_start = random(2, #s-1-size_max)
        A_end = A_start + random(size_min, size_max)

        B_start = random(2, #s-1-size_max)
        B_end = B_start + random(size_min, size_max)


        A_start = info.rnd[info.rnd_index] + 1
        info.rnd_index = info.rnd_index +1
        A_end = A_start + info.rnd[info.rnd_index]
        info.rnd_index = info.rnd_index +1

        B_start = info.rnd[info.rnd_index] + 1
        info.rnd_index = info.rnd_index +1
        B_end = B_start + info.rnd[info.rnd_index]
        info.rnd_index = info.rnd_index +1
    end

    if A_start < B_start then
        reinsert(s, B_start, B_end-1, A_end)
        reinsert(s, A_start, A_end-1, B_end)
    else
        reinsert(s, A_start, A_end-1, B_end)
        reinsert(s, B_start, B_end-1, A_end)
    end

    return s
end

function GILS_RVND(Imax, Iils, R, info, verbose)

    local seq_n = 3 * (info.dimension + 1) * (info.dimension + 1)

    local solut_partial = {
        s = {},
        seq = new_numeric_array(seq_n),
    }

    local solut_crnt = {
        s = {},
        seq = new_numeric_array(seq_n),
    }

    local solut_best = {
        s = {},
        seq = new_numeric_array(seq_n),
    }

    subseq_fill(solut_partial.seq, info)
    solut_partial.cost = 0

    subseq_fill(solut_crnt.seq, info)
    solut_crnt.cost = 0

    subseq_fill(solut_best.seq, info)
    solut_best.cost = 0

    --local solut_crnt = solut_clone(solut_partial)
    --local solut_best = solut_clone(solut_crnt)

    solut_best.cost = huge

    for i=1,Imax do
        local Rsz = #R
        local alpha = R[random(1, Rsz)]
        alpha = R[info.rnd[info.rnd_index] + 1]
        --print(info.rnd[info.rnd_index] + 1)
        info.rnd_index = info.rnd_index + 1


        if verbose then
            print("[+] Local Search", i)
            print("\t[+] Constructing Inital Solution..")
        end
        solut_crnt.s = construction(alpha, info)
        subseq_load(solut_crnt, info)
        if verbose then
            s_print(solut_crnt)
            print("\tConstruction cost  ", solut_crnt.cost)
        end

      --print(solut_crnt.seq[(3 * ((info.dimension+1) * (1 - 1) + (info.dimension+1 - 1)) + info.C)])

      --seq_print(solut_crnt, info)
      --os.exit()
        solut_clone(solut_crnt, solut_partial)
        --solut_partial = solut_clone(solut_crnt)

        --print(solut_partial.seq[1][info.dimension+1][info.C])
        --seq_print(solut_crnt)
        --local rnvd_cost_best = solut_crnt.seq[1][info.dimension+1][info.C]

        if verbose then
            print("\t[+] Looking for the best Neighbor..")
        end
        local iterILS = 0
        while iterILS < Iils do
            RVND(solut_crnt, info)
            --os.exit(0)
            --local rnvd_cost_crnt = solut_crnt.seq[1][info.dimension+1][info.C]
            --if rnvd_cost_crnt < rvnd_cost_best then
               --rnvd_cost_best = rnvd_cost_crnt - info.EPSILON
            if solut_crnt.cost < solut_partial.cost - info.EPSILON then
               solut_partial.cost = solut_crnt.cost - info.EPSILON
               solut_partial.s = table.clone(solut_crnt.s)
               iterILS = 0
            end
            
            solut_crnt.s = perturb(solut_partial, info)
            subseq_load(solut_crnt, info)

            iterILS = iterILS + 1
        end

      --print("partial cost", solut_partial.cost)
      --print("ITER LAST", info.rnd[info.rnd_index])
      --s_print(solut_partial)

        subseq_load(solut_partial, info)
        --local cost_partial = solut_partial.seq[1][info.dimension+1][info.C]

        --if cost_partial < cost_best then
            --cost_best  = cost_partial
        if solut_partial.cost < solut_best.cost then
            solut_clone(solut_partial, solut_best)
            --solut_best = solut_clone(solut_partial)
        end
        if verbose then
            print("\tCurrent best solution cost: ", solut_best.cost)
        end
    end

    print("COST: ", solut_best.cost)
end

function protect(tbl)
    return setmetatable({}, {
        __index = tbl,
        __newindex = function(t, key, value)
            error("attempting to change constant " ..
            tostring(key) .. " to " .. tostring(value), 2)
        end
    })
end

function main() 
    local info = {
        c = {}, 
        T = 1,
        C = 2, W = 3, 
        EPSILON = 1e-15,
        rnd = {},
        rnd_index = 1
    }


    local a =0;
    info.dimension, a = readData(info.c, info.rnd)
    print(info.rnd[a])
    randomseed(os.time()) 

    local Imax = 10
    local Iils = min(100, info.dimension)
    local R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
               0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}


    local verbose = false
    for i = 1, #arg do
        if arg[i] == "-v" or arg[i] == "--verbose" then
            verbose = true
        end
    end

    --info.c = protect(info.c)
    local start = clock()
    GILS_RVND(Imax, Iils, R, info, verbose)

    print("TIME: ", clock()-start)
end

-- jit.off(table.clone)

local do_profile = false
for i = 1, #arg do
    if arg[i] == "--profile" then
        do_profile = true
    end
end

if do_profile then
    local profiler = require("profiler")
    profiler.start()
    main()
    profiler.stop()
    profiler.report("profiler.log")
else
    main()
end

--[[
collectgarbage("stop")
ProFi = require 'ProFi'
ProFi:start()
main()
--coroutine.resume( main )
ProFi:stop()
ProFi:writeReport( 'MyProfilingReport.txt' )
]]--
