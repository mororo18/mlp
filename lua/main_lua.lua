--require("candran").setup()
dofile("Data.lua")

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
            io.write(solut.seq[i][j].T , " ") --= {0, 0, 0}
            io.write(solut.seq[i][j].W , " ") --= {0, 0, 0}
            io.write(solut.seq[i][j].C , " ") --= {0, 0, 0}
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

--#define("to_1D(x, y, z, d)", "(3*((d)*(x-1) + (y-1)) + z)")

function subseq_fill(seq, info)
    for i=1,info.dimension+1 do
        seq[i] = {}
        for j=i,info.dimension+1 do
            seq[i][j] = {T=0, C=0, W=0}
        end
    end
end

function subseq_load(solut, info)
    local dimen = info.dimension
    local seq = solut.seq

    local c = info.c
    local s = solut.s
    local s_size = dimen+1

    for i=1,s_size do
        local k = 1 - i - (i ~= 1 and 0 or 1)

        seq[i][i].T = 0.0
        seq[i][i].C = 0.0
        seq[i][i].W = (i ~= 1 and 1 or 0)

        local some_T = 0.0
        local some_C = 0.0

        local s_j_prev = s[i]

        for j=i+1,s_size do
            local j_prev = j-1

            local s_j = s[j]

            -- set T
            local cost_T= c[s_j_prev][s_j] + some_T
            seq[i][j].T = cost_T

            -- set C
            local cost_C = cost_T + some_C
            seq[i][j].C  = cost_C

            -- set W
            seq[i][j].W = j + k


            s_j_prev = s_j
            some_T = cost_T
            some_C = cost_C
        end
    end

    solut.cost = seq[1][dimen+1].C - info.EPSILON
    --print(solut.cost)
end

function sort(arr, r, info) 
    for i = 1, #arr do
        for j = 1, #arr-i do
            if info.c[r][arr[j]] > info.c[r][arr[j+1]] then
                local tmp = arr[j]
                arr[j] = arr[j+1]
                arr[j+1] = tmp
            end
        end
    end
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
        local range = math.floor(#cList * alpha) + 1
        local i = math.random(1, range)
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
    for k = i,math.floor((j+i)/2) do
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
    local cost_best = math.huge
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

    local EP = info.EPSILON

    for i = 2, dimen-1 do
        local i_prev = i - 1
        local i_next = i + 1

        local s_i = s[i]
        local s_i_next = s[i_next]

        cost_concat_1 =                 seq[1][i_prev].T + c[s[i_prev]][s_i_next]
        cost_concat_2 = cost_concat_1 + seq[i][i_next].T + c[s_i][     s[i_next+1]]

        cost_new = seq[1][i_prev].C                                                    +           --       1st subseq
        seq[i][i_next].W               * (cost_concat_1) + c[s_i_next][s_i]  +           -- concat 2nd subseq
        seq[i_next+1][s_size].W   * (cost_concat_2) + seq[i_next+1][s_size].C   -- concat 3rd subseq

        if cost_new < cost_best then
            cost_best = cost_new - EP
            I = i
            J = i_next
        end

        for j = i_next+1,dimen do
            local j_prev = j-1
            local j_next = j+1
            local s_j = s[j]


            cost_concat_1 =                 seq[1][i_prev].T       + c[s[i_prev]][s_j]
            cost_concat_2 = cost_concat_1                     + c[s_j][s_i_next]
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev].T  + c[s[j_prev]][s_i]
            cost_concat_4 = cost_concat_3                           + c[s_i][s[j_next]]


            cost_new = seq[1][i_prev].C                                                 +      -- 1st subseq
            cost_concat_1 +                                                             -- concat 2nd subseq (single node)
            seq[i_next][j_prev].W      * cost_concat_2 + seq[i_next][j_prev].C +      -- concat 3rd subseq
            cost_concat_3 +                                                             -- concat 4th subseq (single node)
            seq[j_next][s_size].W * cost_concat_4 + seq[j_next][s_size].C   -- concat 5th subseq

            if(cost_new < cost_best) then
                cost_best = cost_new - EP
                I = i;
                J = j;
            end

        end
    end

    if cost_best < solut.seq[1][s_size].C - info.EPSILON then
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
    local cost_best = math.huge
    local I = -1
    local J = -1
    local dimen = info.dimension

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_new = 0.0

    local s = solut.s
    local seq = solut.seq
    local c = info.c

    local EP = info.EPSILON

    local s_size = dimen+1

    for i = 2,dimen-1 do
        local i_prev = i - 1
        local rev_seq_cost = seq[i][i+1].T

        local s_j_prev = s[i+1]

        for j = i+2,dimen do
            local j_next = j+1
            local s_j = s[j]

            rev_seq_cost = rev_seq_cost + c[s_j_prev][s_j] * (seq[i][j].W-1.0)

            cost_concat_1 =                 seq[1][i_prev].T   + c[s_j][     s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[i][j].T        + c[s[j_next]][s[i]]

            cost_new = seq[1][i_prev].C                                                        +   --   1st subseq
                    seq[i][j].W                * cost_concat_1 + rev_seq_cost                  +   -- concat 2nd subseq (reversed seq)
                    seq[j_next][s_size].W * cost_concat_2 + seq[j_next][s_size].C       -- concat 3rd subseq

            if cost_new < cost_best then
                cost_best = cost_new - EP
                I = i
                J = j
            end

            s_j_prev = s_j
        end
    end

    if cost_best < solut.seq[1][s_size].C - EP then
        reverse(solut.s, I, J)
        subseq_load(solut, info)
        return true
    end

    return false
end

function search_reinsertion(solut, info, opt)
    local cost_best = math.huge
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

            cost_concat_1 =                 seq[1][k].T            + c[s[k]][s_i];
            cost_concat_2 = cost_concat_1 + seq[i][j].T            + c[s_j][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev].T  + c[s_i_prev][s_j_next];

            cost_new = seq[1][k].C                                                             +   --       1st subseq
            seq[i][j].W                * cost_concat_1 + seq[i][j].C                  +   -- concat 2nd subseq (reinserted seq)
            seq[k_next][i_prev]. W      * cost_concat_2 + seq[k_next][i_prev]. C        +   -- concat 3rd subseq
            seq[j_next][s_size].W * cost_concat_3     + seq[j_next][s_size].C;       -- concat 4th subseq

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

            cost_concat_1 =                 seq[1][i_prev].T   + c[s_i_prev][s_j_next];
            cost_concat_2 = cost_concat_1 + seq[j_next][k].T   + c[s[k]][     s_i];
            cost_concat_3 = cost_concat_2 + seq[i][j].T        + c[s_j][     s[k_next]];

            cost_new = seq[1][i_prev].C                                                        +   --       1st subseq
                    seq[j_next][k].W           * cost_concat_1 + seq[j_next][k].C             +   -- concat 2nd subseq
                    seq[i][j].W                * cost_concat_2 + seq[i][j].     C                  +   -- concat 3rd subseq (reinserted seq)
                    seq[k_next][s_size].W * cost_concat_3     + seq[k_next][s_size].C;       -- concat 4th subseq

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
        local index = math.random(1, #neighbd_list)

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

    local size_max = math.floor(#s/10)
    size_max = size_max >= 2 and size_max or 2
    local size_min = 2

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do
        A_start = math.random(2, #s-1-size_max)
        A_end = A_start + math.random(size_min, size_max)

        B_start = math.random(2, #s-1-size_max)
        B_end = B_start + math.random(size_min, size_max)


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

function GILS_RVND(Imax, Iils, R, info)

    local solut_partial = {
        s = {},
        seq = {}, 
    }

    local solut_crnt = {
        s = {},
        seq = {}, 
    }

    local solut_best = {
        s = {},
        seq = {}, 
    }

    subseq_fill(solut_partial.seq, info)
    solut_partial.cost = 0

    subseq_fill(solut_crnt.seq, info)
    solut_crnt.cost = 0

    subseq_fill(solut_best.seq, info)
    solut_best.cost = 0

    --local solut_crnt = solut_clone(solut_partial)
    --local solut_best = solut_clone(solut_crnt)

    solut_best.cost = math.huge

    for i=1,Imax do
        local Rsz = #R
        local alpha = R[math.random(1, Rsz)]
        alpha = R[info.rnd[info.rnd_index] + 1]
        --print(info.rnd[info.rnd_index] + 1)
        info.rnd_index = info.rnd_index + 1


        print("[+] Local Search", i)
        print("\t[+] Constructing Inital Solution..")
        solut_crnt.s = construction(alpha, info)
        subseq_load(solut_crnt, info)
        s_print(solut_crnt)
        print("\tConstruction cost  ", solut_crnt.cost)

      --print(solut_crnt.seq[to_1D(1, info.dimension+1, info.C, info.dimension+1)])

      --seq_print(solut_crnt, info)
      --os.exit()
        solut_clone(solut_crnt, solut_partial)
        --solut_partial = solut_clone(solut_crnt)

        --print(solut_partial.seq[1][info.dimension+1][info.C])
        --seq_print(solut_crnt)
        --local rnvd_cost_best = solut_crnt.seq[1][info.dimension+1][info.C]

        print("\t[+] Looking for the best Neighbor..")
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
        print("\tCurrent best solution cost: ", solut_best.cost)
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
    math.randomseed(os.time()) 

    local Imax = 10
    local Iils = math.min(100, info.dimension)
    local R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
               0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}


    --info.c = protect(info.c)
    local start = os.clock()
    GILS_RVND(Imax, Iils, R, info)

    print("TIME: ", os.clock()-start)
end

-- jit.off(table.clone)

local profiler = require("profiler")
profiler.start()
-- Code block and/or called functions to profile --
main()
profiler.stop()
profiler.report("profiler.log")

--[[
collectgarbage("stop")
ProFi = require 'ProFi'
ProFi:start()
main()
--coroutine.resume( main )
ProFi:stop()
ProFi:writeReport( 'MyProfilingReport.txt' )
]]--
