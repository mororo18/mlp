dofile("Data.lua")

function s_print(solut)
    for i=1,#solut.s do
        io.write(solut.s[i], " ")
    end
    print()
end

function seq_print(solut)
    for i=1,#solut.seq do
        for j=1,#solut.seq do
            for s=1,3 do
                io.write(solut.seq[i][j][s] , " ") --= {0, 0, 0}
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

function solut_clone(solut)
    local cpy = {}
    --seq_print(solut)    
    cpy.seq = table.clone(solut.seq)
    cpy.s = table.clone(solut.s)
    cpy.cost = solut.cost

    return cpy
end

function subseq_fill(seq, info)
    for i=1,info.dimension+1 do
        seq[i] = {}
        for j=1,info.dimension+1 do
            seq[i][j] = {0, 0, 0}
        end
    end
end

function subseq_load(solut, info)
    for i=1,info.dimension+1 do
        local k = 1 - i - (i ~= 1 and 0 or 1)

        solut.seq[i][i][info.T] = 0.0
        solut.seq[i][i][info.C] = 0.0
        solut.seq[i][i][info.W] = (i ~= 1 and 1 or 0)
        for j=i+1,info.dimension+1 do
            local j_prev = j-1

            solut.seq[i][j][info.T] = info.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][info.T]
            solut.seq[i][j][info.C] = solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C]
            solut.seq[i][j][info.W] = j + k
        end
    end

    solut.cost = solut.seq[1][info.dimension+1][info.C] - info.EPSILON
    --print(solut.cost)
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

    while #cList > 0 do
    --while table.getn(cList) > 0 do
        table.sort(cList, function(i, j) return info.c[r][i] < info.c[r][j] end)

        local range = math.floor(#cList * alpha) + 1
        local i = math.random(1, range)

        local c = table.remove(cList, i)
        table.insert(s, c)
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

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_concat_4 = 0.0
    local cost_new = 0.0

    for i = 2, info.dimension-1 do
        local i_prev = i - 1
        local i_next = i + 1

        cost_concat_1 =                 solut.seq[1][i_prev][ info.T] + info.c[solut.s[i_prev]][solut.s[i_next]]
        cost_concat_2 = cost_concat_1 + solut.seq[i][i_next][info.T] + info.c[solut.s[i]][solut.s[i_next+1]]

        cost_new = solut.seq[1][i_prev][info.C]                                                    +           --       1st subseq
        solut.seq[i][i_next][info.W]               * (cost_concat_1) + info.c[solut.s[i_next]][solut.s[i]]  +           -- concat 2nd subseq
        solut.seq[i_next+1][info.dimension+1][info.W]   * (cost_concat_2) + solut.seq[i_next+1][info.dimension+1][info.C]   -- concat 3rd subseq

        if cost_new < cost_best then
            cost_best = cost_new - info.EPSILON
            I = i
            J = i_next
        end

        for j = i_next+1,info.dimension do
            local j_prev = j-1
            local j_next = j+1


            cost_concat_1 =                 solut.seq[1][i_prev][info.T]       + info.c[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1                           + info.c[solut.s[j]][solut.s[i_next]]
            cost_concat_3 = cost_concat_2 + solut.seq[i_next][j_prev][info.T]  + info.c[solut.s[j_prev]][solut.s[i]]
            cost_concat_4 = cost_concat_3                           + info.c[solut.s[i]][solut.s[j_next]]


            cost_new = solut.seq[1][i_prev][info.C]                                                 +      -- 1st subseq
            cost_concat_1 +                                                             -- concat 2nd subseq (single node)
            solut.seq[i_next][j_prev][info.W]      * cost_concat_2 + solut.seq[i_next][j_prev][info.C] +      -- concat 3rd subseq
            cost_concat_3 +                                                             -- concat 4th subseq (single node)
            solut.seq[j_next][info.dimension+1][info.W] * cost_concat_4 + solut.seq[j_next][info.dimension+1][info.C]   -- concat 5th subseq

            if(cost_new < cost_best) then
                cost_best = cost_new - info.EPSILON;
                I = i;
                J = j;
            end

        end
    end

    if cost_best < solut.seq[1][info.dimension+1][info.C] - info.EPSILON then
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

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_new = 0.0

    for i = 2,info.dimension-1 do
        local i_prev = i - 1
        local rev_seq_cost = solut.seq[i][i+1][info.T]
        for j = i+2,info.dimension do
            local j_next = j+1

            rev_seq_cost = rev_seq_cost + info.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][j][info.W]-1.0)

            cost_concat_1 =                 solut.seq[1][i_prev][info.T]   + info.c[solut.s[j]][solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T]        + info.c[solut.s[j_next]][solut.s[i]]

            cost_new = solut.seq[1][i_prev][info.C]                                                        +   --   1st subseq
                    solut.seq[i][j][info.W]                * cost_concat_1 + rev_seq_cost                  +   -- concat 2nd subseq (reversed seq)
                    solut.seq[j_next][info.dimension+1][info.W] * cost_concat_2 + solut.seq[j_next][info.dimension+1][info.C]       -- concat 3rd subseq

            if cost_new < cost_best then
                cost_best = cost_new -info.EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.seq[1][info.dimension+1][info.C] - info.EPSILON then
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

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_new = 0.0

    for i = 2, info.dimension-opt+1 do
        local j = opt+i-1
        local i_prev = i-1
        local j_next = j+1

        for k = 1, i_prev-1 do
            local k_next = k+1

            cost_concat_1 =                 solut.seq[1][k][info.T]            + info.c[solut.s[k]][solut.s[i]];
            cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T]            + info.c[solut.s[j]][solut.s[k_next]];
            cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][info.T]  + info.c[solut.s[i_prev]][solut.s[j_next]];

            cost_new = solut.seq[1][k][info.C]                                                             +   --       1st subseq
            solut.seq[i][j][info.W]                * cost_concat_1 + solut.seq[i][j][info.C]                  +   -- concat 2nd subseq (reinserted seq)
            solut.seq[k_next][i_prev][info.W]      * cost_concat_2 + solut.seq[k_next][i_prev][info.C]        +   -- concat 3rd subseq
            solut.seq[j_next][info.dimension+1][info.W] * cost_concat_3 + solut.seq[j_next][info.dimension+1][info.C];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - info.EPSILON
                I = i
                J = j
                POS = k
                --test = table.clone(test_a)
            end
        end

        for k = i+opt,info.dimension-opt-1 do
            local k_next = k+1

            cost_concat_1 =                 solut.seq[1][i_prev][info.T]   + info.c[solut.s[i_prev]][solut.s[j_next]];
            cost_concat_2 = cost_concat_1 + solut.seq[j_next][k][info.T]   + info.c[solut.s[k]][solut.s[i]];
            cost_concat_3 = cost_concat_2 + solut.seq[i][j][info.T]        + info.c[solut.s[j]][solut.s[k_next]];

            cost_new = solut.seq[1][i_prev][info.C]                                                        +   --       1st subseq
                    solut.seq[j_next][k][info.W]           * cost_concat_1 + solut.seq[j_next][k][info.C]             +   -- concat 2nd subseq
                    solut.seq[i][j][info.W]                * cost_concat_2 + solut.seq[i][j][info.C]                  +   -- concat 3rd subseq (reinserted seq)
                    solut.seq[k_next][info.dimension+1][info.W] * cost_concat_3 + solut.seq[k_next][info.dimension+1][info.C];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - info.EPSILON
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

    local neighbd_list = {SWAP, REINSERTION, OR_OPT_2, OR_OPT_3, TWO_OPT}

    --s_print(solut)
    --reinsert(solut.s, 2, 8, 13)
    --s_print(solut)

    --os.exit()

    while #neighbd_list > 0 do
        local index = math.random(1, #neighbd_list)
        local neighbd = neighbd_list[index]

        local improve = false

        if neighbd == SWAP then
            improve = search_swap(solut, info)
        elseif neighbd == REINSERTION then
            improve = search_reinsertion(solut, info, REINSERTION)
        elseif neighbd == OR_OPT_2 then
            improve = search_reinsertion(solut, info, OR_OPT_2)
        elseif neighbd == OR_OPT_3 then
            improve = search_reinsertion(solut, info, OR_OPT_3)
        elseif neighbd == TWO_OPT then
            improve = search_two_opt(solut, info)
        end


        if improve == true then
            local neighbd_list = {SWAP, REINSERTION, OR_OPT_2, OR_OPT_3, TWO_OPT}
        else 
            table.remove(neighbd_list, index)
        end
        
    end
end

function perturb(solut)
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

    subseq_fill(solut_partial.seq, info)
    solut_partial.cost = 0

    local solut_crnt = solut_clone(solut_partial)
    local solut_best = solut_clone(solut_crnt)

    solut_best.cost = math.huge

    for i=1,Imax do
        local Rsz = #R
        local alpha = R[math.random(1, Rsz)]

        solut_crnt.s = construction(alpha, info)
        subseq_load(solut_crnt, info)
        solut_partial = solut_clone(solut_crnt)

        --print(solut_partial.seq[1][info.dimension+1][info.C])
        --seq_print(solut_crnt)
        --local rnvd_cost_best = solut_crnt.seq[1][info.dimension+1][info.C]

        local iterILS = 0
        while iterILS < Iils do
            RVND(solut_crnt, info)
            --local rnvd_cost_crnt = solut_crnt.seq[1][info.dimension+1][info.C]
            --if rnvd_cost_crnt < rvnd_cost_best then
               --rnvd_cost_best = rnvd_cost_crnt - info.EPSILON
            if solut_crnt.cost < solut_partial.cost - info.EPSILON then
               solut_partial.cost = solut_crnt.cost - info.EPSILON
               solut_partial = solut_clone(solut_crnt)
               iterILS = 0
            end
            
            solut_crnt.s = perturb(solut_partial)
            subseq_load(solut_crnt, info)

            iterILS = iterILS + 1
        end

        subseq_load(solut_partial, info)
        --local cost_partial = solut_partial.seq[1][info.dimension+1][info.C]

        --if cost_partial < cost_best then
            --cost_best  = cost_partial
        if solut_partial.cost < solut_best.cost then
            solut_best = solut_clone(solut_partial)
        end
        --print("\tCurrent: ", solut_best.cost)
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
        EPSILON = 1e-15
    }


    info.dimension = readData(info.c)
    math.randomseed(os.time()) 

    local Imax = 10
    local Iils = math.min(100, info.dimension)
    local R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
               0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}


    info = protect(info)
    local start = os.clock()
    GILS_RVND(Imax, Iils, R, info)

    print("TIME: ", os.clock()-start)
end

main()

--[[
collectgarbage("stop")
ProFi = require 'ProFi'
ProFi:start()
main()
--coroutine.resume( main )
ProFi:stop()
ProFi:writeReport( 'MyProfilingReport.txt' )
]]--
