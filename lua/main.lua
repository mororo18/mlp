dofile("Data.lua")

local T = 1;
local C = 2;
local W = 3;

local c = {};
local dimension;
local EPSILON = 1e-15;

function s_print(s)
    for i=1,#s do
        io.write(s[i], " ")
    end
    print()
end

function table_print(tbl)
    for i=1,#tbl do
        io.write(tbl[i], " ")
    end
    print()
end

function matrix_print(c)
    for i = 1, #c do
        for j = 1, #c do
            io.write(c[i][j], " ")
        end
        print()
    end

end

function seq_print(solut)
    for i=1,#seq do
        for j=1,#seq do
            for s=1,3 do
                io.write(seq[i][j][s] , " ") --= {0, 0, 0}
            end
            io.write("| ")
        end
        print()
    end
end

function table.clone(org)
    if type(jit) == 'table' then
        return {unpack(org)}
    else
        return {table.unpack(org)}
    end

end

function solut_clone(solut)
    local cpy = {}
    cpy.seq = table.clone(solut.seq)
    cpy.s = table.clone(solut.s)
    cpy.cost = solut.cost

    return cpy
end

function subseq_fill(seq)
    for i=1,dimension+1 do
        seq[i] = {}
        for j=1,dimension+1 do
            seq[i][j] = {0, 0, 0}
        end
    end
end

function subseq_load(s, seq)
    for i=1,dimension+1 do
        local k = 1 - i - (i ~= 1 and 0 or 1)

        seq[i][i][T] = 0.0
        seq[i][i][C] = 0.0
        seq[i][i][W] = (i ~= 1 and 1 or 0)
        for j=i+1,dimension+1 do
            local j_prev = j-1

            seq[i][j][T] = c[s[j_prev]][s[j]] + seq[i][j_prev][T]
            seq[i][j][C] = seq[i][j][T] + seq[i][j_prev][C]
            seq[i][j][W] = j + k
        end
    end

end

function sort(arr, r)
    quicksort(arr, 1, #arr, r)
end

function quicksort(arr, left, right, r)
    if left < right then
        local pivotIndex = partition(arr, left, right, r)
        quicksort(arr, left, pivotIndex - 1, r)
        quicksort(arr, pivotIndex + 1, right, r)
    end
end

function partition(arr, left, right, r)
    local pivotValue = arr[right]
    local i = left - 1

    for j = left, right - 1 do
        if c[r][arr[j]] < c[r][pivotValue] then
            i = i + 1
            arr[i], arr[j] = arr[j], arr[i] -- Troca os valores
        end
    end

    arr[i + 1], arr[right] = arr[right], arr[i + 1] -- Troca o pivô
    return i + 1
end


function construction(alpha, rnd) 
    local s = {1}

    local cList = {}
    for i=2,dimension do
        cList[#cList + 1] = i
    end

    local r = 1

    while #cList > 0 do
        sort(cList, r)

        local i = rnd.rnd[rnd.rnd_index] + 1
        rnd.rnd_index = rnd.rnd_index + 1

        local cN = table.remove(cList, i)
        table.insert(s, cN)
        r = cN
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

function search_swap(s, seq)
    local cost_best = math.huge
    local I = -1
    local J = -1

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_concat_4 = 0.0
    local cost_new = 0.0

    for i = 2, dimension-1 do
        local i_prev = i - 1
        local i_next = i + 1

        cost_concat_1 =                 seq[1][i_prev][ T] + c[s[i_prev]][s[i_next]]
        cost_concat_2 = cost_concat_1 + seq[i][i_next][T] + c[s[i]][s[i_next+1]]

        cost_new = seq[1][i_prev][C]                                                    +           --       1st subseq
        seq[i][i_next][W]               * (cost_concat_1) + c[s[i_next]][s[i]]  +           -- concat 2nd subseq
        seq[i_next+1][dimension+1][W]   * (cost_concat_2) + seq[i_next+1][dimension+1][C]   -- concat 3rd subseq

        if cost_new < cost_best then
            cost_best = cost_new - EPSILON
            I = i
            J = i_next
        end

        for j = i_next+1,dimension do
            local j_prev = j-1
            local j_next = j+1


            cost_concat_1 =                 seq[1][i_prev][T]       + c[s[i_prev]][s[j]]
            cost_concat_2 = cost_concat_1                           + c[s[j]][s[i_next]]
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev][T]  + c[s[j_prev]][s[i]]
            cost_concat_4 = cost_concat_3                           + c[s[i]][s[j_next]]


            cost_new = seq[1][i_prev][C]                                                 +      -- 1st subseq
            cost_concat_1 +                                                             -- concat 2nd subseq (single node)
            seq[i_next][j_prev][W]      * cost_concat_2 + seq[i_next][j_prev][C] +      -- concat 3rd subseq
            cost_concat_3 +                                                             -- concat 4th subseq (single node)
            seq[j_next][dimension+1][W] * cost_concat_4 + seq[j_next][dimension+1][C]   -- concat 5th subseq

            if(cost_new < cost_best) then
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
            end

        end
    end

    if cost_best < seq[1][dimension+1][C] - EPSILON then
        swap(s, I, J)
        subseq_load(s, seq)
        return true
    end

    return false
end

function search_two_opt(s, seq)
    local cost_best = math.huge
    local I = -1
    local J = -1

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_new = 0.0

    for i = 2,dimension-1 do
        local i_prev = i - 1
        local rev_seq_cost = seq[i][i+1][T]
        for j = i+2,dimension do
            local j_next = j+1

            rev_seq_cost = rev_seq_cost + c[s[j-1]][s[j]] * (seq[i][j][W]-1.0)

            cost_concat_1 =                 seq[1][i_prev][T]   + c[s[j]][s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[i][j][T]        + c[s[j_next]][s[i]]

            cost_new = seq[1][i_prev][C]                                                        +   --   1st subseq
                    seq[i][j][W]                * cost_concat_1 + rev_seq_cost                  +   -- concat 2nd subseq (reversed seq)
                    seq[j_next][dimension+1][W] * cost_concat_2 + seq[j_next][dimension+1][C]       -- concat 3rd subseq

            if cost_new < cost_best then
                cost_best = cost_new -EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < seq[1][dimension+1][C] - EPSILON then
        reverse(s, I, J)
        subseq_load(s, seq)
        return true
    end

    return false
end

function search_reinsertion(s, seq, opt)
    local cost_best = math.huge
    local I = -1
    local J = -1
    local POS = -1

    local cost_concat_1 = 0.0
    local cost_concat_2 = 0.0
    local cost_concat_3 = 0.0
    local cost_new = 0.0

    for i = 2, dimension-opt+1 do
        local j = opt+i-1
        local i_prev = i-1
        local j_next = j+1

        for k = 1, i_prev-1 do
            local k_next = k+1

            cost_concat_1 =                 seq[1][k][T]            + c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[i][j][T]            + c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][T]  + c[s[i_prev]][s[j_next]];

            cost_new = seq[1][k][C]                                                             +   --       1st subseq
            seq[i][j][W]                * cost_concat_1 + seq[i][j][C]                  +   -- concat 2nd subseq (reinserted seq)
            seq[k_next][i_prev][W]      * cost_concat_2 + seq[k_next][i_prev][C]        +   -- concat 3rd subseq
            seq[j_next][dimension+1][W] * cost_concat_3 + seq[j_next][dimension+1][C];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - EPSILON
                I = i
                J = j
                POS = k
            end
        end

        for k = i+opt,dimension do
            local k_next = k+1

            cost_concat_1 =                 seq[1][i_prev][T]   + c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][k][T]   + c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[i][j][T]        + c[s[j]][s[k_next]];

            cost_new = seq[1][i_prev][C]                                                        +   --       1st subseq
                    seq[j_next][k][W]           * cost_concat_1 + seq[j_next][k][C]             +   -- concat 2nd subseq
                    seq[i][j][W]                * cost_concat_2 + seq[i][j][C]                  +   -- concat 3rd subseq (reinserted seq)
                    seq[k_next][dimension+1][W] * cost_concat_3 + seq[k_next][dimension+1][C];       -- concat 4th subseq

            if cost_new < cost_best then
                cost_best = cost_new - EPSILON
                I = i
                J = j
                POS = k
            end

        end
    end

    if cost_best < seq[1][dimension+1][C] - EPSILON then
        reinsert(s, I, J, POS+1)
        subseq_load(s, seq)

        if cost_best ~= seq[1][dimension+1][C] then
            print("ERROR")
            os.exit(1)
        end
        return true
    end

    return false
end

function RVND(s, seq, rnd)
    local SWAP        = 0  
    local REINSERTION = 1
    local OR_OPT_2    = 2
    local OR_OPT_3    = 3
    local TWO_OPT     = 4

    local neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}

    while #neighbd_list > 0 do
        local index = rnd.rnd[rnd.rnd_index] + 1
        rnd.rnd_index = rnd.rnd_index + 1

        local neighbd = neighbd_list[index]

        local improve = false

        if neighbd == SWAP then
            improve = search_swap(s, seq)
        elseif neighbd == REINSERTION then
            improve = search_reinsertion(s, seq, REINSERTION)
        elseif neighbd == OR_OPT_2 then
            improve = search_reinsertion(s, seq, OR_OPT_2)
        elseif neighbd == OR_OPT_3 then
            improve = search_reinsertion(s, seq, OR_OPT_3)
        elseif neighbd == TWO_OPT then
            improve = search_two_opt(s, seq)
        end

        if improve == true then
            neighbd_list = {SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3}
        else 
            table.remove(neighbd_list, index)
        end
        
    end
end

function perturb(sl, rnd)
    local s = table.clone(sl)

    local A_start = 1
    local A_end = 1
    local B_start = 1
    local B_end = 1

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do

        A_start = rnd.rnd[rnd.rnd_index] + 1
        rnd.rnd_index = rnd.rnd_index +1
        A_end = A_start + rnd.rnd[rnd.rnd_index]
        rnd.rnd_index = rnd.rnd_index +1

        B_start = rnd.rnd[rnd.rnd_index] + 1
        rnd.rnd_index = rnd.rnd_index +1
        B_end = B_start + rnd.rnd[rnd.rnd_index]
        rnd.rnd_index = rnd.rnd_index +1
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

function GILS_RVND(Imax, Iils, R, rnd)

    local solut_partial = {};
    local solut_crnt = {};
    local solut_best = {};

    local cost_partial = 0.0;
    local cost_crnt = 0.0;
    local cost_best = 0.0;

    local seq = {};

    subseq_fill(seq)

    cost_best = math.huge

    for i=1,Imax do
        local alpha = R[rnd.rnd[rnd.rnd_index] + 1]
        rnd.rnd_index = rnd.rnd_index + 1

        print("[+] Local Search", i)
        print("\t[+] Constructing Inital Solution..")
        solut_crnt = construction(alpha, rnd)

        subseq_load(solut_crnt, seq)
        cost_crnt = seq[1][dimension+1][C]
        s_print(solut_crnt)
        print("\tConstruction cost  ", cost_crnt)
        cost_partial = cost_crnt
        solut_partial = table.clone(solut_crnt)

        print("\t[+] Looking for the best Neighbor..")
        local iterILS = 0
        while iterILS < Iils do
            RVND(solut_crnt, seq, rnd)
            cost_crnt = seq[1][dimension+1][C]
            if cost_crnt < cost_partial - EPSILON then
               cost_partial = cost_crnt - EPSILON
               solut_partial = table.clone(solut_crnt)
               iterILS = 0
            end
            
            solut_crnt = perturb(solut_partial, rnd)
            subseq_load(solut_crnt, seq)

            iterILS = iterILS + 1
        end

        subseq_load(solut_partial, seq)
        if cost_partial < cost_best then
            solut_best = table.clone(solut_partial)
            cost_best = cost_partial
        end
        print("\tCurrent best solution cost: ", cost_best)
    end

    print("COST: ", cost_best)
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


    local rnd = {
        rnd = {}, 
        rnd_index = 1
    }
    dimension = readData(c, rnd.rnd)

    local Imax = 10
    local Iils = math.min(100, dimension)
    local R = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
               0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25}

    local start = os.clock()
    GILS_RVND(Imax, Iils, R, rnd)

    print("TIME: ", os.clock()-start)
end

main()
