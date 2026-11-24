-- main.lua
-- Translated from main.cpp
--
-- The C++ code indexes cities and solution positions from 0, and its
-- tSeqInfo array is a flat array addressed via to_1D(i,j,size). This
-- translation keeps that same 0-based indexing scheme in Lua tables
-- (rather than shifting everything to 1-based) so the algorithm's logic
-- maps back to the original as directly as possible.

local Data = require("Data")

local EPS = 2.220446049250313e-16 -- DBL_EPSILON

local REINSERTION, OR_OPT_2, OR_OPT_3, SWAP, TWO_OPT = 1, 2, 3, 4, 5

local function to_1D(x, y, size)
    return x * (size + 1) + y
end

local function Solution_init(data)
    local solut = {}
    solut.s = {}
    for i = 0, data.dimen do solut.s[i] = 0 end
    solut.s_n = data.dimen + 1

    solut.seq = {}
    for idx = 0, (data.dimen + 1) * (data.dimen + 1) - 1 do
        solut.seq[idx] = { T = 0.0, C = 0.0, W = 0.0 }
    end

    solut.cost = math.huge
    return solut
end

local function Solution_cpy(src, tgt)
    tgt.s = {}
    for i = 0, src.s_n - 1 do tgt.s[i] = src.s[i] end
    tgt.s_n = src.s_n
    tgt.cost = src.cost
end

local R_TABLE = {
    0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
    0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25,
}

local function R_table(i)
    return R_TABLE[i + 1]
end

local function print_s(s, s_n)
    local out = {}
    for i = 0, s_n - 1 do
        out[#out + 1] = tostring(s[i] + 1)
    end
    print(table.concat(out, " "))
end

local function partition(arr, left, right, data, r)
    local pivot = arr[right]
    local i = left - 1
    for j = left, right - 1 do
        if data.cost[r][arr[j]] < data.cost[r][pivot] then
            i = i + 1
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1
end

local function quicksort(arr, left, right, data, r)
    if left < right then
        local p = partition(arr, left, right, data, r)
        quicksort(arr, left, p - 1, data, r)
        quicksort(arr, p + 1, right, data, r)
    end
end

local function sort_cl(arr, r, data)
    quicksort(arr, 1, #arr, data, r)
end

local function construct(alpha, data)
    local s = {}
    local s_n = 0
    s[s_n] = 0
    s_n = s_n + 1

    local cL = {}
    for i = 1, data.dimen - 1 do
        cL[i] = i
    end

    local r = 0
    while #cL > 0 do
        sort_cl(cL, r, data)

        local index = data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1

        local c = cL[index + 1]
        s[s_n] = c
        s_n = s_n + 1
        r = c
        table.remove(cL, index + 1)
    end

    s[s_n] = 0
    s_n = s_n + 1

    return s, s_n
end

local function vec_swap(vec, i, j)
    vec[i], vec[j] = vec[j], vec[i]
end

local function vec_reverse(vec, i, j)
    while i < j do
        vec[i], vec[j] = vec[j], vec[i]
        i = i + 1
        j = j - 1
    end
end

local function vec_reinsert(vec, n, i, j, pos)
    local seq = {}
    local seq_n = 0
    for k = i, j do
        seq[seq_n] = vec[k]
        seq_n = seq_n + 1
    end

    local new_vec = {}
    local idx = 0

    if pos < i then
        for k = 0, pos - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
        for k = 0, seq_n - 1 do new_vec[idx] = seq[k]; idx = idx + 1 end
        for k = pos, i - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
        for k = j + 1, n - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
    else
        for k = 0, i - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
        for k = j + 1, pos - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
        for k = 0, seq_n - 1 do new_vec[idx] = seq[k]; idx = idx + 1 end
        for k = pos, n - 1 do new_vec[idx] = vec[k]; idx = idx + 1 end
    end

    for k = 0, n - 1 do vec[k] = new_vec[k] end
end

local function update_subseq_info_matrix(solut, data, index)
    index = index or 0
    local dimen = data.dimen
    local from = index

    for i = 0, dimen do
        local k = 1 - i - ((i == 0) and 1 or 0)
        local t = (i == from)

        local ii = to_1D(i, i, dimen)
        solut.seq[ii].T = 0.0
        solut.seq[ii].C = 0.0
        solut.seq[ii].W = (i == 0) and 0.0 or 1.0

        for j = i + 1, dimen do
            local j_prev = j - 1
            local idx = to_1D(i, j, dimen)
            local idx_prev = to_1D(i, j_prev, dimen)
            solut.seq[idx].T = data.cost[solut.s[j_prev]][solut.s[j]] + solut.seq[idx_prev].T
            solut.seq[idx].C = solut.seq[idx].T + solut.seq[idx_prev].C
            solut.seq[idx].W = j + k
        end

        if t then from = from + 1 end
    end

    solut.cost = solut.seq[to_1D(0, dimen, dimen)].C
end

local function search_swap(solut, data)
    local dimen = data.dimen
    local cost_best = math.huge
    local I, J

    for i = 1, dimen - 2 do
        local i_prev = i - 1
        local i_next = i + 1

        local cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[i_next]]
        local cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, dimen)].T + data.cost[solut.s[i]][solut.s[i_next + 1]]

        local cost_new = solut.seq[to_1D(0, i_prev, dimen)].C
            + solut.seq[to_1D(i, i_next, dimen)].W * cost_concat_1 + data.cost[solut.s[i_next]][solut.s[i]]
            + solut.seq[to_1D(i_next + 1, dimen, dimen)].W * cost_concat_2 + solut.seq[to_1D(i_next + 1, dimen, dimen)].C

        if cost_new < cost_best then
            cost_best = cost_new - EPS
            I = i
            J = i_next
        end

        for j = i_next + 1, dimen - 1 do
            local j_next = j + 1
            local j_prev = j - 1

            cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j]]
            cost_concat_2 = cost_concat_1 + data.cost[solut.s[j]][solut.s[i_next]]
            local cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, dimen)].T + data.cost[solut.s[j_prev]][solut.s[i]]
            local cost_concat_4 = cost_concat_3 + data.cost[solut.s[i]][solut.s[j_next]]

            cost_new = solut.seq[to_1D(0, i_prev, dimen)].C
                + cost_concat_1
                + solut.seq[to_1D(i_next, j_prev, dimen)].W * cost_concat_2 + solut.seq[to_1D(i_next, j_prev, dimen)].C
                + cost_concat_3
                + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_4 + solut.seq[to_1D(j_next, dimen, dimen)].C

            if cost_new < cost_best then
                cost_best = cost_new - EPS
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - EPS then
        vec_swap(solut.s, I, J)
        update_subseq_info_matrix(solut, data, I)
        return true
    end

    return false
end

local function search_two_opt(solut, data)
    local dimen = data.dimen
    local cost_best = math.huge
    local I, J

    for i = 1, dimen - 2 do
        local i_prev = i - 1
        local rev_seq_cost = solut.seq[to_1D(i, i + 1, dimen)].T

        for j = i + 2, dimen - 1 do
            local j_next = j + 1

            rev_seq_cost = rev_seq_cost + data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seq[to_1D(i, j, dimen)].W - 1.0)

            local cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[j]][solut.s[i_prev]]
            local cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j_next]][solut.s[i]]

            local cost_new = solut.seq[to_1D(0, i_prev, dimen)].C
                + solut.seq[to_1D(i, j, dimen)].W * cost_concat_1 + rev_seq_cost
                + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_2 + solut.seq[to_1D(j_next, dimen, dimen)].C

            if cost_new < cost_best then
                cost_best = cost_new - EPS
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - EPS then
        vec_reverse(solut.s, I, J)
        update_subseq_info_matrix(solut, data)
        return true
    end

    return false
end

local function search_reinsertion(solut, data, opt)
    local dimen = data.dimen
    local cost_best = math.huge
    local I, J, POS

    local i = 1
    local j = opt + i - 1
    while i < dimen - opt + 1 do
        local j_next = j + 1
        local i_prev = i - 1

        for k = 0, i_prev - 1 do
            local k_next = k + 1

            local cost_concat_1 = solut.seq[to_1D(0, k, dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            local cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]
            local cost_concat_3 = cost_concat_2 + solut.seq[to_1D(k_next, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]

            local cost_new = solut.seq[to_1D(0, k, dimen)].C
                + solut.seq[to_1D(i, j, dimen)].W * cost_concat_1 + solut.seq[to_1D(i, j, dimen)].C
                + solut.seq[to_1D(k_next, i_prev, dimen)].W * cost_concat_2 + solut.seq[to_1D(k_next, i_prev, dimen)].C
                + solut.seq[to_1D(j_next, dimen, dimen)].W * cost_concat_3 + solut.seq[to_1D(j_next, dimen, dimen)].C

            if cost_new < cost_best then
                cost_best = cost_new - EPS
                I = i
                J = j
                POS = k
            end
        end

        for k = i + opt, dimen - 1 do
            local k_next = k + 1

            local cost_concat_1 = solut.seq[to_1D(0, i_prev, dimen)].T + data.cost[solut.s[i_prev]][solut.s[j_next]]
            local cost_concat_2 = cost_concat_1 + solut.seq[to_1D(j_next, k, dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            local cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i, j, dimen)].T + data.cost[solut.s[j]][solut.s[k_next]]

            local cost_new = solut.seq[to_1D(0, i_prev, dimen)].C
                + solut.seq[to_1D(j_next, k, dimen)].W * cost_concat_1 + solut.seq[to_1D(j_next, k, dimen)].C
                + solut.seq[to_1D(i, j, dimen)].W * cost_concat_2 + solut.seq[to_1D(i, j, dimen)].C
                + solut.seq[to_1D(k_next, dimen, dimen)].W * cost_concat_3 + solut.seq[to_1D(k_next, dimen, dimen)].C

            if cost_new < cost_best then
                cost_best = cost_new - EPS
                I = i
                J = j
                POS = k
            end
        end

        i = i + 1
        j = j + 1
    end

    if cost_best < solut.cost - EPS then
        local start_idx = (I < POS + 1) and I or (POS + 1)
        vec_reinsert(solut.s, solut.s_n, I, J, POS + 1)
        update_subseq_info_matrix(solut, data, start_idx)
        return true
    end

    return false
end

local function RVND(solut, data)
    local neighbd_list = { SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3 }

    while #neighbd_list > 0 do
        local index = data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1

        local neighbd = neighbd_list[index + 1]
        local improve_flag = false

        if neighbd == REINSERTION then
            improve_flag = search_reinsertion(solut, data, REINSERTION)
        elseif neighbd == OR_OPT_2 then
            improve_flag = search_reinsertion(solut, data, OR_OPT_2)
        elseif neighbd == OR_OPT_3 then
            improve_flag = search_reinsertion(solut, data, OR_OPT_3)
        elseif neighbd == SWAP then
            improve_flag = search_swap(solut, data)
        elseif neighbd == TWO_OPT then
            improve_flag = search_two_opt(solut, data)
        end

        if improve_flag then
            neighbd_list = { SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3 }
        else
            table.remove(neighbd_list, index + 1)
        end
    end
end

local function perturb(solut, data)
    local s = {}
    for i = 0, solut.s_n - 1 do s[i] = solut.s[i] end

    local A_start, A_end, B_start, B_end = 1, 1, 1, 1

    local size_max = math.floor((data.dimen + 1) / 10)
    size_max = (size_max >= 2) and size_max or 2
    local size_min = 2 -- luacheck: ignore (unused, kept to mirror the C++ source)

    while (A_start <= B_start and B_start <= A_end) or (B_start <= A_start and A_start <= B_end) do
        A_start = data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1
        A_end = A_start + data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1

        B_start = data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1
        B_end = B_start + data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1
    end

    if A_start < B_start then
        vec_reinsert(s, solut.s_n, B_start, B_end - 1, A_end)
        vec_reinsert(s, solut.s_n, A_start, A_end - 1, B_end)
    else
        vec_reinsert(s, solut.s_n, A_start, A_end - 1, B_end)
        vec_reinsert(s, solut.s_n, B_start, B_end - 1, A_end)
    end

    return s
end

local function GILS_RVND(Imax, Iils, data, verbose)
    local solut_partial = Solution_init(data)
    local solut_crnt = Solution_init(data)
    local solut_best = Solution_init(data)

    for i = 0, Imax - 1 do
        local aux = data.rnd[data.rnd_index]
        data.rnd_index = data.rnd_index + 1

        local alpha = R_table(aux)

        if verbose then
            print(string.format("[+] Search %d", i + 1))
            print("\t[+] Constructing..")
        end

        local s, s_n = construct(alpha, data)
        solut_crnt.s = s
        solut_crnt.s_n = s_n
        update_subseq_info_matrix(solut_crnt, data)

        Solution_cpy(solut_crnt, solut_partial)
        if verbose then
            print("\t[+] Looking for the best Neighbor..")
            print(string.format("\t    Construction Cost: %.3f", solut_partial.cost))
        end

        local iterILS = 0
        while iterILS < Iils do
            RVND(solut_crnt, data)
            if solut_crnt.cost < solut_partial.cost - EPS then
                Solution_cpy(solut_crnt, solut_partial)
                iterILS = 0
            end

            solut_crnt.s = perturb(solut_partial, data)
            update_subseq_info_matrix(solut_crnt, data)
            iterILS = iterILS + 1
        end

        if solut_partial.cost < solut_best.cost - EPS then
            Solution_cpy(solut_partial, solut_best)
        end

        if verbose then
            print("\tCurrent best cost: " .. tostring(solut_best.cost))

            io.write("SOLUCAO: ")
            for k = 0, solut_best.s_n - 1 do
                io.write(tostring(solut_best.s[k]), " ")
            end
            io.write("\n")
        end
    end

    print(string.format("COST: %.2f", solut_best.cost))
end

local function main(argv)
    local Imax = 10
    local Iils
    local verbose = false

    for i = 1, #argv do
        if argv[i] == "-v" or argv[i] == "--verbose" then
            verbose = true
        end
    end

    local data = {}
    local dimen, cost, rnd, rnd_n = Data.loadData()
    data.dimen = dimen
    data.cost = cost
    data.rnd = rnd
    data.rnd_index = 0

    local solut = Solution_init(data) -- unused past this point, kept to mirror main.cpp

    math.randomseed(os.clock())

    Iils = (data.dimen < 100) and data.dimen or 100

    local t1 = os.clock()
    local cpu_start = os.clock()
    GILS_RVND(Imax, Iils, data, verbose)
    -- Lua's standard library only exposes CPU time via os.clock(); there is
    -- no separate high-resolution wall-clock source without a C extension,
    -- so both timings below are derived from it.
    local cpu_time = os.clock() - cpu_start
    local t2 = os.clock()

    local res = t2 - t1

    print("TIME: " .. tostring(cpu_time))
    print("wall clock (s): " .. tostring(res))

    print("Tamanho RND " .. tostring(rnd_n))

    os.exit(0)
end

main(arg)
