#
# Minimum Latency Problem solver -- Julia port of main.cpp
#
# All node indices follow the original C++ code (0-based). Accessor helpers
# (`dist`, `node`, `seq_at`) do the +1 translation into Julia's 1-based arrays.
#

using Printf

include("Data.jl")

const TABLE_SZ = 26

const DBL_MAX = floatmax(Float64)
const DBL_EPSILON = eps(Float64)

# Neighborhoods
const REINSERTION = 1
const OR_OPT_2    = 2
const OR_OPT_3    = 3
const SWAP        = 4
const TWO_OPT     = 5

struct tSeqInfo
    T::Float64
    C::Float64
    W::Float64
end

mutable struct tData
    cost::Matrix{Float64}
    dimen::Int
    rnd::Vector{Int}
    rnd_index::Int      # 1-based cursor over rnd
end

mutable struct tSolution
    s::Vector{Int}
    seq::Vector{tSeqInfo}   # flat (dimen+1)x(dimen+1)
    cost::Float64
end

@inline to_1D(x::Int, y::Int, size::Int) = x * (size + 1) + y + 1

@inline dist(data::tData, a::Int, b::Int) = @inbounds data.cost[a+1, b+1]
@inline node(solut::tSolution, i::Int) = @inbounds solut.s[i+1]
@inline seq_at(solut::tSolution, i::Int, j::Int, size::Int) = @inbounds solut.seq[to_1D(i, j, size)]

function Solution_init(data::tData)
    s = zeros(Int, data.dimen + 1)
    seq = fill(tSeqInfo(0.0, 0.0, 0.0), (data.dimen + 1) * (data.dimen + 1))
    return tSolution(s, seq, DBL_MAX)
end

function Solution_cpy!(src::tSolution, tgt::tSolution)
    resize!(tgt.s, length(src.s))
    copyto!(tgt.s, src.s)
    tgt.cost = src.cost
    return nothing
end

function R_table(i::Int)
    table = (0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
             0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25)
    return table[i+1]
end

function print_s(s::Vector{Int})
    for i in eachindex(s)
        print(s[i] + 1, " ")
    end
    println()
end

function partition!(arr::Vector{Int}, left::Int, right::Int, data::tData, r::Int)
    pivot = arr[right]
    i = left - 1
    for j in left:right-1
        if dist(data, r, arr[j]) < dist(data, r, pivot)
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
    arr[i+1], arr[right] = arr[right], arr[i+1]
    return i + 1
end

function quicksort!(arr::Vector{Int}, left::Int, right::Int, data::tData, r::Int)
    if left < right
        pivot = partition!(arr, left, right, data, r)
        quicksort!(arr, left, pivot - 1, data, r)
        quicksort!(arr, pivot + 1, right, data, r)
    end
    return nothing
end

function sort_cl!(arr::Vector{Int}, r::Int, data::tData)
    quicksort!(arr, 1, length(arr), data, r)
    return nothing
end

function construct(alpha::Float64, data::tData)
    s = Int[0]
    sizehint!(s, data.dimen + 1)
    cL = collect(1:data.dimen-1)

    r = 0
    while !isempty(cL)
        sort_cl!(cL, r, data)

        index = data.rnd[data.rnd_index]
        data.rnd_index += 1

        c = cL[index+1]
        push!(s, c)
        r = c
        deleteat!(cL, index + 1)
    end

    push!(s, 0)

    return s
end

function swap!(vec::Vector{Int}, i::Int, j::Int)
    vec[i+1], vec[j+1] = vec[j+1], vec[i+1]
    return nothing
end

function reverse_seq!(vec::Vector{Int}, i::Int, j::Int)
    reverse!(vec, i + 1, j + 1)
    return nothing
end

function reinsert!(vec::Vector{Int}, i::Int, j::Int, pos::Int)
    seq = vec[i+1:j+1]
    if pos < i
        deleteat!(vec, i+1:j+1)
        splice!(vec, pos+1:pos, seq)
    else
        splice!(vec, pos+1:pos, seq)
        deleteat!(vec, i+1:j+1)
    end
    return nothing
end

function update_subseq_info_matrix!(solut::tSolution, data::tData, index::Int = 0)
    dimen = data.dimen
    from = index

    for i in 0:dimen
        k = 1 - i - (i == 0 ? 1 : 0)
        t = (i == from)

        @inbounds solut.seq[to_1D(i, i, dimen)] = tSeqInfo(0.0, 0.0, i == 0 ? 0.0 : 1.0)

        for j in i+1:dimen
            j_prev = j - 1

            prev = seq_at(solut, i, j_prev, dimen)
            T = dist(data, node(solut, j_prev), node(solut, j)) + prev.T
            C = T + prev.C
            W = Float64(j + k)

            @inbounds solut.seq[to_1D(i, j, dimen)] = tSeqInfo(T, C, W)
        end

        from += t ? 1 : 0
    end

    solut.cost = seq_at(solut, 0, dimen, dimen).C
    return nothing
end

function search_swap(solut::tSolution, data::tData)
    dimen = data.dimen
    cost_best = DBL_MAX
    I = 0
    J = 0

    for i in 1:dimen-2
        i_prev = i - 1
        i_next = i + 1

        # consecutive nodes
        cost_concat_1 =                 seq_at(solut, 0, i_prev, dimen).T + dist(data, node(solut, i_prev), node(solut, i_next))
        cost_concat_2 = cost_concat_1 + seq_at(solut, i, i_next, dimen).T + dist(data, node(solut, i), node(solut, i_next+1))

        cost_new = seq_at(solut, 0, i_prev, dimen).C                                                                       +   #       1st subseq
            seq_at(solut, i, i_next, dimen).W        * cost_concat_1 + dist(data, node(solut, i_next), node(solut, i))      +   # concat 2nd subseq
            seq_at(solut, i_next+1, dimen, dimen).W  * cost_concat_2 + seq_at(solut, i_next+1, dimen, dimen).C                  # concat 3rd subseq

        if cost_new < cost_best
            cost_best = cost_new - DBL_EPSILON
            I = i
            J = i_next
        end

        for j in i_next+1:dimen-1
            j_next = j + 1
            j_prev = j - 1

            cost_concat_1 =                 seq_at(solut, 0, i_prev, dimen).T      + dist(data, node(solut, i_prev), node(solut, j))
            cost_concat_2 = cost_concat_1                                          + dist(data, node(solut, j), node(solut, i_next))
            cost_concat_3 = cost_concat_2 + seq_at(solut, i_next, j_prev, dimen).T + dist(data, node(solut, j_prev), node(solut, i))
            cost_concat_4 = cost_concat_3                                          + dist(data, node(solut, i), node(solut, j_next))

            cost_new = seq_at(solut, 0, i_prev, dimen).C                                                    +   # 1st subseq
                cost_concat_1                                                                               +   # concat 2nd subseq (single node)
                seq_at(solut, i_next, j_prev, dimen).W * cost_concat_2 + seq_at(solut, i_next, j_prev, dimen).C +   # concat 3rd subseq
                cost_concat_3                                                                               +   # concat 4th subseq (single node)
                seq_at(solut, j_next, dimen, dimen).W  * cost_concat_4 + seq_at(solut, j_next, dimen, dimen).C     # concat 5th subseq

            if cost_new < cost_best
                cost_best = cost_new - DBL_EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - DBL_EPSILON
        swap!(solut.s, I, J)
        update_subseq_info_matrix!(solut, data, I)
        return true
    end

    return false
end

function search_two_opt(solut::tSolution, data::tData)
    dimen = data.dimen
    cost_best = DBL_MAX
    I = 0
    J = 0

    for i in 1:dimen-2
        i_prev = i - 1

        rev_seq_cost = seq_at(solut, i, i+1, dimen).T
        for j in i+2:dimen-1
            j_next = j + 1

            rev_seq_cost += dist(data, node(solut, j-1), node(solut, j)) * (seq_at(solut, i, j, dimen).W - 1.0)

            cost_concat_1 =                 seq_at(solut, 0, i_prev, dimen).T + dist(data, node(solut, j), node(solut, i_prev))
            cost_concat_2 = cost_concat_1 + seq_at(solut, i, j, dimen).T      + dist(data, node(solut, j_next), node(solut, i))

            cost_new = seq_at(solut, 0, i_prev, dimen).C                                                  +   #  1st subseq
                seq_at(solut, i, j, dimen).W          * cost_concat_1 + rev_seq_cost                      +   # concat 2nd subseq (reversed seq)
                seq_at(solut, j_next, dimen, dimen).W * cost_concat_2 + seq_at(solut, j_next, dimen, dimen).C # concat 3rd subseq

            if cost_new < cost_best
                cost_best = cost_new - DBL_EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - DBL_EPSILON
        reverse_seq!(solut.s, I, J)
        update_subseq_info_matrix!(solut, data)
        return true
    end

    return false
end

function search_reinsertion(solut::tSolution, data::tData, opt::Int)
    dimen = data.dimen
    cost_best = DBL_MAX
    I = 0
    J = 0
    POS = 0

    for i in 1:dimen-opt
        j = opt + i - 1
        j_next = j + 1
        i_prev = i - 1

        # k -> edges
        for k in 0:i_prev-1
            k_next = k + 1

            cost_concat_1 =                 seq_at(solut, 0, k, dimen).T             + dist(data, node(solut, k), node(solut, i))
            cost_concat_2 = cost_concat_1 + seq_at(solut, i, j, dimen).T             + dist(data, node(solut, j), node(solut, k_next))
            cost_concat_3 = cost_concat_2 + seq_at(solut, k_next, i_prev, dimen).T   + dist(data, node(solut, i_prev), node(solut, j_next))

            cost_new = seq_at(solut, 0, k, dimen).C                                                                  +   #       1st subseq
                seq_at(solut, i, j, dimen).W           * cost_concat_1 + seq_at(solut, i, j, dimen).C                 +   #  concat 2nd subseq (reinserted seq)
                seq_at(solut, k_next, i_prev, dimen).W * cost_concat_2 + seq_at(solut, k_next, i_prev, dimen).C       +   #  concat 3rd subseq
                seq_at(solut, j_next, dimen, dimen).W  * cost_concat_3 + seq_at(solut, j_next, dimen, dimen).C            # concat 4th subseq

            if cost_new < cost_best
                cost_best = cost_new - DBL_EPSILON
                I = i
                J = j
                POS = k
            end
        end

        for k in i+opt:dimen-1
            k_next = k + 1

            cost_concat_1 =                 seq_at(solut, 0, i_prev, dimen).T   + dist(data, node(solut, i_prev), node(solut, j_next))
            cost_concat_2 = cost_concat_1 + seq_at(solut, j_next, k, dimen).T   + dist(data, node(solut, k), node(solut, i))
            cost_concat_3 = cost_concat_2 + seq_at(solut, i, j, dimen).T        + dist(data, node(solut, j), node(solut, k_next))

            cost_new = seq_at(solut, 0, i_prev, dimen).C                                                             +   #       1st subseq
                seq_at(solut, j_next, k, dimen).W     * cost_concat_1 + seq_at(solut, j_next, k, dimen).C             +   # concat 2nd subseq
                seq_at(solut, i, j, dimen).W          * cost_concat_2 + seq_at(solut, i, j, dimen).C                  +   # concat 3rd subseq (reinserted seq)
                seq_at(solut, k_next, dimen, dimen).W * cost_concat_3 + seq_at(solut, k_next, dimen, dimen).C             # concat 4th subseq

            if cost_new < cost_best
                cost_best = cost_new - DBL_EPSILON
                I = i
                J = j
                POS = k
            end
        end
    end

    if cost_best < solut.cost - DBL_EPSILON
        reinsert!(solut.s, I, J, POS + 1)
        update_subseq_info_matrix!(solut, data, I < POS + 1 ? I : POS + 1)
        return true
    end

    return false
end

function RVND!(solut::tSolution, data::tData)
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]

    while !isempty(neighbd_list)

        index = data.rnd[data.rnd_index]
        data.rnd_index += 1
        neighbd = neighbd_list[index+1]
        improve_flag = false

        if neighbd == REINSERTION
            improve_flag = search_reinsertion(solut, data, REINSERTION)
        elseif neighbd == OR_OPT_2
            improve_flag = search_reinsertion(solut, data, OR_OPT_2)
        elseif neighbd == OR_OPT_3
            improve_flag = search_reinsertion(solut, data, OR_OPT_3)
        elseif neighbd == SWAP
            improve_flag = search_swap(solut, data)
        elseif neighbd == TWO_OPT
            improve_flag = search_two_opt(solut, data)
        end

        if improve_flag
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3]
        else
            deleteat!(neighbd_list, index + 1)
        end

    end

    return nothing
end

function perturb(solut::tSolution, data::tData)
    s = copy(solut.s)
    A_start = 1
    A_end = 1
    B_start = 1
    B_end = 1

    size_max = floor(Int, (data.dimen + 1) / 10)
    size_max = size_max >= 2 ? size_max : 2
    size_min = 2
    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)

        A_start = data.rnd[data.rnd_index]; data.rnd_index += 1
        A_end = A_start + data.rnd[data.rnd_index]; data.rnd_index += 1

        B_start = data.rnd[data.rnd_index]; data.rnd_index += 1
        B_end = B_start + data.rnd[data.rnd_index]; data.rnd_index += 1
    end

    if A_start < B_start
        reinsert!(s, B_start, B_end - 1, A_end)
        reinsert!(s, A_start, A_end - 1, B_end)
    else
        reinsert!(s, A_start, A_end - 1, B_end)
        reinsert!(s, B_start, B_end - 1, A_end)
    end

    return s
end

function GILS_RVND(Imax::Int, Iils::Int, data::tData, verbose::Bool)

    solut_partial = Solution_init(data)
    solut_crnt = Solution_init(data)
    solut_best = Solution_init(data)

    for i in 0:Imax-1
        aux = data.rnd[data.rnd_index]
        data.rnd_index += 1

        alpha = R_table(aux)

        if verbose
            @printf("[+] Search %d\n", i + 1)
            @printf("\t[+] Constructing..\n")
        end

        solut_crnt.s = construct(alpha, data)
        update_subseq_info_matrix!(solut_crnt, data)

        Solution_cpy!(solut_crnt, solut_partial)
        if verbose
            @printf("\t[+] Looking for the best Neighbor..\n")
            @printf("\t    Construction Cost: %.3f\n", solut_partial.cost)
        end

        iterILS = 0
        while iterILS < Iils
            RVND!(solut_crnt, data)
            if solut_crnt.cost < solut_partial.cost - DBL_EPSILON
                Solution_cpy!(solut_crnt, solut_partial)
                iterILS = 0
            end

            solut_crnt.s = perturb(solut_partial, data)
            update_subseq_info_matrix!(solut_crnt, data)
            iterILS += 1
        end

        if solut_partial.cost < solut_best.cost - DBL_EPSILON
            Solution_cpy!(solut_partial, solut_best)
        end

        if verbose
            @printf("\tCurrent best cost: %g\n", solut_best.cost)

            print("SOLUCAO: ")
            for k in eachindex(solut_best.s)
                print(solut_best.s[k], " ")
            end
            println()
        end

    end

    @printf("COST: %.2f\n", solut_best.cost)
    return nothing
end

function main(args::Vector{String})
    Imax = 10
    verbose = false

    for a in args
        if a == "-v" || a == "--verbose"
            verbose = true
        end
    end

    dimen, cost, rnd = load_data()

    data = tData(cost, dimen, rnd, 1)

    solut = Solution_init(data)

    Iils = data.dimen < 100 ? data.dimen : 100

    t1 = time_ns()
    GILS_RVND(Imax, Iils, data, verbose)
    t2 = time_ns()

    elapsed = (t2 - t1) / 1e9
    @printf("TIME: %g\n", elapsed)
    @printf("wall clock (s): %g\n", elapsed)

    println("Tamanho RND ", length(rnd))

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
