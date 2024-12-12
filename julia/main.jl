#! /usr/bin/julia
using Printf
include("Data.jl")

mutable struct tData
    cost::Matrix{Float64}
    dimen::Int

    T::Int
    W::Int
    C::Int

    REINSERTION::Int
    OR_OPT2::Int
    OR_OPT3::Int
    TWO_OPT::Int
    SWAP::Int

    EPSILON::Float64

    reinsert_count::Int

    rnd::Array{Int, 1}
    rnd_index::Int
end

mutable struct tSeqInfo 
    T::Float64
    C::Float64
    W::Float64
end

mutable struct tSolution
    s::Array{Int, 1}
    #seq::Array{tSeqInfo, 2}
    seq::Array{Float64, 3}
    cost::Float64
end

function subseq_load(solut::tSolution, data::tData)
    for i in 1:data.dimen+1
        k::Int = 1 - i - (i==1)

         solut.seq[data.T, i, i] = 0.0
         solut.seq[data.C, i, i] = 0.0
         solut.seq[data.W, i, i] = i!=1

        for j in i+1:data.dimen+1
            j_prev = j-1

            solut.seq[data.T, j, i] = data.cost[solut.s[j_prev], solut.s[j]] + solut.seq[data.T, j_prev, i]
            solut.seq[data.C, j, i] = solut.seq[data.T, j, i] + solut.seq[data.C, j_prev, i]
            solut.seq[data.W, j, i] = j + k

        end
    end
    solut.cost = solut.seq[data.C, data.dimen+1, 1]

    return
end

function sort(arr::Array{Int}, r::Int, data::tData)
    quicksort(arr, 1, length(arr), data, r)
end

function quicksort(arr::Array{Int}, left::Int, right::Int, data::tData, r::Int)
    if left < right
        pivot = partition(arr, left, right, data, r)
        quicksort(arr, left, pivot - 1, data, r)
        quicksort(arr, pivot + 1, right, data, r)
    end
end

function partition(arr::Array{Int}, left::Int, right::Int, data::tData, r::Int)
    pivot = arr[right]
    i = left - 1
    for j in left:right-1
        if data.cost[r, arr[j]] < data.cost[r, pivot]
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1
end

function construction(alpha::Float64, data::tData)::Array{Int, 1}
    s = [1]
    cList = [2:data.dimen;]

    r = 1
    while length(cList) > 0
        sort(cList, r, data)

        r_value = data.rnd[data.rnd_index] + 1
        data.rnd_index += 1

        cN = cList[r_value]

        push!(s, cN)
        r = cN
        deleteat!(cList, findfirst(x-> x==cN, cList))

    end

    push!(s, 1)
    return s
end

function swap(s::Array{Int,1}, i::Int, j::Int)
    s[i], s[j] = s[j], s[i]
end

# func from https://github.com/JuliaLang/julia/blob/master/base/array.jl
_deleteat!(a::Vector, i::Integer, delta::Integer) =
    ccall(:jl_array_del_at, Cvoid, (Any, Int, UInt), a, i - 1, delta)

function reinsert(s::Array{Int,1}, i::Int, j::Int, pos::Int)
    sz = j - i + 1
    if i < pos
        splice!(s, pos:pos-1, s[i:j]) 
        _deleteat!(s, i, sz) # substituir por splice
    else
        splice!(s, pos:pos-1, s[i:j]) 
        _deleteat!(s, i+sz, sz)
    end
end

function search_swap(solut::tSolution, data::tData)::Bool

    cost_best = Inf
    I = -1
    J = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_concat_4 = 0.0
    cost_new = 0.0

    for i in 2:data.dimen-1
            i_prev::Int = i - 1
            i_next::Int = i + 1

            cost_concat_1 =                 solut.seq[data.T, i_prev, 1] + data.cost[solut.s[i_prev], solut.s[i_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[data.T, i_next, i] + data.cost[solut.s[i], solut.s[i_next+1]]

            cost_new = solut.seq[data.C, i_prev, 1]                                                    +           #       1st subseq
            solut.seq[data.W, i_next, i]               * (cost_concat_1) + data.cost[solut.s[i_next], solut.s[i]]  +           # concat 2nd subseq
            solut.seq[data.W, data.dimen+1, i_next+1]   * (cost_concat_2) + solut.seq[data.C, data.dimen+1, i_next+1]   # concat 3rd subseq

            if cost_new < cost_best
                cost_best = cost_new - data.EPSILON
                I = i
                J = i_next
            end

            for j in i_next+1:data.dimen
                j_prev = j-1
                j_next = j+1


                cost_concat_1 =                 solut.seq[data.T, i_prev, 1]       + data.cost[solut.s[i_prev], solut.s[j]]
                cost_concat_2 = cost_concat_1                           + data.cost[solut.s[j], solut.s[i_next]]
                cost_concat_3 = cost_concat_2 + solut.seq[data.T, j_prev, i_next]  + data.cost[solut.s[j_prev], solut.s[i]]
                cost_concat_4 = cost_concat_3                           + data.cost[solut.s[i], solut.s[j_next]]


                cost_new = solut.seq[data.C, i_prev, 1]                                                 +      # 1st subseq
                        cost_concat_1 +                                                             # concat 2nd subseq (single node)
                        solut.seq[data.W, j_prev, i_next]      * cost_concat_2 + solut.seq[data.C, j_prev, i_next] +      # concat 3rd subseq
                        cost_concat_3 +                                                             # concat 4th subseq (single node)
                        solut.seq[data.W, data.dimen+1, j_next] * cost_concat_4 + solut.seq[data.C, data.dimen+1, j_next]   # concat 5th subseq

                if(cost_new < cost_best)
                    cost_best = cost_new - data.EPSILON;
                    I = i;
                    J = j;
                end

        end
    end

    if cost_best < solut.cost - data.EPSILON
        swap(solut.s, I, J)
        subseq_load(solut, data)
        return true
    end

    return false
end

function search_two_opt(solut::tSolution, data::tData)::Bool
    cost_best = Inf
    I = -1
    J = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_new = 0.0

    for i in 2:data.dimen-1
        i_prev = i - 1
        rev_seq_cost = solut.seq[data.T, i+1, i]
        for j in i+2:data.dimen
            j_next = j+1

            rev_seq_cost += data.cost[solut.s[j-1], solut.s[j]] * (solut.seq[data.W, j, i]-1.0)

            cost_concat_1 =                 solut.seq[data.T, i_prev, 1]   + data.cost[solut.s[j], solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[data.T, j, i]        + data.cost[solut.s[j_next], solut.s[i]]

            cost_new = solut.seq[data.C, i_prev, 1]                                                        +   #   1st subseq
                    solut.seq[data.W, j, i]                * cost_concat_1 + rev_seq_cost                  +   # concat 2nd subseq (reversed seq)
                    solut.seq[data.W, data.dimen+1, j_next] * cost_concat_2 + solut.seq[data.C, data.dimen+1, j_next]       # concat 3rd subseq

            if cost_new < cost_best
                cost_best = cost_new - data.EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - data.EPSILON
        reverse!(solut.s, I, J)
        subseq_load(solut, data)
        return true
    end

    return false
end

function search_reinsertion(solut::tSolution, data::tData, opt::Int)::Bool
    data.reinsert_count += 1
    cost_best = Inf
    I = -1
    J = -1
    POS = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_new = 0.0

    for i in 2:data.dimen-opt+1
        j = opt+i-1
        i_prev = i-1
        j_next = j+1

        for k in 1:i_prev-1
                k_next = k+1

                cost_concat_1 =                 solut.seq[data.T, k, 1]            + data.cost[solut.s[k], solut.s[i]]
                cost_concat_2 = cost_concat_1 + solut.seq[data.T, j, i]            + data.cost[solut.s[j], solut.s[k_next]]
                cost_concat_3 = cost_concat_2 + solut.seq[data.T, i_prev, k_next]  + data.cost[solut.s[i_prev], solut.s[j_next]]

                cost_new = solut.seq[data.C, k, 1]                                                             +   #       1st subseq
                solut.seq[data.W, j, i]                * cost_concat_1 + solut.seq[data.C, j, i]                  +   # concat 2nd subseq (reinserted seq)
                solut.seq[data.W, i_prev, k_next]      * cost_concat_2 + solut.seq[data.C, i_prev, k_next]        +   # concat 3rd subseq
                solut.seq[data.W, data.dimen+1, j_next] * cost_concat_3 + solut.seq[data.C, data.dimen+1, j_next]       # concat 4th subseq

                if cost_new < cost_best
                    cost_best = cost_new - data.EPSILON
                    I = i
                    J = j
                    POS = k
                end

        end

        for k in i+opt:data.dimen
                k_next = k+1

                cost_concat_1 =                 solut.seq[data.T, i_prev, 1]   + data.cost[solut.s[i_prev], solut.s[j_next]]
                cost_concat_2 = cost_concat_1 + solut.seq[data.T, k, j_next]   + data.cost[solut.s[k], solut.s[i]]
                cost_concat_3 = cost_concat_2 + solut.seq[data.T, j, i]        + data.cost[solut.s[j], solut.s[k_next]]

                #println(i_prev, "  ", solut.seq[i_prev, 1].C)
                cost_new = solut.seq[data.C, i_prev, 1]                                                        +   #       1st subseq
                solut.seq[data.W, k, j_next]           * cost_concat_1 + solut.seq[data.C, k, j_next]             +   # concat 2nd subseq
                solut.seq[data.W, j, i]                * cost_concat_2 + solut.seq[data.C, j, i]                  +   # concat 3rd subseq (reinserted seq)
                solut.seq[data.W, data.dimen+1, k_next] * cost_concat_3 + solut.seq[data.C, data.dimen+1, k_next]       # concat 4th subseq

                if cost_new < cost_best
                    cost_best = cost_new - data.EPSILON
                    I = i
                    J = j
                    POS = k
                end

        end
    end

    if cost_best < solut.cost - data.EPSILON
        reinsert(solut.s, I, J, POS+1)
        subseq_load(solut, data)
        return true
    end

    return false
end


function RVND(solut::tSolution, data::tData)
    neighbd_list = [data.SWAP, data.TWO_OPT, data.REINSERTION, data.OR_OPT2, data.OR_OPT3]
    """
    t_reinsertion_local = 0
    t_or_opt2_local = 0
    t_or_opt3_local = 0
    t_two_opt_local = 0
    t_swap_local = 0
    """

    while length(neighbd_list) > 0
        
        i = data.rnd[data.rnd_index] + 1
        data.rnd_index += 1
        
        neighbd = neighbd_list[i]

        improve = false

        if neighbd == data.REINSERTION
            improve = search_reinsertion(solut, data, data.REINSERTION)
        elseif neighbd == data.OR_OPT2
            improve = search_reinsertion(solut, data, data.OR_OPT2)
        elseif neighbd == data.OR_OPT3
            improve = search_reinsertion(solut, data, data.OR_OPT3)
        elseif neighbd == data.SWAP
            improve = search_swap(solut, data)
        elseif neighbd == data.TWO_OPT
            improve = search_two_opt(solut, data)
        end

        if improve
            neighbd_list = [data.SWAP, data.TWO_OPT, data.REINSERTION, data.OR_OPT2, data.OR_OPT3]
        else
            deleteat!(neighbd_list, i)
        end

    end

    return
end

function perturb(solut_partial::tSolution, data::tData)::Array{Int64, 1}
    s = copy(solut_partial.s)

    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    size_max = Int(floor(length(s)/10))
    size_max = (size_max >= 2 ? size_max : 2)
    size_min = 2

    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)

        A_start = data.rnd[data.rnd_index] + 1
        data.rnd_index += 1
        A_end = A_start + data.rnd[data.rnd_index]
        data.rnd_index += 1

        B_start = data.rnd[data.rnd_index] + 1
        data.rnd_index += 1
        B_end = B_start + data.rnd[data.rnd_index]
        data.rnd_index += 1
    end

    if A_start < B_start
        reinsert(s, B_start, B_end-1, A_end)
        reinsert(s, A_start, A_end-1, B_end)
    else
        reinsert(s, A_start, A_end-1, B_end)
        reinsert(s, B_start, B_end-1, A_end)
    end

    return s
end

function seq_init(seq::Array{tSeqInfo, 2}, data::tData)
    for i in 1:data.dimen+1
        for j in i:data.dimen+1
            seq[j, i] = tSeqInfo(0.0, 0.0, 0.0)
        end
    end
end

function GILS_RVND(Imax::Int, Iils::Int, R::Vector{Float64}, data::tData)

    solut_best::tSolution = tSolution(zeros(Int, data.dimen+1), Array{Float64, 3}(undef, 3, data.dimen+1, data.dimen+1), Inf)
    solut_partial::tSolution = tSolution(zeros(Int, data.dimen+1), Array{Float64, 3}(undef, 3, data.dimen+1, data.dimen+1), 0)
    solut_crnt::tSolution = tSolution(zeros(Int, data.dimen+1), Array{Float64, 3}(undef, 3, data.dimen+1, data.dimen+1), 0)

    for i in 1:Imax
        r_value = data.rnd[data.rnd_index] + 1
        data.rnd_index += 1

        alpha = R[r_value]

        @printf "[+] Local Search %d\n" i
        solut_crnt.s = construction(alpha, data)
        subseq_load(solut_crnt, data)
        @printf "\t[+] Constructing Inital Solution.. %.2lf\n" solut_crnt.cost

        solut_partial.cost = solut_crnt.cost
        solut_partial.s = copy(solut_crnt.s)

        @printf "\t[+] Looking for the best Neighbor..\n"
        iterILS = 0
        while iterILS < Iils
            RVND(solut_crnt, data)
            if solut_crnt.cost < solut_partial.cost - data.EPSILON
                solut_partial.cost = solut_crnt.cost
                solut_partial.s = copy(solut_crnt.s)
                iterILS = 0
            end

            solut_crnt.s = perturb(solut_partial, data)
            subseq_load(solut_crnt, data)

            iterILS += 1
        end

        subseq_load(solut_partial, data)

        if solut_partial.cost < solut_best.cost
            solut_best.cost = solut_partial.cost
            solut_best.s = copy(solut_partial.s)
        end

        @printf "\tCurrent best solution cost: %.2lf\n" solut_best.cost
    end
    @printf "COST: %.2lf\n" solut_best.cost
    print("SOLUTION: ")
    println(solut_best.s)
end

function main()

    dimension, cost, rand_values = get_instance_info()


    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
         0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25] 

    Imax = 10
    Iils = min(dimension, 100)

    data::tData = tData(cost, dimension, 1, 2, 3, 1, 2, 3, 4, 5, 1e-15, 0, rand_values, 1)

    time = @elapsed GILS_RVND(Imax, Iils, R, data)

    @printf "TIME: %.6lf\n" time
    println("reinsertion calls ", data.reinsert_count)

end

main()
