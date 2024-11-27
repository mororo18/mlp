#! /usr/bin/julia
using Printf
#using Profile, PProf
#using BenchmarkTools
#using InteractiveUtils
include("Data.jl")

mutable struct tInfo
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

function subseq_load(solut::tSolution, info::tInfo) # dimension::Int, c::Matrix{Float64})
    for i in 1:info.dimen+1
        k::Int = 1 - i - (i==1)#convert(Int, i==0)

         solut.seq[info.T, i, i] = 0.0
         solut.seq[info.C, i, i] = 0.0
         solut.seq[info.W, i, i] = i!=1#convert(Float64, i != 0)

        for j in i+1:info.dimen+1
            j_prev = j-1

            solut.seq[info.T, j, i] = info.cost[solut.s[j_prev], solut.s[j]] + solut.seq[info.T, j_prev, i]
            solut.seq[info.C, j, i] = solut.seq[info.T, j, i] + solut.seq[info.C, j_prev, i]
            solut.seq[info.W, j, i] = j + k

        end
    end
    solut.cost = solut.seq[info.C, info.dimen+1, 1]

    return
end

function sort(arr::Array{Int}, r::Int, info::tInfo)
    quicksort(arr, 1, length(arr), info, r)
end

function quicksort(arr::Array{Int}, left::Int, right::Int, info::tInfo, r::Int)
    if left < right
        pivot = partition(arr, left, right, info, r)
        quicksort(arr, left, pivot - 1, info, r)
        quicksort(arr, pivot + 1, right, info, r)
    end
end

function partition(arr::Array{Int}, left::Int, right::Int, info::tInfo, r::Int)
    pivot = arr[right]
    i = left - 1
    for j in left:right-1
        if info.cost[r, arr[j]] < info.cost[r, pivot]
            i += 1
            arr[i], arr[j] = arr[j], arr[i]
        end
    end
    arr[i + 1], arr[right] = arr[right], arr[i + 1]
    return i + 1
end

function construction(alpha::Float64, info::tInfo)::Array{Int, 1}
    s = [1]
    cList = [2:info.dimen;]

    r = 1
    while length(cList) > 0
        sort(cList, r, info)

        r_value = info.rnd[info.rnd_index] + 1
        info.rnd_index += 1

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

function search_swap(solut::tSolution, info::tInfo)::Bool

    cost_best = Inf
    I = -1
    J = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_concat_4 = 0.0
    cost_new = 0.0

    for i in 2:info.dimen-1
            i_prev::Int = i - 1
            i_next::Int = i + 1

            cost_concat_1 =                 solut.seq[info.T, i_prev, 1] + info.cost[solut.s[i_prev], solut.s[i_next]]
            cost_concat_2 = cost_concat_1 + solut.seq[info.T, i_next, i] + info.cost[solut.s[i], solut.s[i_next+1]]

            cost_new = solut.seq[info.C, i_prev, 1]                                                    +           #       1st subseq
            solut.seq[info.W, i_next, i]               * (cost_concat_1) + info.cost[solut.s[i_next], solut.s[i]]  +           # concat 2nd subseq
            solut.seq[info.W, info.dimen+1, i_next+1]   * (cost_concat_2) + solut.seq[info.C, info.dimen+1, i_next+1]   # concat 3rd subseq

            if cost_new < cost_best
                cost_best = cost_new - info.EPSILON
                I = i
                J = i_next
            end

            for j in i_next+1:info.dimen
                j_prev = j-1
                j_next = j+1


                cost_concat_1 =                 solut.seq[info.T, i_prev, 1]       + info.cost[solut.s[i_prev], solut.s[j]]
                cost_concat_2 = cost_concat_1                           + info.cost[solut.s[j], solut.s[i_next]]
                cost_concat_3 = cost_concat_2 + solut.seq[info.T, j_prev, i_next]  + info.cost[solut.s[j_prev], solut.s[i]]
                cost_concat_4 = cost_concat_3                           + info.cost[solut.s[i], solut.s[j_next]]


                cost_new = solut.seq[info.C, i_prev, 1]                                                 +      # 1st subseq
                        cost_concat_1 +                                                             # concat 2nd subseq (single node)
                        solut.seq[info.W, j_prev, i_next]      * cost_concat_2 + solut.seq[info.C, j_prev, i_next] +      # concat 3rd subseq
                        cost_concat_3 +                                                             # concat 4th subseq (single node)
                        solut.seq[info.W, info.dimen+1, j_next] * cost_concat_4 + solut.seq[info.C, info.dimen+1, j_next]   # concat 5th subseq

                if(cost_new < cost_best)
                    cost_best = cost_new - info.EPSILON;
                    I = i;
                    J = j;
                end

        end
    end

    if cost_best < solut.cost - info.EPSILON
        swap(solut.s, I, J)
        subseq_load(solut, info)
        return true
    end

    return false
end

function search_two_opt(solut::tSolution, info::tInfo)::Bool
    cost_best = Inf
    I = -1
    J = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_new = 0.0

    for i in 2:info.dimen-1
        i_prev = i - 1
        rev_seq_cost = solut.seq[info.T, i+1, i]
        for j in i+2:info.dimen
            j_next = j+1

            rev_seq_cost += info.cost[solut.s[j-1], solut.s[j]] * (solut.seq[info.W, j, i]-1.0)

            cost_concat_1 =                 solut.seq[info.T, i_prev, 1]   + info.cost[solut.s[j], solut.s[i_prev]]
            cost_concat_2 = cost_concat_1 + solut.seq[info.T, j, i]        + info.cost[solut.s[j_next], solut.s[i]]

            cost_new = solut.seq[info.C, i_prev, 1]                                                        +   #   1st subseq
                    solut.seq[info.W, j, i]                * cost_concat_1 + rev_seq_cost                  +   # concat 2nd subseq (reversed seq)
                    solut.seq[info.W, info.dimen+1, j_next] * cost_concat_2 + solut.seq[info.C, info.dimen+1, j_next]       # concat 3rd subseq

            if cost_new < cost_best
                cost_best = cost_new - info.EPSILON
                I = i
                J = j
            end
        end
    end

    if cost_best < solut.cost - info.EPSILON
        reverse!(solut.s, I, J)
        subseq_load(solut, info)
        return true
    end

    return false
end

function search_reinsertion(solut::tSolution, info::tInfo, opt::Int)::Bool
    info.reinsert_count += 1
    cost_best = Inf
    I = -1
    J = -1
    POS = -1

    cost_concat_1 = 0.0
    cost_concat_2 = 0.0
    cost_concat_3 = 0.0
    cost_new = 0.0

    for i in 2:info.dimen-opt+1
        j = opt+i-1
        i_prev = i-1
        j_next = j+1

        for k in 1:i_prev-1
                k_next = k+1

                cost_concat_1 =                 solut.seq[info.T, k, 1]            + info.cost[solut.s[k], solut.s[i]]
                cost_concat_2 = cost_concat_1 + solut.seq[info.T, j, i]            + info.cost[solut.s[j], solut.s[k_next]]
                cost_concat_3 = cost_concat_2 + solut.seq[info.T, i_prev, k_next]  + info.cost[solut.s[i_prev], solut.s[j_next]]

                cost_new = solut.seq[info.C, k, 1]                                                             +   #       1st subseq
                solut.seq[info.W, j, i]                * cost_concat_1 + solut.seq[info.C, j, i]                  +   # concat 2nd subseq (reinserted seq)
                solut.seq[info.W, i_prev, k_next]      * cost_concat_2 + solut.seq[info.C, i_prev, k_next]        +   # concat 3rd subseq
                solut.seq[info.W, info.dimen+1, j_next] * cost_concat_3 + solut.seq[info.C, info.dimen+1, j_next]       # concat 4th subseq

                if cost_new < cost_best
                    cost_best = cost_new - info.EPSILON
                    I = i
                    J = j
                    POS = k
                end

        end

        for k in i+opt:info.dimen
                k_next = k+1

                cost_concat_1 =                 solut.seq[info.T, i_prev, 1]   + info.cost[solut.s[i_prev], solut.s[j_next]]
                cost_concat_2 = cost_concat_1 + solut.seq[info.T, k, j_next]   + info.cost[solut.s[k], solut.s[i]]
                cost_concat_3 = cost_concat_2 + solut.seq[info.T, j, i]        + info.cost[solut.s[j], solut.s[k_next]]

                #println(i_prev, "  ", solut.seq[i_prev, 1].C)
                cost_new = solut.seq[info.C, i_prev, 1]                                                        +   #       1st subseq
                solut.seq[info.W, k, j_next]           * cost_concat_1 + solut.seq[info.C, k, j_next]             +   # concat 2nd subseq
                solut.seq[info.W, j, i]                * cost_concat_2 + solut.seq[info.C, j, i]                  +   # concat 3rd subseq (reinserted seq)
                solut.seq[info.W, info.dimen+1, k_next] * cost_concat_3 + solut.seq[info.C, info.dimen+1, k_next]       # concat 4th subseq

                if cost_new < cost_best
                    cost_best = cost_new - info.EPSILON
                    I = i
                    J = j
                    POS = k
                end

        end
    end

    if cost_best < solut.cost - info.EPSILON
        reinsert(solut.s, I, J, POS+1)
        subseq_load(solut, info)
        return true
    end

    return false
end


function RVND(solut::tSolution, info::tInfo)
    neighbd_list = [info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT2, info.OR_OPT3]
    """
    t_reinsertion_local = 0
    t_or_opt2_local = 0
    t_or_opt3_local = 0
    t_two_opt_local = 0
    t_swap_local = 0
    """

    while length(neighbd_list) > 0
        
        i = info.rnd[info.rnd_index] + 1
        info.rnd_index += 1
        
        neighbd = neighbd_list[i]

        improve = false

        if neighbd == info.REINSERTION
            improve = search_reinsertion(solut, info, info.REINSERTION)
        elseif neighbd == info.OR_OPT2
            improve = search_reinsertion(solut, info, info.OR_OPT2)
        elseif neighbd == info.OR_OPT3
            improve = search_reinsertion(solut, info, info.OR_OPT3)
        elseif neighbd == info.SWAP
            improve = search_swap(solut, info)
        elseif neighbd == info.TWO_OPT
            improve = search_two_opt(solut, info)
        end

        if improve
            neighbd_list = [info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT2, info.OR_OPT3]
        else
            deleteat!(neighbd_list, i)
        end

    end

    return
end

function perturb(solut_partial::tSolution, info::tInfo)::Array{Int64, 1}
    s = copy(solut_partial.s)

    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    size_max = Int(floor(length(s)/10))
    size_max = (size_max >= 2 ? size_max : 2)
    size_min = 2

    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)

        A_start = info.rnd[info.rnd_index] + 1
        info.rnd_index += 1
        A_end = A_start + info.rnd[info.rnd_index]
        info.rnd_index += 1

        B_start = info.rnd[info.rnd_index] + 1
        info.rnd_index += 1
        B_end = B_start + info.rnd[info.rnd_index]
        info.rnd_index += 1
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

function seq_init(seq::Array{tSeqInfo, 2}, info::tInfo)
    for i in 1:info.dimen+1
        for j in i:info.dimen+1
            seq[j, i] = tSeqInfo(0.0, 0.0, 0.0)
        end
    end
end

function GILS_RVND(Imax::Int, Iils::Int, R::Vector{Float64}, info::tInfo)

    solut_best::tSolution = tSolution(zeros(Int, info.dimen+1), Array{Float64, 3}(undef, 3, info.dimen+1, info.dimen+1), Inf)
    solut_partial::tSolution = tSolution(zeros(Int, info.dimen+1), Array{Float64, 3}(undef, 3, info.dimen+1, info.dimen+1), 0)
    solut_crnt::tSolution = tSolution(zeros(Int, info.dimen+1), Array{Float64, 3}(undef, 3, info.dimen+1, info.dimen+1), 0)

    for i in 1:Imax
        r_value = info.rnd[info.rnd_index] + 1
        info.rnd_index += 1

        alpha = R[r_value]

        @printf "[+] Local Search %d\n" i
        solut_crnt.s = construction(alpha, info)
        subseq_load(solut_crnt, info)
        @printf "\t[+] Constructing Inital Solution.. %.2lf\n" solut_crnt.cost

        solut_partial.cost = solut_crnt.cost
        solut_partial.s = copy(solut_crnt.s)

        @printf "\t[+] Looking for the best Neighbor..\n"
        iterILS = 0
        while iterILS < Iils
            RVND(solut_crnt, info)
            if solut_crnt.cost < solut_partial.cost - info.EPSILON
                solut_partial.cost = solut_crnt.cost
                solut_partial.s = copy(solut_crnt.s)
                iterILS = 0
            end

            solut_crnt.s = perturb(solut_partial, info)
            subseq_load(solut_crnt, info)

            iterILS += 1
        end

        subseq_load(solut_partial, info)

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

    info::tInfo = tInfo(cost, dimension, 1, 2, 3, 1, 2, 3, 4, 5, 1e-15, 0, rand_values, 1)

    time = @elapsed GILS_RVND(Imax, Iils, R, info)

    @printf "TIME: %.6lf\n" time
    println("reinsertion calls ", info.reinsert_count)

end

main()
