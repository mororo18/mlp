#! /usr/bin/julia
include("Data.jl")

const T = 1
const C = 2
const W = 3

const EPSILON = 1e-15

const REINSERTION = 1
const OR_OPT2     = 2
const OR_OPT3     = 3
const SWAP        = 4
const TWO_OPT     = 5

dimension, c = get_instance_info()

improv_flag = true

function subseq_load(s::Array{Int64, 1}, seq::Array{Float64, 3})

    for i in 1:dimension+1
        k = 1 - i -(i==0)#convert(Int64, i==0)

        seq[i,i,T] = 0.0
        seq[i,i,C] = 0.0
        seq[i,i,W] = i!=0#convert(Float64, i != 0)

        for j in i+1:dimension+1
            j_prev = j-1
            
            seq[i,j,T] = c[s[j_prev], s[j]] + seq[i,j_prev,T]
            seq[i,j,C] = seq[i,j,T] + seq[i,j_prev,C]
            seq[i,j,W] = j+k
        end
    end

end

function construction(alpha::Float64)
    s = [1]
    cList = [2:dimension;]

    r = 1
    while length(cList) > 0
        sort!(cList, by= i -> c[i,r])

        i = convert(Int64, floor(length(cList)*alpha + 1))
        cN = cList[rand(1:i)]
        push!(s, cN)
        r = cN
        deleteat!(cList, findfirst(x->x==cN, cList))

        println(cList)
    end

    push!(s, 1)
    return s
end

function swap(s, i, j)
    s[i], s[j] = s[j], s[i]
end

function reinsert(s, i, j, pos)
    # func removed from https://github.com/JuliaLang/julia/blob/master/base/array.jl
    _deleteat!(a::Vector, i::Integer, delta::Integer) =
        ccall(:jl_array_del_at, Cvoid, (Any, Int, UInt), a, i - 1, delta)

    sz = j - i + 1
    if i < pos
        #insert!(s, pos, s[i:j])
        splice!(s, pos:pos-1, s[i:j]) 
        _deleteat!(s, i, sz)
    else
        #insert!(s, pos, s[i:j])
        splice!(s, pos:pos-1, s[i:j]) 
        _deleteat!(s, i+sz, sz)
    end
end

function search_swap(s::Array{Int64,1}, seq::Array{Float64,3})
    cost_best = Inf
    I = -1
    J = -1

    for i in 2:dimension-1
        i_prev = i - 1
        i_next = i + 1

        cost_concat_1 = seq[1, i_prev, T] + c[s[i_prev], s[i_next]]
        cost_concat_2 = cost_concat_1 + seq[i, i_next, T] + c[s[i], s[i_next]]


        cost_new = seq[1, i_prev, C]
                + seq[i, i_next, W]             * (cost_concat_1) + c[s[i_next], s[i]]
                + seq[i_next+1, dimension+1, W]   * (cost_concat_2) + seq[i_next+1, dimension+1, C]

        if cost_new < cost_best
            cost_best = cost_new - EPSILON
            I = i
            J = i_next
        end

        for j in i_next+1:dimension
            j_prev = j-1
            j_next = j+1


            cost_concat_1 = seq[1, i_prev, T] + c[s[i_prev], s[j]]
            cost_concat_2 = cost_concat_1 + c[s[j], s[i_next]]
            cost_concat_3 = cost_concat_2 + seq[i_next, j_prev, T] + c[s[j_prev], s[i]]
            cost_concat_4 = cost_concat_3 + c[s[i], s[j_next]]


            cost_new = seq[1, i_prev, C]
                    + cost_concat_1
                    + seq[i_next, j_prev, W] * cost_concat_2 + seq[i_next, j_prev, C]
                    + cost_concat_3
                    + seq[j_next, dimension+1, W] * cost_concat_4 + seq[j_next, dimension+1, C]

            if(cost_new < cost_best)
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
            end
        end
    end

    if cost_best < seq[1,dimension+1,C] - EPSILON
        swap(s, I, J)
        subseq_load(s, seq)
        improv_flag = true
    end
end

function search_two_opt(s, seq)
    cost_best = Inf
    I = -1
    J = -1

    for i in 2:dimension-1
        i_prev = i - 1
        rev_seq_cost = seq[i, i+1, T]
        for j in i+2:dimension
            j_next = j+1

            rev_seq_cost += c[s[j-1], s[j]] * (seq[i, j, W]-1.0)

            cost_concat_1 = seq[1, i_prev, T] + c[s[j], s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[i, j, T] + c[s[j_next], s[i]]

            cost_new = seq[i, i_prev, C]
                    + seq[i, j, W]              * cost_concat_1 + rev_seq_cost
                    + seq[j_next, dimension+1, W] * cost_concat_2 + seq[j_next, dimension+1, C]

            if cost_new < cost_best
                cost_best = cost_new - EPSILON
                I = i
                J = j
            end
        end

    end

    if cost_best < seq[1, dimension+1, C] - EPSILON
        reverse!(s, I, J)
        subseq_load(s, seq)
        improv_flag = true
        println(s)
    end
end

function search_reinsertion(s, seq, opt)
    cost_best = Inf
    I = -1
    J = -1
    POS = -1

    for i in 2:dimension-opt+1
        j = opt+1-1
        i_prev = i-1
        j_next = j+1

        for k in 1:i_prev
            k_next = k+1

            cost_concat_1 = seq[1, k, T] + c[s[k], s[i]]
            cost_concat_2 = cost_concat_1 + seq[i, j, T] + c[s[j], s[k_next]]
            cost_concat_3 = cost_concat_2 + seq[k_next, i_prev, T] + c[s[i_prev], s[j_next]]

            cost_new = seq[1, k, C]
                    + seq[i, j, W]              * cost_concat_1 + seq[i, j, C]
                    + seq[k_next, i_prev, W]    * cost_concat_2 + seq[k_next, i_prev, C]
                    + seq[j_next, dimension+1, W] * cost_concat_3 + seq[j_next, dimension+1, C]

            if cost_new < cost_best
                cost_best = cost_new - EPSILON
                I = i
                J = j
                POS = k
            end

        end

        for k in i+opt+1:dimension-opt-1
            k_next = k+1

            cost_concat_1 = seq[1, i_prev, T] + c[s[i_prev], s[j_next]]
            cost_concat_2 = cost_concat_1 + seq[j_next, k, T] + c[s[k], s[i]]
            cost_concat_3 = cost_concat_2 + seq[i, j, T] + c[s[j], s[k_next]]

            cost_new = seq[1, i_prev, C]
                    + seq[j_next, k, W]         * cost_concat_1 + seq[j_next, k, C]
                    + seq[i, j, W]              * cost_concat_2 + seq[i, j, C]
                    + seq[k_next, dimension+1, W] * cost_concat_3 + seq[k_next, dimension+1, C]

            if cost_new < cost_best
                cost_best = cost_new - EPSILON
                I = i
                J = j
                POS = k
            end

        end

    end

    if cost_best < seq[1, dimension+1, C]
        reinsert(s, I, J, POS+1)
        subseq_load(s, seq)
        improv_flag = true
        println(s)
    end
end

function RVND(s::Array{Int64, 1}, seq::Array{Float64, 3})
    #neighbd_list = [SWAP]
    #neighbd_list = [TWO_OPT]
    #neighbd_list = [OR_OPT2]
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3]
    println(s)

    while length(neighbd_list) > 0
        i = rand(1:length(neighbd_list))
        neighbd = neighbd_list[i]

        improv_flag = false

        if neighbd == REINSERTION
            search_reinsertion(s, seq, REINSERTION)
        elseif neighbd == OR_OPT2
            search_reinsertion(s, seq, OR_OPT2)
        elseif neighbd == OR_OPT3
            search_reinsertion(s, seq, OR_OPT3)
        elseif neighbd == SWAP
            search_swap(s, seq)
        elseif neighbd == TWO_OPT
            search_two_opt(s, seq)
        end

        if improv_flag
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3]
        else
            deleteat!(neighbd_list, i)
        end

        println(neighbd_list)
        #exit(0)
    end
end

function perturb(sl::Array{Int64, 1})
    s = copy(sl)

    A_start, A_end = 1, 1
    B_start, B_end = 1, 1

    size_max = convert(Int64, floor(length(s)/10))
    size_max = (size_max >= 2 ? size_max : 2)
    size_min = 2

    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)
        A_start = rand(2:length(s)-1-size_max)
        A_end = A_start + rand(size_min:size_max)

        B_start = rand(2:length(s)-1-size_max)
        B_end = B_start + rand(size_min:size_max)
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

function GILS_RVND(Imax::Int64, Iils::Int64, R)
    cost_best = Inf
    s_best = []

    subseq = zeros(dimension+1, dimension+1, 3)

    """
    s = Array{Int64, 1}()
    for i in 1:dimension
        push!(s, i)
    end
    push!(s, 1)
    println(s)

    subseq_load(s, subseq)

    println(subseq)
    println(subseq[1, dimension+1, C])
    """
    for i in 1:Imax
        alpha = R[rand(1:26)]
        s = construction(alpha)
        sl = copy(s)
        subseq_load(s, subseq)

        rvnd_cost_best = subseq[1,dimension+1,C] - EPSILON

        iterILS = 0
        while iterILS < Iils
            RVND(s, subseq)
            rvnd_cost_crnt = subseq[1,dimension+1,C] - EPSILON
            if rvnd_cost_crnt < rvnd_cost_best
                rvnd_cost_best = rvnd_cost_crnt
                sl = copy(s)
                iterILS = 0
            end

            s = perturb(sl)
            subseq_load(s, subseq)

            iterILS += 1
        end

        subseq_load(sl, subseq)
        sl_cost = subseq[1,dimension+1,C] - EPSILON

        if sl_cost < cost_best
            s_best = sl
            cost_best = sl_cost
        end
    end
    println(cost_best)
end

function main()
    println(dimension)
    println(c)

    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
         0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25] 
    
    Imax = 10
    Iils = min(dimension, 100)

    GILS_RVND(Imax, Iils, R)

end

main()
