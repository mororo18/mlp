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


# func removed from https://github.com/JuliaLang/julia/blob/master/base/array.jl
_deleteat!(a::Vector, i::Integer, delta::Integer) =
    ccall(:jl_array_del_at, Cvoid, (Any, Int, UInt), a, i - 1, delta)

function reinsert(s, i, j, pos)
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

function search_swap(s, seq)
end

function search_two_opt(s, seq)

    reverse!(s, 2, 6)
    println(s)
end

function search_reinsertion(s, seq, opt)
    reinsert(s, 2,4,9)
    println(s)
end

function RVND(s::Array{Int64, 1}, seq::Array{Float64, 3})
    #neighbd_list = [SWAP]
    #neighbd_list = [TWO_OPT]
    neighbd_list = [OR_OPT2]
    #neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3]
    println(s)

    while length(neighbd_list) > 0
        i = rand(1:length(neighbd_list))
        neighbd = neighbd_list[i]

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
        exit(0)
    end
end

function perturb(sl::Array{Int64, 1})
    s = copy(sl)
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
