module GILS_RVND

#using Printf

const T = 1
const C = 2
const W = 3

const EPSILON = 1e-15

const REINSERTION = 1
const OR_OPT2     = 2
const OR_OPT3     = 3
const SWAP        = 4
const TWO_OPT     = 5

dimension = nothing
matrix = nothing

function get_instance_info()

    f = open("../distance_matrix")

    line = readline(f)
    dimension = parse(Int32, line)
    #matrix = Matrix{Float64}(undef, dimension,dimension)
    matrix = zeros(dimension, dimension)

    i = 1
    for line in eachline(f)
        #println(line)
        if line == "EOF"
            break
        end


        j = i+1
        while length(line) > 0
            index = findfirst(isequal(' '), line)
            d = parse(Float64, line[begin:index])
            #setindex!(matrix, d, i, j)
            matrix[i,j] = d
            matrix[j,i] = d

            line = line[index+1:end]

            j += 1
        end

        i += 1
    end

    close(f)

    return dimension, matrix
end

dimension, c = get_instance_info()

improv_flag = true
it = 0

t_swap = 0
t_reinsertion = 0
t_or_opt2 = 0
t_or_opt3 = 0
t_two_opt = 0
t_seq = 0

function subseq_load(s::Array{Int64, 1}, seq::Array{Float64, 3})

    global t_seq -= time_ns()
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

            """
            seq[j,i,T] = seq[i,j,T]
            seq[j,i,C] = seq[i,j,C]
            seq[j,i,W] = seq[i,j,W]
            """
        end
    end
    t_seq += time_ns()

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

    end

    push!(s, 1)
    return s
end

function swap(s::Array{Int64,1}, i::Int64, j::Int64)
    s[i], s[j] = s[j], s[i]
end

function reinsert(s::Array{Int64,1}, i::Int64, j::Int64, pos::Int64)
    # func from https://github.com/JuliaLang/julia/blob/master/base/array.jl
    _deleteat!(a::Vector, i::Integer, delta::Integer) =
        ccall(:jl_array_del_at, Cvoid, (Any, Int, UInt), a, i - 1, delta)

    sz = j - i + 1
    if i < pos
        splice!(s, pos:pos-1, s[i:j]) 
        _deleteat!(s, i, sz)
    else
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

        cost_concat_1 =                 seq[1, i_prev, T] + c[s[i_prev], s[i_next]]
        cost_concat_2 = cost_concat_1 + seq[i, i_next, T] + c[s[i], s[i_next+1]]

        cost_new = seq[1, i_prev, C]                                                    +           #       1st subseq
                seq[i, i_next, W]               * (cost_concat_1) + c[s[i_next], s[i]]  +           # concat 2nd subseq
                seq[i_next+1, dimension+1, W]   * (cost_concat_2) + seq[i_next+1, dimension+1, C]   # concat 3rd subseq

        if cost_new < cost_best
            cost_best = cost_new - EPSILON
            I = i
            J = i_next
        end

        for j in i_next+1:dimension
            j_prev = j-1
            j_next = j+1


            cost_concat_1 =                 seq[1, i_prev, T]       + c[s[i_prev], s[j]]
            cost_concat_2 = cost_concat_1                           + c[s[j], s[i_next]]
            cost_concat_3 = cost_concat_2 + seq[i_next, j_prev, T]  + c[s[j_prev], s[i]]
            cost_concat_4 = cost_concat_3                           + c[s[i], s[j_next]]


            cost_new = seq[1, i_prev, C]                                                 +      # 1st subseq
                    cost_concat_1 +                                                             # concat 2nd subseq (single node)
                    seq[i_next, j_prev, W]      * cost_concat_2 + seq[i_next, j_prev, C] +      # concat 3rd subseq
                    cost_concat_3 +                                                             # concat 4th subseq (single node)
                    seq[j_next, dimension+1, W] * cost_concat_4 + seq[j_next, dimension+1, C]   # concat 5th subseq

            if(cost_new < cost_best)
                """
                @printf "%.2lf %d %d\n" cost_new i j
                @printf "\tconcat 1 %.2lf\n" cost_concat_1
                @printf "\tconcat 2 %.2lf\n" cost_concat_2
                @printf "\tconcat 3 %.2lf\n" cost_concat_3
                @printf "\tconcat 4 %.2lf\n" cost_concat_4
                """
                cost_best = cost_new - EPSILON;
                I = i;
                J = j;
            end
        end
    end

    if cost_best < seq[1,dimension+1,C] - EPSILON
        swap(s, I, J)
        subseq_load(s, seq)
        global improv_flag = true
    end
end

function search_two_opt(s::Array{Int64,1}, seq::Array{Float64,3})
    cost_best = Inf
    I = -1
    J = -1

    for i in 2:dimension-1
        i_prev = i - 1
        rev_seq_cost = seq[i, i+1, T]
        for j in i+2:dimension
            j_next = j+1

            rev_seq_cost += c[s[j-1], s[j]] * (seq[i, j, W]-1.0)

            cost_concat_1 =                 seq[1, i_prev, T]   + c[s[j], s[i_prev]]
            cost_concat_2 = cost_concat_1 + seq[i, j, T]        + c[s[j_next], s[i]]

            cost_new = seq[1, i_prev, C]                                                        +   #   1st subseq
            seq[i, j, W]                        * cost_concat_1 + rev_seq_cost                  +   # concat 2nd subseq (reversed seq)
                    seq[j_next, dimension+1, W] * cost_concat_2 + seq[j_next, dimension+1, C]       # concat 3rd subseq

            if cost_new < cost_best
                """
                @printf "%.2lf %d %d\n" cost_new i j
                @printf "\treve_seq %.2lf\n" rev_seq_cost
                @printf "\tconcat 1 %.2lf\n" cost_concat_1
                @printf "\tconcat 2 %.2lf\n" cost_concat_2
                @printf "\ttest %.2lf\n" seq[i, i_prev, C]
                """
                cost_best = cost_new - EPSILON
                I = i
                J = j
            end
        end

    end

    if cost_best < seq[1, dimension+1, C] - EPSILON
        reverse!(s, I, J)
        subseq_load(s, seq)
        global improv_flag = true
    end
end

function search_reinsertion(s::Array{Int64,1}, seq::Array{Float64,3}, opt::Int64)
    cost_best = Inf
    I = -1
    J = -1
    POS = -1

    for i in 2:dimension-opt+1
        j = opt+i-1
        i_prev = i-1
        j_next = j+1

        for k in 1:i_prev-1
            k_next = k+1

            cost_concat_1 =                 seq[1, k, T]            + c[s[k], s[i]]
            cost_concat_2 = cost_concat_1 + seq[i, j, T]            + c[s[j], s[k_next]]
            cost_concat_3 = cost_concat_2 + seq[k_next, i_prev, T]  + c[s[i_prev], s[j_next]]

            cost_new = seq[1, k, C]                                                             +   #       1st subseq
                    seq[i, j, W]                * cost_concat_1 + seq[i, j, C]                  +   # concat 2nd subseq (reinserted seq)
                    seq[k_next, i_prev, W]      * cost_concat_2 + seq[k_next, i_prev, C]        +   # concat 3rd subseq
                    seq[j_next, dimension+1, W] * cost_concat_3 + seq[j_next, dimension+1, C]       # concat 4th subseq

            if cost_new < cost_best
                cost_best = cost_new - EPSILON
                I = i
                J = j
                POS = k
            end
        end

        for k in i+opt:dimension-opt-1
            k_next = k+1

            cost_concat_1 = seq[1, i_prev, T] + c[s[i_prev], s[j_next]]
            cost_concat_2 = cost_concat_1 + seq[j_next, k, T] + c[s[k], s[i]]
            cost_concat_3 = cost_concat_2 + seq[i, j, T] + c[s[j], s[k_next]]

            cost_new = seq[1, i_prev, C]                                                        +   #       1st subseq
                    seq[j_next, k, W]           * cost_concat_1 + seq[j_next, k, C]             +   # concat 2nd subseq
                    seq[i, j, W]                * cost_concat_2 + seq[i, j, C]                  +   # concat 3rd subseq (reinserted seq)
                    seq[k_next, dimension+1, W] * cost_concat_3 + seq[k_next, dimension+1, C]       # concat 4th subseq

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
        global improv_flag = true
    end
end


function RVND(s::Array{Int64, 1}, seq::Array{Float64, 3})
    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3]
    t_reinsertion_local = 0
    t_or_opt2_local = 0
    t_or_opt3_local = 0
    t_two_opt_local = 0
    t_swap_local = 0

    while length(neighbd_list) > 0
        global it += 1
        i = rand(1:length(neighbd_list))
        neighbd = neighbd_list[i]

        global improv_flag = false

        if neighbd == REINSERTION
            t_reinsertion_local    += @elapsed search_reinsertion(s, seq, REINSERTION)
        elseif neighbd == OR_OPT2
            t_or_opt2_local        += @elapsed search_reinsertion(s, seq, OR_OPT2)
        elseif neighbd == OR_OPT3
            t_or_opt3_local        += @elapsed search_reinsertion(s, seq, OR_OPT3)
        elseif neighbd == SWAP
            t_swap_local           += @elapsed search_swap(s, seq)
        elseif neighbd == TWO_OPT
            t_two_opt_local        += @elapsed search_two_opt(s, seq)
        end

        if improv_flag
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT2, OR_OPT3]
        else
            deleteat!(neighbd_list, i)
        end

    end

    global t_swap +=t_swap_local 
    global t_reinsertion +=t_reinsertion_local
    global t_or_opt2 += t_or_opt2_local
    global t_or_opt3 += t_or_opt3_local
    global t_two_opt += t_two_opt_local
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

function main(Imax::Int64, Iils::Int64, R)
    cost_best = Inf
    s_best = []

    subseq = zeros(dimension+1, dimension+1, 3)

    for i in 1:Imax
        alpha = R[rand(1:26)]
        #@printf "[+] Local Search %d\n" i
        #@printf "\t[+] Constructing Inital Solution..\n"
        s = construction(alpha)
        sl = copy(s)
        subseq_load(s, subseq)

        rvnd_cost_best = subseq[1,dimension+1,C] - EPSILON

        #@printf "\t[+] Looking for the best Neighbor..\n"
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

        #@printf "\tCurrent best solution cost: %.2lf\n" cost_best
    end
    #@printf "COST: %.2lf\n" cost_best
    print("SOLUTION: ")
    println(s_best)
    println(cost_best)
    #@printf "Swap time: %.6lf\n" t_swap
    #@printf "Reinsertion time: %.6lf\n" t_reinsertion
    #@printf "Or_opt2 time: %.6lf\n" t_or_opt2
    #@printf "Or_opt3 time: %.6lf\n" t_or_opt3
    #@printf "Two_opt time: %.6lf\n" t_two_opt
    #@printf "subseq time: %.6lf\n" t_seq/1e9
end

function solve()
    R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
         0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25] 
    
    Imax = 10
    Iils = min(dimension, 100)

    time = @elapsed main(Imax, Iils, R)

    #@printf "TIME %.6lf\n" time
    #@printf "RVND iteracoes %d\n" it

end

precompile(solve, (Any,))

#main()

end # module
