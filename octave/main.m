1;
mainn()

global dimen;
global cost;
global fmax;

global T;
global C;
global W;

global REINSERTION;
global OR_OPT_2;
global OR_OPT_3;
global SWAP;
global TWO_OPT;
    
function print_s(s)
    sz = size(s);
    sz = sz(2);

    for i = 1:sz
        fprintf('%d ', s(i))
    end
    fprintf('\n')
end

function ret = subseq_load (s)
global dimen;
global cost;
global T;
global C;
global W;
    tern1 = [1, 0];
    tern2 = [0, 1];

    seq = zeros(dimen+1, dimen+1, 3);
    for i = 1:dimen+1
        k = 1 - i - tern1((i ~= 1)  + 1);

        seq(i, i, T) = 0.0;
        seq(i, i, C) = 0.0;
        seq(i, i, W) = tern2((i ~= 1)  + 1);
        for j = i+1:dimen+1
            j_prev = j-1;
            s;
            seq(i, j, T) = cost(s(j_prev), s(j)) + seq(i, j_prev, T);
            seq(i, j, C) = seq(i, j, T) + seq(i, j_prev, C);
            seq(i, j, W) = j+k;

        end
    end

    ret = seq;
end

function ret = sort_by(arr, r)
    global cost
    sz = size(arr);
    sz = sz(2);

    for i = 1:sz
        for j = 1:sz-i
            if (cost(r, arr(j)) > cost(r, arr(j+1)))
                tmp = arr(j);
                arr(j) = arr(j+1);
                arr(j+1) = tmp;
            end
        end
    end

    ret = arr;
end

function [ret, index_new] = construction(alpha, rnd)
    global dimen
    s(1) = 1;
    for i = 1:dimen-1
        cL(i) = i+1;
    end

    r = 1;
    sz_cL = size(cL);
    sz_cL = sz_cL(2);
    while (sz_cL > 0) 
        cL = sort_by(cL, r );

        rng = ceil(sz_cL * alpha);
        RND = rand(1);
        %RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        index = ceil(rng * RND);
        
        index = rnd.rnd(rnd.rnd_index) + 1;
        %info.rnd(info.rnd_index);
        rnd.rnd_index = rnd.rnd_index+1;

        cN = cL(index);

        cL(index) = [];
        sz_s = size(s);
        sz_s = sz_s(2);
        s(sz_s+1) = cN;
        r = cN;
        sz_cL = sz_cL - 1;
    end

    sz_s = size(s);
    sz_s = sz_s(2);
    s(sz_s+1) = 1;
    ret = s;
    index_new = rnd.rnd_index;
end

function ret = swap(s, i, j)
    tmp = s(i);
    s(i) = s(j);
    s(j) = tmp;

    ret = s;
    end

function ret = reverse(s, i, j)
    s(i:j) = flip(s(i:j));
    ret = s;
end

function ret = reinsert(s, i, j, pos)
    if (i < pos)
        tmp = s(i:j);
        inter = s(j+1:pos-1);
        s(i:i+length(inter)-1) = inter;
        s(i+length(inter):pos-1) = tmp;
    else
        tmp = s(i:j);
        inter = s(pos:i-1);
        s(pos:pos+length(tmp)-1) = tmp;
        s(pos+length(tmp):j) = inter;
    end

    ret = s;
end

function [solut_new, seq_new, ret] = search_swap(s, seq)
    global fmax;
    global dimen;
    global cost;
    global T;
    global C;
    global W;

    cost_best = fmax;
    cost_new = 0.0;
    cost_concat_1 = 0.0;
    cost_concat_2 = 0.0;
    cost_concat_3 = 0.0;
    cost_concat_4 = 0.0;

    for i = 2:dimen-1
        i_prev = i-1;
        i_next = i+1;

        cost_concat_1 =                 seq(1, i_prev, T) + cost(s(i_prev), s(i_next));
        cost_concat_2 = cost_concat_1 + seq(i, i_next, T) + cost(s(i), s(i_next+1));

        cost_new = seq(1, i_prev, C)                                                    + ...
                seq(i, i_next, W)               * (cost_concat_1) + cost(s(i_next), s(i))  + ...
                seq(i_next+1, dimen+1, W)   * (cost_concat_2) + seq(i_next+1, dimen+1, C);

        if (cost_new < cost_best)
            cost_best = cost_new - eps;
            I_best = i;
            J_best = i_next;
        end

        for j = i_next+1:dimen

            j_prev = j-1;
            j_next = j+1;


            cost_concat_1 =                 seq(1, i_prev, T)       + cost(s(i_prev), s(j));
            cost_concat_2 = cost_concat_1                           + cost(s(j), s(i_next));
            cost_concat_3 = cost_concat_2 + seq(i_next, j_prev, T)  + cost(s(j_prev), s(i));
            cost_concat_4 = cost_concat_3                           + cost(s(i), s(j_next));


            cost_new = seq(1, i_prev, C)                                                 + ...     % 1st subseq
                    cost_concat_1 + ...                                                           % concat 2nd subseq (single node)
                    seq(i_next, j_prev, W)      * cost_concat_2 + seq(i_next, j_prev, C) + ...    % concat 3rd subseq
                    cost_concat_3 + ...                                                           % concat 4th subseq (single node)
                    seq(j_next, dimen+1, W) * cost_concat_4 + seq(j_next, dimen+1, C);   % concat 5th subseq

            if (cost_new < cost_best)
                cost_best = cost_new - eps;
                I_best = i;
                J_best = j;
            end

        end
        
    end

    if (cost_best < seq(1, dimen+1, C) - eps)
        %cost_best
        s = swap(s, I_best, J_best);
        seq = subseq_load(s);
        %solut.cost
        solut_new = s;
        seq_new = seq;
        ret = true;
    else
        solut_new = s;
        seq_new = seq;
        ret = false;
    end

end

function [solut_new, seq_new, ret] = search_two_opt(s, seq)
    global fmax;
    global dimen;
    global cost;

    global T;
    global C;
    global W;

    cost_best = fmax;
    cost_new = 0.0;
    cost_concat_1 = 0.0;
    cost_concat_2 = 0.0;

    for i = 2:dimen-1
        i_prev = i-1;
        rev_seq_cost = seq(i, i+1, T);

        for j = i+2:dimen
            j_next = j+1;

            rev_seq_cost = rev_seq_cost + cost(s(j-1), s(j)) * (seq(i, j, W)-1.0);

            cost_concat_1 =                 seq(1, i_prev, T)   + cost(s(j), s(i_prev));
            cost_concat_2 = cost_concat_1 + seq(i, j, T)        + cost(s(j_next), s(i));

            cost_new = seq(1, i_prev, C)                                                        + ... %   1st subseq
                    seq(i, j, W)                * cost_concat_1 + rev_seq_cost                  + ... % concat 2nd subseq (reversed seq)
                    seq(j_next, dimen+1, W) * cost_concat_2 + seq(j_next, dimen+1, C);       % concat 3rd subseq

            if (cost_new < cost_best)
                cost_best = cost_new - eps;
                I_best = i;
                J_best = j;
            end
        end
    end

    if (cost_best < seq(1, dimen+1, C) - eps)
        %cost_best
        s = reverse(s, I_best, J_best);
        seq = subseq_load(s);
        %solut.cost
        solut_new = s;
        seq_new = seq;
        ret = true;
    else
        solut_new = s;
        seq_new = seq;
        ret = false;
    end

end

function [solut_new, seq_new, ret] = search_reinsertion(s, seq, opt)
    global fmax;
    global dimen;
    global cost;

    global T;
    global C;
    global W;

    cost_best = fmax;
    cost_new = 0.0;
    cost_concat_1 = 0.0;
    cost_concat_2 = 0.0;
    cost_concat_3 = 0.0;

    for i = 2:dimen-opt+1
        j = opt+i-1;
        i_prev = i-1;
        j_next = j+1;

        for k = 1:i_prev-1
            k_next = k+1;

            cost_new =  seq(1, dimen+1, C);
            cost_concat_1 =                 seq(1, k, T)            + cost(s(k), s(i));
            cost_concat_2 = cost_concat_1 + seq(i, j, T)            + cost(s(j), s(k_next));
            cost_concat_3 = cost_concat_2 + seq(k_next, i_prev, T)  + cost(s(i_prev), s(j_next));

            cost_new = seq(1, k, C)                                                             + ... %       1st subseq
                    seq(i, j, W)                * cost_concat_1 + seq(i, j, C)                  + ... % concat 2nd subseq (reinserted seq)
                    seq(k_next, i_prev, W)      * cost_concat_2 + seq(k_next, i_prev, C)        + ... % concat 3rd subseq
                    seq(j_next, dimen+1, W) * cost_concat_3 + seq(j_next, dimen+1, C);      % concat 4th subseq

            if (cost_new < cost_best)
                cost_best = cost_new - eps;
                I_best = i;
                J_best = j;
                POS_best = k;
            end

        end

        for k = i+opt:dimen
            k_next = k+1;

            cost_concat_1 =                 seq(1, i_prev, T)   + cost(s(i_prev), s(j_next));
            cost_concat_2 = cost_concat_1 + seq(j_next, k, T)   + cost(s(k), s(i));
            cost_concat_3 = cost_concat_2 + seq(i, j, T)        + cost(s(j), s(k_next));

            cost_new = seq(1, i_prev, C)                                                        + ... %       1st subseq
                    seq(j_next, k, W)           * cost_concat_1 + seq(j_next, k, C)             + ... % concat 2nd subseq
                    seq(i, j, W)                * cost_concat_2 + seq(i, j, C)                  + ... % concat 3rd subseq (reinserted seq)
                    seq(k_next, dimen+1, W) * cost_concat_3 + seq(k_next, dimen+1, C);      % concat 4th subseq

            if (cost_new < cost_best)
                cost_best = cost_new - eps;
                I_best = i;
                J_best = j;
                POS_best = k;
            end
        end
    end

    if (cost_best < seq(1, dimen+1, C))
        %cost_best
        s = reinsert(s, I_best, J_best, POS_best+1);
        seq = subseq_load(s);
        %solut.cost

        solut_new = s;
        seq_new = seq;
        ret = true;
    else
        solut_new = s;
        seq_new = seq;
        ret = false;
    end

end

function [s, cost, index_new] = RVND(solut, seq, rnd)
    global SWAP;
    global TWO_OPT;
    global REINSERTION;
    global OR_OPT_2;
    global OR_OPT_3;
    global C;
    global dimen;

    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
    while (isempty(neighbd_list) == false)
        RND = rand(1);
        %RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        %index = ceil(RND*columns(neighbd_list));
        %info.rnd(info.rnd_index);
        index = rnd.rnd(rnd.rnd_index) + 1;
        rnd.rnd_index = rnd.rnd_index + 1;

        neighbd = neighbd_list(index);

        improve_flag = false;
        switch (neighbd)
            case REINSERTION
                [solut, seq, improve_flag] = search_reinsertion(solut, seq, REINSERTION);
            case OR_OPT_2
                [solut, seq, improve_flag] = search_reinsertion(solut, seq, OR_OPT_2);
            case OR_OPT_3
                [solut, seq, improve_flag] = search_reinsertion(solut, seq, OR_OPT_3);
            case SWAP
                [solut, seq, improve_flag] = search_swap(solut, seq);
            case TWO_OPT
                [solut, seq, improve_flag] = search_two_opt(solut, seq);
        end

        if (improve_flag)
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
            index;
        else
            neighbd_list(index) = [];
        end
    end

    s = solut;
    cost = seq(1, dimen+1, C);
    index_new = rnd.rnd_index;
end

function ret = notnull_rnd()
    RND = rand(1);
    %RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
    ret = RND;
end

function [ret, index_new] = perturb(s, rnd)
    global dimen;
    A_start = 1;
    A_end = 1;
    B_start = 1;
    B_end = 1;

    size_max = ceil((dimen+1) / 10);
    %size_max = [2, size_max]((size_max >= 2) + 1);
    size_min = 2;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end))
        max_ = (dimen+1) -2 -size_max;
        RND = notnull_rnd();
        A_start = ceil(max_ * RND) + 1;
        RND = notnull_rnd();
        A_end = A_start + ceil(((size_max-size_min) * RND) + size_min);

        RND = notnull_rnd();
        B_start = ceil(max_ * RND) + 1;
        RND = notnull_rnd();
        B_end = B_start + ceil(((size_max-size_min) * RND) + size_min);


        A_start = rnd.rnd(rnd.rnd_index) + 1;
        %info.rnd(info.rnd_index)
        rnd.rnd_index = rnd.rnd_index + 1;
        A_end = A_start + rnd.rnd(rnd.rnd_index);
        %info.rnd(info.rnd_index)
        rnd.rnd_index = rnd.rnd_index + 1;

        B_start = rnd.rnd(rnd.rnd_index) + 1;
        %info.rnd(info.rnd_index)
        rnd.rnd_index = rnd.rnd_index + 1;
        B_end = B_start + rnd.rnd(rnd.rnd_index);
        %info.rnd(info.rnd_index)
        rnd.rnd_index = rnd.rnd_index + 1;
    end

    if (A_start < B_start)
        s = reinsert(s, B_start, B_end - 1, A_end);
        s = reinsert(s, A_start, A_end - 1, B_end);
    else
        s = reinsert(s, A_start, A_end - 1, B_end);
        s = reinsert(s, B_start, B_end - 1, A_end);
    end

    ret = s;
    index_new = rnd.rnd_index;
end

function ret = solut_init(info)
    s.s = zeros(info.dimen+1);
    s.seq = zeros(info.dimen+1, info.dimen+1, 3);
    s.cost = Inf(1);

    ret = s;
end

function ret = GILS_RVND(Imax, Iils, R, rnd)
    global dimen;
    global C;
    solut_crnt = zeros(dimen+1);
    solut_partial = zeros(dimen+1);
    solut_best = zeros(dimen+1);

    cost_crnt = 0;
    cost_partial = 0;
    cost_best = Inf(1);

    seq = zeros(dimen+1, dimen+1, 3);

    for i = 1:Imax
        RND = rand(1);
        %RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        %index = ceil(columns(R) * RND);
        index = rnd.rnd(rnd.rnd_index) + 1;
        rnd.rnd_index = rnd.rnd_index + 1;
        alpha = R(index);

        fprintf("[+] Search %d\n", i)
        fprintf("\t[+] Constructing..\n");
        %"ITER ", i

        [solut_crnt, rnd.rnd_index] = construction(alpha, rnd);
        seq = subseq_load(solut_crnt);

        cost_crnt = seq(1, dimen+1, C);
        cost_partial = cost_crnt;

        solut_partial = solut_crnt;
        fprintf("\t[+] Looking for the best Neighbor..\n")
        fprintf("\t    Construction Cost: %.2f\n", cost_partial)

        iterILS = 0;
        while (iterILS < Iils)
            [solut_crnt, cost_crnt, rnd.rnd_index] = RVND(solut_crnt, seq, rnd);

            if (cost_crnt < cost_partial - eps)
                solut_partial = solut_crnt;
                cost_partial = cost_crnt;
                iterILS = 0;
            end

            [solut_crnt, rnd.rnd_index] = perturb(solut_partial, rnd);
            seq = subseq_load(solut_crnt);
            iterILS = iterILS + 1;
        end

        if (cost_partial < cost_best)
            solut_best = solut_partial;
            cost_best = cost_partial;
        end

        fprintf("\tCurrent best cost: %.2f\n", cost_best)
        fprintf('SOLUCAO: ')
        print_s(solut_best)
    end

    fprintf('COST: %.2f\n', cost_best)

    ret = solut_best;
end

function mainn
    global dimen
    global cost
    global fmax;
    global T
    global C
    global W
    global REINSERTION;
    global OR_OPT_2;
    global OR_OPT_3;
    global SWAP;
    global TWO_OPT;

    REINSERTION = 1;
    OR_OPT_2 = 2;
    OR_OPT_3 = 3;
    SWAP = 4;
    TWO_OPT = 5;
    
    

    fmax = Inf(1);
    [dimen, cost, rnd.rnd] = Data();
    rnd.rnd_index = 1;
    T = 1;
    C = 2;
    W = 3;


   %for i = 1:info.dimen
   %    sol.s(i) = i;
   %end
   %sol.s(info.dimen+1) = 1;
   %sol.cost = 0;
   %sol.seq = zeros(info.dimen+1, info.dimen+1, 3);
   %sol = subseq_load(sol, info);

   %sol.seq(1, info.dimen+1, info.C);
   %s = construction(0.1, info);
   %%s(2:7)
   %%s(2:7) = flip(s(2:7))
   %s = reinsert(s, 6, 11, 2);

    Imax = 10;

    Iils = dimen;
    if (dimen > 100)
        Iils = 100;
    end
    R = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26];
    t0 = clock ();
    s = GILS_RVND(Imax, Iils, R, rnd);
    elapsed_time = etime (clock (), t0);
    fprintf('TIME: %.2f\n', elapsed_time)
    exit
end
