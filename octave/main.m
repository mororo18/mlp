1;

function ret = subseq_load (sol, info)
    tern1 = [1, 0];
    tern2 = [0, 1];
    for i = 1:info.dimen+1
        k = 1 - i - tern1((i != 1) + 1);

        sol.seq(i, i, info.T) = 0.0;
        sol.seq(i, i, info.C) = 0.0;
        sol.seq(i, i, info.W) = tern2((i != 1) + 1);
        for j = i+1:info.dimen+1
            j_prev = j-1;

            sol.seq(i, j, info.T) = info.cost(sol.s(j_prev), sol.s(j)) + sol.seq(i, j_prev, info.T);
            sol.seq(i, j, info.C) = sol.seq(i, j, info.T) + sol.seq(i, j_prev, info.C);
            sol.seq(i, j, info.W) = j+k;

        end
    end

    sol.cost = sol.seq(1, info.dimen+1, info.C);
    ret = sol;
end

function ret = sort_by(arr, r, info)
    sz = columns(arr);

    for i = 1:sz
        for j = 1:sz-i
            if (info.cost(r, arr(j)) > info.cost(r, arr(j+1)))
                tmp = arr(j);
                arr(j) = arr(j+1);
                arr(j+1) = tmp;
            end
        end
    end

    ret = arr;
end

function ret = construction(alpha, info)
    s(1) = 1;
    for i = 1:info.dimen-1
        cL(i) = i+1;
    end

    r = 1;
    while (columns(cL) > 0) 
        cL = sort_by(cL, r, info);

        rng = ceil(columns(cL) * alpha);
        RND = rand(1);
        RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        index = ceil(rng * RND);

        cN = cL(index);

        cL(index) = [];
        s(columns(s)+1) = cN;
        r = cN;
    end

    s(columns(s)+1) = 1;
    ret = s;
end

function [solut_new, flag] = search_swap(solut, info)
end

function [solut_new, flag] = search_two_opt(solut, info)
end

function [solut_new, flag] = search_reinsertion(solut, info, opt)
end

function ret = RVND(solut, info)

    neighbd_list = [info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3, info.TWO_OPT, info.SWAP]

    while (columns(neighbd_list) > 0)
        RND = rand(1);
        RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        index = ceil(RND*columns(neighbd_list))
        neighbd = neighbd_list(index)

        improve_flag = false
        switch (neighbd)
            case info.REINSERTION
            case info.OR_OPT_2
            case info.OR_OPT_3
            case info.SWAP
            case info.TWO_OPT
        end

        if (improve_flag)
            neighbd_list = [info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3, info.TWO_OPT, info.SWAP]
        else
            neighbd_list(index) = []
        end
    end

    ret = solut
end

function ret = solut_init(info)
    s.s = zeros(info.dimen+1);
    s.seq = zeros(info.dimen+1, info.dimen+1, 3);
    s.cost = Inf(1);

    ret = s;
end

function ret = GILS_RVND(Imax, Iils, R, info)
    solut_crnt = solut_init(info);
    solut_partial = solut_init(info);
    solut_best = solut_init(info);

    for i = 1:Imax
        RND = rand(1);
        RND = [RND, RND+0.0000000001]((RND < 0.0000000001) + 1);
        index = ceil(columns(R) * RND);
        alpha = R(index);

        solut_crnt.s = construction(alpha, info);
        solut_crnt = subseq_load(solut_crnt, info);

        solut_partial = solut_crnt;
        iterILS = 0;
        while (iterILS < Iils)
            solut_crnt = RVND(solut_crnt, info)

            if (solut_crnt.cost < solut_partial.cost - eps)
                solut_partial = solut_crnt
                iterILS = 0
            end
            iterILS++;
        end

        if (solut_partial.cost < solut_best.cost)
            solut_best = solut_partial
        end
    end
end

function main
    [info.dimen, info.cost] = Data();
    info;
    info.T = 1;
    info.C = 2;
    info.W = 3;
    info.REINSERTION = 1
    info.OR_OPT_2 = 2
    info.OR_OPT_3 = 3
    info.SWAP = 4
    info.TWO_OPT = 5

    for i = 1:info.dimen
        sol.s(i) = i;
    end
    sol.s(info.dimen+1) = 1;
    sol.cost = 0;
    sol.seq = zeros(info.dimen+1, info.dimen+1, 3);
    sol = subseq_load(sol, info);

    s = construction(0.1, info);
end

main();
