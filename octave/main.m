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

function main
    [info.dimen, info.cost] = Data();
    info;
    info.T = 1;
    info.C = 2;
    info.W = 3;

    for i = 1:info.dimen
        sol.s(i) = i;
    end
    sol.s(info.dimen+1) = 1;
    sol.cost = 0;
    sol.seq = zeros(info.dimen+1, info.dimen+1, 3);
    sol = subseq_load(sol, info);

    s = construction(0.1, info)
end

main();
