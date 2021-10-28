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
    sol.cost
    sol;
    ret = sol;
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
end

main();

