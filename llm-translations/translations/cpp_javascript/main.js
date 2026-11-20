'use strict';

const { loadData } = require('./data.js');

const TABLE_SZ = 26;

const DBL_MAX = Number.MAX_VALUE;
const DBL_EPSILON = Number.EPSILON;

const Neighborhood = {
    REINSERTION: 1,
    OR_OPT_2: 2,
    OR_OPT_3: 3,
    SWAP: 4,
    TWO_OPT: 5,
};

function to_1D(x, y, size) {
    const a = x * (size + 1) + y;

    return a;
}

function Solution_init(data) {
    const solut = {};
    solut.s = new Array(data.dimen + 1).fill(0);

    const sz = (data.dimen + 1) * (data.dimen + 1);
    solut.seq_T = new Float64Array(sz);
    solut.seq_C = new Float64Array(sz);
    solut.seq_W = new Float64Array(sz);

    solut.cost = DBL_MAX;

    return solut;
}

function Solution_cpy(src, tgt, data) {

    tgt.s = src.s.slice();
    tgt.cost = src.cost;

}

function R_table(i) {
    const table = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                   0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];
    return table[i];
}

function print_s(s) {

    let out = '';
    for (let i = 0; i < s.length; i++)
        out += (s[i] + 1) + ' ';
    console.log(out);
}

function partition(arr, left, right, data, r) {
    const pivot = arr[right];
    let i = left - 1;
    for (let j = left; j < right; j++) {
        if (data.cost[r][arr[j]] < data.cost[r][pivot]) {
            i++;
            const tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
        }
    }
    const tmp = arr[i + 1]; arr[i + 1] = arr[right]; arr[right] = tmp;
    return i + 1;
}

function quicksort(arr, left, right, data, r) {
    if (left < right) {
        const pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

function sort(arr, r, data) {
    quicksort(arr, 0, arr.length - 1, data, r);
}

function construct(alpha, data) {

    const s = [0];
    const cL = new Array(data.dimen - 1);
    for (let i = 1; i < data.dimen; ++i) {
        cL[i - 1] = i;
    }

    let r = 0;
    while (cL.length !== 0) {
        sort(cL, r, data);

        const index = data.rnd[data.rnd_index++];

        const c = cL[index];
        s.push(c);
        r = c;
        cL.splice(index, 1);
    }

    s.push(0);

    return s;
}

function swap(vec, i, j) {
    const tmp = vec[i]; vec[i] = vec[j]; vec[j] = tmp;
}

function reverse(vec, i, j) {
    while (i < j) {
        const tmp = vec[i]; vec[i] = vec[j]; vec[j] = tmp;
        i++;
        j--;
    }
}

function reinsert(vec, i, j, pos) {
    const seq = vec.slice(i, j + 1);
    if (pos < i) {
        vec.splice(i, j + 1 - i);
        vec.splice(pos, 0, ...seq);
    } else {
        vec.splice(pos, 0, ...seq);
        vec.splice(i, j + 1 - i);
    }

}

function update_subseq_info_matrix(solut, data, index = 0) {
    let i, j, j_prev, k;
    let from = index;
    let t;

    const s = solut.s;
    const cost = data.cost;
    const dimen = data.dimen;
    const T = solut.seq_T, C = solut.seq_C, W = solut.seq_W;

    for (i = 0; i < dimen + 1; i++) {
        k = 1 - i - (i === 0 ? 1 : 0);
        t = (i === from);

        T[to_1D(i, i, dimen)] = 0.0;
        C[to_1D(i, i, dimen)] = 0.0;
        W[to_1D(i, i, dimen)] = (i === 0 ? 0.0 : 1.0);

        for (j = i + 1; j < dimen + 1; j++) {
            j_prev = j - 1;

            T[to_1D(i, j, dimen)] = cost[s[j_prev]][s[j]] + T[to_1D(i, j_prev, dimen)];
            C[to_1D(i, j, dimen)] = T[to_1D(i, j, dimen)] + C[to_1D(i, j_prev, dimen)];
            W[to_1D(i, j, dimen)] = j + k;
        }
        from += (t ? 1 : 0);
    }

    solut.cost = C[to_1D(0, dimen, dimen)];
}

function search_swap(solut, data) {
    let cost_new,
        cost_concat_1, cost_concat_2, cost_concat_3, cost_concat_4;
    let cost_best = DBL_MAX;
    let i, j, j_prev, j_next, i_prev, i_next;
    let I;
    let J;

    const s = solut.s;
    const cost = data.cost;
    const dimen = data.dimen;
    const T = solut.seq_T, C = solut.seq_C, W = solut.seq_W;

    for (i = 1; i < dimen - 1; ++i) {
        i_prev = i - 1;
        i_next = i + 1;

        //consecutive nodes
        cost_concat_1 =                 T[to_1D(0, i_prev, dimen)] + cost[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + T[to_1D(i, i_next, dimen)] + cost[s[i]][s[i_next + 1]];

        cost_new = C[to_1D(0, i_prev, dimen)]                                                      +   //       1st subseq
        W[to_1D(i, i_next, dimen)]             * (cost_concat_1) + cost[s[i_next]][s[i]]           +   // concat 2nd subseq
        W[to_1D(i_next + 1, dimen, dimen)]     * (cost_concat_2) + C[to_1D(i_next + 1, dimen, dimen)];  // concat 3rd subseq

        if (cost_new < cost_best) {
            cost_best = cost_new - DBL_EPSILON;
            I = i;
            J = i_next;
        }

        for (j = i_next + 1; j < dimen; ++j) {
            j_next = j + 1;
            j_prev = j - 1;

            cost_concat_1 =                 T[to_1D(0, i_prev, dimen)]      + cost[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1                                   + cost[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + T[to_1D(i_next, j_prev, dimen)] + cost[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3                                   + cost[s[i]][s[j_next]];

            cost_new = C[to_1D(0, i_prev, dimen)]                                                  +   // 1st subseq
            cost_concat_1 +                                                                            // concat 2nd subseq (single node)
            W[to_1D(i_next, j_prev, dimen)] * cost_concat_2 + C[to_1D(i_next, j_prev, dimen)]      +   // concat 3rd subseq
            cost_concat_3 +                                                                            // concat 4th subseq (single node)
            W[to_1D(j_next, dimen, dimen)]  * cost_concat_4 + C[to_1D(j_next, dimen, dimen)];          // concat 5th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        swap(solut.s, I, J);
        update_subseq_info_matrix(solut, data, I);
        return true;
    }

    return false;
}

function search_two_opt(solut, data) {
    let cost_new,
        cost_concat_1, cost_concat_2;
    let cost_best = DBL_MAX; // cost_l1, cost_l2;
    let rev_seq_cost;
    let i, j, i_prev, j_next;
    let I;
    let J;

    const s = solut.s;
    const cost = data.cost;
    const dimen = data.dimen;
    const T = solut.seq_T, C = solut.seq_C, W = solut.seq_W;

    for (i = 1; i < dimen - 1; ++i) {
        i_prev = i - 1;

        rev_seq_cost = T[to_1D(i, i + 1, dimen)];
        for (j = i + 2; j < dimen; ++j) {
            j_next = j + 1;

            rev_seq_cost += cost[s[j - 1]][s[j]] * (W[to_1D(i, j, dimen)] - 1.0);

            cost_concat_1 =                 T[to_1D(0, i_prev, dimen)] + cost[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + T[to_1D(i, j, dimen)]      + cost[s[j_next]][s[i]];

            cost_new = C[to_1D(0, i_prev, dimen)]                                                  +   //  1st subseq
                    W[to_1D(i, j, dimen)]          * cost_concat_1 + rev_seq_cost                  +   // concat 2nd subseq (reversed seq)
                    W[to_1D(j_next, dimen, dimen)] * cost_concat_2 + C[to_1D(j_next, dimen, dimen)];   // concat 3rd subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reverse(solut.s, I, J);
        update_subseq_info_matrix(solut, data);
        return true;
    }

    return false;
}

function search_reinsertion(solut, data, opt) {
    let cost_new, cost_concat_1, cost_concat_2, cost_concat_3;
    let cost_best = DBL_MAX; //, cost_l1, cost_l2, cost_l3;
    let i, j, k, k_next, i_prev, j_next;
    let I;
    let J;
    let POS;

    const s = solut.s;
    const cost = data.cost;
    const dimen = data.dimen;
    const T = solut.seq_T, C = solut.seq_C, W = solut.seq_W;

    for (i = 1, j = opt + i - 1; i < dimen - opt + 1; ++i, ++j) {
        j_next = j + 1;
        i_prev = i - 1;

        //k -> edges
        for (k = 0; k < i_prev; ++k) {
            k_next = k + 1;

            cost_concat_1 =                 T[to_1D(0, k, dimen)]             + cost[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + T[to_1D(i, j, dimen)]             + cost[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + T[to_1D(k_next, i_prev, dimen)]   + cost[s[i_prev]][s[j_next]];

            cost_new = C[to_1D(0, k, dimen)]                                                          +   //       1st subseq
                W[to_1D(i, j, dimen)]           * cost_concat_1 + C[to_1D(i, j, dimen)]               +   //  concat 2nd subseq (reinserted seq)
                W[to_1D(k_next, i_prev, dimen)] * cost_concat_2 + C[to_1D(k_next, i_prev, dimen)]     +   //  concat 3rd subseq
                W[to_1D(j_next, dimen, dimen)]  * cost_concat_3 + C[to_1D(j_next, dimen, dimen)];         // concat 4th subseq

            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (k = i + opt; k < dimen; ++k) {
            k_next = k + 1;

            cost_concat_1 =                 T[to_1D(0, i_prev, dimen)] + cost[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + T[to_1D(j_next, k, dimen)] + cost[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + T[to_1D(i, j, dimen)]      + cost[s[j]][s[k_next]];

            cost_new = C[to_1D(0, i_prev, dimen)]                                                     +   //       1st subseq
                W[to_1D(j_next, k, dimen)]      * cost_concat_1 + C[to_1D(j_next, k, dimen)]          +   // concat 2nd subseq
                W[to_1D(i, j, dimen)]           * cost_concat_2 + C[to_1D(i, j, dimen)]               +   // concat 3rd subseq (reinserted seq)
                W[to_1D(k_next, dimen, dimen)]  * cost_concat_3 + C[to_1D(k_next, dimen, dimen)];         // concat 4th subseq


            if (cost_new < cost_best) {
                cost_best = cost_new - DBL_EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if (cost_best < solut.cost - DBL_EPSILON) {
        reinsert(solut.s, I, J, POS + 1);
        update_subseq_info_matrix(solut, data, I < POS + 1 ? I : POS + 1);
        return true;
    }

    return false;
}

function RVND(solut, data) {

    let neighbd_list = [Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3];
    let index;
    let neighbd;
    let improve_flag;

    while (neighbd_list.length !== 0) {

        index = data.rnd[data.rnd_index++];
        neighbd = neighbd_list[index];
        improve_flag = false;

        switch (neighbd) {
            case Neighborhood.REINSERTION:
                improve_flag = search_reinsertion(solut, data, Neighborhood.REINSERTION);
                break;
            case Neighborhood.OR_OPT_2:
                improve_flag = search_reinsertion(solut, data, Neighborhood.OR_OPT_2);
                break;
            case Neighborhood.OR_OPT_3:
                improve_flag = search_reinsertion(solut, data, Neighborhood.OR_OPT_3);
                break;
            case Neighborhood.SWAP:
                improve_flag = search_swap(solut, data);
                break;
            case Neighborhood.TWO_OPT:
                improve_flag = search_two_opt(solut, data);
                break;
        }
        if (improve_flag) {
            neighbd_list = [Neighborhood.SWAP, Neighborhood.TWO_OPT, Neighborhood.REINSERTION, Neighborhood.OR_OPT_2, Neighborhood.OR_OPT_3];
        } else {
            neighbd_list.splice(index, 1);
        }

    }

}

function perturb(solut, data) {
    const s = solut.s.slice();
    let A_start = 1;
    let A_end = 1;
    let B_start = 1;
    let B_end = 1;

    let size_max = Math.floor((data.dimen + 1) / 10);
    size_max = size_max >= 2 ? size_max : 2;
    const size_min = 2;
    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data.rnd[data.rnd_index++];
        A_end = A_start + data.rnd[data.rnd_index++];

        B_start = data.rnd[data.rnd_index++];
        B_end = B_start + data.rnd[data.rnd_index++];
    }

    if (A_start < B_start) {
        reinsert(s, B_start, B_end - 1, A_end);
        reinsert(s, A_start, A_end - 1, B_end);
    } else {
        reinsert(s, A_start, A_end - 1, B_end);
        reinsert(s, B_start, B_end - 1, A_end);
    }

    return s;
}


function GILS_RVND(Imax, Iils, data, verbose) {

    const solut_partial = Solution_init(data);
    const solut_crnt = Solution_init(data);
    const solut_best = Solution_init(data);

    for (let i = 0; i < Imax; ++i) {
        const aux = data.rnd[data.rnd_index++];

        const alpha = R_table(aux);

        if (verbose) {
            process.stdout.write(`[+] Search ${i + 1}\n`);
            process.stdout.write('\t[+] Constructing..\n');
        }


        solut_crnt.s = construct(alpha, data);
        update_subseq_info_matrix(solut_crnt, data);

        Solution_cpy(solut_crnt, solut_partial, data);
        if (verbose) {
            process.stdout.write('\t[+] Looking for the best Neighbor..\n');
            process.stdout.write(`\t    Construction Cost: ${solut_partial.cost.toFixed(3)}\n`);
        }

        let iterILS = 0;
        while (iterILS < Iils) {
            RVND(solut_crnt, data);
            if (solut_crnt.cost < solut_partial.cost - DBL_EPSILON) {
                Solution_cpy(solut_crnt, solut_partial, data);
                iterILS = 0;
            }

            solut_crnt.s = perturb(solut_partial, data);
            update_subseq_info_matrix(solut_crnt, data);
            iterILS++;
        }

        if (solut_partial.cost < solut_best.cost - DBL_EPSILON) {
            Solution_cpy(solut_partial, solut_best, data);
        }

        if (verbose) {
            console.log('\tCurrent best cost: ' + ostream_double(solut_best.cost));

            let out = 'SOLUCAO: ';
            for (let i = 0; i < solut_best.s.length; i++) {
                out += solut_best.s[i] + ' ';
            }
            console.log(out);
        }

    }
    process.stdout.write(`COST: ${solut_best.cost.toFixed(2)}\n`);
}

// mimics the default formatting of std::cout for doubles (6 significant digits)
function ostream_double(v) {
    if (!isFinite(v)) return String(v);
    if (v === 0) return '0';

    const exp = Math.floor(Math.log10(Math.abs(v)));
    if (exp < -5 || exp >= 6) {
        let str = v.toExponential(5);
        let [mant, e] = str.split('e');
        if (mant.indexOf('.') !== -1) mant = mant.replace(/0+$/, '').replace(/\.$/, '');
        const sign = e[0];
        let digits = e.slice(1);
        if (digits.length < 2) digits = '0' + digits;
        return mant + 'e' + sign + digits;
    }

    let str = v.toPrecision(6);
    if (str.indexOf('.') !== -1) str = str.replace(/0+$/, '').replace(/\.$/, '');
    return str;
}

function main(argv) {
    const Imax = 10;
    let Iils;
    let verbose = false;

    for (let i = 0; i < argv.length; i++) {
        if (argv[i] === '-v' || argv[i] === '--verbose') {
            verbose = true;
        }
    }

    const data = {
        cost: null,
        dimen: 0,
        rnd: null,
        rnd_index: 0,
    };


    const rnd = [];

    const loaded = loadData(rnd);
    data.cost = loaded.matrix;
    data.dimen = loaded.dimension;
    data.rnd = rnd;
    data.rnd_index = 0;

    const solut = Solution_init(data);

    Iils = data.dimen < 100 ? data.dimen : 100;
    const t1 = process.hrtime.bigint();
    const start = process.cpuUsage();
    GILS_RVND(Imax, Iils, data, verbose);
    const cpu_usage = process.cpuUsage(start);
    const cpu_time = (cpu_usage.user + cpu_usage.system) / 1e6;
    const t2 = process.hrtime.bigint();
    const duration = Number((t2 - t1) / 1000000n);

    const res = duration / 10e2;
    console.log('TIME: ' + ostream_double(cpu_time));
    console.log('wall clock (s): ' + ostream_double(res));

    console.log('Tamanho RND ' + rnd.length);

    process.exit(0);
}

main(process.argv.slice(2));
