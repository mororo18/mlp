#! /usr/bin/node

var to_1D = function (x, y, z, n) {
    return (x*n + y)*3 + z;
}

function subseq_fill(dimension) {
    seq = Array((dimension+1)*(dimension+1)*3);
    //seq = [(dimension+1)*(dimension+1)*3];
    for (var i = 0; i < dimension+1; i++) {
        //seq[i] = Array(dimension+1);
        for (var j = 0; j < dimension+1; j++) {
            //seq[i][j] = Array(3);
        }
    }

    return seq;
}

function update_subseq_info_matrix(s, seq, info, data) {

    for (var i = 0; i < data.dimension+1; i++) {
        var k = 1 - i - (i != 0 ? 0 : 1);

        seq[to_1D(i, i, info.T, data.dimension)] = 0.0;
        seq[to_1D(i, i, info.C, data.dimension)] = 0.0;
        seq[to_1D(i, i, info.W, data.dimension)] = (i != 0 ? 1.0 : 0.0);

        for (var j = i+1; j < data.dimension+1; j++) {
            let j_prev = j-1;
            seq[to_1D(i, j, info.T, data.dimension)] = data.c[s[j_prev]][s[j]] + seq[to_1D(i, j_prev, info.T, data.dimension)];
            seq[to_1D(i, j, info.C, data.dimension)] = seq[to_1D(i, j, info.T, data.dimension)] + seq[to_1D(i, j_prev, info.C, data.dimension)];
            seq[to_1D(i, j, info.W, data.dimension)] = j + k;
        }
    }
}

function sort(arr, r, data) {
    quicksort(arr, 0, arr.length - 1, r, data);
}

function quicksort(arr, left, right, r, data) {
    if (left < right) {
        let pivotIndex = partition(arr, left, right, r, data);
        quicksort(arr, left, pivotIndex - 1, r, data);
        quicksort(arr, pivotIndex + 1, right, r, data);
    }
}

function partition(arr, left, right, r, data) {
    let pivotValue = arr[right];
    let i = left - 1;
    for (let j = left; j < right; j++) {
        if (data.c[r][arr[j]] < data.c[r][pivotValue]) {
            i++;
            [arr[i], arr[j]] = [arr[j], arr[i]]; // Swap elements
        }
    }
    [arr[i + 1], arr[right]] = [arr[right], arr[i + 1]]; // Swap pivot
    return i + 1;
}

function construction(alpha, data) {
    s = [0];
    var cList = Array.from({length: data.dimension-1}, (_, i) => i + 1)

    var r = 0;
    while (cList.length > 0) {
        sort(cList, r, data);

        var i = data.rnd[data.rnd_index++];
        var c = cList.splice(i, 1);
        c = c[0];
        s.push(c);
        r = c;
    }

    s[data.dimension] = 0;

    return s;
}

function swap(s, i, j) {
    [s[i], s[j]] = [s[j], s[i]];
}

function reverse(s, i, j) {
    for (var f = i, b = j; f <= (i+j)/2; f++, b--) {
        swap(s, f, b);
    }
}

function reinsert(s, i, j, pos) {
    var sz = j - i + 1;
    if (i < pos) {
        let sub = s.splice(i, sz);
        s.splice(pos - sz, 0, ...sub);
    } else {
        let sub = s.splice(i, sz);
        s.splice(pos, 0, ...sub);
    }
}

function search_swap(s, seq, info, data) {
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;
    var cost_concat_3;
    var cost_concat_4;
    var I = -1;
    var J = -1;

    for (var i = 1; i < data.dimension-1; i++) {
        var i_prev = i - 1;
        var i_next = i + 1;

        // immediate nodes case

        cost_concat_1 =                 seq[to_1D(0, i_prev, info.T, data.dimension)] + data.c[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[to_1D(i, i_next, info.T, data.dimension)] + data.c[s[i]][s[i_next+1]];

        cost_new = seq[to_1D(0, i_prev, info.C, data.dimension)]                                                            // 1st subseq 
            + seq[to_1D(i, i_next, info.W, data.dimension)]           * (cost_concat_1) + data.c[s[i_next]][s[i]]              // concat 2nd subseq
            + seq[to_1D(i_next+1, data.dimension, info.W, data.dimension)] * (cost_concat_2) + seq[to_1D(i_next+1, data.dimension, info.C, data.dimension)];  // concat 3rd subseq

        if(cost_new < cost_best){
            cost_best = cost_new - Number.EPSILON;
            I = i;
            J = i_next;
        }

        for(var j = i_next+1; j < data.dimension; j++){
            var j_next = j+1;
            var j_prev = j-1;

            cost_concat_1 = seq[to_1D(0, i_prev, info.T, data.dimension)] + data.c[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + data.c[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(i_next, j_prev, info.T, data.dimension)] + data.c[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3 + data.c[s[i]][s[j_next]];

            cost_new = seq[to_1D(0, i_prev, info.C, data.dimension)]                                                        /* first subseq */
                + cost_concat_1                                                             /* concatenate second subseq (single node) */
                + seq[to_1D(i_next, j_prev, info.W, data.dimension)] * cost_concat_2 + seq[to_1D(i_next, j_prev, info.C, data.dimension)]           /* concatenate third subseq */
                + cost_concat_3                                                             /* concatenate fourth subseq (single node) */
                + seq[to_1D(j_next, data.dimension, info.W, data.dimension)] * cost_concat_4 + seq[to_1D(j_next, data.dimension, info.C, data.dimension)];    /* concatenate fifth subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
            }

        }

    }


    if(cost_best < seq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON){
        swap(s, I, J);
        /*
        console.log("swap");
        console.log(cost_best);
        */
        update_subseq_info_matrix(s, seq, info, data);
        /*
        console.log(seq[0][data.dimension][info.C]);
        console.log();
        */

        return true;
    }

    return false;
}

function search_two_opt(s, seq, info, data) {
    var I = -1;
    var J = -1;
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;

    for(var i = 1; i < data.dimension-1; i++){
        var i_prev = i -1;
        var rev_seq_cost = seq[to_1D(i, i+1, info.T, data.dimension)];

        for(var j = i+2; j < data.dimension; j++){
            var j_next = j+1;
            var j_prev = j-1;

            rev_seq_cost += data.c[s[j_prev]][s[j]] * (seq[to_1D(i, j, info.W, data.dimension)]-1);

            cost_concat_1 =                 seq[to_1D(0, i_prev, info.T, data.dimension)] + data.c[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(i, j, info.T, data.dimension)] + data.c[s[j_next]][s[i]];

            cost_new = seq[to_1D(0, i_prev, info.C, data.dimension)]                                                        /*        1st subseq */
                + seq[to_1D(i, j, info.W, data.dimension)]              * cost_concat_1 + rev_seq_cost                  /* concat 2nd subseq (reversed seq) */
                + seq[to_1D(j_next, data.dimension, info.W, data.dimension)] * cost_concat_2 + seq[to_1D(j_next, data.dimension, info.C, data.dimension)];    /* concat 3rd subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if(cost_best < seq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON){
        reverse(s, I, J);
        /*
        console.log("two opt");
        console.log(cost_best);
        */
        update_subseq_info_matrix(s, seq, info, data);
        /*
        console.log(seq[0][data.dimension][info.C]);
        console.log();
        */

        return true;
    }
    return false;
}

function search_reinsertion(s, seq, info, data, opt) {
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;
    var cost_concat_3;
    var I = -1;
    var J = -1;
    var POS = -1;

    for (var i = 1; i < data.dimension - opt + 1; i++) { 
        var j = opt + i - 1;
        var i_prev = i-1;
        var j_next = j+1;

        // k -> reinsertion places
        for (var k = 0; k < i_prev; k++) {
            var k_next = k+1;

            cost_concat_1 = seq[to_1D(0, k, info.T, data.dimension)] + data.c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(i, j, info.T, data.dimension)] + data.c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(k_next, i_prev, info.T, data.dimension)] + data.c[s[i_prev]][s[j_next]];

            cost_new = seq[to_1D(0, k, info.C, data.dimension)]                                                             /*        1st subseq */
                + seq[to_1D(i, j, info.W, data.dimension)]              * cost_concat_1 + seq[to_1D(i, j, info.C, data.dimension)]                  /* concat 2nd subseq (reinserted seq) */
                + seq[to_1D(k_next, i_prev, info.W, data.dimension)]    * cost_concat_2 + seq[to_1D(k_next, i_prev, info.C, data.dimension)]        /* concat 3rd subseq */
                + seq[to_1D(j_next, data.dimension, info.W, data.dimension)] * cost_concat_3 + seq[to_1D(j_next, data.dimension, info.C, data.dimension)];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (var k = i+opt; k < data.dimension; k++) {
            var k_next = k+1;

            cost_concat_1 = seq[to_1D(0, i_prev, info.T, data.dimension)] + data.c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[to_1D(j_next, k, info.T, data.dimension)] + data.c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[to_1D(i, j, info.T, data.dimension)] + data.c[s[j]][s[k_next]];

            cost_new = seq[to_1D(0, i_prev, info.C, data.dimension)]                                                        /*      1st subseq */
                + seq[to_1D(j_next, k, info.W, data.dimension)]         * cost_concat_1 + seq[to_1D(j_next, k, info.C, data.dimension)]             /* concat 2nd subseq */
                + seq[to_1D(i, j, info.W, data.dimension)]              * cost_concat_2 + seq[to_1D(i, j, info.C, data.dimension)]                  /* concat 3rd subseq (reinserted seq) */
                + seq[to_1D(k_next, data.dimension, info.W, data.dimension)] * cost_concat_3 + seq[to_1D(k_next, data.dimension, info.C, data.dimension)];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if(cost_best < seq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON){
        reinsert(s, I, J, POS+1);
        update_subseq_info_matrix(s, seq, info, data);

        return true;
    }

    return false;
}

function RVND(s, subseq, info, data) {
    const SWAP        = 0;
    const REINSERTION = 1;
    const OR_OPT_2    = 2;
    const OR_OPT_3    = 3;
    const TWO_OPT     = 4;

    neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];

    var improve = false;

    while (neighbd_list.length > 0) {
        i = data.rnd[data.rnd_index++];
        let neighbd = neighbd_list[i];

        switch (neighbd) {
            case SWAP:
                improve = search_swap(s, subseq, info, data);
                break;
            case REINSERTION:
                improve = search_reinsertion(s, subseq, info, data, REINSERTION);
                break;
            case OR_OPT_2:
                improve = search_reinsertion(s, subseq, info, data, OR_OPT_2);
                break;
            case OR_OPT_3:
                improve = search_reinsertion(s, subseq, info, data, OR_OPT_3);
                break;
            case TWO_OPT:
                improve = search_two_opt(s, subseq, info, data);
                break;
        }

        if (improve) {
            neighbd_list = [SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
        } else {
            neighbd_list.splice(i, 1);
        }
    }

}

function perturb(sl, data) {
    var s = [...sl];

    var A_start = 1, A_end = 1;
    var B_start = 1, B_end = 1;

    var size_max = parseInt(Math.floor(sl.length/10));
    size_max = (size_max >= 2 ? size_max : 2);
    var size_min = 2;

    while ((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)) {

        A_start = data.rnd[data.rnd_index++];
        A_end = A_start + data.rnd[data.rnd_index++];

        B_start = data.rnd[data.rnd_index++];
        B_end = B_start + data.rnd[data.rnd_index++];
    }

    if(A_start < B_start){
        reinsert(s, B_start, B_end-1, A_end);
        reinsert(s, A_start, A_end-1, B_end);
    }else{
        reinsert(s, A_start, A_end-1, B_end);
        reinsert(s, B_start, B_end-1, A_end);
    }

    return s;

}

function GILS_RVND(Iils, Imax, R, info, data) {

    subseq = subseq_fill(data.dimension);
    var s_best = [];
    var cost_best = Number.MAX_VALUE;

    for (var i = 0; i < Imax; i++) {
        var alpha = R[data.rnd[data.rnd_index++]];
        console.log("[+] Local Search ", i+1);
        console.log("\t[+] Constructing Inital Solution..");
        var s = construction(alpha, data);
        var sl = [...s];

        update_subseq_info_matrix(s, subseq, info, data);
        var rvnd_cost_best = subseq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON;
        console.log("Construction cost", rvnd_cost_best);
        var iterILS = 0;

        console.log("\t[+] Looking for the best Neighbor..");
        while (iterILS < Iils) {
            RVND(s, subseq, info, data);
            var rvnd_cost_crnt = subseq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON;
            if (rvnd_cost_crnt < rvnd_cost_best) {
                rvnd_cost_best = rvnd_cost_crnt;
                sl = [...s];
                iterILS = 0;
            }

            s = perturb(sl, data);
            update_subseq_info_matrix(s, subseq, info, data);
            iterILS++;
        }


        update_subseq_info_matrix(sl, subseq, info, data);
        var sl_cost = subseq[to_1D(0, data.dimension, info.C, data.dimension)] - Number.EPSILON;

        if (sl_cost < cost_best) {
            s_best = sl;
            cost_best = sl_cost;
        }
        console.log("\tcurrent best solution cost",  cost_best);
        console.log();
    }

    console.log(s_best);
    console.log("COST: ", cost_best);
}

function main() {
    var dimension;
    var c = [];
    var rnd = [];
    var Data = require("./Data"); 

    var t = Data.data_load(c);
    dimension = t.dimension;
    rnd = t.rnd;
    Iils = Math.min(dimension, 100);
    const Imax = 10;
    const R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21,0.22, 0.23, 0.24, 0.25];

    var data = {
        c : c,
        rnd : rnd,
        dimension : dimension, 
        rnd_index : 0
    }

    var info = {
        T : 0,
        C : 1, 
        W : 2, 
    };

    var start = new Date();
    GILS_RVND(Iils, Imax, R, info, data);
    var end = new Date();

    console.log("TIME: ", (end-start)/1000);
}

main();
