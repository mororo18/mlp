#! /usr/bin/node

function subseq_fill(dimension) {
    seq = [];
    for (var i = 0; i < dimension+1; i++) {
        seq[i] = Array(dimension+1);
        for (var j = 0; j < dimension+1; j++) {
            seq[i][j] = Array(3);
        }
    }

    return seq;
}

function subseq_load(s, seq, info) {
    //console.log(info);
    for (var i = 0; i < info.dimension+1; i++) {
        let k = 1 - i - (i != 0) ? 0 : 1;

        seq[i][i][info.T] = 0.0;
        seq[i][i][info.C] = 0.0;
        seq[i][i][info.W] = (i != 0) ? 1.0 : 0.0;

        for (var j = i+1; j < info.dimension+1; j++) {
            let j_prev = j-1;
            seq[i][j][info.T] = info.c[s[j_prev]][s[j]] + seq[i][j_prev][info.T];
            seq[i][j][info.C] = seq[i][j][info.T] + seq[i][j_prev][info.C];
            seq[i][j][info.W] = j + k;
        }
    }
    //console.log(seq);
}

function construction(alpha, info) {
    s = [0];
    //var cList = [...Array(info.dimension).keys()];
    var cList = Array.from({length: info.dimension-1}, (_, i) => i + 1)

    var r = 0;
    while (cList.length > 0) {
        cList.sort((i, j) => info.c[i][r] - info.c[j][r]);

        var a = Math.random()*(cList.length)*alpha;
        var i = parseInt(a);
        var c = cList.splice(i, 1);
        c = c[0];
        console.log(c);
        s.push(c);
        r = c;
    }

    s[info.dimension] = 0;

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
        var sub = s.splice(i, sz);
        s.splice(pos - sz, 0, ...sub);
    } else {
        var sub = s.splice(i, sz);
        s.splice(pos, 0, ...sub);
    }
}

function search_swap(s, subseq, info) {
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;
    var cost_concat_3;
    var cost_concat_4;
    var I = -1;
    var J = -1;

    for (var i = 1; i < info.dimension-1; i++) {
        let i_prev = i - 1;
        let i_next = i + 1;

        // immediate nodes case

        cost_concat_1 =                 seq[0][i_prev][info.T] + info.c[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][i_next][info.T] + info.c[s[i]][s[i_next+1]];

        cost_new = seq[0][i_prev][info.C]                                                            // 1st subseq 
            + seq[i][i_next][info.W]           * (cost_concat_1) + info.c[s[i_next]][s[i]]              // concat 2nd subseq
            + seq[i_next+1][info.dimension][info.W] * (cost_concat_2) + seq[i_next+1][info.dimension][info.C];  // concat 3rd subseq

        if(cost_new < cost_best){
            cost_best = cost_new - Number.EPSILON;
            I = i;
            J = i_next;
        }

        for(var j = i_next+1; j < info.dimension; j++){
            let j_next = j+1;
            let j_prev = j-1;

            cost_concat_1 = seq[0][i_prev][info.T] + info.c[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + info.c[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev][info.T] + info.c[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3 + info.c[s[i]][s[j_next]];

            cost_new = seq[0][i_prev][info.C]                                                        /* first subseq */
                + cost_concat_1                                                             /* concatenate second subseq (single node) */
                + seq[i_next][j_prev][info.W] * cost_concat_2 + seq[i_next][j_prev][info.C]           /* concatenate third subseq */
                + cost_concat_3                                                             /* concatenate fourth subseq (single node) */
                + seq[j_next][info.dimension][info.W] * cost_concat_4 + seq[j_next][info.dimension][info.C];    /* concatenate fifth subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
            }

        }

    }

    if(cost_best < seq[0][info.dimension][info.C] - Number.EPSILON){
        swap(s, I, J);
        console.log(seq[0][info.dimension][info.C]);
        subseq_load(s, seq, info);
        console.log(seq[0][info.dimension][info.C]);

        return true;
    }

    return false;
}

function search_two_opt(s, subseq, info) {
    var I = -1;
    var J = -1;
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;

    for(var i = 1; i < info.dimension-1; i++){
        let i_prev = i -1;
        let rev_seq_cost = seq[i][i+1][info.T];

        for(var j = i+2; j < info.dimension; j++){
            let j_next = j+1;

            rev_seq_cost += info.c[s[j_prev]][s[j]] * (seq[i][j][info.W]-1);

            cost_concat_1 =                 seq[0][i_prev][info.T] + info.c[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[i][j][info.T] + info.c[s[j_next]][s[i]];

            cost_new = seq[0][i_prev][info.C]                                                        /*        1st subseq */
                + seq[i][j][info.W]              * cost_concat_1 + rev_seq_cost                  /* concat 2nd subseq (reversed seq) */
                + seq[j_next][info.dimension][info.W] * cost_concat_2 + seq[j_next][info.dimension][info.C];    /* concat 3rd subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if(cost_best < seq[0][info.dimension][info.C] - Number.EPSILON){
        reverse(s, I, J);
        console.log(seq[0][info.dimension][info.C]);
        subseq_load(s, seq, info);
        console.log(seq[0][info.dimension][info.C]);

        return true;
    }
    return false;
}

function search_reinsertion(s, subseq, info, opt) {
    var cost_best = Number.MAX_VALUE;
    var cost_new;
    var cost_concat_1;
    var cost_concat_2;
    var cost_concat_3;
    var I = -1;
    var J = -1;
    var POS = -1;

    for (var i = 1; i < info.dimension - opt + 1; i++) { 
        let j = opt + i - 1;
        let i_prev = i-1;
        let j_next = j+1;

        // k -> reinsertion places
        for (var k = 0; k < i_prev; k++) {
            let k_next = k+1;

            cost_concat_1 = seq[0][k][info.T] + info.c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[i][j][info.T] + info.c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][info.T] + info.c[s[i_prev]][s[j_next]];

            cost_new = seq[0][k][info.C]                                                             /*        1st subseq */
                + seq[i][j][info.W]              * cost_concat_1 + seq[i][j][info.C]                  /* concat 2nd subseq (reinserted seq) */
                + seq[k_next][i_prev][info.W]    * cost_concat_2 + seq[k_next][i_prev][info.C]        /* concat 3rd subseq */
                + seq[j_next][info.dimension][info.W] * cost_concat_3 + seq[j_next][info.dimension][info.C];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for (var k = i+opt; k < info.dimension - opt - 1; k++) {
            let k_next = k+1;

            cost_concat_1 = seq[0][i_prev][info.T] + info.c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][k][info.T] + info.c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[i][j][info.T] + info.c[s[j]][s[k_next]];

            cost_new = seq[0][i_prev][info.C]                                                        /*      1st subseq */
                + seq[j_next][k][info.W]         * cost_concat_1 + seq[j_next][k][info.C]             /* concat 2nd subseq */
                + seq[i][j][info.W]              * cost_concat_2 + seq[i][j][info.C]                  /* concat 3rd subseq (reinserted seq) */
                + seq[k_next][info.dimension][info.W] * cost_concat_3 + seq[k_next][info.dimension][info.C];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - Number.EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if(cost_best < seq[0][info.dimension][info.C] - Number.EPSILON){
        reinsert(s, I, J, POS+1);
        console.log(seq[0][info.dimension][info.C]);
        subseq_load(s, seq, info);
        console.log(seq[0][info.dimension][info.C]);
        return true;
    }

    return false;
}

function RVND(s, subseq, info) {
    const SWAP = 0;
    const REINSERTION = 1;
    const OR_OPT_2 = 2;
    const OR_OPT_3 = 3;
    const TWO_OPT = 4;

    neighbd_list = [SWAP, REINSERTION, OR_OPT_2, OR_OPT_3, TWO_OPT];

    var improve = false;

    console.log("opa");
    while (neighbd_list.length > 0) {
        let i = parseInt(Math.random() * neighbd_list.length);
        let neighbd = neighbd_list[i];
        console.log(i, neighbd);

        switch (neighbd) {
            case SWAP:
                improv = search_swap(s, subseq, info);
                break;
            case REINSERTION:
                improv = search_reinsertion(s, subseq, info, REINSERTION);
                break;
            case OR_OPT_2:
                improv = search_reinsertion(s, subseq, info, OR_OPT_2);
                break;
            case OR_OPT_3:
                improv = search_reinsertion(s, subseq, info, OR_OPT_3);
                break;
            case TWO_OPT:
                improv = search_swap(s, subseq, info);
                break;
        }

        if (improv) {
            neighbd_list = [SWAP, REINSERTION, OR_OPT_2, OR_OPT_3, TWO_OPT];
        } else {
            neighbd_list.splice(i, 1);
        }
    }
    console.log(s, subseq[0][info.dimension][info.C]);
    process.exit();

}

function perturb(sl) {
}

function GILS_RVND(Iils, Imax, R, info) {

    subseq = subseq_fill(info.dimension);

    for (var i = 0; i < info.dimension; i++) {
        var alpha = R[parseInt(Math.random() * 26)];
        var s = construction(alpha, info);
        console.log(s);
        var sl = [...s];

        subseq_load(s, subseq, info);
        var rvnd_cost_best = subseq[0][info.dimension][info.C] - Number.EPSILON;
        var iterILS = 0;

        while (iterILS < Iils) {
            RVND(s, subseq, info);
            var rvnd_cost_crnt = subseq[0][info.dimension][info.C] - Number.EPSILON;
            if (rvnd_cost_crnt < rvnd_cost_best) {
                rvnd_cost_best = rvnd_cost_crnt -  Number.EPSILON;
                sl = [...s];
                iterILS = 0;
            }

            s = perturb(sl);
            subseq_load(s, subseq, info);
            iterILS++;
        }

        subseq_load(s, subseq, info);
        var sl_cost = subseq[0][info.dimension][info.C] - Number.EPSILON;
    }
}

function main() {
    var dimension;
    var c = [];
    var Data = require("./Data"); 

    dimension = Data.info_load(c);
    Iils = Math.min(dimension, 100);
    const Imax = 10;
    const R = [0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.20,0.20,0.20,0.20,0.20];

    /*
    s = [];
    for (var i = 0; i < dimension+1; i++) {
        s[i] = i;
    }
    s[dimension] = 0;
    */

    var info = Object.freeze({
        c : c,
        dimension : dimension, 
        T : 0,
        C : 1, 
        W : 2
    });

    GILS_RVND(Iils, Imax, R, info);

}

main();
