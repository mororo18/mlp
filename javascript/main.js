const T = 0;
const C = 1;
const W = 2;

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

        seq[i][i][T] = 0.0;
        seq[i][i][C] = 0.0;
        seq[i][i][W] = (i != 0) ? 1.0 : 0.0;

        for (var j = i+1; j < info.dimension+1; j++) {
            let j_prev = j-1;
            seq[i][j][T] = info.c[s[j_prev]][s[j]] + seq[i][j_prev][T];
            seq[i][j][C] = seq[i][j][T] + seq[i][j_prev][C];
            seq[i][j][W] = j + k;
        }
    }
    console.log(seq);
}

function GILS_RVND(Iils, Imax, R) {
}

function main() {
    var dimension;
    var c = [];
    var Data = require("./Data"); 

    dimension = Data.info_load(c);
    Iils = Math.min(dimension, 100);
    Imax = 10;

    subseq = subseq_fill(dimension);

    s = [];
    for (var i = 0; i < dimension+1; i++) {
        s[i] = i;
    }
    s[dimension] = 0;

    info = {c : c, dimension : dimension};

    subseq_load(s, subseq, info);
    //console.log(info.c);
}

main();
