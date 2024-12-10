mod data;
extern crate rand;

use rand::Rng;
use std::cmp::Ordering;
use std::process::exit;
use std::time::Instant;

static SIZE: usize = 319;

#[derive(Debug, Clone)]
struct tRnd {
    rnd: Vec<usize>,
    rnd_index: usize,
}

#[derive(Debug, Copy, Clone)]
struct tSeqInfo {
    T: f64,
    C: f64,
    W: f64,
}

#[derive(Debug, Clone)]
struct tSolution {
    seq: Box<[[tSeqInfo; 319]; 319]>,
    s: Vec<usize>,
    cost: f64,
}

#[inline(always)]
unsafe fn to_1D(x: usize, y: usize, size: usize) -> usize {
    return x * (size + 1) + y;
}

fn update_subseq_info_matrix(
    s: &Vec<usize>,
    seq: &mut Vec<Vec<tSeqInfo>>,
    c: &Vec<Vec<f64>>,
    dimen: usize,
) {
    for i in 0..dimen + 1 {
        let k: i32 = 1 - (i as i32) - if i == 0 { 1 } else { 0 };

        seq[i][i].T = 0.0;
        seq[i][i].C = 0.0;
        seq[i][i].W = if i != 0 { 1.0 } else { 0.0 };

        for j in (i + 1)..(dimen + 1) {
            let j_prev: usize = j - 1;

            seq[i][j].T = c[s[j_prev]][s[j]] + seq[i][j_prev].T;
            seq[i][j].C = seq[i][j].T + seq[i][j_prev].C;
            seq[i][j].W = (j as i32 + k) as f64;
        }
    }
}

fn sort(arr: &mut Vec<usize>, r: usize, cost: &Vec<Vec<f64>>) {
    quicksort(arr, 0, arr.len() as isize - 1, cost, r);
}

fn quicksort(arr: &mut Vec<usize>, left: isize, right: isize, cost: &Vec<Vec<f64>>, r: usize) {
    if left < right {
        let pivot = partition(arr, left, right, cost, r);
        quicksort(arr, left, pivot - 1, cost, r);
        quicksort(arr, pivot + 1, right, cost, r);
    }
}

fn partition(
    arr: &mut Vec<usize>,
    left: isize,
    right: isize,
    cost: &Vec<Vec<f64>>,
    r: usize,
) -> isize {
    let pivot = arr[right as usize];
    let mut i = left - 1;
    for j in left..right {
        if cost[r][arr[j as usize]] < cost[r][pivot] {
            i += 1;
            arr.swap(i as usize, j as usize);
        }
    }
    arr.swap((i + 1) as usize, right as usize);
    i + 1
}

fn construction(alpha: f64, c: &Vec<Vec<f64>>, dimen: usize, rnd: &mut tRnd) -> Vec<usize> {
    let mut s = vec![0; 1];

    let mut c_list = vec![];
    for i in 1..dimen {
        c_list.push(i);
    }

    let mut r: usize = 0;
    while c_list.is_empty() == false {
        sort(&mut c_list, r, c);

        let index = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
        let c = c_list[index];
        r = c;
        c_list.remove(index);
        s.push(c);
    }
    s.push(0);
    return s;
}

fn swap(s: &mut Vec<usize>, i: usize, j: usize) {
    s.swap(i, j);
}

fn reverse(s: &mut Vec<usize>, i: usize, j: usize) {
    s[i..=j].reverse();
}

fn reinsert(s: &mut Vec<usize>, i: usize, j: usize, pos: usize) {
    if i < pos {
        for _k in i..=j {
            s.insert(pos, s[i]);
            s.remove(i);
        }
    } else {
        for _k in i..=j {
            let tmp = s[j];
            s.remove(j);
            s.insert(pos, tmp);
        }
    }
}

fn search_swap(
    s: &mut Vec<usize>,
    seq: &mut Vec<Vec<tSeqInfo>>,
    c: &Vec<Vec<f64>>,
    dimen: usize,
) -> bool {
    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;
    let mut cost_concat_3: f64;
    let mut cost_concat_4: f64;

    let mut cost_best: f64 = f64::MAX;
    let mut cost_new: f64; // = std::f64::MAX;
    let mut I: usize = 0;
    let mut J: usize = 0;

    for i in 1..dimen - 1 {
        let i_prev: usize = i - 1;
        let i_next: usize = i + 1;

        cost_concat_1 = seq[0][i_prev].T + c[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][i_next].T + c[s[i]][s[i_next + 1]];

        cost_new = seq[0][i_prev].C
            + seq[i][i_next].W * cost_concat_1
            + c[s[i_next]][s[i]]
            + seq[i_next + 1][dimen].W * cost_concat_2
            + seq[i_next + 1][dimen].C;

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next + 1..dimen {
            let j_next = j + 1;
            let j_prev = j - 1;

            cost_concat_1 = seq[0][i_prev].T + c[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + c[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev].T + c[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3 + c[s[i]][s[j_next]];

            cost_new = seq[0][i_prev].C
                + cost_concat_1
                + seq[i_next][j_prev].W * cost_concat_2
                + seq[i_next][j_prev].C
                + cost_concat_3
                + seq[j_next][dimen].W * cost_concat_4
                + seq[j_next][dimen].C;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < seq[0][dimen].C - f64::EPSILON {
        //println!("swap \n{}", cost_best);
        swap(s, I, J);
        update_subseq_info_matrix(s, seq, c, dimen);
        //update_subseq_info_matrix(s, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn search_two_opt(
    s: &mut Vec<usize>,
    seq: &mut Vec<Vec<tSeqInfo>>,
    c: &Vec<Vec<f64>>,
    dimen: usize,
) -> bool {
    let mut cost_new: f64;
    let mut cost_best: f64 = f64::MAX;

    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;

    let mut I: usize = 0;
    let mut J: usize = 0;

    for i in 1..dimen - 1 {
        let i_prev: usize = i - 1;
        let mut rev_seq_cost: f64 = seq[i][i + 1].T;
        //let mut rev_seq_cost : f64 = solut.seq.get_unchecked(i).get_unchecked( i+1).T;

        for j in i + 2..dimen {
            let j_next = j + 1;

            rev_seq_cost += c[s[j - 1]][s[j]] * (seq[i][j].W - 1.0);

            cost_concat_1 = seq[0][i_prev].T + c[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[i][j].T + c[s[j_next]][s[i]];

            cost_new = seq[0][i_prev].C
                + seq[i][j].W * cost_concat_1
                + rev_seq_cost
                + seq[j_next][dimen].W * cost_concat_2
                + seq[j_next][dimen].C;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < seq[0][dimen].C - f64::EPSILON {
        //println!("two_opt \n{}", cost_best);
        reverse(s, I, J);
        update_subseq_info_matrix(s, seq, c, dimen);
        //update_subseq_info_matrix(s, seq, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn search_reinsertion(
    s: &mut Vec<usize>,
    seq: &mut Vec<Vec<tSeqInfo>>,
    opt: usize,
    c: &Vec<Vec<f64>>,
    dimen: usize,
) -> bool {
    let mut cost_best: f64 = f64::MAX;
    let mut cost_new: f64;

    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;
    let mut cost_concat_3: f64;

    let mut I: usize = 0;
    let mut J: usize = 0;
    let mut POS: usize = 0;

    for i in 1..dimen - opt + 1 {
        let j: usize = opt + i - 1;
        let i_prev: usize = i - 1;
        let j_next: usize = j + 1;

        for k in 0..i_prev {
            let k_next: usize = k + 1;

            cost_concat_1 = seq[0][k].T + c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[i][j].T + c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev].T + c[s[i_prev]][s[j_next]];

            cost_new = seq[0][k].C                                                             /*        1st subseq */
                + seq[i][j].W              * cost_concat_1 + seq[i][j].C                  /* concat 2nd subseq (reinserted seq) */
                + seq[k_next][i_prev].W   * cost_concat_2 + seq[k_next][ i_prev].C        /* concat 3rd subseq */
                + seq[j_next][ dimen].W * cost_concat_3 + seq[j_next][ dimen].C; /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for k in i + opt..dimen {
            let k_next: usize = k + 1;
            cost_concat_1 = seq[0][i_prev].T + c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][k].T + c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[i][j].T + c[s[j]][s[k_next]];

            cost_new = seq[0][ i_prev].C                                                        /*      1st subseq */
                + seq[j_next][ k].W         * cost_concat_1 + seq[j_next][ k].C             /* concat 2nd subseq */
                + seq[i][ j].W              * cost_concat_2 + seq[i][ j].C                  /* concat 3rd subseq (reinserted seq) */
                + seq[k_next][ dimen].W * cost_concat_3 + seq[k_next][dimen].C; /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if cost_best < seq[0][dimen].C - f64::EPSILON {
        //println!("reinsertion {}   {} {} {}\n{}", opt,I,J, 1+POS, cost_best);
        reinsert(s, I, J, POS + 1);
        update_subseq_info_matrix(s, seq, c, dimen);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn RVND(
    solut: &mut Vec<usize>,
    seq: &mut Vec<Vec<tSeqInfo>>,
    c: &Vec<Vec<f64>>,
    dimen: usize,
    rnd: &mut tRnd,
) {
    const SWAP: usize = 0;
    const REINSERTION: usize = 1;
    const OR_OPT_2: usize = 2;
    const OR_OPT_3: usize = 3;
    const TWO_OPT: usize = 4;

    let mut neighbd_list = vec![SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
    let mut improv_flag = false;

    while neighbd_list.is_empty() == false {
        let index = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;

        let neighbd: usize = neighbd_list[index];

        improv_flag = false;

        if neighbd == SWAP {
            improv_flag = search_swap(solut, seq, c, dimen);
        } else if neighbd == TWO_OPT {
            improv_flag = search_two_opt(solut, seq, c, dimen);
        } else if neighbd == REINSERTION {
            improv_flag = search_reinsertion(solut, seq, REINSERTION, c, dimen);
        } else if neighbd == OR_OPT_2 {
            improv_flag = search_reinsertion(solut, seq, OR_OPT_2, c, dimen);
        } else if neighbd == OR_OPT_3 {
            improv_flag = search_reinsertion(solut, seq, OR_OPT_3, c, dimen);
        }

        if improv_flag {
            neighbd_list = vec![SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
        } else {
            neighbd_list.remove(index);
        }
    }
}

fn perturb(sl: &Vec<usize>, dimen: usize, rnd: &mut tRnd) -> Vec<usize> {
    let mut s = sl.clone();
    let mut A_start: usize = 1;
    let mut A_end: usize = 1;
    let mut B_start: usize = 1;
    let mut B_end: usize = 1;

    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {
        A_start = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
        A_end = A_start + rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;

        B_start = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
        B_end = B_start + rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
    }

    if A_start < B_start {
        reinsert(&mut s, B_start, B_end - 1, A_end);
        reinsert(&mut s, A_start, A_end - 1, B_end);
    } else {
        reinsert(&mut s, A_start, A_end - 1, B_end);
        reinsert(&mut s, B_start, B_end - 1, A_end);
    }

    return s;
}

fn GILS_RVND(
    Imax: usize,
    Iils: usize,
    R: [f64; 26],
    c: &Vec<Vec<f64>>,
    dimen: usize,
    rnd: &mut tRnd,
) {
    let mut seq: Vec<Vec<tSeqInfo>> = vec![
        vec![
            tSeqInfo {
                C: 0.0,
                W: 0.0,
                T: 0.0
            };
            dimen + 1
        ];
        dimen + 1
    ];

    let mut solut_best: Vec<usize> = vec![0; dimen];
    let mut solut_partial: Vec<usize> = vec![0; dimen];
    let mut solut_crnt: Vec<usize> = vec![0; dimen];

    let mut cost_best = f64::MAX;
    let mut cost_partial = f64::MAX;
    let mut cost_crnt = f64::MAX;

    for _i in 0..Imax {
        let r_value = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
        let alpha = R[r_value];

        println!("[+] Local Search {}", _i);
        solut_crnt = construction(alpha, &c, dimen, rnd);
        update_subseq_info_matrix(&solut_crnt, &mut seq, &c, dimen);
        cost_crnt = seq[0][dimen].C;
        println!("{:?}", solut_crnt);

        println!("\t[+] Constructing Inital Solution.. {}", cost_crnt);
        solut_partial = solut_crnt.clone();
        cost_partial = cost_crnt;

        println!("\t[+] Looking for the best Neighbor..");
        let mut iterILS: usize = 0;
        while iterILS < Iils {
            RVND(&mut solut_crnt, &mut seq, &c, dimen, rnd);
            cost_crnt = seq[0][dimen].C;
            if cost_crnt < cost_partial {
                solut_partial = solut_crnt.clone();
                cost_partial = cost_crnt;
                iterILS = 0;
            }

            solut_crnt = perturb(&solut_partial, dimen, rnd);
            update_subseq_info_matrix(&mut solut_crnt, &mut seq, &c, dimen);
            iterILS += 1;
        }

        if cost_partial < cost_best {
            solut_best = solut_partial.clone();
            cost_best = cost_partial;
        }

        println!("\tCurrent Best Cost {}", cost_best);
    }

    println!("{:?}", solut_best);
    println!("COST: {}", cost_best);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn main() {
    let mut dimension: usize = 0;
    let mut c: Vec<Vec<f64>> = vec![vec![]; 0];

    let mut rnd_vec: Vec<usize> = vec![];

    data::load(&mut dimension, &mut c, &mut rnd_vec);

    let mut rnd = tRnd {
        rnd: rnd_vec,
        rnd_index: 0,
    };

    let Imax = 10;
    let Iils = if dimension < 100 { dimension } else { 100 };

    let R = [
        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14,
        0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25,
    ];

    let now = Instant::now();

    GILS_RVND(Imax, Iils, R, &c, dimension, &mut rnd);

    let new_now = Instant::now();
    println!("TIME: {}", new_now.duration_since(now).as_secs_f64());
}
