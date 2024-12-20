mod data;
extern crate rand;

use rand::Rng;
use std::time::Instant;

#[derive(Debug, Clone)]
enum Info {
    T = 0,
    C = 1,
    W = 2,
}

#[derive(Debug, Clone)]
struct Data {
    c: Vec<Vec<f64>>,
    dimen: usize,
}

#[derive(Debug, Clone)]
struct Rnd {
    rnd: Vec<usize>,
    rnd_index: usize,
}

#[derive(Debug, Clone)]
enum Moves {
    SWAP(usize),
    REINSERTION(usize),
    OR_OPT_2(usize),
    OR_OPT_3(usize),
    TWO_OPT(usize),
}

#[derive(Debug, Clone)]
struct tSolution {
    seq: Vec<Vec<[f64; 3]>>,
    s: Vec<usize>,
    cost: f64,
}

fn update_subseq_info_matrix(solut: &mut tSolution, data: &Data) {
    for i in 0..data.dimen + 1 {
        let k: i32 = 1 - (i as i32) - if i == 0 { 1 } else { 0 };

        solut.seq[i][i][Info::T as usize] = 0.0;
        solut.seq[i][i][Info::C as usize] = 0.0;
        solut.seq[i][i][Info::W as usize] = if i != 0 { 1.0 } else { 0.0 };

        for j in (i + 1)..(data.dimen + 1) {
            let j_prev: usize = j - 1;

            solut.seq[i][j][Info::T as usize] =
                data.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][Info::T as usize];
            solut.seq[i][j][Info::C as usize] =
                solut.seq[i][j][Info::T as usize] + solut.seq[i][j_prev][Info::C as usize];
            solut.seq[i][j][Info::W as usize] = (j as i32 + k) as f64;
        }
    }

    solut.cost = solut.seq[0][data.dimen][Info::C as usize];
}

fn sort(arr: &mut Vec<usize>, r: usize, data: &Data) {
    quicksort(arr, 0, arr.len() as isize - 1, r, data);
}

fn quicksort(arr: &mut Vec<usize>, left: isize, right: isize, r: usize, data: &Data) {
    if left < right {
        let pivot = partition(arr, left, right, r, data);
        quicksort(arr, left, pivot - 1, r, data);
        quicksort(arr, pivot + 1, right, r, data);
    }
}

fn partition(arr: &mut Vec<usize>, left: isize, right: isize, r: usize, data: &Data) -> isize {
    let pivot = arr[right as usize];
    let mut i = left - 1;
    for j in left..right {
        if data.c[r][arr[j as usize]] < data.c[r][pivot] {
            i += 1;
            arr.swap(i as usize, j as usize);
        }
    }
    arr.swap((i + 1) as usize, right as usize);
    i + 1
}

fn construction(alpha: f64, rnd: &mut Rnd, data: &Data) -> Vec<usize> {
    let mut s = vec![0; 1];
    let mut c_list = vec![];

    for i in 1..data.dimen {
        c_list.push(i);
    }

    let mut r: usize = 0;
    while c_list.is_empty() == false {
        sort(&mut c_list, r, data);

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

fn search_swap(solut: &mut tSolution, data: &Data) -> bool {
    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;
    let mut cost_concat_3: f64;
    let mut cost_concat_4: f64;

    let mut cost_best: f64 = f64::MAX;
    let mut cost_new: f64; // = std::f64::MAX;
    let mut I: usize = 0;
    let mut J: usize = 0;

    for i in 1..data.dimen - 1 {
        let i_prev: usize = i - 1;
        let i_next: usize = i + 1;

        cost_concat_1 =
            solut.seq[0][i_prev][Info::T as usize] + data.c[solut.s[i_prev]][solut.s[i_next]];
        cost_concat_2 = cost_concat_1
            + solut.seq[i][i_next][Info::T as usize]
            + data.c[solut.s[i]][solut.s[i_next + 1]];

        cost_new = solut.seq[0][i_prev][Info::C as usize]
            + solut.seq[i][i_next][Info::W as usize] * cost_concat_1
            + data.c[solut.s[i_next]][solut.s[i]]
            + solut.seq[i_next + 1][data.dimen][Info::W as usize] * cost_concat_2
            + solut.seq[i_next + 1][data.dimen][Info::C as usize];

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next + 1..data.dimen {
            let j_next = j + 1;
            let j_prev = j - 1;

            cost_concat_1 =
                solut.seq[0][i_prev][Info::T as usize] + data.c[solut.s[i_prev]][solut.s[j]];
            cost_concat_2 = cost_concat_1 + data.c[solut.s[j]][solut.s[i_next]];
            cost_concat_3 = cost_concat_2
                + solut.seq[i_next][j_prev][Info::T as usize]
                + data.c[solut.s[j_prev]][solut.s[i]];
            cost_concat_4 = cost_concat_3 + data.c[solut.s[i]][solut.s[j_next]];

            cost_new = solut.seq[0][i_prev][Info::C as usize]
                + cost_concat_1
                + solut.seq[i_next][j_prev][Info::W as usize] * cost_concat_2
                + solut.seq[i_next][j_prev][Info::C as usize]
                + cost_concat_3
                + solut.seq[j_next][data.dimen][Info::W as usize] * cost_concat_4
                + solut.seq[j_next][data.dimen][Info::C as usize];

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        swap(&mut solut.s, I, J);
        update_subseq_info_matrix(solut, data);
        return true;
    } else {
        return false;
    }
}

fn search_two_opt(solut: &mut tSolution, data: &Data) -> bool {
    let mut cost_new: f64;
    let mut cost_best: f64 = f64::MAX;

    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;

    let mut I: usize = 0;
    let mut J: usize = 0;

    for i in 1..data.dimen - 1 {
        let i_prev: usize = i - 1;
        let mut rev_seq_cost: f64 = solut.seq[i][i + 1][Info::T as usize];

        for j in i + 2..data.dimen {
            let j_next = j + 1;

            rev_seq_cost +=
                data.c[solut.s[j - 1]][solut.s[j]] * (solut.seq[i][j][Info::W as usize] - 1.0);

            cost_concat_1 =
                solut.seq[0][i_prev][Info::T as usize] + data.c[solut.s[j]][solut.s[i_prev]];
            cost_concat_2 = cost_concat_1
                + solut.seq[i][j][Info::T as usize]
                + data.c[solut.s[j_next]][solut.s[i]];

            cost_new = solut.seq[0][i_prev][Info::C as usize]
                + solut.seq[i][j][Info::W as usize] * cost_concat_1
                + rev_seq_cost
                + solut.seq[j_next][data.dimen][Info::W as usize] * cost_concat_2
                + solut.seq[j_next][data.dimen][Info::C as usize];

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        reverse(&mut solut.s, I, J);
        update_subseq_info_matrix(solut, data);
        return true;
    } else {
        return false;
    }
}

fn search_reinsertion(solut: &mut tSolution, opt: usize, data: &Data) -> bool {
    let mut cost_best: f64 = f64::MAX;
    let mut cost_new: f64;

    let mut cost_concat_1: f64;
    let mut cost_concat_2: f64;
    let mut cost_concat_3: f64;

    let mut I: usize = 0;
    let mut J: usize = 0;
    let mut POS: usize = 0;

    for i in 1..data.dimen - opt + 1 {
        let j: usize = opt + i - 1;
        let i_prev: usize = i - 1;
        let j_next: usize = j + 1;

        for k in 0..i_prev {
            let k_next: usize = k + 1;

            cost_concat_1 = solut.seq[0][k][Info::T as usize] + data.c[solut.s[k]][solut.s[i]];
            cost_concat_2 = cost_concat_1
                + solut.seq[i][j][Info::T as usize]
                + data.c[solut.s[j]][solut.s[k_next]];
            cost_concat_3 = cost_concat_2
                + solut.seq[k_next][i_prev][Info::T as usize]
                + data.c[solut.s[i_prev]][solut.s[j_next]];

            cost_new = solut.seq[0][k][Info::C as usize]                                                             /*        1st subseq */
                + solut.seq[i][j][Info::W as usize]              * cost_concat_1 + solut.seq[i][j][Info::C as usize]                  /* concat 2nd subseq (reinserted seq) */
                + solut.seq[k_next][i_prev][Info::W as usize]   * cost_concat_2 + solut.seq[k_next][ i_prev][Info::C as usize]        /* concat 3rd subseq */
                + solut.seq[j_next][ data.dimen][Info::W as usize] * cost_concat_3 + solut.seq[j_next][ data.dimen][Info::C as usize]; /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for k in i + opt..data.dimen {
            let k_next: usize = k + 1;

            cost_concat_1 =
                solut.seq[0][i_prev][Info::T as usize] + data.c[solut.s[i_prev]][solut.s[j_next]];
            cost_concat_2 = cost_concat_1
                + solut.seq[j_next][k][Info::T as usize]
                + data.c[solut.s[k]][solut.s[i]];
            cost_concat_3 = cost_concat_2
                + solut.seq[i][j][Info::T as usize]
                + data.c[solut.s[j]][solut.s[k_next]];

            cost_new = solut.seq[0][ i_prev][Info::C as usize]                                                        /*      1st subseq */
                + solut.seq[j_next][ k][Info::W as usize]         * cost_concat_1 + solut.seq[j_next][ k][Info::C as usize]             /* concat 2nd subseq */
                + solut.seq[i][ j][Info::W as usize]              * cost_concat_2 + solut.seq[i][ j][Info::C as usize]                  /* concat 3rd subseq (reinserted seq) */
                + solut.seq[k_next][ data.dimen][Info::W as usize] * cost_concat_3 + solut.seq[k_next][ data.dimen][Info::C as usize]; /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        reinsert(&mut solut.s, I, J, POS + 1);
        update_subseq_info_matrix(solut, data);
        return true;
    } else {
        return false;
    }
}

fn RVND(solut: &mut tSolution, rnd: &mut Rnd, data: &Data) {
    let mut neighbd_list: Vec<Moves> = vec![
        Moves::SWAP(0),
        Moves::TWO_OPT(4),
        Moves::REINSERTION(1),
        Moves::OR_OPT_2(2),
        Moves::OR_OPT_3(3),
    ];
    let mut improv_flag: bool;

    while neighbd_list.is_empty() == false {
        let index: usize = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;

        let neighbd: &Moves = &neighbd_list[index];

        improv_flag = match *neighbd {
            Moves::SWAP(_) => search_swap(solut, data),
            Moves::TWO_OPT(_) => search_two_opt(solut, data),
            Moves::OR_OPT_2(opt) => search_reinsertion(solut, opt, data),
            Moves::OR_OPT_3(opt) => search_reinsertion(solut, opt, data),
            Moves::REINSERTION(opt) => search_reinsertion(solut, opt, data),
        };

        if improv_flag {
            neighbd_list = vec![
                Moves::SWAP(0),
                Moves::TWO_OPT(4),
                Moves::REINSERTION(1),
                Moves::OR_OPT_2(2),
                Moves::OR_OPT_3(3),
            ];
        } else {
            neighbd_list.remove(index);
        }

    }
}

fn perturb(sl: &Vec<usize>, rnd: &mut Rnd) -> Vec<usize> {
    let mut rng = rand::thread_rng();
    let mut s = sl.clone();
    let mut A_start: usize = 1;
    let mut A_end: usize = 1;
    let mut B_start: usize = 1;
    let mut B_end: usize = 1;

    let size_max = if (s.len() as f64 / 10.0) as usize >= 2 {
        (s.len() as f64 / 10.0) as usize
    } else {
        2
    };
    let size_min = 2;

    let mut range = size_min..size_max;

    if size_max == 2 {
        range = 0..1;
    }

    while (A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {
        A_start = rng.gen_range(1..s.len() - 1 - size_max);
        A_end = A_start + rng.gen_range(range.clone());

        B_start = rng.gen_range(1..s.len() - 1 - size_max);
        B_end = B_start + rng.gen_range(range.clone());

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

fn GILS_RVND(Imax: usize, Iils: usize, R: [f64; 26], rnd: &mut Rnd, data: &Data) {
    let mut solut_best = tSolution {
        seq: vec![vec![[0.0; 3]; data.dimen + 1]; data.dimen + 1],
        s: vec![0; data.dimen],
        cost: f64::MAX,
    };

    let mut solut_partial = tSolution {
        seq: vec![vec![[0.0; 3]; data.dimen + 1]; data.dimen + 1],
        s: vec![0; data.dimen],
        cost: 0.0,
    };

    let mut solut_crnt = tSolution {
        seq: vec![vec![[0.0; 3]; data.dimen + 1]; data.dimen + 1],
        s: vec![0; data.dimen],
        cost: 0.0,
    };

    for _i in 0..Imax {
        let r_value = rnd.rnd[rnd.rnd_index];
        rnd.rnd_index += 1;
        let alpha = R[r_value];

        println!("[+] Local Search {}", _i);
        solut_crnt.s = construction(alpha, rnd, data);
        println!("{:?}", solut_crnt.s);
        update_subseq_info_matrix(&mut solut_crnt, data);
        println!("\t[+] Constructing Inital Solution.. {}", solut_crnt.cost);
        solut_partial.s = solut_crnt.s.clone();
        solut_partial.cost = solut_crnt.cost;

        println!("\t[+] Looking for the best Neighbor..");
        let mut iterILS: usize = 0;
        while iterILS < Iils {
            RVND(&mut solut_crnt, rnd, data);
            if solut_crnt.cost < solut_partial.cost {
                solut_partial.s = solut_crnt.s.clone();
                solut_partial.cost = solut_crnt.cost;
                iterILS = 0;

            }

            solut_crnt.s = perturb(&solut_partial.s, rnd);
            update_subseq_info_matrix(&mut solut_crnt, data);
            iterILS += 1;
        }


        if solut_partial.cost < solut_best.cost {
            solut_best.s = solut_partial.s.clone();
            solut_best.cost = solut_partial.cost;
        }

        println!("\tCurrent Best Cost {}", solut_best.cost);
    }

    println!("{:?}", solut_best.s);
    println!("COST: {}", solut_best.cost);
}

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>())
}

fn main() {
    let mut dimension: usize = 0;
    let mut c: Vec<Vec<f64>> = vec![]; //vec![0.0; 0]];

    let mut rnd: Vec<usize> = vec![];

    data::load(&mut dimension, &mut c, &mut rnd);

    let data = Data {
        dimen: dimension,
        c: c.clone(),
    };

    let mut rnd = Rnd {
        rnd: rnd,
        rnd_index: 0,
    };

    //let mut a  = array![0; 500];

    //println!(" {:#?} , {}", a);

    //exit(0);

    //let mut temperature = Array3::<f64>::zeros((3, 4, 5));
    //temperature[[2, 2, 2]] += 0.5;

    //println!("{}", temperature);

    /*
    let mut a : Vec<f64> = vec![];

    for i in 0..(info.dimen+1) * (info.dimen+1) * 3 {
        a.push(0.0);
    }

    let mut solut_test = tSolution {
        seq : a.into_boxed_slice(),
        //seq : vec![vec![[0.0; 3]; info.dimen+1]; info.dimen+1],
        s : vec![0; info.dimen],
        cost : f64::MAX,
    };

    for i in 0..info.dimen {
        solut_test.s[i] = i;
    }
    solut_test.s.push(0);

    update_subseq_info_matrix(&mut solut_test);


    println!(" {:?}\nCost {}", solut_test.s, solut_test.cost);
    //exit(0);
    */

    println!("TEST");

    let Imax = 10;
    let Iils = if dimension < 100 { dimension } else { 100 };

    let R = [
        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14,
        0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25,
    ];

    let now = Instant::now();

    //test!();
    GILS_RVND(Imax, Iils, R, &mut rnd, &data);

    let new_now = Instant::now();
    println!("TIME: {}", new_now.duration_since(now).as_secs_f64());

    /*
    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", Info::C as usize[i][j]);
        }
        println!();
    }
    */
}
