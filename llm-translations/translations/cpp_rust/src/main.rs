mod data;
mod read_data;

use data::load_data;
use std::time::Instant;

#[derive(Clone, Copy)]
struct SeqInfo {
    t: f64,
    c: f64,
    w: f64,
}

impl Default for SeqInfo {
    fn default() -> Self {
        SeqInfo {
            t: 0.0,
            c: 0.0,
            w: 0.0,
        }
    }
}

struct AlgoData {
    cost: Vec<Vec<f64>>,
    dimen: i32,
    rnd: Vec<i32>,
    rnd_index: usize,
}

#[derive(Clone, Copy, PartialEq)]
enum Neighborhood {
    Reinsertion = 1,
    OrOpt2 = 2,
    OrOpt3 = 3,
    Swap = 4,
    TwoOpt = 5,
}

struct Solution {
    s: Vec<i32>,
    seq: Vec<SeqInfo>,
    cost: f64,
}

fn to_1d(x: i32, y: i32, size: i32) -> usize {
    (x * (size + 1) + y) as usize
}

fn solution_init(data: &AlgoData) -> Solution {
    let n = (data.dimen + 1) as usize;
    Solution {
        s: vec![0; n],
        seq: vec![SeqInfo::default(); n * n],
        cost: f64::MAX,
    }
}

fn solution_cpy(src: &Solution, tgt: &mut Solution) {
    tgt.s = src.s.clone();
    tgt.cost = src.cost;
}

fn r_table(i: usize) -> f64 {
    const TABLE: [f64; 26] = [
        0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14,
        0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25,
    ];
    TABLE[i]
}

#[allow(dead_code)]
fn print_s(s: &[i32]) {
    for v in s {
        print!("{} ", v + 1);
    }
    println!();
}

fn partition(arr: &mut Vec<i32>, left: i32, right: i32, data: &AlgoData, r: i32) -> i32 {
    let pivot = arr[right as usize];
    let mut i = left - 1;
    for j in left..right {
        if data.cost[r as usize][arr[j as usize] as usize] < data.cost[r as usize][pivot as usize]
        {
            i += 1;
            arr.swap(i as usize, j as usize);
        }
    }
    arr.swap((i + 1) as usize, right as usize);
    i + 1
}

fn quicksort(arr: &mut Vec<i32>, left: i32, right: i32, data: &AlgoData, r: i32) {
    if left < right {
        let pivot = partition(arr, left, right, data, r);
        quicksort(arr, left, pivot - 1, data, r);
        quicksort(arr, pivot + 1, right, data, r);
    }
}

fn sort_cl(arr: &mut Vec<i32>, r: i32, data: &AlgoData) {
    let right = arr.len() as i32 - 1;
    quicksort(arr, 0, right, data, r);
}

fn construct(_alpha: f64, data: &mut AlgoData) -> Vec<i32> {
    let mut s: Vec<i32> = vec![0];
    s.reserve(data.dimen as usize + 1);
    let mut c_l: Vec<i32> = (1..data.dimen).collect();

    let mut r: i32 = 0;
    while !c_l.is_empty() {
        sort_cl(&mut c_l, r, data);

        let index = data.rnd[data.rnd_index] as usize;
        data.rnd_index += 1;

        let c = c_l[index];
        s.push(c);
        r = c;
        c_l.remove(index);
    }

    s.push(0);
    s
}

fn vec_swap(vec: &mut [i32], i: usize, j: usize) {
    vec.swap(i, j);
}

fn vec_reverse(vec: &mut [i32], i: usize, j: usize) {
    vec[i..=j].reverse();
}

fn reinsert(vec: &mut Vec<i32>, i: usize, j: usize, pos: usize) {
    let seq: Vec<i32> = vec[i..=j].to_vec();
    if pos < i {
        vec.drain(i..=j);
        for (offset, v) in seq.into_iter().enumerate() {
            vec.insert(pos + offset, v);
        }
    } else {
        for (offset, v) in seq.iter().enumerate() {
            vec.insert(pos + offset, *v);
        }
        vec.drain(i..=j);
    }
}

fn update_subseq_info_matrix(solut: &mut Solution, data: &AlgoData, index: i32) {
    let dimen = data.dimen;
    let mut from = index;

    for i in 0..=dimen {
        let indicator = if i == 0 { 1 } else { 0 };
        let k = 1 - i - indicator;
        let t = i == from;

        let idx_ii = to_1d(i, i, dimen);
        solut.seq[idx_ii].t = 0.0;
        solut.seq[idx_ii].c = 0.0;
        solut.seq[idx_ii].w = if i == 0 { 0.0 } else { 1.0 };

        for j in (i + 1)..=dimen {
            let j_prev = j - 1;
            let idx_ij = to_1d(i, j, dimen);
            let idx_ijprev = to_1d(i, j_prev, dimen);

            let t_val = data.cost[solut.s[j_prev as usize] as usize][solut.s[j as usize] as usize]
                + solut.seq[idx_ijprev].t;
            solut.seq[idx_ij].t = t_val;
            solut.seq[idx_ij].c = t_val + solut.seq[idx_ijprev].c;
            solut.seq[idx_ij].w = (j + k) as f64;
        }

        from += if t { 1 } else { 0 };
    }

    solut.cost = solut.seq[to_1d(0, dimen, dimen)].c;
}

fn search_swap(solut: &mut Solution, data: &AlgoData) -> bool {
    let dimen = data.dimen;
    let mut cost_best = f64::MAX;
    let mut best_i = 0i32;
    let mut best_j = 0i32;

    for i in 1..(dimen - 1) {
        let i_prev = i - 1;
        let i_next = i + 1;

        let cost_concat_1 = solut.seq[to_1d(0, i_prev, dimen)].t
            + data.cost[solut.s[i_prev as usize] as usize][solut.s[i_next as usize] as usize];
        let cost_concat_2 = cost_concat_1
            + solut.seq[to_1d(i, i_next, dimen)].t
            + data.cost[solut.s[i as usize] as usize][solut.s[(i_next + 1) as usize] as usize];

        let mut cost_new = solut.seq[to_1d(0, i_prev, dimen)].c
            + solut.seq[to_1d(i, i_next, dimen)].w * cost_concat_1
            + data.cost[solut.s[i_next as usize] as usize][solut.s[i as usize] as usize]
            + solut.seq[to_1d(i_next + 1, dimen, dimen)].w * cost_concat_2
            + solut.seq[to_1d(i_next + 1, dimen, dimen)].c;

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            best_i = i;
            best_j = i_next;
        }

        for j in (i_next + 1)..dimen {
            let j_next = j + 1;
            let j_prev = j - 1;

            let cc1 = solut.seq[to_1d(0, i_prev, dimen)].t
                + data.cost[solut.s[i_prev as usize] as usize][solut.s[j as usize] as usize];
            let cc2 =
                cc1 + data.cost[solut.s[j as usize] as usize][solut.s[i_next as usize] as usize];
            let cc3 = cc2
                + solut.seq[to_1d(i_next, j_prev, dimen)].t
                + data.cost[solut.s[j_prev as usize] as usize][solut.s[i as usize] as usize];
            let cc4 =
                cc3 + data.cost[solut.s[i as usize] as usize][solut.s[j_next as usize] as usize];

            cost_new = solut.seq[to_1d(0, i_prev, dimen)].c
                + cc1
                + solut.seq[to_1d(i_next, j_prev, dimen)].w * cc2
                + solut.seq[to_1d(i_next, j_prev, dimen)].c
                + cc3
                + solut.seq[to_1d(j_next, dimen, dimen)].w * cc4
                + solut.seq[to_1d(j_next, dimen, dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                best_i = i;
                best_j = j;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        vec_swap(&mut solut.s, best_i as usize, best_j as usize);
        update_subseq_info_matrix(solut, data, best_i);
        true
    } else {
        false
    }
}

fn search_two_opt(solut: &mut Solution, data: &AlgoData) -> bool {
    let dimen = data.dimen;
    let mut cost_best = f64::MAX;
    let mut best_i = 0i32;
    let mut best_j = 0i32;

    for i in 1..(dimen - 1) {
        let i_prev = i - 1;
        let mut rev_seq_cost = solut.seq[to_1d(i, i + 1, dimen)].t;

        for j in (i + 2)..dimen {
            let j_next = j + 1;

            rev_seq_cost += data.cost[solut.s[(j - 1) as usize] as usize]
                [solut.s[j as usize] as usize]
                * (solut.seq[to_1d(i, j, dimen)].w - 1.0);

            let cc1 = solut.seq[to_1d(0, i_prev, dimen)].t
                + data.cost[solut.s[j as usize] as usize][solut.s[i_prev as usize] as usize];
            let cc2 = cc1
                + solut.seq[to_1d(i, j, dimen)].t
                + data.cost[solut.s[j_next as usize] as usize][solut.s[i as usize] as usize];

            let cost_new = solut.seq[to_1d(0, i_prev, dimen)].c
                + solut.seq[to_1d(i, j, dimen)].w * cc1
                + rev_seq_cost
                + solut.seq[to_1d(j_next, dimen, dimen)].w * cc2
                + solut.seq[to_1d(j_next, dimen, dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                best_i = i;
                best_j = j;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        vec_reverse(&mut solut.s, best_i as usize, best_j as usize);
        update_subseq_info_matrix(solut, data, 0);
        true
    } else {
        false
    }
}

fn search_reinsertion(solut: &mut Solution, data: &AlgoData, opt: i32) -> bool {
    let dimen = data.dimen;
    let mut cost_best = f64::MAX;
    let mut best_i = 0i32;
    let mut best_j = 0i32;
    let mut best_pos = 0i32;

    let mut i = 1;
    let mut j = opt + i - 1;
    while i < dimen - opt + 1 {
        let j_next = j + 1;
        let i_prev = i - 1;

        let mut k = 0;
        while k < i_prev {
            let k_next = k + 1;

            let cc1 = solut.seq[to_1d(0, k, dimen)].t
                + data.cost[solut.s[k as usize] as usize][solut.s[i as usize] as usize];
            let cc2 = cc1
                + solut.seq[to_1d(i, j, dimen)].t
                + data.cost[solut.s[j as usize] as usize][solut.s[k_next as usize] as usize];
            let cc3 = cc2
                + solut.seq[to_1d(k_next, i_prev, dimen)].t
                + data.cost[solut.s[i_prev as usize] as usize][solut.s[j_next as usize] as usize];

            let cost_new = solut.seq[to_1d(0, k, dimen)].c
                + solut.seq[to_1d(i, j, dimen)].w * cc1
                + solut.seq[to_1d(i, j, dimen)].c
                + solut.seq[to_1d(k_next, i_prev, dimen)].w * cc2
                + solut.seq[to_1d(k_next, i_prev, dimen)].c
                + solut.seq[to_1d(j_next, dimen, dimen)].w * cc3
                + solut.seq[to_1d(j_next, dimen, dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                best_i = i;
                best_j = j;
                best_pos = k;
            }
            k += 1;
        }

        let mut k = i + opt;
        while k < dimen {
            let k_next = k + 1;

            let cc1 = solut.seq[to_1d(0, i_prev, dimen)].t
                + data.cost[solut.s[i_prev as usize] as usize][solut.s[j_next as usize] as usize];
            let cc2 = cc1
                + solut.seq[to_1d(j_next, k, dimen)].t
                + data.cost[solut.s[k as usize] as usize][solut.s[i as usize] as usize];
            let cc3 = cc2
                + solut.seq[to_1d(i, j, dimen)].t
                + data.cost[solut.s[j as usize] as usize][solut.s[k_next as usize] as usize];

            let cost_new = solut.seq[to_1d(0, i_prev, dimen)].c
                + solut.seq[to_1d(j_next, k, dimen)].w * cc1
                + solut.seq[to_1d(j_next, k, dimen)].c
                + solut.seq[to_1d(i, j, dimen)].w * cc2
                + solut.seq[to_1d(i, j, dimen)].c
                + solut.seq[to_1d(k_next, dimen, dimen)].w * cc3
                + solut.seq[to_1d(k_next, dimen, dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                best_i = i;
                best_j = j;
                best_pos = k;
            }
            k += 1;
        }

        i += 1;
        j += 1;
    }

    if cost_best < solut.cost - f64::EPSILON {
        reinsert(&mut solut.s, best_i as usize, best_j as usize, (best_pos + 1) as usize);
        let idx = if best_i < best_pos + 1 { best_i } else { best_pos + 1 };
        update_subseq_info_matrix(solut, data, idx);
        true
    } else {
        false
    }
}

fn rvnd(solut: &mut Solution, data: &mut AlgoData) {
    let mut neighbd_list = vec![
        Neighborhood::Swap,
        Neighborhood::TwoOpt,
        Neighborhood::Reinsertion,
        Neighborhood::OrOpt2,
        Neighborhood::OrOpt3,
    ];

    while !neighbd_list.is_empty() {
        let index = data.rnd[data.rnd_index] as usize;
        data.rnd_index += 1;

        let neighbd = neighbd_list[index];
        let improve_flag = match neighbd {
            Neighborhood::Reinsertion => {
                search_reinsertion(solut, data, Neighborhood::Reinsertion as i32)
            }
            Neighborhood::OrOpt2 => search_reinsertion(solut, data, Neighborhood::OrOpt2 as i32),
            Neighborhood::OrOpt3 => search_reinsertion(solut, data, Neighborhood::OrOpt3 as i32),
            Neighborhood::Swap => search_swap(solut, data),
            Neighborhood::TwoOpt => search_two_opt(solut, data),
        };

        if improve_flag {
            neighbd_list = vec![
                Neighborhood::Swap,
                Neighborhood::TwoOpt,
                Neighborhood::Reinsertion,
                Neighborhood::OrOpt2,
                Neighborhood::OrOpt3,
            ];
        } else {
            neighbd_list.remove(index);
        }
    }
}

fn perturb(solut: &Solution, data: &mut AlgoData) -> Vec<i32> {
    let mut s = solut.s.clone();
    let mut a_start = 1i32;
    let mut a_end = 1i32;
    let mut b_start = 1i32;
    let mut b_end = 1i32;

    let size_max_raw = (data.dimen + 1) / 10;
    let _size_max = if size_max_raw >= 2 { size_max_raw } else { 2 };
    let _size_min = 2;

    while (a_start <= b_start && b_start <= a_end) || (b_start <= a_start && a_start <= b_end) {
        a_start = data.rnd[data.rnd_index];
        data.rnd_index += 1;
        a_end = a_start + data.rnd[data.rnd_index];
        data.rnd_index += 1;

        b_start = data.rnd[data.rnd_index];
        data.rnd_index += 1;
        b_end = b_start + data.rnd[data.rnd_index];
        data.rnd_index += 1;
    }

    if a_start < b_start {
        reinsert(&mut s, b_start as usize, (b_end - 1) as usize, a_end as usize);
        reinsert(&mut s, a_start as usize, (a_end - 1) as usize, b_end as usize);
    } else {
        reinsert(&mut s, a_start as usize, (a_end - 1) as usize, b_end as usize);
        reinsert(&mut s, b_start as usize, (b_end - 1) as usize, a_end as usize);
    }

    s
}

fn gils_rvnd(imax: i32, iils: i32, data: &mut AlgoData, verbose: bool) {
    let mut solut_partial = solution_init(data);
    let mut solut_crnt = solution_init(data);
    let mut solut_best = solution_init(data);

    for i in 0..imax {
        let aux = data.rnd[data.rnd_index] as usize;
        data.rnd_index += 1;

        let alpha = r_table(aux);

        if verbose {
            println!("[+] Search {}", i + 1);
            println!("\t[+] Constructing..");
        }

        solut_crnt.s = construct(alpha, data);
        update_subseq_info_matrix(&mut solut_crnt, data, 0);

        solution_cpy(&solut_crnt, &mut solut_partial);
        if verbose {
            println!("\t[+] Looking for the best Neighbor..");
            println!("\t    Construction Cost: {:.3}", solut_partial.cost);
        }

        let mut iter_ils = 0;
        while iter_ils < iils {
            rvnd(&mut solut_crnt, data);
            if solut_crnt.cost < solut_partial.cost - f64::EPSILON {
                solution_cpy(&solut_crnt, &mut solut_partial);
                iter_ils = 0;
            }

            solut_crnt.s = perturb(&solut_partial, data);
            update_subseq_info_matrix(&mut solut_crnt, data, 0);
            iter_ils += 1;
        }

        if solut_partial.cost < solut_best.cost - f64::EPSILON {
            solution_cpy(&solut_partial, &mut solut_best);
        }

        if verbose {
            println!("\tCurrent best cost: {}", solut_best.cost);

            print!("SOLUCAO: ");
            for v in &solut_best.s {
                print!("{} ", v);
            }
            println!();
        }
    }

    println!("COST: {:.2}", solut_best.cost);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let imax = 10;
    let mut verbose = false;

    for a in &args[1..] {
        if a == "-v" || a == "--verbose" {
            verbose = true;
        }
    }

    let mut rnd: Vec<i32> = Vec::new();
    let (cost, dimen) = load_data(&mut rnd);

    let mut data = AlgoData {
        cost,
        dimen,
        rnd: rnd.clone(),
        rnd_index: 0,
    };

    let iils = if data.dimen < 100 { data.dimen } else { 100 };

    let t1 = Instant::now();
    let start_cpu = Instant::now();
    gils_rvnd(imax, iils, &mut data, verbose);
    let cpu_time = start_cpu.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let duration_ms = t2.duration_since(t1).as_millis();

    let res = duration_ms as f64 / 10e2;
    println!("TIME: {}", cpu_time);
    println!("wall clock (s): {}", res);

    println!("Tamanho RND {}", rnd.len());
}
