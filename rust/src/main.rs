mod subseq;
mod data;

use subseq::*;

use std::time::Duration;
use cpu_time::ProcessTime;

fn update_subseq_info_matrix(solut: &mut Solution, info: &Data) {
    for i in 0..=info.dimen {
        let k: i32 = 1 - (i as i32) - if i == 0 { 1 } else { 0 };

        solut.seq[(i, i)].t = 0.0;
        solut.seq[(i, i)].c = 0.0;
        solut.seq[(i, i)].w = if i != 0 { 1.0 } else { 0.0 };

        for j in (i + 1)..=info.dimen {
            let j_prev: usize = j - 1;

            // Cálculo de 't' usando acesso direto à matriz de custos
            let t = info.c[(solut[j_prev], solut[j])] + solut.seq[(i, j_prev)].t;
            solut.seq[(i, j)].t = t;

            // Cálculo de 'c'
            let c = solut.seq[(i, j)].t + solut.seq[(i, j_prev)].c;
            solut.seq[(i, j)].c = c;

            // Cálculo de 'w'
            let w = (j as i32 + k) as f64;
            solut.seq[(i, j)].w = w;
        }
    }

    // Cálculo do custo total diretamente
    solut.cost = solut.seq[(0, info.dimen)].c;
}


fn sort(arr: &mut Vec<usize>, r: usize, info: &Data) {
    quicksort(arr, 0, arr.len() as isize - 1, info, r);
}

fn quicksort(arr: &mut Vec<usize>, left: isize, right: isize, info: &Data, r: usize) {
    if left < right {
        let pivot = partition(arr, left, right, info, r);
        quicksort(arr, left, pivot - 1, info, r);
        quicksort(arr, pivot + 1, right, info, r);
    }
}

fn partition(arr: &mut Vec<usize>, left: isize, right: isize, info: &Data, r: usize) -> isize {
    let pivot = arr[right as usize];
    let mut i = left - 1;
    for j in left..right {
        if info.c[(r, arr[j as usize])] < info.c[(r, pivot)] {
            i += 1;
            arr.swap(i as usize, j as usize);
        }
    }
    arr.swap((i + 1) as usize, right as usize);
    i + 1
}

fn construction(alpha : f64, info : &mut Data) -> Vec<usize> {
    let mut s = vec![0; 1];

    let mut c_list = vec![];
    for i in 1..info.dimen {
        c_list.push(i);
    }

    let mut r : usize = 0;
    while !c_list.is_empty() {
        sort(&mut c_list, r, & info);

        let index = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        let c = c_list[index];
        r = c;
        c_list.remove(index);
        s.push(c);
    }
    s.push(0);

    return s;
}

fn swap(s : &mut Vec<usize>, i : usize, j : usize){
    s.swap(i, j);
}

fn reverse(s : &mut Vec<usize>, i : usize, j : usize){
    s[i..=j].reverse();
}

fn reinsert(s : &mut Vec<usize>, i : usize, j : usize, pos : usize){
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

fn search_swap(solut : &mut Solution, info : & Data) -> bool {
    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;
    let mut cost_concat_3 : f64;
    let mut cost_concat_4 : f64;

    let mut cost_best : f64 = f64::MAX;
    let mut cost_new : f64;
    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimen - 1 {
        let i_prev: usize = i - 1;
        let i_next: usize = i + 1;

        cost_concat_1 = solut.seq[(0, i_prev)].t + info.c[(solut[i_prev], solut[i_next])];
        cost_concat_2 = cost_concat_1 + solut.seq[(i, i_next)].t + info.c[(solut[i], solut[i_next + 1])];

        cost_new = solut.seq[(0, i_prev)].c
            + solut.seq[(i, i_next)].w * cost_concat_1
            + info.c[(solut[i_next], solut[i])]
            + solut.seq[(i_next + 1, info.dimen)].w * cost_concat_2
            + solut.seq[(i_next + 1, info.dimen)].c;

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next + 1..info.dimen {
            let j_next = j + 1;
            let j_prev = j - 1;

            cost_concat_1 = solut.seq[(0, i_prev)].t + info.c[(solut[i_prev], solut[j])];
            cost_concat_2 = cost_concat_1 + info.c[(solut[j], solut[i_next])];
            cost_concat_3 = cost_concat_2 + solut.seq[(i_next, j_prev)].t + info.c[(solut[j_prev], solut[i])];
            cost_concat_4 = cost_concat_3 + info.c[(solut[i], solut[j_next])];

            cost_new = solut.seq[(0, i_prev)].c
                + cost_concat_1
                + solut.seq[(i_next, j_prev)].w * cost_concat_2
                + solut.seq[(i_next, j_prev)].c
                + cost_concat_3
                + solut.seq[(j_next, info.dimen)].w * cost_concat_4
                + solut.seq[(j_next, info.dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }


    if cost_best < solut.cost - f64::EPSILON {
        //println!("swap \n{}", cost_best);
        swap(&mut solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    } else {
        return false;
    }
}

fn search_two_opt(solut: &mut Solution, info: & Data) -> bool {
    let mut cost_new : f64;
    let mut cost_best : f64 = f64::MAX;

    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;

    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimen - 1 {
        let i_prev: usize = i - 1;
        let mut rev_seq_cost: f64 = solut.seq[(i, i + 1)].t;

        for j in i + 2..info.dimen {
            let j_next = j + 1;

            rev_seq_cost += info.c[(solut[j - 1], solut[j])] * (solut.seq[(i, j)].w - 1.0);

            cost_concat_1 = solut.seq[(0, i_prev)].t + info.c[(solut[j], solut[i_prev])];
            cost_concat_2 = cost_concat_1 + solut.seq[(i, j)].t + info.c[(solut[j_next], solut[i])];

            cost_new = solut.seq[(0, i_prev)].c
                + solut.seq[(i, j)].w * cost_concat_1
                + rev_seq_cost
                + solut.seq[(j_next, info.dimen)].w * cost_concat_2
                + solut.seq[(j_next, info.dimen)].c;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }


    if cost_best < solut.cost - f64::EPSILON {
        reverse(&mut solut.s, I, J);
        update_subseq_info_matrix(solut, info);
        return true;
    } else {
        return false;
    }
}

fn search_reinsertion(solut : &mut Solution, opt: usize, info: & Data) -> bool {
    let mut cost_best : f64 = f64::MAX;
    let mut cost_new : f64;

    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;
    let mut cost_concat_3 : f64;

    let mut I : usize = 0;
    let mut J : usize = 0;
    let mut POS : usize = 0;

    for i in 1..info.dimen-opt+1 {
        let j : usize = opt+i-1;
        let i_prev : usize = i-1;
        let j_next : usize = j+1;

        let seq_i_j_w = solut.seq[(i, j)].w;
        let seq_i_j_c = solut.seq[(i, j)].c;
        let seq_jnext_dimen_w = solut.seq[(j_next, info.dimen)].w;
        let seq_jnext_dimen_c = solut.seq[(j_next, info.dimen)].c;

        for k in 0..i_prev {
            let k_next: usize = k + 1;

           // Cálculo dos custos concatenados
            cost_concat_1 = solut.seq[(0, k)].t + info.c[(solut[k], solut[i])];
            cost_concat_2 = cost_concat_1 + solut.seq[(i, j)].t + info.c[(solut[j], solut[k_next])];
            cost_concat_3 = cost_concat_2 + solut.seq[(k_next, i_prev)].t + info.c[(solut[i_prev], solut[j_next])];

            // Cálculo do novo custo
            cost_new = solut.seq[(0, k)].c
                + seq_i_j_w * cost_concat_1
                + seq_i_j_c
                + solut.seq[(k_next, i_prev)].w * cost_concat_2
                + solut.seq[(k_next, i_prev)].c
                + seq_jnext_dimen_w * cost_concat_3
                + seq_jnext_dimen_c;

            // Verifica se o custo atual é melhor que o melhor encontrado
            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }


        for k in (i + opt)..info.dimen {
            let k_next: usize = k + 1;

            // Cálculo dos custos concatenados
            cost_concat_1 = solut.seq[(0, i_prev)].t + info.c[(solut[i_prev], solut[j_next])];
            cost_concat_2 = cost_concat_1 + solut.seq[(j_next, k)].t + info.c[(solut[k], solut[i])];
            cost_concat_3 = cost_concat_2 + solut.seq[(i, j)].t + info.c[(solut[j], solut[k_next])];

            // Cálculo do novo custo
            cost_new = solut.seq[(0, i_prev)].c
                + solut.seq[(j_next, k)].w * cost_concat_1
                + solut.seq[(j_next, k)].c
                + solut.seq[(i, j)].w * cost_concat_2
                + solut.seq[(i, j)].c
                + solut.seq[(k_next, info.dimen)].w * cost_concat_3
                + solut.seq[(k_next, info.dimen)].c;

            // Verifica se o custo atual é melhor que o melhor encontrado
            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

    }

    if cost_best < solut.cost - f64::EPSILON {
        //println!("reinsertion {}   {} {} {}\n{}", opt,I,J, 1+POS, cost_best);
        reinsert(&mut solut.s, I, J, POS+1);
        update_subseq_info_matrix(solut, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn rvnd(solut: &mut Solution, info: &mut Data) {

    let mut neighbd_list = vec![
        Operator::Swap,
        Operator::TwoOpt,
        Operator::Reinsertion,
        Operator::OrOpt2,
        Operator::OrOpt3,
    ];


    while neighbd_list.is_empty() == false {

        let index = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        let neighbd = neighbd_list[index];

        let improv_flag = match neighbd {
            Operator::Swap          => search_swap(solut, info),
            Operator::TwoOpt        => search_two_opt(solut, info),
            Operator::Reinsertion   => search_reinsertion(solut, Operator::Reinsertion as usize, info),
            Operator::OrOpt2        => search_reinsertion(solut, Operator::OrOpt2 as usize, info),
            Operator::OrOpt3        => search_reinsertion(solut, Operator::OrOpt3 as usize, info),
        };


        if improv_flag {
            neighbd_list.resize(5, Operator::Swap);
            neighbd_list.copy_from_slice(&[
                Operator::Swap,
                Operator::TwoOpt,
                Operator::Reinsertion,
                Operator::OrOpt2,
                Operator::OrOpt3,
            ]);

        } else {
            neighbd_list.remove(index);
        }

    }
}

fn perturb(sl: & Vec<usize>, info: &mut Data) -> Vec<usize> {
    let mut s = sl.clone();
    let mut a_start : usize = 1;
    let mut a_end : usize = 1;
    let mut b_start : usize = 1;
    let mut b_end : usize = 1;

    while (a_start <= b_start &&  b_start <= a_end) || (b_start <= a_start && a_start <= b_end) {

        a_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        a_end = a_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

        b_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        b_end = b_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

    }

    if a_start < b_start {
        reinsert(&mut s, b_start, b_end-1, a_end);
        reinsert(&mut s, a_start, a_end-1, b_end);
    } else {
        reinsert(&mut s, a_start, a_end-1, b_end);
        reinsert(&mut s, b_start, b_end-1, a_end);
    }

    return s;
}

fn gils_rvnd(imax : usize, iils : usize, r : [f64; 26], info: &mut Data) {

    let mut solut_best = Solution {
        seq : SubseqMatrix::new(info.dimen),
        s : vec![0; info.dimen],
        cost : f64::MAX,
    };

    let mut solut_partial = Solution {
        seq : SubseqMatrix::new(info.dimen),
        s : vec![0; info.dimen],
        cost : 0.0,
    };

    let mut solut_crnt = Solution {
        seq : SubseqMatrix::new(info.dimen),
        s : vec![0; info.dimen],
        cost : 0.0,
    };


    for _i in 0..imax {

        let r_value = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        let alpha = r[r_value];

        println!("[+] Local Search {}", _i);
        solut_crnt.s = construction(alpha, info);
        println!("{:?}", solut_crnt.s);

        update_subseq_info_matrix(&mut solut_crnt, info);
        println!("\t[+] Constructing Inital Solution.. {}", solut_crnt.cost);
        solut_partial.s = solut_crnt.s.clone();
        solut_partial.cost = solut_crnt.cost;

        println!("\t[+] Looking for the best Neighbor..");
        let mut iter_ils : usize = 0;
        while iter_ils < iils {
            rvnd(&mut solut_crnt, info);
            if solut_crnt.cost < solut_partial.cost {
                solut_partial.s = solut_crnt.s.clone();
                solut_partial.cost = solut_crnt.cost;
                iter_ils = 0;

            }

            solut_crnt.s = perturb(&solut_partial.s, info);
            update_subseq_info_matrix(&mut solut_crnt, info);
            iter_ils += 1;
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

#[allow(non_snake_case)]
fn main() {
    let filename = "../distance_matrix";

    let (dimension, c, rnd) = data::load(filename);

    let mut info = Data {
        dimen : dimension,
        c,
        rnd,
        rnd_index : 0,
    };

    let imax = 10;
    let iils = if dimension < 100 {dimension} else {100};

    let r = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];

    let start = ProcessTime::try_now().expect("Getting process time failed");

    gils_rvnd(imax, iils, r, &mut info);

    let cpu_time: Duration = start.try_elapsed().expect("Getting process time failed");
    println!("TIME: {:?}", cpu_time.as_secs_f64());

}
