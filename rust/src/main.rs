mod gils_rvnd;

use subseq::*;

use std::time::Instant;
use std::process::exit;
use std::cmp::Ordering;
use rand::Rng;

use std::time::Duration;
use cpu_time::ProcessTime;



fn subseq_load(solut : &mut tSolution, info : & tInfo) {
    for i in 0..info.dimen + 1 {
        let k : i32 = 1 - (i as i32) - if i == 0 {1} else {0};

        // Flat[x + HEIGHT* (y + WIDTH* z)] = Original[x, y, z]

        solut.seq.set_T(i, i, 0.0);
        solut.seq.set_C(i, i, 0.0);
        solut.seq.set_W(i, i, if i != 0 {1.0} else {0.0});
        //solut.seq[i][i].W = if i != 0 {1.0} else {0.0};

      //solut.seq[i][i][info.T] = 0.0;
      //solut.seq[i][i][info.C] = 0.0;
      //solut.seq[i][i][info.W] = if i != 0 {1.0} else {0.0};

        for j in (i+1)..(info.dimen + 1) {
            let j_prev : usize = j - 1;

            let T = info.c.get(solut.s.get(j_prev), solut.s.get(j)) + solut.seq.get_T(i, j_prev); 
            solut.seq.set_T(i, j, T);

            let C = solut.seq.get_T(i, j) + solut.seq.get_C(i, j_prev); 
            solut.seq.set_C(i, j, C);
            //solut.seq.get_C_mut(i, j) = C;

            let W = (j as i32 + k) as f64;
            solut.seq.set_W(i, j, W);

        }
    }

    solut.cost = solut.seq.get_C(0, info.dimen);
}


fn sort(arr : &mut Vec<usize>, r : usize, info : & tInfo) {
    for i in 0..arr.len() {
        for j in 0..(arr.len()-1) {
            if info.c.get(r, arr.get(j)) > info.c.get(r, arr.get(j+1)) {
                arr.swap(j, j+1);
            }
        }
    }
}

fn construction(alpha : f64, info : &mut tInfo) -> Vec<usize> {
    let mut s = vec![0; 1];

    //let mut c_list = vec![0; info.dimen -1];
    let mut c_list = vec![];
    for i in 1..info.dimen {
        c_list.push(i);
    }

    let mut rng = rand::thread_rng();
    let mut r : usize = 0;
    while c_list.is_empty() == false {
        sort(&mut c_list, r, & info);

        let range = (c_list.len() as f64 * alpha + 1.0) as usize;
        let mut index = rng.gen::<usize>() % range;
        index = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        let c = c_list[index];
        r = c;
        c_list.remove(index);
        s.push(c);
    }
    s.push(0);
    //println!("{:?}", s);

    return s;
}

fn swap(s : &mut Vec<usize>, i : usize, j : usize){
    //unchecked
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

fn search_swap(solut : &mut tSolution, info : & tInfo) -> bool {
    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;
    let mut cost_concat_3 : f64;
    let mut cost_concat_4 : f64;

    let mut cost_best : f64 = f64::MAX;
    let mut cost_new : f64;// = std::f64::MAX;
    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        let i_next : usize = i + 1;

        cost_concat_1 = solut.seq.get_T(0, i_prev) + info.c.get(solut.s.get(i_prev), solut.s.get(i_next));
        cost_concat_2 = cost_concat_1 + solut.seq.get_T(i, i_next) + info.c.get(solut.s.get(i), solut.s.get(i_next+1));

        cost_new = solut.seq.get_C(0, i_prev)
            + solut.seq.get_W(i, i_next)     * cost_concat_1 + info.c.get(solut.s.get(i_next), solut.s.get(i))
            + solut.seq.get_W(i_next+1, info.dimen) * cost_concat_2 + solut.seq.get_C(i_next+1, info.dimen);


        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next+1..info.dimen {
            let j_next = j+1;
            let j_prev = j-1;

            cost_concat_1 = solut.seq.get_T(0, i_prev) + info.c.get(solut.s.get(i_prev), solut.s.get(j));
            cost_concat_2 = cost_concat_1 + info.c.get(solut.s.get(j), solut.s.get(i_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_T(i_next, j_prev) + info.c.get(solut.s.get(j_prev), solut.s.get(i));
            cost_concat_4 = cost_concat_3  + info.c.get(solut.s.get(i), solut.s.get(j_next));

            cost_new = solut.seq.get_C(0, i_prev)
                    + cost_concat_1
                    + solut.seq.get_W(i_next, j_prev) * cost_concat_2 + solut.seq.get_C(i_next, j_prev)
                    + cost_concat_3
                    + solut.seq.get_W(j_next, info.dimen) * cost_concat_4 + solut.seq.get_C(j_next, info.dimen);



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
        subseq_load(solut, info);
        return true;
    } else {
        return false;
    }
}

fn search_two_opt(solut : &mut tSolution, info : & tInfo) -> bool {
    let mut cost_new : f64;
    let mut cost_best : f64 = f64::MAX;

    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;

    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        let mut rev_seq_cost : f64 = solut.seq.get_T(i, i+1);

        for j in i+2..info.dimen {
            let j_next = j + 1;

            rev_seq_cost += info.c.get(solut.s.get(j-1), solut.s.get(j)) * (solut.seq.get_W(i, j)-1.0);

            cost_concat_1 =  solut.seq.get_T(0, i_prev) + info.c.get(solut.s.get(j), solut.s.get(i_prev));
            cost_concat_2 = cost_concat_1 + solut.seq.get_T(i, j) + info.c.get(solut.s.get(j_next), solut.s.get(i));

            cost_new = solut.seq.get_C(0, i_prev)
                    + solut.seq.get_W(i, j)      * cost_concat_1 + rev_seq_cost
                    + solut.seq.get_W(j_next, info.dimen) * cost_concat_2 + solut.seq.get_C(j_next, info.dimen);

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < solut.cost - f64::EPSILON {
        //println!("two_opt \n{}", cost_best);
        reverse(&mut solut.s, I, J);
        subseq_load(solut, info);
        return true;
    } else {
        return false;
    }
}

fn search_reinsertion(solut : &mut tSolution, opt : usize, info : & tInfo) -> bool {
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

        let seq_i_j_W = solut.seq.get_W(i, j);
        let seq_i_j_C = solut.seq.get_C(i, j);
        let seq_jnext_dimen_W = solut.seq.get_W(j_next, info.dimen);
        let seq_jnext_dimen_C = solut.seq.get_C(j_next, info.dimen);

        for k in 0..i_prev {
            let k_next : usize = k+1;

            cost_concat_1 = solut.seq.get_T(0, k) + info.c.get(solut.s.get(k), solut.s.get(i));
            cost_concat_2 = cost_concat_1 + solut.seq.get_T(i, j) + info.c.get(solut.s.get(j), solut.s.get(k_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_T(k_next, i_prev) + info.c.get(solut.s.get(i_prev), solut.s.get(j_next));

            cost_new = solut.seq.get_C(0, k)                                                             /*        1st subseq */
                + seq_i_j_W              * cost_concat_1 + seq_i_j_C                  /* concat 2nd subseq (reinserted seq) */
                + solut.seq.get_W(k_next, i_prev)   * cost_concat_2 + solut.seq.get_C(k_next, i_prev)        /* concat 3rd subseq */
                + seq_jnext_dimen_W * cost_concat_3 + seq_jnext_dimen_C;    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for k in i+opt..info.dimen {
            let k_next : usize = k+1;

            cost_concat_1 = solut.seq.get_T(0, i_prev) + info.c.get(solut.s.get(i_prev), solut.s.get(j_next));
            cost_concat_2 = cost_concat_1 + solut.seq.get_T(j_next, k) + info.c.get(solut.s.get(k), solut.s.get(i));
            cost_concat_3 = cost_concat_2 + solut.seq.get_T(i, j) + info.c.get(solut.s.get(j), solut.s.get(k_next));

            cost_new = solut.seq.get_C(0, i_prev)                                                        /*      1st subseq */
                + solut.seq.get_W(j_next, k)         * cost_concat_1 + solut.seq.get_C(j_next, k)             /* concat 2nd subseq */
                + solut.seq.get_W(i, j)              * cost_concat_2 + solut.seq.get_C(i, j)                  /* concat 3rd subseq (reinserted seq) */
                + solut.seq.get_W(k_next, info.dimen) * cost_concat_3 + solut.seq.get_C(k_next, info.dimen);    /* concat 4th subseq */

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
        subseq_load(solut, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn RVND(solut : &mut tSolution, info : &mut tInfo) {

    //let mut neighbd_list = vec![REINSERTION, OR_OPT_2, OR_OPT_3];
    let mut neighbd_list = vec![info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3];
    let mut improv_flag = false;

    let mut rng = rand::thread_rng();
    while neighbd_list.is_empty() == false {
        let mut index : usize = rng.gen::<usize>() % neighbd_list.len();

        index = info.rnd[info.rnd_index];
        info.rnd_index += 1;

        //println!("{} {}", info.rnd_index, index);

        let neighbd : usize = neighbd_list[index];

        improv_flag = false;

        if neighbd == info.SWAP {
            improv_flag = search_swap(solut, info);
        } else if neighbd == info.TWO_OPT {
            improv_flag = search_two_opt(solut, info);
        } else if neighbd == info.REINSERTION {
            improv_flag = search_reinsertion(solut, info.REINSERTION, info);
        } else if neighbd == info.OR_OPT_2 {
            improv_flag = search_reinsertion(solut,  info.OR_OPT_2, info);
        } else if neighbd == info.OR_OPT_3 {
            improv_flag = search_reinsertion(solut,  info.OR_OPT_3, info);
        }

        if improv_flag {
            neighbd_list = vec![info.SWAP, info.TWO_OPT, info.REINSERTION, info.OR_OPT_2, info.OR_OPT_3];
        } else {
            neighbd_list.remove(index);
        }

    }
}

fn perturb(sl : & Vec<usize>, info : &mut tInfo) -> Vec<usize> {
    let mut rng = rand::thread_rng();
    let mut s = sl.clone();
    let mut A_start : usize = 1;
    let mut A_end : usize = 1;
    let mut B_start : usize = 1;
    let mut B_end : usize = 1;

    let size_max = if (s.len() as f64 / 10.0) as usize >= 2 {(s.len() as f64 / 10.0) as usize} else {2};
    let size_min = 2;

    let mut range = size_min..size_max;

    if size_max == 2 {
        range = 0..1;
    }

    while (A_start <= B_start &&  B_start <= A_end) || (B_start <= A_start && A_start <= B_end) {
        A_start = rng.gen_range(1.. s.len() - 1 - size_max);
        A_end = A_start + rng.gen_range(range.clone());

        B_start = rng.gen_range(1.. s.len() - 1 - size_max);
        B_end = B_start + rng.gen_range(range.clone());

        A_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        A_end = A_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

        B_start = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        B_end = B_start + info.rnd[info.rnd_index];
        info.rnd_index += 1;

    }

    if A_start < B_start {
        reinsert(&mut s, B_start, B_end-1, A_end);
        reinsert(&mut s, A_start, A_end-1, B_end);
    } else {
        reinsert(&mut s, A_start, A_end-1, B_end);
        reinsert(&mut s, B_start, B_end-1, A_end);
    }

    return s;
}

fn GILS_RVND(Imax : usize, Iils : usize, R : [f64; 26], info : &mut tInfo) {

    let mut solut_best = tSolution {
        seq : tSeqData::New(),
        s : vec![0; info.dimen],
        cost : f64::MAX,
    };

    let mut solut_partial = tSolution {
        seq : tSeqData::New(),
        s : vec![0; info.dimen],
        cost : 0.0,
    };

    let mut solut_crnt = tSolution {
        seq : tSeqData::New(),
        s : vec![0; info.dimen],
        cost : 0.0,
    };


    let mut rng = rand::thread_rng();
    for _i in 0..Imax {
        let mut alpha : f64 = R[rng.gen::<usize>() % 26];

        let r_value = info.rnd[info.rnd_index];
        info.rnd_index += 1;
        alpha = R[r_value];

        println!("[+] Local Search {}", _i);
        solut_crnt.s = construction(alpha, info);
        println!("{:?}", solut_crnt.s);

        subseq_load(&mut solut_crnt, info);
        println!("\t[+] Constructing Inital Solution.. {}", solut_crnt.cost);
        solut_partial.s = solut_crnt.s.clone();
        solut_partial.cost = solut_crnt.cost;

        println!("\t[+] Looking for the best Neighbor..");
        let mut iterILS : usize = 0;
        while iterILS < Iils {
            RVND(&mut solut_crnt, info);
            if solut_crnt.cost < solut_partial.cost {
                solut_partial.s = solut_crnt.s.clone();
                solut_partial.cost = solut_crnt.cost;
                iterILS = 0;

            }

            solut_crnt.s = perturb(&solut_partial.s, info);
            subseq_load(&mut solut_crnt, info);
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

#[allow(non_snake_case)]
fn main() {
    let filename = "../distance_matrix";

    let (dimension, c, rnd) = data::load(filename);

    let mut info = tInfo {
        dimen : dimension,
        c : c.clone(),
      //T : 0,
      //C : 1,
      //W : 2,
        SWAP : 0,
        REINSERTION : 1,
        OR_OPT_2 : 2,
        OR_OPT_3 : 3,
        TWO_OPT : 4,
        rnd : rnd,
        rnd_index : 0,
    };

    let Imax = 10;
    let Iils = if dimension < 100 {dimension} else {100};

    let R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];

    let start = ProcessTime::try_now().expect("Getting process time failed");

    GILS_RVND(Imax, Iils, R, &mut info);

    let cpu_time: Duration = start.try_elapsed().expect("Getting process time failed");
    println!("TIME: {:?}", cpu_time.as_secs_f64());

}
