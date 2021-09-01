mod Data;
extern crate rand;
use std::time::Instant;
use std::process::exit;
use rand::Rng;
use std::cmp::Ordering;

static T : usize = 0;
static C : usize = 1;
static W : usize = 2;

static SWAP         : usize = 0;
static REINSERTION  : usize = 1;
static OR_OPT_2     : usize = 2;
static OR_OPT_3     : usize = 3;
static TWO_OPT      : usize = 4;

// 3d matrix
type SeqStruct = Vec<Vec<Vec<f64>>>;

#[derive(Debug)]
struct InstanceInfo {
    dimension : usize,
    c : Vec<Vec<f64>>,
}

fn subseq_load(s : & Vec<usize>, seq : &mut SeqStruct, info : & InstanceInfo) {
    for i in 0..info.dimension + 1 as usize{
        let k : i32 = 1 - (i as i32) - if i == 0 {1} else {0};

        seq[i][i][T] = 0.0;
        seq[i][i][C] = 0.0;
        seq[i][i][W] = if i != 0 {1.0} else {0.0};

        for j in (i+1)..(info.dimension + 1) as usize {
            let j_prev : usize = j - 1;

            seq[i][j][T] = info.c[s[j_prev]][s[j]] + seq[i][j_prev][T];
            seq[i][j][C] = seq[i][j][T] + seq[i][j_prev][C];
            seq[i][j][W] = (j as i32 + k) as f64;
        }
    }
}

fn construction(alpha : f64, info : & InstanceInfo) -> Vec<usize> {
    let mut s = vec![0; 1];

    //let mut cList = vec![0; info.dimension -1];
    let mut cList = vec![];
    for i in 1..info.dimension {
        cList.push(i);
    }

    let mut rng = rand::thread_rng();
    let mut r : usize = 0;
    while cList.is_empty() == false {
        cList.sort_by(|i, j| 
                    if info.c[*i as usize][r] < info.c[*j as usize][r] 
                        {Ordering::Less} 
                    else 
                        {Ordering::Greater});

        let range = (cList.len() as f64 * alpha + 1.0) as usize;
        let index = rng.gen::<usize>() % range;
        let c = cList[index];
        r = c;
        cList.remove(index);
        s.push(c);
    }
    s.push(0);
    //println!("{:?}", s);

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
        for k in i..=j {
            s.insert(pos, s[i]);
            s.remove(i);
        }
    } else {
        for k in i..=j {
            let tmp = s[j];
            s.remove(j);
            s.insert(pos, tmp);
        }
    }
}

fn search_swap(s : &mut Vec<usize>, seq : &mut SeqStruct, info : & InstanceInfo) -> bool {
    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;
    let mut cost_concat_3 : f64;
    let mut cost_concat_4 : f64;

    let mut cost_best : f64 = f64::MAX;
    let mut cost_new : f64;// = std::f64::MAX;
    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimension-1 {
        let mut i_prev : usize = i - 1;
        let mut i_next : usize = i + 1;

        cost_concat_1 = seq[0][i_prev][T] + info.c[s[i_prev]][s[i_next]];
        cost_concat_2 = cost_concat_1 + seq[i][i_next][T] + info.c[s[i]][s[i_next+1]];

        cost_new = seq[0][i_prev][C]
            + seq[i][i_next][W]     * cost_concat_1 + info.c[s[i_next]][s[i]]
            + seq[i_next+1][info.dimension][W] * cost_concat_2 + seq[i_next+1][info.dimension][C];

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next+1..info.dimension {
            let mut j_next = j+1;
            let mut j_prev = j-1;

            cost_concat_1 = seq[0][i_prev][T] + info.c[s[i_prev]][s[j]];
            cost_concat_2 = cost_concat_1 + info.c[s[j]][s[i_next]];
            cost_concat_3 = cost_concat_2 + seq[i_next][j_prev][T] + info.c[s[j_prev]][s[i]];
            cost_concat_4 = cost_concat_3  + info.c[s[i]][s[j_next]];

            cost_new = seq[0][i_prev][C]
                    + cost_concat_1
                    + seq[i_next][j_prev][W] * cost_concat_2 + seq[i_next][j_prev][C]
                    + cost_concat_3
                    + seq[j_next][info.dimension][W] * cost_concat_4 + seq[j_next][info.dimension][C];

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < seq[0][info.dimension][C] - f64::EPSILON {
        //println!("swap \n{}", cost_best);
        swap(s, I, J);
        subseq_load(s, seq, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn search_two_opt(s : &mut Vec<usize>, seq : &mut SeqStruct, info : & InstanceInfo) -> bool {
    let mut cost_new : f64;
    let mut cost_best : f64 = f64::MAX;

    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;

    let mut I : usize = 0;
    let mut J : usize = 0;

    for i in 1..info.dimension-1 {
        let i_prev : usize = i - 1;
        let mut rev_seq_cost : f64 = seq[i][i+1][T];

        for j in i+2..info.dimension {
            let j_next = j + 1;

            rev_seq_cost += info.c[s[j-1]][s[j]] * (seq[i][j][W]-1.0);

            cost_concat_1 =  seq[0][i_prev][T] + info.c[s[j]][s[i_prev]];
            cost_concat_2 = cost_concat_1 + seq[i][j][T] + info.c[s[j_next]][s[i]];

            cost_new = seq[0][i_prev][C]
                    + seq[i][j][W]      * cost_concat_1 + rev_seq_cost
                    + seq[j_next][info.dimension][W] * cost_concat_2 + seq[j_next][info.dimension][C];

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }

    if cost_best < seq[0][info.dimension][C] - f64::EPSILON {
        //println!("two_opt \n{}", cost_best);
        reverse(s, I, J);
        subseq_load(s, seq, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}


fn search_reinsertion(s : &mut Vec<usize>, seq : &mut SeqStruct, opt : usize, info : & InstanceInfo) -> bool {
    let mut cost_best : f64 = f64::MAX;
    let mut cost_new : f64 = f64::MAX;

    let mut cost_concat_1 : f64;
    let mut cost_concat_2 : f64;
    let mut cost_concat_3 : f64;

    let mut I : usize = 0;
    let mut J : usize = 0;
    let mut POS : usize = 0;

    for i in 1..info.dimension-opt+1 {
        let j : usize = opt+i-1;
        let i_prev : usize = i-1;
        let j_next : usize = j+1;

        for k in 0..i_prev {
            let k_next : usize = k+1;

            cost_concat_1 = seq[0][k][T] + info.c[s[k]][s[i]];
            cost_concat_2 = cost_concat_1 + seq[i][j][T] + info.c[s[j]][s[k_next]];
            cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][T] + info.c[s[i_prev]][s[j_next]];

            cost_new = seq[0][k][C]                                                             /*        1st subseq */
                + seq[i][j][W]              * cost_concat_1 + seq[i][j][C]                  /* concat 2nd subseq (reinserted seq) */
                + seq[k_next][i_prev][W]    * cost_concat_2 + seq[k_next][i_prev][C]        /* concat 3rd subseq */
                + seq[j_next][info.dimension][W] * cost_concat_3 + seq[j_next][info.dimension][C];    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for k in i+opt..info.dimension-opt-1 {
            let k_next : usize = k+1;

            cost_concat_1 = seq[0][i_prev][T] + info.c[s[i_prev]][s[j_next]];
            cost_concat_2 = cost_concat_1 + seq[j_next][k][T] + info.c[s[k]][s[i]];
            cost_concat_3 = cost_concat_2 + seq[i][j][T] + info.c[s[j]][s[k_next]];

            cost_new = seq[0][i_prev][C]                                                        /*      1st subseq */
                + seq[j_next][k][W]         * cost_concat_1 + seq[j_next][k][C]             /* concat 2nd subseq */
                + seq[i][j][W]              * cost_concat_2 + seq[i][j][C]                  /* concat 3rd subseq (reinserted seq) */
                + seq[k_next][info.dimension][W] * cost_concat_3 + seq[k_next][info.dimension][C];    /* concat 4th subseq */

            if(cost_new < cost_best){
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }

    if cost_best < seq[0][info.dimension][C] - f64::EPSILON {
        //println!("reinsertion {}   {} {} {}\n{}", opt,I,J, 1+POS, cost_best);
        reinsert(s, I, J, POS+1);
        subseq_load(s, seq, info);
        //println!("{}", seq[0][info.dimension][C]);
        return true;
    } else {
        return false;
    }
}

fn RVND(s : &mut Vec<usize>, seq : &mut SeqStruct, info : & InstanceInfo) {

    //let mut neighbd_list = vec![REINSERTION, OR_OPT_2, OR_OPT_3];
    let mut neighbd_list = vec![SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
    let mut improv_flag = false;

    /*
    let mut s = &mut vec![];
    for i in 0..info.dimension {
        s.push(i);
    }
    s.push(0);
    subseq_load(s, seq, info);
    */

    let mut rng = rand::thread_rng();
    while neighbd_list.is_empty() == false {
        //println!("is empty {} \n", neighbd_list.is_empty());
        let index : usize = rng.gen::<usize>() % neighbd_list.len();
        let neighbd : usize = neighbd_list[index];

        //println!("{:?} \n", s);

        if neighbd == SWAP {
            improv_flag = search_swap(s, seq, info);
        } else if neighbd == TWO_OPT {
            improv_flag = search_two_opt(s, seq, info);
        } else if neighbd == REINSERTION {
            improv_flag = search_reinsertion(s, seq, REINSERTION, info);
        } else if neighbd == OR_OPT_2 {
            improv_flag = search_reinsertion(s, seq, OR_OPT_2, info);
        } else if neighbd == OR_OPT_3 {
            improv_flag = search_reinsertion(s, seq, OR_OPT_3, info);
        }

        if improv_flag {
            neighbd_list = vec![SWAP, TWO_OPT, REINSERTION, OR_OPT_2, OR_OPT_3];
        } else {
            neighbd_list.remove(index);
        }

        //println!("{:?} \n{}\n", s, seq[0][info.dimension][C]);
        //println!("{:?} \n", s);
        //exit(1);
    }
}

fn perturb(sl : & Vec<usize>) -> Vec<usize> {
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

fn GILS_RVND(Imax : usize, Iils : usize, R : [f64; 26], info : & InstanceInfo) {

    let mut subseq : SeqStruct = 
        vec![
            vec![
                vec![0.0; 3]; 
            info.dimension+1]; 
        info.dimension+1];

    let mut sl_cost;
    let mut cost_best = f64::MAX;
    let mut s_best = vec![];

    //println!("{:#?} \n", subseq);
    let mut rng = rand::thread_rng();
    for i in 0..Imax {
        let alpha : f64 = R[rng.gen::<usize>() % 26];

        println!("Construction");
        let mut s : Vec<usize> = construction(alpha, info);
        let mut sl = s.clone();

        subseq_load(&s, &mut subseq, info);
        let mut rvnd_cost_best : f64 = subseq[0][info.dimension][C];
        //println!("Construction cost {} \n", rvnd_cost_best);

        println!("\tLocal search");
        let mut iterILS : usize = 0;
        while iterILS < Iils {
            RVND(&mut s, &mut subseq, info);
            //exit(0);
            let rvnd_cost_crnt = subseq[0][info.dimension][C];
            if rvnd_cost_crnt < rvnd_cost_best {
                rvnd_cost_best = rvnd_cost_crnt;
                sl = s.clone();
                iterILS = 0;
            }

            s = perturb(&sl);
            subseq_load(&s, &mut subseq, info);
            iterILS += 1;
        }

        subseq_load(&sl, &mut subseq, info);
        sl_cost = subseq[0][info.dimension][C] - f64::EPSILON;

        if sl_cost < cost_best {
            s_best = sl.clone();
            cost_best = sl_cost;
        }

        println!("\tCurrent Best Cost {}", cost_best);
    }


    println!("{:?}", s_best);
    println!("COST {}", cost_best);
}

fn main() {
    let mut dimension : usize = 0;
    let mut c : Vec<Vec<f64>> = vec![vec![];0];

    Data::load(&mut dimension, &mut c);

    let mut info = InstanceInfo {
        dimension : dimension,
        c : c,
    };

    let Imax = 10;
    let Iils = if dimension < 100 {dimension} else {dimension};;

    let R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];

    let now = Instant::now();

    GILS_RVND(Imax, Iils, R, &info);

    let new_now = Instant::now();
    println!("{:?}", new_now.duration_since(now));

    /*
    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", info.c[i][j]);
        }
        println!();
    }
    */

}
