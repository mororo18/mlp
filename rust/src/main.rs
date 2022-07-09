//mod sz;
mod data;
extern crate rand;

use data::sz;
use ndarray::array;
use std::time::Instant;
use std::process::exit;
use rand::Rng;
use std::cmp::Ordering;

#[derive(Debug, Clone)]
struct tInfo {
    c : Box<[[f64; sz::SIZE]; sz::SIZE]>,
    //c : Vec<Vec<f64>>,
    dimen : usize,
  //T : usize,
  //C : usize,
  //W : usize,
    SWAP         : usize ,
    REINSERTION  : usize ,
    OR_OPT_2     : usize ,
    OR_OPT_3     : usize ,
    TWO_OPT      : usize ,
    rnd : Vec<usize>,
    rnd_index : usize,
}

#[derive(Debug, Copy, Clone)]
struct tSeqInfo {
    T : f64,
    C : f64,
    W : f64,
}

#[derive(Debug, Clone)]
struct tSolution {
    //seq : Box<[f64]>,
    //seq : Vec<f64>,
    //seq : Vec<tSeqInfo>,
    seq : Box<[[tSeqInfo; sz::SIZE+1]; sz::SIZE+1]>,
    //seq : Vec<Vec<tSeqInfo>>,
    s : Vec<usize>,
    cost : f64,
}

#[inline(always)]
unsafe
fn to_1D(x : usize, y : usize, size : usize) -> usize {

    return x * (size+1) + y;
    //return x + (size+1) * (y + (size+1) * z);
}

fn subseq_load(solut : &mut tSolution, info : & tInfo) {
    for i in 0..info.dimen + 1 {
        let k : i32 = 1 - (i as i32) - if i == 0 {1} else {0};

        // Flat[x + HEIGHT* (y + WIDTH* z)] = Original[x, y, z]

        unsafe {
        solut.seq[i][i].T = 0.0;
        solut.seq[i][i].C = 0.0;
        solut.seq[i][i].W = if i != 0 {1.0} else {0.0};
        }

      //solut.seq[i][i][info.T] = 0.0;
      //solut.seq[i][i][info.C] = 0.0;
      //solut.seq[i][i][info.W] = if i != 0 {1.0} else {0.0};

      //solut.seq[i][i][info.T] = 0.0;
      //solut.seq[i][i][info.C] = 0.0;
      //solut.seq[i][i][info.W] = if i != 0 {1.0} else {0.0};

        //println!(" i = {}", i);
        //
        for j in (i+1)..(info.dimen + 1) {
            let j_prev : usize = j - 1;

            //println!(" index[{}, {}, {}] = {}", i, j, info.T, to_1D(i, j, info.T, info));
            unsafe {
          //solut.seq[to_1D(i, j, info.T, info)] = info.c[solut.s[j_prev]][solut.s[j]] + solut.seq[to_1D(i, j_prev, info.T, info)];
          //solut.seq[to_1D(i, j, info.C, info)] = solut.seq[to_1D(i, j, info.T, info)] + solut.seq[to_1D(i, j_prev, info.C, info)];
          //solut.seq[to_1D(i, j, info.W, info)] = (j as i32 + k) as f64;
            solut.seq.get_unchecked_mut(i).get_unchecked_mut(j).T = info.c.get_unchecked(*solut.s.get_unchecked(j_prev)).get_unchecked(*solut.s.get_unchecked(j)) + solut.seq.get_unchecked(i).get_unchecked( j_prev).T ; 
            solut.seq.get_unchecked_mut(i).get_unchecked_mut( j).C = solut.seq.get_unchecked(i).get_unchecked(j).T + solut.seq.get_unchecked(i).get_unchecked( j_prev).C; 
            solut.seq.get_unchecked_mut(i).get_unchecked_mut( j).W = (j as i32 + k) as f64;
            }

            /*
            unsafe{
            *solut.seq.get_unchecked_mut(i).get_unchecked_mut(j).get_unchecked_mut(info.T) = info.c.get_unchecked(*solut.s.get_unchecked(j_prev)).get_unchecked(*solut.s.get_unchecked(j)) + solut.seq.get_unchecked(i).get_unchecked(j_prev).get_unchecked(info.T);
            *solut.seq.get_unchecked_mut(i).get_unchecked_mut(j).get_unchecked_mut(info.C) = solut.seq.get_unchecked(i).get_unchecked(j).get_unchecked(info.T) + solut.seq.get_unchecked(i).get_unchecked(j_prev).get_unchecked(info.C);
            *solut.seq.get_unchecked_mut(i).get_unchecked_mut(j).get_unchecked_mut(info.W) = (j as i32 + k) as f64;
            }
            */

          //solut.seq[i][j][info.T] = info.c[solut.s[j_prev]][solut.s[j]] + solut.seq[i][j_prev][info.T];
          //solut.seq[i][j][info.C] = solut.seq[i][j][info.T] + solut.seq[i][j_prev][info.C];
          //solut.seq[i][j][info.W] = (j as i32 + k) as f64;
        }
    }

    //solut.cost = solut.seq[0][info.dimen][info.C];
    unsafe{solut.cost = solut.seq.get_unchecked(0).get_unchecked( info.dimen).C;}
}


fn sort(arr : &mut Vec<usize>, r : usize, info : & tInfo) {
    for i in 0..arr.len() {
        for j in 0..(arr.len()-1) {
            if info.c[r][arr[j]] > info.c[r][arr[j+1]] {
                let tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;
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
        /*c_list.sort_by(|i, j| 
                    if info.c[*i as usize][r] < info.c[*j as usize][r] 
                        {Ordering::Less} 
                    else 
                        {Ordering::Greater});
        */

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

    unsafe {
    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        let i_next : usize = i + 1;

        cost_concat_1 = solut.seq.get_unchecked(0).get_unchecked( i_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(i_next));
        cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(i).get_unchecked( i_next).T + info.c.get_unchecked(*solut.s.get_unchecked(i)).get_unchecked(*solut.s.get_unchecked(i_next+1));

        cost_new = solut.seq.get_unchecked(0).get_unchecked( i_prev).C
            + solut.seq.get_unchecked(i).get_unchecked( i_next).W     * cost_concat_1 + info.c.get_unchecked(*solut.s.get_unchecked(i_next)).get_unchecked(*solut.s.get_unchecked(i))
            + solut.seq.get_unchecked(i_next+1).get_unchecked( info.dimen).W * cost_concat_2 + solut.seq.get_unchecked(i_next+1).get_unchecked(info.dimen ).C;


      //cost_concat_1 = solut.seq[0][ i_prev][ info.T] + info.c[solut.s[i_prev]][solut.s[i_next]];
      //cost_concat_2 = cost_concat_1 + solut.seq[i][ i_next][ info.T] + info.c[solut.s[i]][solut.s[i_next+1]];

      //cost_new = solut.seq[0][ i_prev][ info.C]
      //    + solut.seq[i][ i_next][ info.W]     * cost_concat_1 + info.c[solut.s[i_next]][solut.s[i]]
      //    + solut.seq[i_next+1][ info.dimen][ info.W] * cost_concat_2 + solut.seq[i_next+1][ info.dimen ][info.C];

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next+1..info.dimen {
            let j_next = j+1;
            let j_prev = j-1;

            cost_concat_1 = solut.seq.get_unchecked(0).get_unchecked( i_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j));
            cost_concat_2 = cost_concat_1 + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(i_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(i_next).get_unchecked( j_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(j_prev)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_4 = cost_concat_3  + info.c.get_unchecked(*solut.s.get_unchecked(i)).get_unchecked(*solut.s.get_unchecked(j_next));

            cost_new = solut.seq.get_unchecked(0).get_unchecked( i_prev).C
                    + cost_concat_1
                    + solut.seq.get_unchecked(i_next).get_unchecked( j_prev).W * cost_concat_2 + solut.seq.get_unchecked(i_next).get_unchecked( j_prev).C
                    + cost_concat_3
                    + solut.seq.get_unchecked(j_next).get_unchecked( info.dimen).W * cost_concat_4 + solut.seq.get_unchecked(j_next).get_unchecked( info.dimen).C;


          //cost_concat_1 = solut.seq[0][ i_prev][ info.T] + info.c[solut.s[i_prev]][solut.s[j]];
          //cost_concat_2 = cost_concat_1 + info.c[solut.s[j]][solut.s[i_next]];
          //cost_concat_3 = cost_concat_2 + solut.seq[i_next][ j_prev][ info.T] + info.c[solut.s[j_prev]][solut.s[i]];
          //cost_concat_4 = cost_concat_3  + info.c[solut.s[i]][solut.s[j_next]];

          ////cost_new = *solut.seq.get_unchecked(0).get_unchecked(i_prev).get_unchecked( info.C)
          //cost_new = solut.seq[0][ i_prev][ info.C]
          //        + cost_concat_1
          //        + solut.seq[i_next][ j_prev][ info.W] * cost_concat_2 + solut.seq[i_next][ j_prev][ info.C]
          //        + cost_concat_3
          //        + solut.seq[j_next][ info.dimen][ info.W] * cost_concat_4 + solut.seq[j_next][ info.dimen][ info.C];


            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }
    }

    /*
    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        let i_next : usize = i + 1;

      //cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, &info)] + info.c[solut.s[i_prev]][solut.s[i_next]];
      //cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, i_next, info.T, &info)] + info.c[solut.s[i]][solut.s[i_next+1]];

      //cost_new = solut.seq[to_1D(0, i_prev, info.C, &info)]
      //    + solut.seq[to_1D(i, i_next, info.W, &info)]     * cost_concat_1 + info.c[solut.s[i_next]][solut.s[i]]
      //    + solut.seq[to_1D(i_next+1, info.dimen, info.W, &info)] * cost_concat_2 + solut.seq[to_1D(i_next+1, info.dimen ,info.C, &info)];

        unsafe {
        cost_concat_1 = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(i_next));
        cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(to_1D(i, i_next, info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(i)).get_unchecked(*solut.s.get_unchecked(i_next+1));

        cost_new = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).C
            + solut.seq.get_unchecked(to_1D(i, i_next,  info.dimen)).W     * cost_concat_1 + info.c.get_unchecked(*solut.s.get_unchecked(i_next)).get_unchecked(*solut.s.get_unchecked(i))
            + solut.seq.get_unchecked(to_1D(i_next+1, info.dimen,  info.dimen)).W * cost_concat_2 + solut.seq.get_unchecked(to_1D(i_next+1, info.dimen , info.dimen)).C;
        }

        if cost_new < cost_best {
            cost_best = cost_new - f64::EPSILON;
            I = i;
            J = i_next;
        }

        for j in i_next+1..info.dimen {
            let j_next = j+1;
            let j_prev = j-1;

          //cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, &info)] + info.c[solut.s[i_prev]][solut.s[j]];
          //cost_concat_2 = cost_concat_1 + info.c[solut.s[j]][solut.s[i_next]];
          //cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i_next, j_prev, info.T, &info)] + info.c[solut.s[j_prev]][solut.s[i]];
          //cost_concat_4 = cost_concat_3  + info.c[solut.s[i]][solut.s[j_next]];

          //cost_new = solut.seq[to_1D(0, i_prev, info.C, &info)]
          //        + cost_concat_1
          //        + solut.seq[to_1D(i_next, j_prev, info.W, &info)] * cost_concat_2 + solut.seq[to_1D(i_next, j_prev, info.C, &info)]
          //        + cost_concat_3
          //        + solut.seq[to_1D(j_next, info.dimen, info.W, &info)] * cost_concat_4 + solut.seq[to_1D(j_next, info.dimen, info.C, &info)];

            unsafe {
            cost_concat_1 = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j)) ;
            cost_concat_2 = cost_concat_1 + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(i_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(to_1D(i_next, j_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(j_prev)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_4 = cost_concat_3  + info.c.get_unchecked(*solut.s.get_unchecked(i)).get_unchecked(*solut.s.get_unchecked(j_next));

            cost_new = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).C
                    + cost_concat_1
                    + solut.seq.get_unchecked(to_1D(i_next, j_prev,  info.dimen)).W * cost_concat_2 + solut.seq.get_unchecked(to_1D(i_next, j_prev,  info.dimen)).C
                    + cost_concat_3
                    + solut.seq.get_unchecked(to_1D(j_next, info.dimen,  info.dimen)).W * cost_concat_4 + solut.seq.get_unchecked(to_1D(j_next, info.dimen,  info.dimen)).C;
            }

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }
    */

    if cost_best < solut.cost - f64::EPSILON {
        //println!("swap \n{}", cost_best);
        swap(&mut solut.s, I, J);
        subseq_load(solut, info);
        //subseq_load(s, info);
        //println!("{}", seq[0][info.dimension][C]);
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

    unsafe {
    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        //let mut rev_seq_cost : f64 = solut.seq[i][ i+1][ info.T];
        let mut rev_seq_cost : f64 = solut.seq.get_unchecked(i).get_unchecked( i+1).T;

        for j in i+2..info.dimen {
            let j_next = j + 1;

            rev_seq_cost += info.c.get_unchecked(*solut.s.get_unchecked(j-1)).get_unchecked(*solut.s.get_unchecked(j)) * (solut.seq.get_unchecked(i).get_unchecked( j).W-1.0);

            cost_concat_1 =  solut.seq.get_unchecked(0).get_unchecked( i_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(i_prev));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(i).get_unchecked( j).T + info.c.get_unchecked(*solut.s.get_unchecked(j_next)).get_unchecked(*solut.s.get_unchecked(i));

            cost_new = solut.seq.get_unchecked(0).get_unchecked( i_prev).C
                    + solut.seq.get_unchecked(i).get_unchecked( j).W      * cost_concat_1 + rev_seq_cost
                    + solut.seq.get_unchecked(j_next).get_unchecked( info.dimen).W * cost_concat_2 + solut.seq.get_unchecked(j_next).get_unchecked( info.dimen).C;

          //rev_seq_cost += info.c[solut.s[j-1]][solut.s[j]] * (solut.seq[i][ j][ info.W]-1.0);

          //cost_concat_1 =  solut.seq[0][ i_prev][ info.T] + info.c[solut.s[j]][solut.s[i_prev]];
          //cost_concat_2 = cost_concat_1 + solut.seq[i][ j][ info.T] + info.c[solut.s[j_next]][solut.s[i]];

          //cost_new = solut.seq[0][ i_prev][ info.C]
          //        + solut.seq[i][ j][ info.W]      * cost_concat_1 + rev_seq_cost
          //        + solut.seq[j_next][ info.dimen][ info.W] * cost_concat_2 + solut.seq[j_next][ info.dimen][ info.C];

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
    }
    }

    /*
    for i in 1..info.dimen-1 {
        let i_prev : usize = i - 1;
        unsafe { 
        let mut rev_seq_cost : f64 = solut.seq.get_unchecked(to_1D(i, i+1,  info.dimen)).T;

        for j in i+2..info.dimen {
            let j_next = j + 1;

          //rev_seq_cost += info.c[solut.s[j-1]][solut.s[j]] * (solut.seq[to_1D(i, j, info.W, &info)]-1.0);

          //cost_concat_1 =  solut.seq[to_1D(0, i_prev, info.T, &info)] + info.c[solut.s[j]][solut.s[i_prev]];
          //cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, info.T, &info)] + info.c[solut.s[j_next]][solut.s[i]];

          //cost_new = solut.seq[to_1D(0, i_prev, info.C, &info)]
          //        + solut.seq[to_1D(i, j, info.W, &info)]      * cost_concat_1 + rev_seq_cost
          //        + solut.seq[to_1D(j_next, info.dimen, info.W, &info)] * cost_concat_2 + solut.seq[to_1D(j_next, info.dimen, info.C, &info)];

            rev_seq_cost += info.c.get_unchecked(*solut.s.get_unchecked(j-1)).get_unchecked(*solut.s.get_unchecked(j)) * (solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).W-1.0);

            cost_concat_1 =  solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(i_prev));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(j_next)).get_unchecked(*solut.s.get_unchecked(i));

            cost_new = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).C
                    + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).W      * cost_concat_1 + rev_seq_cost
                    + solut.seq.get_unchecked(to_1D(j_next, info.dimen,  info.dimen)).W * cost_concat_2 + solut.seq.get_unchecked(to_1D(j_next, info.dimen, info.dimen)).C;

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
            }
        }
        }
    }
    */


    if cost_best < solut.cost - f64::EPSILON {
        //println!("two_opt \n{}", cost_best);
        reverse(&mut solut.s, I, J);
        subseq_load(solut, info);
        //subseq_load(s, seq, info);
        //println!("{}", seq[0][info.dimension][C]);
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

    unsafe {
    for i in 1..info.dimen-opt+1 {
        let j : usize = opt+i-1;
        let i_prev : usize = i-1;
        let j_next : usize = j+1;

        let seq_i_j_W = solut.seq.get_unchecked(i).get_unchecked(j).W;
        let seq_i_j_C = solut.seq.get_unchecked(i).get_unchecked(j).C;
        let seq_jnext_dimen_W = solut.seq.get_unchecked(j_next).get_unchecked(info.dimen).W;
        let seq_jnext_dimen_C = solut.seq.get_unchecked(j_next).get_unchecked(info.dimen).C;

        for k in 0..i_prev {
            let k_next : usize = k+1;

            cost_concat_1 = solut.seq.get_unchecked(0).get_unchecked(k).T + info.c.get_unchecked(*solut.s.get_unchecked(k)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(i).get_unchecked(j).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(k_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(k_next).get_unchecked(i_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j_next));

          //cost_concat_1 = solut.seq[0][k][info.T] + info.c[solut.s[k]][solut.s[i]];
          //cost_concat_2 = cost_concat_1 + solut.seq[i][j][info.T] + info.c[solut.s[j]][solut.s[k_next]];
          //cost_concat_3 = cost_concat_2 + solut.seq[k_next][i_prev][info.T] + info.c[solut.s[i_prev]][solut.s[j_next]];







            //cost_new = solut.seq[0][k][info.C]                                                             /*        1st subseq */
            //  + solut.seq[i][j][info.W]              * cost_concat_1 + solut.seq[i][j][info.C]                  /* concat 2nd subseq (reinserted seq) */
            //  + solut.seq[k_next][i_prev][info.W]   * cost_concat_2 + solut.seq[k_next][ i_prev][ info.C]        /* concat 3rd subseq */
            //  + solut.seq[j_next][ info.dimen][ info.W] * cost_concat_3 + solut.seq[j_next][ info.dimen][ info.C];    /* concat 4th subseq */

            cost_new = solut.seq.get_unchecked(0).get_unchecked(k).C                                                             /*        1st subseq */
                + seq_i_j_W              * cost_concat_1 + seq_i_j_C                  /* concat 2nd subseq (reinserted seq) */
                + solut.seq.get_unchecked(k_next).get_unchecked(i_prev).W   * cost_concat_2 + solut.seq.get_unchecked(k_next).get_unchecked(i_prev).C        /* concat 3rd subseq */
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

            cost_concat_1 = solut.seq.get_unchecked(0).get_unchecked(i_prev).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j_next));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(j_next).get_unchecked(k).T + info.c.get_unchecked(*solut.s.get_unchecked(k)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(i).get_unchecked(j).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(k_next));

            cost_new = solut.seq.get_unchecked(0).get_unchecked( i_prev).C                                                        /*      1st subseq */
                + solut.seq.get_unchecked(j_next).get_unchecked( k).W         * cost_concat_1 + solut.seq.get_unchecked(j_next).get_unchecked( k).C             /* concat 2nd subseq */
                + solut.seq.get_unchecked(i).get_unchecked( j).W              * cost_concat_2 + solut.seq.get_unchecked(i).get_unchecked( j).C                  /* concat 3rd subseq (reinserted seq) */
                + solut.seq.get_unchecked(k_next).get_unchecked( info.dimen).W * cost_concat_3 + solut.seq.get_unchecked(k_next).get_unchecked( info.dimen).C;    /* concat 4th subseq */

          //cost_concat_1 = solut.seq[0][ i_prev][ info.T] + info.c[solut.s[i_prev]][solut.s[j_next]];
          //cost_concat_2 = cost_concat_1 + solut.seq[j_next][ k][ info.T] + info.c[solut.s[k]][solut.s[i]];
          //cost_concat_3 = cost_concat_2 + solut.seq[i][ j][ info.T] + info.c[solut.s[j]][solut.s[k_next]];

          //cost_new = solut.seq[0][ i_prev][ info.C]                                                        /*      1st subseq */
          //    + solut.seq[j_next][ k][ info.W]         * cost_concat_1 + solut.seq[j_next][ k][ info.C]             /* concat 2nd subseq */
          //    + solut.seq[i][ j][ info.W]              * cost_concat_2 + solut.seq[i][ j][ info.C]                  /* concat 3rd subseq (reinserted seq) */
          //    + solut.seq[k_next][ info.dimen][ info.W] * cost_concat_3 + solut.seq[k_next][ info.dimen][ info.C];    /* concat 4th subseq */

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }
    }

        /*
    for i in 1..info.dimen-opt+1 {
        let j : usize = opt+i-1;
        let i_prev : usize = i-1;
        let j_next : usize = j+1;

        for k in 0..i_prev {
            let k_next : usize = k+1;

          //cost_concat_1 = solut.seq[to_1D(0, k, info.T, &info)] + info.c[solut.s[k]][solut.s[i]];
          //cost_concat_2 = cost_concat_1 + solut.seq[to_1D(i, j, info.T, &info)] + info.c[solut.s[j]][solut.s[k_next]];
          //cost_concat_3 = cost_concat_2 + solut.seq[to_1D(k_next, i_prev, info.T, &info)] + info.c[solut.s[i_prev]][solut.s[j_next]];

          //cost_new = solut.seq[to_1D(0, k, info.C, &info)]                                                             /*        1st subseq */
          //    + solut.seq[to_1D(i, j, info.W, &info)]              * cost_concat_1 + solut.seq[to_1D(i, j, info.C, &info)]                  /* concat 2nd subseq (reinserted seq) */
          //    + solut.seq[to_1D(k_next, i_prev, info.W, &info)]    * cost_concat_2 + solut.seq[to_1D(k_next, i_prev, info.C, &info)]        /* concat 3rd subseq */
          //    + solut.seq[to_1D(j_next, info.dimen, info.W, &info)] * cost_concat_3 + solut.seq[to_1D(j_next, info.dimen, info.C, &info)];    /* concat 4th subseq */

            unsafe {
            cost_concat_1 = solut.seq.get_unchecked(to_1D(0, k,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(k)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(k_next));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(to_1D(k_next, i_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j_next));

            cost_new = solut.seq.get_unchecked(to_1D(0, k,  info.dimen)).C                                                             /*        1st subseq */
                + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).W              * cost_concat_1 + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).C                  /* concat 2nd subseq (reinserted seq) */
                + solut.seq.get_unchecked(to_1D(k_next, i_prev,  info.dimen)).W    * cost_concat_2 + solut.seq.get_unchecked(to_1D(k_next, i_prev,  info.dimen)).C
                + solut.seq.get_unchecked(to_1D(j_next, info.dimen, info.dimen)).W * cost_concat_3 + solut.seq.get_unchecked(to_1D(j_next, info.dimen, info.dimen)).C;    /* concat 4th subseq */
            }

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }

        for k in i+opt..info.dimen-opt-1 {
            let k_next : usize = k+1;

          //cost_concat_1 = solut.seq[to_1D(0, i_prev, info.T, &info)] + info.c[solut.s[i_prev]][solut.s[j_next]];
          //cost_concat_2 = cost_concat_1 + solut.seq[to_1D(j_next, k, info.T, &info)] + info.c[solut.s[k]][solut.s[i]];
          //cost_concat_3 = cost_concat_2 + solut.seq[to_1D(i, j, info.T, &info)] + info.c[solut.s[j]][solut.s[k_next]];

          //cost_new = solut.seq[to_1D(0, i_prev, info.C, &info)]                                                        /*      1st subseq */
          //    + solut.seq[to_1D(j_next, k, info.W, &info)]         * cost_concat_1 + solut.seq[to_1D(j_next, k, info.C, &info)]             /* concat 2nd subseq */
          //    + solut.seq[to_1D(i, j, info.W, &info)]              * cost_concat_2 + solut.seq[to_1D(i, j, info.C, &info)]                  /* concat 3rd subseq (reinserted seq) */
          //    + solut.seq[to_1D(k_next, info.dimen, info.W, &info)] * cost_concat_3 + solut.seq[to_1D(k_next, info.dimen, info.C, &info)];    /* concat 4th subseq */
            unsafe {

            cost_concat_1 = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(i_prev)).get_unchecked(*solut.s.get_unchecked(j_next));
            cost_concat_2 = cost_concat_1 + solut.seq.get_unchecked(to_1D(j_next, k,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(k)).get_unchecked(*solut.s.get_unchecked(i));
            cost_concat_3 = cost_concat_2 + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).T + info.c.get_unchecked(*solut.s.get_unchecked(j)).get_unchecked(*solut.s.get_unchecked(k_next));

            cost_new = solut.seq.get_unchecked(to_1D(0, i_prev,  info.dimen)).C                                                        /*      1st subseq */
                + solut.seq.get_unchecked(to_1D(j_next, k,  info.dimen)).W         * cost_concat_1 + solut.seq.get_unchecked(to_1D(j_next, k,  info.dimen)).C             /* concat 2nd subseq */
                + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).W              * cost_concat_2 + solut.seq.get_unchecked(to_1D(i, j,  info.dimen)).C                  /* concat 3rd subseq (reinserted seq) */
                + solut.seq.get_unchecked(to_1D(k_next, info.dimen,  info.dimen)).W * cost_concat_3 + solut.seq.get_unchecked(to_1D(k_next, info.dimen,  info.dimen)).C;    /* concat 4th subseq */
            }

            if cost_new < cost_best {
                cost_best = cost_new - f64::EPSILON;
                I = i;
                J = j;
                POS = k;
            }
        }
    }
*/

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
        //println!("is empty {} \n", neighbd_list.is_empty());
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

        //println!("{:?} \n{}\n", s, seq[0][info.dimension][C]);
        //println!("{:?} \n", s);
        //exit(1);
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


        //println!("perturbaa");
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
        //seq : vec![tSeqInfo {C:0.0, W:0.0, T:0.0} ; (info.dimen+1)*(info.dimen+1)],
        //seq : vec![0.0; (info.dimen+1)*(info.dimen+1)*3].into_boxed_slice(),
        seq : Box::new([[tSeqInfo {C:0.0, W:0.0, T:0.0}; sz::SIZE+1]; sz::SIZE+1]),
        //seq : vec![vec![tSeqInfo {C:0.0, W:0.0, T:0.0}; info.dimen+1]; info.dimen+1],
        s : vec![0; info.dimen],
        cost : f64::MAX,
    };

    let mut solut_partial = tSolution {
        //seq : vec![tSeqInfo {C:0.0, W:0.0, T:0.0}; (info.dimen+1)*(info.dimen+1)],
        //seq : vec![0.0; (info.dimen+1)*(info.dimen+1)*3].into_boxed_slice(),
        seq : Box::new([[tSeqInfo {C:0.0, W:0.0, T:0.0}; sz::SIZE+1]; sz::SIZE+1]),
        //seq : vec![vec![tSeqInfo {C:0.0, W:0.0, T:0.0}; info.dimen+1]; info.dimen+1],
        s : vec![0; info.dimen],
        cost : 0.0,
    };

    let mut solut_crnt = tSolution {
        //seq : vec![tSeqInfo {C:0.0, W:0.0, T:0.0}; (info.dimen+1)*(info.dimen+1)],
        //seq : vec![0.0; (info.dimen+1)*(info.dimen+1)*3].into_boxed_slice(),
        seq : Box::new([[tSeqInfo {C:0.0, W:0.0, T:0.0}; sz::SIZE+1]; sz::SIZE+1]),
        //seq : vec![vec![tSeqInfo {C:0.0, W:0.0, T:0.0}; info.dimen+1]; info.dimen+1],
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
        //exit(0);
        subseq_load(&mut solut_crnt, info);
        println!("\t[+] Constructing Inital Solution.. {}", solut_crnt.cost);
        solut_partial.s = solut_crnt.s.clone();
        solut_partial.cost = solut_crnt.cost;
        //let mut s : Vec<usize> = construction(alpha, info);
        //let mut sl = s.clone();

        //subseq_load(&s, &mut subseq, info);
        //println!("Construction cost {} \n", rvnd_cost_best);

        println!("\t[+] Looking for the best Neighbor..");
        let mut iterILS : usize = 0;
        while iterILS < Iils {
            RVND(&mut solut_crnt, info);
            if solut_crnt.cost < solut_partial.cost {
                solut_partial.s = solut_crnt.s.clone();
                solut_partial.cost = solut_crnt.cost;
                iterILS = 0;

                //println!("{}  {:?}", solut_partial.cost, solut_partial.s);
            }

            solut_crnt.s = perturb(&solut_partial.s, info);
            subseq_load(&mut solut_crnt, info);
            iterILS += 1;
        }

        //subseq_load(&sl, &mut subseq, info);

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
    let mut dimension : usize = 0;
    let mut c = Box::new([[0.0; sz::SIZE]; sz::SIZE]);
    //let mut c : [[f64; 500]; 500] = [[0.0; 500]; 500];
    //let mut c : Vec<Vec<f64>> = vec![vec![];0];

    let mut rnd : Vec<usize> = vec![];

    data::load(&mut dimension, &mut c, &mut rnd);

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

    subseq_load(&mut solut_test, &info);


    println!(" {:?}\nCost {}", solut_test.s, solut_test.cost);
    //exit(0);
    */


    println!("TEST");

    let Imax = 10;
    let Iils = if dimension < 100 {dimension} else {100};

    let R = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25];

    let now = Instant::now();


    //test!();
    GILS_RVND(Imax, Iils, R, &mut info);

    let new_now = Instant::now();
    println!("TIME: {}", new_now.duration_since(now).as_secs_f64());

    /*
    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", info.c[i][j]);
        }
        println!();
    }
    */

}
