mod Data;
extern crate opa;
extern crate rand;
use rand::Rng;
use std::cmp::Ordering;

static T : usize = 0;
static C : usize = 1;
static W : usize = 2;

#[derive(Debug)]
struct Instance_info {
    dimension : usize,
    c : Vec<Vec<f64>>,
}

fn subseq_load(s : Vec<usize>, seq : &mut Vec<Vec<Vec<f64>>>, info : &mut Instance_info) {
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

fn construction(alpha : f64, info : Instance_info) -> Vec<usize> {
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
        let index = rng.gen::<usize>() % i;
        let c = cList[index];
        cList.remove(index);
        s.push(c);
    }
    println!("{:?}", s);

    return s;
}

fn main() {
    let mut dimension : usize = 0;
    let mut c : Vec<Vec<f64>> = vec![vec![];0];

    Data::load(&mut dimension, &mut c);

    let mut info = Instance_info {
        dimension : dimension,
        c : c,
    };

    //println!("{:#?}\n", info);
    let mut s : Vec<usize> = vec![];
    for i in 0..dimension {
        s.push(i);
    }
    s.push(0);

    let mut subseq : Vec<Vec<Vec<f64>>> = 
        vec![
            vec![
                vec![0.0; 3]; 
            dimension+1]; 
        dimension+1];

    subseq_load(s, &mut subseq, &mut info);
    println!("{} \n", subseq[0][info.dimension][C]);

    construction(0.2, info);
    /*
    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", info.c[i][j]);
        }
        println!();
    }
    */

}
