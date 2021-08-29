mod Data;

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

fn main() {
    let mut dimension : usize = 0;
    let mut c : Vec<Vec<f64>> = vec![vec![];0];

    Data::load(&mut dimension, &mut c);
    //println!("{} \n", dimension);

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

    /*
    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", info.c[i][j]);
        }
        println!();
    }
    */

}
