mod Data;

fn main() {
    let mut dimension : usize = 0;
    let mut c : Vec<Vec<f64>> = vec![vec![];0];

    Data::load(&mut dimension, &mut c);
    println!("{} \n", dimension);

    for i in 0..dimension {
        for j in 0..dimension {
            print!("{} ", c[i][j]);
        }
        println!();
    }

}
