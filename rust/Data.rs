use std::process::exit;
use std::fs;
use std::vec;

fn print_type_of<T>(_: &T) {
        println!("{}", std::any::type_name::<T>())
}

pub fn load(dimension : &mut usize, c : &mut Vec<Vec<f64>>) {

    let filename = "../distance_matrix";
    //println!("In file {}", filename);

    let contents = fs::read_to_string(filename)
        .expect("Something went wrong reading the file");

    let mut l_i = -1;
    print_type_of(&contents);
    for l in contents.lines() {
        let mut iter = l.split_ascii_whitespace();
        let mut n = iter.next();

        if l == "EOF" {
            break;
        }

        print_type_of(&n.unwrap().to_string().parse::<u32>().unwrap());

        if l_i == -1 {
            *dimension = n.unwrap().to_string().parse::<usize>().unwrap();
            for k in 0..*dimension {c.push(vec![0.0; *dimension]);}
            /*
            let mut k : usize = (*dimension as usize) - 1; 
            c[k][k] = 999f64;
            println!("{:#?} \n", c[(*dimension as usize) - 1][k]);
            return;
            */
            l_i += 1;
            continue;
        }

        let mut i = l_i as usize;
        let mut j = (l_i as usize) + 1;

        while j < *dimension as usize{
            c[i][j] = n.unwrap().to_string().parse::<f64>().unwrap();
            c[j][i] = c[i][j];
            //print!("{} ", n.unwrap());
            n = iter.next();
            j += 1;
        }
        //print!("\n");
        l_i += 1;
    }

}
