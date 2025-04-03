use std::process::exit;
use std::fs;
use std::vec;
use std::fs::File;
use std::io::prelude::*;

pub fn load(dimension : &mut usize, c : &mut Vec<Vec<f64>>, rnd : &mut Vec<usize>) {

    let filename = "../distance_matrix";

    let contents = fs::read_to_string(filename)
        .expect("Something went wrong reading the file");

    let mut line_i = 0;

    let mut lines : Vec<String> = vec![];
    for l in contents.lines() {
        lines.push(l.to_string());
    }

    *dimension = lines[line_i].parse::<usize>().unwrap();
    for k in 0..*dimension {c.push(vec![0.0; *dimension]);}
    line_i += 1;

    for l_i in 0..*dimension {
        let mut iter = lines[line_i].split_ascii_whitespace();
        line_i += 1;
        let mut n = iter.next();

        let mut i = l_i;
        let mut j = l_i + 1;

        while j < *dimension {
            c[i][j] = n.unwrap().to_string().parse::<f64>().unwrap();
            c[j][i] = c[i][j];
            n = iter.next();
            j += 1;
        }
        
    }

    line_i += 2;
    let rnd_size = lines[line_i].parse::<usize>().unwrap();
    line_i += 1;

    for i in 0..rnd_size {
        rnd.push(lines[line_i].parse::<usize>().unwrap());
        line_i += 1;
    }

}
