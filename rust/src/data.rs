use std::process::exit;
use std::fs;
use std::fs::File;
use std::io::prelude::*;


use crate::subseq::CostMatrix;

pub fn load(filename: &str) -> (usize, CostMatrix, Vec<usize>) {

    let contents = fs::read_to_string(filename)
        .expect("Something went wrong reading the file");

    let mut line_i = 0;

    let mut lines : Vec<String> = vec![];
    for l in contents.lines() {
        lines.push(l.to_string());
    }

    let dimension = lines[line_i].parse::<usize>().unwrap();
    line_i += 1;

    let c = CostMatrix::new(dimension);

    for l_i in 0..dimension {
        let mut iter = lines[line_i].split_ascii_whitespace();
        line_i += 1;
        let mut n = iter.next();

        let mut i = l_i;
        let mut j = l_i + 1;

        while j < dimension {
            c[i][j] = n.unwrap().to_string().parse::<f64>().unwrap();
            c[j][i] = c[i][j];
            n = iter.next();
            j += 1;
        }
        //print!("\n");
        //l_i += 1;
        
    }

    line_i += 2;
    let rnd_size = lines[line_i].parse::<usize>().unwrap();
    line_i += 1;

    let rnd: Vec<usize> = vec![];
    for i in 0..rnd_size {
        rnd.push(lines[line_i].parse::<usize>().unwrap());
        line_i += 1;
    }

    (dimension, c, rnd)
}
