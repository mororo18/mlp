#![allow(dead_code)]

use std::fs::File;
use std::io::{BufRead, BufReader};

pub fn load_data(rnd: &mut Vec<i32>) -> (Vec<Vec<f64>>, i32) {
    let file = File::open("../distance_matrix").expect("could not open ../distance_matrix");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let dimension: i32 = lines
        .next()
        .expect("empty distance_matrix file")
        .unwrap()
        .parse()
        .unwrap();

    let mut matrix = vec![vec![0.0f64; dimension as usize]; dimension as usize];

    let mut flag = false;
    let mut rnd_size: i32 = -1;
    let mut i: usize = 0;

    for line in lines {
        let mut line = line.unwrap();
        let mut j = i + 1;

        while let Some(pos) = line.find(' ') {
            let value: f64 = line[..pos].parse().unwrap();
            line = line[pos + 1..].to_string();
            matrix[i][j] = value;
            matrix[j][i] = value;
            j += 1;
        }

        i += 1;

        if flag {
            if rnd_size < 0 {
                rnd_size = line.parse().unwrap();
                rnd.reserve(rnd_size as usize);
            } else {
                rnd.push(line.parse().unwrap());
            }
        }

        if line == "RND" {
            flag = true;
        }
    }

    (matrix, dimension)
}

pub struct Data {
    dimen: i32,
    cost: Vec<Vec<f64>>,
    rnd: Vec<i32>,
}

impl Data {
    pub fn new() -> Self {
        Data {
            dimen: 0,
            cost: Vec::new(),
            rnd: Vec::new(),
        }
    }

    pub fn load(&mut self) {
        let mut rnd = Vec::new();
        let (cost, dimen) = load_data(&mut rnd);
        self.cost = cost;
        self.dimen = dimen;
        self.rnd = rnd;
    }

    pub fn get_rnd_vec(&self) -> Vec<i32> {
        self.rnd.clone()
    }

    pub fn get_dimen(&self) -> i32 {
        self.dimen
    }

    pub fn get_cost_ptr(&self) -> &Vec<Vec<f64>> {
        &self.cost
    }
}
