use std::fs;
use std::process::exit;
use std::fmt::format;
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;

fn main() {
    let inst_file = "../distance_matrix";

    let content = fs::read_to_string(inst_file).expect("qebro");
    let size = content.lines().nth(0).unwrap().parse::<usize>().unwrap();

    let code = format!("pub const SIZE : usize = {};\n", size); 

    let mut file = File::create("src/data/sz.rs").expect("qebro tmb");
    file.write(code.as_bytes());

}
