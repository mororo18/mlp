use std::env;
use std::fs;

fn print_type_of<T>(_: &T) {
        println!("{}", std::any::type_name::<T>())
}

pub fn opa() {
    println!("opa !!");
    let filename = "../distance_matrix";
    println!("In file {}", filename);

    let contents = fs::read_to_string(filename)
        .expect("Something went wrong reading the file");

    //println!("With text:\n{}", contents);
    print_type_of(&contents);
    for l in contents.lines() {
        println!("{}", l);
    }

}
