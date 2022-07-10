use crate::data::sz;

#[derive(Debug, Clone)]
pub struct tInfo {
    pub c : Box<[[f64; sz::SIZE]; sz::SIZE]>,
    //c : Vec<Vec<f64>>,
    pub dimen : usize,
  //T : usize,
  //C : usize,
  //W : usize,
    pub SWAP         : usize ,
    pub REINSERTION  : usize ,
    pub OR_OPT_2     : usize ,
    pub OR_OPT_3     : usize ,
    pub TWO_OPT      : usize ,
    pub rnd : Vec<usize>,
    pub  rnd_index : usize,
}

#[derive(Debug, Copy, Clone)]
pub struct tSeqInfo {
   pub  T : f64,
   pub  C : f64,
   pub  W : f64,
}


pub type tSeqData = Box<[[tSeqInfo; sz::SIZE+1]; sz::SIZE+1]>;

/*
trait Kba {
    fn new(size : usize) -> Self ;
}

impl Kba for tSeqData {
    fn new(size : usize) -> tSeqData {
        Box::new([[tSeqInfo {C:0.0, W:0.0, T:0.0}; size]; size])
    }
}
*/

#[derive(Debug, Clone)]
pub struct tSolution {
    //seq : Box<[f64]>,
    //seq : Vec<f64>,
    //seq : Vec<tSeqInfo>,
    //seq : Vec<Vec<tSeqInfo>>,
    pub seq : tSeqData,
    pub s : Vec<usize>,
    pub cost : f64,
}

pub trait Test<T> {
    fn get(&self, i : usize) -> T;
}

impl<T> Test<T> for Vec<T> 
where T: Copy, {
    
    fn get(&self, i : usize) -> T {
        unsafe {*self.get_unchecked(i)}
    }
}

pub trait Access {
    fn set_C(&mut self, i : usize, j : usize, value : f64);
    fn set_T(&mut self, i : usize, j : usize, value : f64);
    fn set_W(&mut self, i : usize, j : usize, value : f64);

    fn get_C(&self, i : usize, j : usize) -> f64;
    fn get_T(&self, i : usize, j : usize) -> f64;
    fn get_W(&self, i : usize, j : usize) -> f64;
}

impl Access for Box<[[tSeqInfo; sz::SIZE+1]; sz::SIZE+1]> {

    fn set_C(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).C = value;}
    }
     
    fn set_T(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).T = value;}
    }

    fn set_W(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).W = value;}
    }




    fn get_C(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).C}
    }

    fn get_T(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).T}
    }

    fn get_W(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).W}
    }
}

