use crate::data::sz;
use std::ops::Index;


#[derive(Debug, Copy, Clone)]
pub struct tSeqInfo {
   pub  T : f64,
   pub  C : f64,
   pub  W : f64,
}

pub type tSeqData = Box<[[tSeqInfo; sz::SIZE+1]; sz::SIZE+1]>;
pub type tCostData = Box<[[f64; sz::SIZE]; sz::SIZE]>;


#[derive(Debug, Clone)]
pub struct tInfo {
    pub c : tCostData,
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
