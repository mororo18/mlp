use std::ops::Index;

pub const SWAP        : usize = 0;
pub const REINSERTION : usize = 1;
pub const OR_OPT_2    : usize = 2;
pub const OR_OPT_3    : usize = 3;
pub const TWO_OPT     : usize = 4;

#[derive(Debug, Copy, Clone)]
pub struct SubseqInfo {
   pub  T : f64,
   pub  C : f64,
   pub  W : f64,
}

impl SubseqInfo {
    fn zeored () -> Self {
        Self { T: 0.0, C: 0.0, W: 0.0 }
    }
}

struct SubseqMatrix {
    row_size: usize,
#[cfg(feature="flat")]
    mem: Vec<SubseqInfo>,
#[cfg(not(feature="flat"))]
    mem:  Vec<Vec<SubseqInfo>>,
}

impl SubseqMatrix {
    fn new (n: usize) -> Self {
        Self { 
            size: n+1,
#[cfg(feature="flat")]
            mem: vec![SubseqInfo::zeored(); (n+1) * (n+1)],
#[cfg(not(feature="flat"))]
            mem: vec![ vec![ SubseqInfo::zeored(); n+1]; n+1],
        }
    }
}

impl Index<(usize, usize)> for SubseqMatrix {
    type Output = SubseqInfo;

    fn index(&self,  index: (usize, usize)) -> &Self::Output {
        let (row, col) = row_index;
#[cfg(feature="flat")]
        unsafe {
            self.mem.get_unchecked(self.row_size * row + col)
        }
#[cfg(not(feature="flat"))]
        unsafe {
            self.mem.get_unchecked(row).get_unchecked(col)
        }
    }
}

/*
impl Index<usize> for SubseqMatrixRow {

    fn index(&self,  index : usize) -> SubseqInfo {
        unsafe { self.get_unchecked(index) }
    }
}
*/

pub struct CostMatrix {
    row_size: usize,
#[cfg(feature="flat")]
    data:     Vec<f64>,
#[cfg(not(feature="flat"))]
    data:     Vec< Vec<f64> >,
}

impl CostMatrix {
    fn new(n: usize) -> CostMatrix {
        Self {
            row_size: n,
#[cfg(feature="flat")]
            data: vec![0.0; n*n],
#[cfg(not(feature="flat"))]
            data: vec![ vec![0.0; n]; n],
        }
    }

    fn get(&self, i : usize, j : usize) -> f64 {
        unsafe {*self.get_unchecked(i).get_unchecked(j)}
    }
}

impl Index<(usize, usize)> for CostMatrix {
    type Output = f64;

    fn index (&'a self,  row_index: (usize, usize)) -> Self::Output {
        let (row, col) = row_index;
#[cfg(feature="flat")]
        unsafe {
            self.data.get_unchecked(self.row_size * row + col)
        }
#[cfg(not(feature="flat"))]
        unsafe {
            self.data.get_unchecked(row).get_unchecked(col)
        }
    }
}

#[derive(Debug, Clone)]
pub struct Info {
    pub c : CostMatrix,
    pub dimen : usize,
    pub rnd : Vec<usize>,
    pub rnd_index : usize,
}

#[derive(Debug, Clone)]
pub struct Solution {
    pub seq : SubseqMatrix,
    pub s : Vec<usize>,
    pub cost : f64,
}

pub trait Test<T> {
    fn get(&self, i : usize) -> T;
    fn set(&mut self, i : usize, value : T);
    fn swap(&mut self, i : usize, j : usize);
    fn index(&self, index : usize) -> T;
}


pub trait Access {
    fn set_C(&mut self, i : usize, j : usize, value : f64);
    fn set_T(&mut self, i : usize, j : usize, value : f64);
    fn set_W(&mut self, i : usize, j : usize, value : f64);

    fn get_C(&self, i : usize, j : usize) -> f64;
    fn get_T(&self, i : usize, j : usize) -> f64;
    fn get_W(&self, i : usize, j : usize) -> f64;

    //fn get_C_mut(&mut self, i : usize, j : usize) -> &mut f64;
}

/*
#[inline(always)]
fn to_1D(i : usize, j : usize, size : usize) -> usize {
    return i * size + j;
}

impl Access for tSeqData {

    /*
    fn get_C_mut(&mut self, i : usize, j : usize) -> &mut f64 {
        unsafe {mut self.get_unchecked_mut(i).get_unchecked_mut(j).C}
    }
    */

    #[inline]
    fn set_C(&mut self, i : usize, j : usize, value : f64) {
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked_mut(to_1D(i, j, sz::SIZE+1)).C = value;}}
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).C = value;}}
        //unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).C = value;}
    }
     
    #[inline]
    fn set_T(&mut self, i : usize, j : usize, value : f64) {
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked_mut(to_1D(i, j, sz::SIZE+1)).T = value;}}
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).T = value;}}
        //unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).T = value;}
    }

    #[inline]
    fn set_W(&mut self, i : usize, j : usize, value : f64) {
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked_mut(to_1D(i, j, sz::SIZE+1)).W = value;}}
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).W = value;}}
        //unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).W = value;}
    }

    #[inline]
    fn get_C(&self, i : usize, j : usize) -> f64 {
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked(i).get_unchecked(j).C}}
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked(to_1D(i, j, sz::SIZE+1)).C}}
        //unsafe {self.get_unchecked(i).get_unchecked(j).C}
    }

    #[inline]
    fn get_T(&self, i : usize, j : usize) -> f64 {
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked(i).get_unchecked(j).T}}
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked(to_1D(i, j, sz::SIZE+1)).T}}
        //unsafe {self.get_unchecked(i).get_unchecked(j).T}
    }

    #[inline]
    fn get_W(&self, i : usize, j : usize) -> f64 {
#[cfg(not(feature="flat"))]
        {unsafe {self.get_unchecked(i).get_unchecked(j).W}}
#[cfg(feature="flat")]
        {unsafe {self.get_unchecked(to_1D(i, j, sz::SIZE+1)).W}}
        //unsafe {self.get_unchecked(i).get_unchecked(j).W}
    }
}
*/
