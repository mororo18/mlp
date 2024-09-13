use std::ops::{Index, IndexMut};

pub const SWAP        : usize = 0;
pub const REINSERTION : usize = 1;
pub const OR_OPT_2    : usize = 2;
pub const OR_OPT_3    : usize = 3;
pub const TWO_OPT     : usize = 4;

#[derive(Debug, Copy, Clone)]
pub struct SubseqInfo {
   pub  t : f64,
   pub  c : f64,
   pub  w : f64,
}

impl SubseqInfo {
    fn zeored () -> Self {
        Self { t: 0.0, c: 0.0, w: 0.0 }
    }
}

#[derive(Debug, Clone)]
pub
struct SubseqMatrix {
    row_size: usize,
#[cfg(feature="flat")]
    mem: Vec<SubseqInfo>,
#[cfg(not(feature="flat"))]
    data:  Vec<Vec<SubseqInfo>>,
}

impl SubseqMatrix {
    pub
    fn new (n: usize) -> Self {
        Self { 
            row_size: n+1,
#[cfg(feature="flat")]
            mem: vec![SubseqInfo::zeored(); (n+1) * (n+1)],
#[cfg(not(feature="flat"))]
            data: vec![ vec![ SubseqInfo::zeored(); n+1]; n+1],
        }
    }
}

impl Index<(usize, usize)> for SubseqMatrix {
    type Output = SubseqInfo;

    fn index(&self,  index: (usize, usize)) -> &Self::Output {
        let (row, col) = index;
        unsafe {
            #[cfg(feature = "flat")]
            {
                self.data.get_unchecked(self.row_size * row + col)
            }

            #[cfg(not(feature = "flat"))]
            {
                self.data.get_unchecked(row).get_unchecked(col)
            }
        }
    }
}

impl IndexMut<(usize, usize)> for SubseqMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (row, col) = index;
        unsafe {
            #[cfg(feature = "flat")]
            {
                self.data.get_unchecked_mut(self.row_size * row + col)
            }

            #[cfg(not(feature = "flat"))]
            {
                self.data.get_unchecked_mut(row).get_unchecked_mut(col)
            }
        }
    }
}


#[derive(Debug, Clone)]
pub struct CostMatrix {
    row_size: usize,
#[cfg(feature="flat")]
    data:     Vec<f64>,
#[cfg(not(feature="flat"))]
    data:     Vec< Vec<f64> >,
}

impl CostMatrix {
    pub
    fn new(n: usize) -> CostMatrix {
        Self {
            row_size: n,
#[cfg(feature="flat")]
            data: vec![0.0; n*n],
#[cfg(not(feature="flat"))]
            data: vec![ vec![0.0; n]; n],
        }
    }
}

impl Index<(usize, usize)> for CostMatrix {
    type Output = f64;

    fn index (& self,  index: (usize, usize)) -> &Self::Output {
        let (row, col) = index;
        unsafe {
            #[cfg(feature = "flat")]
            {
                self.data.get_unchecked(self.row_size * row + col)
            }

            #[cfg(not(feature = "flat"))]
            {
                self.data.get_unchecked(row).get_unchecked(col)
            }
        }
    }
}

impl IndexMut<(usize, usize)> for CostMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (row, col) = index;
        unsafe {
            #[cfg(feature = "flat")]
            {
                self.data.get_unchecked_mut(self.row_size * row + col)
            }

            #[cfg(not(feature = "flat"))]
            {
                self.data.get_unchecked_mut(row).get_unchecked_mut(col)
            }
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

impl Index<usize> for Solution {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.s[index]
    }
}

