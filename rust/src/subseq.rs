use std::ops::{Index, IndexMut};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operator {
    Swap = 0,
    Reinsertion = 1,
    OrOpt2 = 2,
    OrOpt3 = 3,
    TwoOpt = 4,
}

#[derive(Debug, Copy, Clone)]
pub struct SubseqInfo {
   pub  t : f64,
   pub  c : f64,
   pub  w : f64,
}

impl SubseqInfo {
    #[inline]
    fn zeored () -> Self {
        Self { t: 0.0, c: 0.0, w: 0.0 }
    }
}

#[derive(Debug, Clone)]
pub
struct SubseqMatrix {
#[cfg(feature="flat")]
    row_size: usize,
#[cfg(feature="flat")]
    data: Vec<SubseqInfo>,
#[cfg(not(feature="flat"))]
    data:  Vec<Vec<SubseqInfo>>,
}

impl SubseqMatrix {
    pub
    fn new (n: usize) -> Self {
        Self { 
#[cfg(feature="flat")]
            row_size: n+1,
#[cfg(feature="flat")]
            data: vec![SubseqInfo::zeored(); (n+1) * (n+1)],
#[cfg(not(feature="flat"))]
            data: vec![ vec![ SubseqInfo::zeored(); n+1]; n+1],
        }
    }
}

impl Index<(usize, usize)> for SubseqMatrix {
    type Output = SubseqInfo;
    #[inline]
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
    #[inline]
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
#[cfg(feature="flat")]
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
#[cfg(feature="flat")]
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

    #[inline]
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
    #[inline]
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

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        unsafe { self.s.get_unchecked(index) }
    }
}

