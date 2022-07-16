use crate::types::{tSeqData, tCostData, tSeqInfo};
use crate::data::sz;

pub trait Kba {
    fn New() -> Self ;
}

impl Kba for tSeqData {
    fn New() -> tSeqData {
        Box::new([[tSeqInfo {C:0.0, W:0.0, T:0.0}; sz::SIZE+1]; sz::SIZE+1])
    }
}

pub trait Q {
    fn get(&self, i : usize, j : usize) -> f64;
    fn New() -> Self;
}

impl Q for tCostData {
    fn get(&self, i : usize, j : usize) -> f64 {
        unsafe {*self.get_unchecked(i).get_unchecked(j)}
    }

    fn New() -> tCostData {
        Box::new([[0.0; sz::SIZE]; sz::SIZE])
    }
}


pub trait Test<T> {
    fn get(&self, i : usize) -> T;
    fn set(&mut self, i : usize, value : T);
    fn swap(&mut self, i : usize, j : usize);
    fn index(&self, index : usize) -> T;
}

impl<T> Test<T> for Vec<T> 
where T: Copy, {

    fn set(&mut self, i : usize, value : T) {
        unsafe {*self.get_unchecked_mut(i) = value;}
    }
    
    fn get(&self, i : usize) -> T {
        unsafe {*self.get_unchecked(i)}
    }

    fn swap(&mut self, i : usize, j : usize) {
        let tmp = self.get(i); 
        self.set(i, self.get(j)); 
        self.set(j, tmp);
    }

    #[inline]
    fn index(&self, index : usize) -> T {
        println!("opa");
        unsafe {*self.get_unchecked(index)}
    }
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

impl Access for tSeqData {

    /*
    fn get_C_mut(&mut self, i : usize, j : usize) -> &mut f64 {
        unsafe {mut self.get_unchecked_mut(i).get_unchecked_mut(j).C}
    }
    */

    #[inline]
    fn set_C(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).C = value;}
    }
     
    #[inline]
    fn set_T(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).T = value;}
    }

    #[inline]
    fn set_W(&mut self, i : usize, j : usize, value : f64) {
        unsafe {self.get_unchecked_mut(i).get_unchecked_mut(j).W = value;}
    }




    #[inline]
    fn get_C(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).C}
    }

    #[inline]
    fn get_T(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).T}
    }

    #[inline]
    fn get_W(&self, i : usize, j : usize) -> f64 {
        unsafe {self.get_unchecked(i).get_unchecked(j).W}
    }
}

/*
impl<T> Index<usize> for Vec<T> {

    fn index(&self,  index : usize) -> T {
        unsafe {self.get_unchecked(index)}
    }
}
*/

