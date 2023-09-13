use std::sync::atomic::{AtomicU64, Ordering};

use crate::types::{BasisIdx, Complex};
use crate::utility;

use super::Table;

#[derive(Debug)]
pub struct DenseStateTable {
    pub array: Vec<AtomicU64>,
}

impl DenseStateTable {
    pub fn new(num_qubits: usize) -> Self {
        let capacity = 1 << num_qubits;
        // TODO: Check if the initialization is performance bottleneck
        Self {
            array: (0..capacity).map(|_| AtomicU64::new(0)).collect(),
        }
    }

    pub fn num_nonzeros(&self) -> usize {
        self.array
            .iter()
            .filter(|v| {
                let (re, im) = utility::unpack_complex(v.load(Ordering::Relaxed));
                utility::is_real_nonzero(re) || utility::is_real_nonzero(im)
            })
            .count()
    }
    pub fn atomic_put(&self, bidx: BasisIdx, weight: Complex) {
        // FIXME: We can use `put` method instead of `atomic_put` method if we
        // change the signature of `put` method from `&mut self` to &self
        let idx = bidx.into_idx();

        atomic_put(&self.array[idx], weight);
    }
}

impl Table for DenseStateTable {
    fn put(&mut self, bidx: BasisIdx, weight: Complex) {
        let idx = bidx.into_idx();

        atomic_put(&self.array[idx], weight);
    }
}

fn atomic_put(to: &AtomicU64, c: Complex) {
    let mut old = to.load(Ordering::Relaxed);
    let (mut old_re, mut old_im) = utility::unpack_complex(old);
    let (mut new_re, mut new_im) = (old_re + c.re, old_im + c.im);
    let mut new = utility::pack_complex(new_re, new_im);
    while let Err(actual) = to.compare_exchange(old, new, Ordering::Relaxed, Ordering::Relaxed) {
        old = actual;
        (old_re, old_im) = utility::unpack_complex(old);
        (new_re, new_im) = (old_re + c.re, old_im + c.im);
        new = utility::pack_complex(new_re, new_im);
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_atomic_put() {
        todo!()
    }
}