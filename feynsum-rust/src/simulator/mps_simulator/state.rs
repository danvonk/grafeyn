use super::{dense_table::DenseStateTable, mps::MPSState, sparse_table::SparseStateTable};
use crate::utility;

use super::super::Compactifiable;
use crate::types::{BasisIdx, Complex};

#[derive(Debug)]
pub enum State<B: BasisIdx> {
    MPS(MPSState),
    Sparse(SparseStateTable<B>),
    Dense(DenseStateTable),
}

pub trait Table<B: BasisIdx> {
    fn put(&mut self, bidx: B, weight: Complex);
}

impl<B: BasisIdx> State<B> {
    pub fn num_nonzeros(&self) -> usize {
        match self {
            State::MPS(mps) => mps.num_nonzeros::<B>(),
            State::Sparse(table) => table.num_nonzeros(),
            State::Dense(table) => table.num_nonzeros(),
        }
    }
}

impl<B: BasisIdx> Compactifiable<B> for State<B> {
    fn compactify(self) -> Box<dyn Iterator<Item = (B, Complex)>> {
        match self {
            State::MPS(mps) => Box::new(mps.nonzeros().into_iter()),
            State::Sparse(table) => Box::new(
                table
                    .table
                    .into_iter()
                    .filter(|(_b, c)| utility::is_nonzero(*c)),
            ),
            State::Dense(table) => {
                Box::new(table.array.into_iter().enumerate().filter_map(|(idx, c)| {
                    if utility::is_nonzero(c) {
                        Some((B::from_idx(idx), c))
                    } else {
                        None
                    }
                }))
            }
        }
    }
}
