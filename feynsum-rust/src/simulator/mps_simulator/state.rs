use super::mps::MPSState;

use super::super::Compactifiable;
use crate::types::{BasisIdx, Complex};

#[derive(Debug)]
pub enum State {
    MPS(MPSState),
}

impl<B: BasisIdx> Compactifiable<B> for State {
    fn compactify(self) -> Box<dyn Iterator<Item = (B, Complex)>> {
        match self {
            State::MPS(mps) => Box::new(mps.nonzeros().into_iter()),
        }
    }
}
