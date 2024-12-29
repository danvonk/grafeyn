use super::state::State;
use crate::circuit::{Gate, GateDefn};
use crate::config::Config;
use crate::types::BasisIdx;

pub struct ExpandResult {
    pub state: State,
    pub num_gate_apps: usize,
    // pub num_nonzeros: usize,
}

pub fn expand<B: BasisIdx>(
    gates: Vec<&Gate<B>>,
    config: &Config,
    _num_qubits: usize,
    state: State,
) -> ExpandResult {
    let mut mps = match state {
        State::MPS(mps_state) => mps_state,
    };

    let mut gate_apps = 0;

    for g in gates {
        mps.apply_gate(config, g);
        gate_apps += 1;
    }

    ExpandResult {
        state: State::MPS(mps),
        num_gate_apps: gate_apps,
    }
}
