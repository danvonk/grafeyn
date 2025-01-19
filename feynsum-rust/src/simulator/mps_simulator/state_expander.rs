use super::state::State;
use crate::circuit::Gate;
use crate::config::Config;
use crate::types::BasisIdx;

pub struct ExpandResult {
    pub state: State,
    pub num_gate_apps: usize,
}

pub fn expand<B: BasisIdx>(
    gates: Vec<&Gate<B>>,
    _config: &Config,
    _num_qubits: usize,
    state: State,
) -> ExpandResult {
    let mut mps = match state {
        State::MPS(mps_state) => mps_state,
    };

    let mut gate_apps = 0;

    for g in gates {
        let mat = g.defn.gate_to_matrix().unwrap();
        match g.defn.affects_qubits() {
            1 => match g.defn {
                crate::circuit::GateDefn::Hadamard(qindex)
                | crate::circuit::GateDefn::X(qindex) => mps.apply_single_qubit_gate(&mat, qindex),
                _ => (),
            },
            _ => (),
        }
        gate_apps += 1;
    }

    ExpandResult {
        state: State::MPS(mps),
        num_gate_apps: gate_apps,
    }
}
