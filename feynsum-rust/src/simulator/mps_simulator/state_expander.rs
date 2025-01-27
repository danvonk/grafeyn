use super::sparse_table::SparseStateTable;
use super::{
    mps::MPSState,
    state::{State, Table},
};
use crate::{
    circuit::{Gate, GateDefn, PullApplyOutput, PushApplicable, PushApplyOutput},
    config::Config,
    simulator::{expected_cost, Compactifiable},
    types::{BasisIdx, Complex},
    utility,
};
use std::fmt::{self, Display, Formatter};

#[derive(Debug)]
pub enum ExpandMethod {
    MPS,
    Sparse,
}

impl Display for ExpandMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ExpandMethod::Sparse => write!(f, "sparse"),
            ExpandMethod::MPS => write!(f, "mps"),
        }
    }
}

pub struct ExpandResult<B: BasisIdx> {
    pub state: State<B>,
    pub num_nonzeros: usize,
    pub num_gate_apps: usize,
    pub method: ExpandMethod,
}

pub fn expand<B: BasisIdx>(
    gates: Vec<&Gate<B>>,
    config: &Config,
    num_qubits: usize,
    prev_num_nonzeros: usize,
    state: State<B>,
) -> ExpandResult<B> {
    // let (expected_density, _) = expected_cost(num_qubits, state.num_nonzeros(), prev_num_nonzeros);

    //if expected_density <= config.dense_threshold {
    //    log::warn!("Switching to sparse");
    //    expand_sparse(gates, state)
    //} else {
    //    log::warn!("Switching to MPS");
    expand_mps(gates, config, num_qubits, state)
    //}
}

fn expand_mps<B: BasisIdx>(
    gates: Vec<&Gate<B>>,
    config: &Config,
    num_qubits: usize,
    state: State<B>,
) -> ExpandResult<B> {
    let mut gate_apps = 0;

    let mut mps = match state {
        State::MPS(mps_state) => mps_state,
        State::Sparse(_) => MPSState::from_nonzeros(config, state, num_qubits),
        _ => todo!("Implement dense table"),
    };

    for g in gates {
        mps.apply_gate(config, g);
        gate_apps += 1;
    }

    ExpandResult {
        state: State::MPS(mps),
        num_gate_apps: gate_apps,
        num_nonzeros: 0,
        method: ExpandMethod::MPS,
    }
}

fn expand_sparse<B: BasisIdx>(gates: Vec<&Gate<B>>, state: State<B>) -> ExpandResult<B> {
    let mut table = SparseStateTable::<B>::new();

    let num_gate_apps = state
        .compactify()
        .map(|(bidx, weight)| apply_gates_sparsely(&gates, &mut table, bidx, weight))
        .sum();

    let num_nonzeros = table.num_nonzeros();

    ExpandResult {
        state: State::Sparse(table),
        num_nonzeros,
        num_gate_apps,
        method: ExpandMethod::Sparse,
    }
}

fn apply_gates_sparsely<B: BasisIdx>(
    gates: &[&Gate<B>],
    table: &mut impl Table<B>,
    bidx: B,
    weight: Complex,
) -> usize {
    if utility::is_zero(weight) {
        return 0;
    }
    if gates.is_empty() {
        table.put(bidx, weight);
        return 0;
    }

    match gates[0].push_apply(bidx, weight) {
        PushApplyOutput::Nonbranching(new_bidx, new_weight) => {
            1 + apply_gates_sparsely(&gates[1..], table, new_bidx, new_weight)
        }
        PushApplyOutput::Branching((new_bidx1, new_weight1), (new_bidx2, new_weight2)) => {
            let num_gate_apps_1 = apply_gates_sparsely(&gates[1..], table, new_bidx1, new_weight1);
            let num_gate_apps_2 = apply_gates_sparsely(&gates[1..], table, new_bidx2, new_weight2);
            1 + num_gate_apps_1 + num_gate_apps_2
        }
    }
}
