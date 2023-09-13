use std::cmp;
use std::fmt::{self, Display, Formatter};

use crate::circuit::{Gate, MaybeBranchingOutput, PushApplicable};
use crate::config::Config;
use crate::types::{BasisIdx, Complex, Real};
use crate::utility::is_zero;

use super::super::{Compactifiable, SimulatorError};
use super::state::{DenseStateTable, SparseStateTable, State, Table};

pub enum ExpandMethod {
    Sparse,
    PushDense,
    #[allow(dead_code)]
    PullDense,
}

impl Display for ExpandMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            ExpandMethod::Sparse => write!(f, "push sparse"),
            ExpandMethod::PushDense => write!(f, "push dense"),
            ExpandMethod::PullDense => write!(f, "pull dense"),
        }
    }
}

pub struct ExpandResult {
    pub state: State,
    pub num_nonzeros: usize,
    pub num_gate_apps: usize,
    pub method: ExpandMethod,
}

pub fn expand(
    gates: Vec<&Gate>,
    config: &Config,
    num_qubits: usize,
    prev_num_nonzeros: usize,
    state: State,
) -> Result<ExpandResult, SimulatorError> {
    let expected_cost = expected_cost(num_qubits, state.num_nonzeros(), prev_num_nonzeros);

    let all_gates_pullable = gates.iter().all(|_| false); // FIXME

    assert!(config.dense_threshold <= config.pull_threshold);

    if expected_cost < config.dense_threshold {
        expand_sparse(gates, state)
    } else if expected_cost >= config.pull_threshold && all_gates_pullable {
        expand_pull_dense(gates, num_qubits, state)
    } else {
        expand_push_dense(gates, num_qubits, state)
    }
}

fn expected_cost(num_qubits: usize, num_nonzeros: usize, prev_num_nonzeros: usize) -> Real {
    let max_num_states = 1 << num_qubits;
    let rate = Real::max(1.0, num_nonzeros as Real / prev_num_nonzeros as Real);
    let expected = cmp::min(max_num_states, (rate * num_nonzeros as Real) as i64);
    let expected_density = expected as Real / max_num_states as Real;
    let current_density = num_nonzeros as Real / max_num_states as Real;
    expected_density.max(current_density)
}

fn expand_sparse(gates: Vec<&Gate>, state: State) -> Result<ExpandResult, SimulatorError> {
    let mut table = SparseStateTable::new();

    let num_gate_apps = state
        .compactify()
        .try_fold(0, |num_gate_apps, (bidx, weight)| {
            Ok::<usize, SimulatorError>(
                num_gate_apps + apply_gates(&gates, &mut table, bidx, weight)?,
            )
        })?;

    let num_nonzeros = table.num_nonzeros();

    Ok(ExpandResult {
        state: State::Sparse(table),
        num_nonzeros,
        num_gate_apps,
        method: ExpandMethod::Sparse,
    })
}

fn expand_push_dense(
    gates: Vec<&Gate>,
    num_qubits: usize,
    state: State,
) -> Result<ExpandResult, SimulatorError> {
    let mut table = DenseStateTable::new(num_qubits);

    let num_gate_apps = state
        .compactify()
        .map(|(bidx, weight)| -> Result<usize, SimulatorError> {
            apply_gates(&gates, &mut table, bidx, weight)
        })
        .sum::<Result<usize, SimulatorError>>()?;

    let num_nonzeros = table.num_nonzeros();

    Ok(ExpandResult {
        state: State::Dense(table),
        num_nonzeros,
        num_gate_apps,
        method: ExpandMethod::PushDense,
    })
}

fn expand_pull_dense(
    _gates: Vec<&Gate>,
    _num_qubits: usize,
    _state: State,
) -> Result<ExpandResult, SimulatorError> {
    unimplemented!()
}

fn apply_gates(
    gates: &[&Gate],
    table: &mut impl Table,
    bidx: BasisIdx,
    weight: Complex,
) -> Result<usize, SimulatorError> {
    if is_zero(weight) {
        return Ok(0);
    }
    if gates.is_empty() {
        table.put(bidx, weight);
        return Ok(0);
    }

    match gates[0].push_apply(bidx, weight)? {
        MaybeBranchingOutput::OuptutOne((new_bidx, new_weight)) => {
            Ok(1 + apply_gates(&gates[1..], table, new_bidx, new_weight)?)
        }
        MaybeBranchingOutput::OutputTwo((new_bidx1, new_weight1), (new_bidx2, new_weight2)) => {
            let num_gate_apps_1 = apply_gates(&gates[1..], table, new_bidx1, new_weight1)?;
            let num_gate_apps_2 = apply_gates(&gates[1..], table, new_bidx2, new_weight2)?;
            Ok(1 + num_gate_apps_1 + num_gate_apps_2)
        }
    }
}