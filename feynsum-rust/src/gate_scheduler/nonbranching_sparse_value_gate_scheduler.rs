use super::{utility, GateScheduler};
use crate::gate_scheduler::utility::next_touch;
use crate::types::{GateIndex, QubitIndex};

pub struct NonBranchingSparseValueGateScheduler<'a> {
    frontier: Vec<GateIndex>,
    num_gates: usize,
    num_qubits: usize,
    gate_touches: Vec<&'a [QubitIndex]>,
    gate_is_branching: Vec<bool>,
    num_of_next_gates: usize,
    sparsity: Vec<usize>,
}

impl<'a> GateScheduler for NonBranchingSparseValueGateScheduler<'a> {
    fn pick_next_gates(&mut self) -> Vec<GateIndex> {
        let mut next_gates = Vec::<GateIndex>::new();
        let mut selection: GateIndex;
        let mut next_gi: GateIndex;
        let mut num_branches: usize = 0;

        'pool: while num_branches < self.num_of_next_gates {
            selection = match self.first_ok_touch() {
                Some(x) => x,
                None => break 'pool,
            };

            for qi in 0..self.num_qubits {
                next_gi = self.frontier[qi];

                if next_gi < self.num_gates
                    && utility::okay_to_visit(
                        self.num_gates,
                        &self.gate_touches,
                        &self.frontier,
                        next_gi,
                    )
                {
                    selection = match (
                        self.sparsity[selection] > self.sparsity[next_gi],
                        self.gate_is_branching[selection],
                        self.gate_is_branching[next_gi],
                    ) {
                        (false, false, false) => next_gi,
                        (false, false, true) => selection, //next_gi if sparsity is valued over branching
                        (false, true, false) => next_gi,
                        (false, true, true) => next_gi,
                        (true, false, false) => selection,
                        (true, false, true) => selection,
                        (true, true, false) => next_gi, //selection if sparsity is valued over branching
                        (true, true, true) => selection,
                    };
                }
            }

            utility::mark_as_visit(
                self.num_gates,
                &self.gate_touches,
                &mut self.frontier,
                selection,
            );

            if self.gate_is_branching[selection] {
                num_branches += 1;
            }

            next_gates.push(selection);
        }

        next_gates
    }
}

impl<'a> NonBranchingSparseValueGateScheduler<'a> {
    pub fn new(
        num_gates: usize,
        num_qubits: usize,
        gate_touches: Vec<&'a [QubitIndex]>,
        gate_is_branching: Vec<bool>,
        sparsity: Vec<usize>,
    ) -> Self {
        Self {
            frontier: (0..num_qubits)
                .map(|qi| next_touch(num_gates, &gate_touches, qi, 0))
                .collect(),
            num_gates,
            num_qubits,
            gate_touches,
            gate_is_branching,
            num_of_next_gates: 2,
            sparsity,
        }
    }

    fn first_ok_touch(&mut self) -> Option<GateIndex> {
        let mut next_gi: GateIndex;
        for qi in 0..self.num_qubits {
            next_gi = self.frontier[qi];
            if next_gi < self.num_gates
                && utility::okay_to_visit(
                    self.num_gates,
                    &self.gate_touches,
                    &self.frontier,
                    next_gi,
                )
            {
                return Some(next_gi);
            }
        }

        None
    }
}
