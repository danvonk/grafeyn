mod mps;
mod state;
mod state_expander;

use crate::circuit::Circuit;
use crate::config::Config;
use crate::gate_scheduler;
use crate::profile;
use crate::types::{BasisIdx, Real};

pub use state::State;
pub use state_expander::{expand, ExpandResult};

pub fn run<B: BasisIdx>(config: &Config, circuit: Circuit<B>) -> State {
    let num_gates = circuit.num_gates();
    let num_qubits = circuit.num_qubits;

    let mut num_gates_visited = 0;
    let mut state = State::MPS(mps::MPSState::singleton(num_qubits));
    let mut num_gate_apps = 0;

    let mut gate_scheduler = gate_scheduler::create_gate_scheduler(config, &circuit);

    let (duration, _) = profile!(loop {
        let these_gates = gate_scheduler
            .pick_next_gates()
            .into_iter()
            .map(|idx| &circuit.gates[idx])
            .collect::<Vec<_>>();

        log::debug!("applying gates: {:?}", these_gates);

        if these_gates.is_empty() {
            break;
        }

        let num_gates_visited_here = these_gates.len();

        let (
            duration,
            ExpandResult {
                state: new_state,
                num_gate_apps: num_gate_apps_here,
            },
        ) = profile!(expand::<B>(these_gates, config, num_qubits, state));

        //let density = {
        //    let max_num_states: u64 = 1 << num_qubits;
        //    num_nonzeros as Real / max_num_states as Real
        //};

        let throughput = (num_gate_apps_here as Real / 1e6) / duration.as_secs_f32();

        println!(
            "gate: {:<3} hop: {:<2} time: {:.4}s throughput: {:.2}M gates/s",
            num_gates_visited,
            num_gates_visited_here,
            duration.as_secs_f32(),
            throughput
        );

        num_gates_visited += num_gates_visited_here;
        num_gate_apps += num_gate_apps_here;
        state = new_state;
    });

    //let final_density = {
    //    let max_num_states: u64 = 1 << num_qubits;
    //    num_nonzeros as f64 / max_num_states as f64
    //};

    println!(
        "gate: {:<2} \ngate app count: {}, time: {}s",
        num_gates_visited,
        //final_density,
        //num_nonzeros,
        num_gate_apps,
        duration.as_secs_f32()
    );

    assert!(num_gates_visited >= num_gates);
    state
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser;
    use crate::types::BasisIdx64;

    #[test]
    fn test_run() {
        let config = Config::default();
        // bv_n30.qasm
        let circuit = Circuit::new(
            parser::parse_program(
                r#"
OPENQASM 2.0;
include "qelib1.inc";
qreg q0[30];
creg c0[30];
h q0[0];
h q0[1];
h q0[2];
h q0[3];
h q0[4];
h q0[5];
h q0[6];
h q0[7];
h q0[8];
h q0[9];
h q0[10];
h q0[11];
h q0[12];
h q0[13];
h q0[14];
h q0[15];
h q0[16];
h q0[17];
h q0[18];
h q0[19];
h q0[20];
h q0[21];
h q0[22];
h q0[23];
h q0[24];
h q0[25];
h q0[26];
h q0[27];
h q0[28];
x q0[29];
h q0[29];
            "#,
            )
            .unwrap(),
        )
        .unwrap();

        let _state = run::<BasisIdx64>(&config, circuit);

        println!("{:?}", _state);
    }
}
