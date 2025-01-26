use crate::config::Config;
use crate::utility::is_zero;
use nalgebra::*;

use crate::{
    circuit::{Gate, GateDefn, Unitary, UnitaryMatrix},
    types::{BasisIdx, Complex},
};

// TODO: we should be able to switch to/from this and the other representations
#[derive(Debug)]
pub struct MPSState {
    // MPS tensors, each represented as a collection of matrices (2 slices per site)
    // One component in the pair for Left |0> and Right |1> respectively.
    pub tensors: Vec<(DMatrix<Complex>, DMatrix<Complex>)>,
    // Bond dimensions between sites: (dL, dR)
    pub bond_dims: Vec<(usize, usize)>,
    // Number of sites (qubits)
    pub n_sites: usize,
}

impl MPSState {
    /// Create a new MPS (0,...,0) state
    pub fn singleton(num_qubits: usize) -> Self {
        let mut tensors = Vec::new();
        let mut bond_dims = Vec::new();

        for _ in 0..num_qubits {
            // Initialize |0> tensor
            let tensor_0 = DMatrix::from_element(1, 1, Complex::new(1.0, 0.0)); // |0>
                                                                                // no contribution from |1> state
            let tensor_1 = DMatrix::from_element(1, 1, Complex::new(0.0, 0.0)); // |1>

            tensors.push((tensor_0, tensor_1));
            bond_dims.push((1, 1)); // Initial bond dimensions are 1
        }

        Self {
            tensors,
            bond_dims,
            n_sites: num_qubits,
        }
    }

    pub fn nonzeros<B: BasisIdx>(&self) -> Vec<(B, Complex)> {
        let n = self.n_sites;
        if n == 0 {
            // No qubits => empty state
            return vec![(B::zeros(), Complex::new(1.0, 0.0))];
        }

        // Maintain a list of partial expansions: (bitstring_so_far, bond_vector)
        let mut partials: Vec<(B, DVector<Complex>)> = Vec::new();

        // Start from site 0
        {
            let (tensor_0, tensor_1) = &self.tensors[0]; // Matrices for |0⟩ and |1⟩
                                                         // let d_out = tensor_0.ncols(); // Right bond dimension

            // Process |0⟩
            let bond_vec_0 = tensor_0.row(0).transpose(); // Extract bond vector for |0⟩
            if !bond_vec_0.iter().all(|&c| is_zero(c)) {
                partials.push((B::from_idx(0), bond_vec_0));
            }

            // Process |1⟩
            let bond_vec_1 = tensor_1.row(0).transpose(); // Extract bond vector for |1⟩
            if !bond_vec_1.iter().all(|&c| is_zero(c)) {
                partials.push((B::from_idx(1), bond_vec_1));
            }
        }

        // Iterate through the rest of the sites
        for site in 1..n {
            let mut next_partials: Vec<(B, DVector<Complex>)> = Vec::new();
            let (tensor_0, tensor_1) = &self.tensors[site]; // Matrices for |0⟩ and |1⟩
            let d_in = tensor_0.nrows(); // Left bond dimension
            let d_out = tensor_0.ncols(); // Right bond dimension

            // Process each partial expansion from the previous site
            for (bits, bond_vec) in partials {
                // Process |0⟩
                let mut new_bond_vec_0 = DVector::zeros(d_out);
                for alpha_in in 0..d_in {
                    for alpha_out in 0..d_out {
                        new_bond_vec_0[alpha_out] +=
                            bond_vec[alpha_in] * tensor_0[(alpha_in, alpha_out)];
                    }
                }
                if !new_bond_vec_0.iter().all(|&c| is_zero(c)) {
                    next_partials.push((bits.clone(), new_bond_vec_0)); // No change to bits for |0⟩
                }

                // Process |1⟩
                let mut new_bond_vec_1 = DVector::zeros(d_out);
                for alpha_in in 0..d_in {
                    for alpha_out in 0..d_out {
                        new_bond_vec_1[alpha_out] +=
                            bond_vec[alpha_in] * tensor_1[(alpha_in, alpha_out)];
                    }
                }
                if !new_bond_vec_1.iter().all(|&c| is_zero(c)) {
                    next_partials.push((bits.set(site), new_bond_vec_1)); // Set the bit for |1⟩
                }
            }

            partials = next_partials;
        }

        // Extract the final non-zero states and their amplitudes
        let mut nonzeros = Vec::new();
        for (bits, bond_vec) in partials {
            // The final bond vector should have size 1 (d_out = 1 at the last site)
            assert_eq!(bond_vec.len(), 1);
            nonzeros.push((bits, bond_vec[0]));
        }
        nonzeros
    }

    pub fn num_nonzeros<B: BasisIdx>(&self) -> usize {
        self.nonzeros::<B>().iter().count()
    }

    /// Apply the unitary matrix of a gate to a chosen site
    fn apply_single_qubit_gate(&mut self, gate: &UnitaryMatrix, site: usize) {
        let (tensor_0, tensor_1) = &mut self.tensors[site];
        let mat = &gate.mat;

        // Apply the gate to the |0> and |1> components of the site
        let new_tensor_0 = tensor_0.clone() * mat[(0, 0)] + tensor_1.clone() * mat[(0, 1)];
        let new_tensor_1 = tensor_0.clone() * mat[(1, 0)] + tensor_1.clone() * mat[(1, 1)];

        *tensor_0 = new_tensor_0;
        *tensor_1 = new_tensor_1;
    }

    /// Apply the unitary matrix of a two-qubit gate to the chosen sites. They must be adjacent.
    fn apply_two_qubit_gate(
        &mut self,
        config: &Config,
        gate: &UnitaryMatrix,
        site1: usize,
        site2: usize,
    ) {
        // Ensure site1 < site2 for simplicity (non-symmetric gates must be switched for their equivalent)
        let (left_site, right_site) = if site1 < site2 {
            (site1, site2)
        } else {
            (site2, site1)
        };

        // Ensure the sites are adjacent
        assert_eq!(
            right_site,
            left_site + 1,
            "apply_two_qubit_gate is only implemented for adjacent sites."
        );

        // Retrieve tensors for the two sites
        let (tensor1_0, tensor1_1) = &self.tensors[left_site];
        let (tensor2_0, tensor2_1) = &self.tensors[right_site];

        let (bond_left, bond_middle) = self.bond_dims[left_site];
        let (_, bond_right) = self.bond_dims[right_site];

        println!(
            "Bond left: {}, Bond middle: {}, Bond right: {}",
            bond_left, bond_middle, bond_right
        );

        let mut combined = DMatrix::from_fn(bond_left * bond_right, 4, |row, col| {
            // row in [0, .. bondL * bondR)
            // col in [0, .. 4)
            // Decompose row into (alphaL, alphaR) pair
            let alpha_l = row / bond_right; // in [0..bond_left)
            let alpha_r = row % bond_right; // in [0..bond_right)

            // Decompose col into (i, j)
            let i = col / 2; // in [0..2)
            let j = col % 2; // in [0..2)

            // Sum over alphaM (the middle bond), to contract the MPS
            let mut val = Complex::new(0.0, 0.0);
            for alpha_middle in 0..bond_middle {
                // Choose correct site1/site2 blocks depending on (i,j)
                // i => {tensor1_0, tensor1_1}, j => {tensor2_0, tensor2_1}
                match (i, j) {
                    (0, 0) => {
                        val +=
                            tensor1_0[(alpha_l, alpha_middle)] * tensor2_0[(alpha_middle, alpha_r)]
                    }
                    (0, 1) => {
                        val +=
                            tensor1_0[(alpha_l, alpha_middle)] * tensor2_1[(alpha_middle, alpha_r)];
                    }
                    (1, 0) => {
                        val +=
                            tensor1_1[(alpha_l, alpha_middle)] * tensor2_0[(alpha_middle, alpha_r)];
                    }
                    (1, 1) => {
                        val +=
                            tensor1_1[(alpha_l, alpha_middle)] * tensor2_1[(alpha_middle, alpha_r)];
                    }
                    _ => unreachable!(),
                }
            }
            val
        });

        println!(
            "shapes: doing combined * gate = {:?} * {:?}",
            combined.shape(),
            gate.mat.shape(),
        );

        combined = combined * gate.mat.clone();

        // We'll create `expanded` with shape (dL*2, 2*dR).
        // row => alphaL*2 + i
        // col => alphaR*2 + j
        let mut expanded =
            DMatrix::from_element(bond_left * 2, 2 * bond_right, Complex::new(0.0, 0.0));

        for new_row in 0..(bond_left * 2) {
            let alpha_left = new_row / 2; // 0..(bond_left-1)
            let i = new_row % 2; // i in {0,1}
            for new_col in 0..(2 * bond_right) {
                let alpha_right = new_col / 2;
                let j = new_col % 2;

                // get from row = alpha_left*dR + alpha_right, col = i*2 + j
                let row = alpha_left * bond_right + alpha_right; // 0..(dL*dR)
                let col = i * 2 + j; // 0..4

                let val = combined[(row, col)];
                expanded[(new_row, new_col)] = val;
            }
        }

        println!(
            "Expanded has shape {:?}, whilst (dL * 2, 2 * dR) = ({},{})",
            expanded.shape(),
            2 * bond_left,
            2 * bond_right
        );

        // Perform SVD on updated tensor
        let svd = expanded.svd(true, true);
        let u = svd.u.unwrap();
        let sigma = svd.singular_values;
        let vt = svd.v_t.unwrap();

        // Truncate if needed
        let chi = sigma.len().min(config.bond_dimension_threshold);

        // Scale symmetrically by the singular values
        let mut u_scaled = u.columns(0, chi).into_owned();
        let mut vt_scaled = vt.rows(0, chi).into_owned();

        for i in 0..chi {
            let sqrt_sigma = sigma[i].sqrt();
            u_scaled.column_mut(i).scale_mut(sqrt_sigma);
            vt_scaled.row_mut(i).scale_mut(sqrt_sigma);
        }

        println!(
            "Shape of U: {:?}, V^T: {:?} (Want (dL*2, Chi) and (Chi, dR*2))",
            u_scaled.shape(),
            vt_scaled.shape()
        );

        // Left site: (tensor0_left, tensor1_left), each have shape (dL, chi)
        let mut left_0 = DMatrix::zeros(bond_left, chi);
        let mut left_1 = DMatrix::zeros(bond_left, chi);

        // Right site: (tensor0_right, tensor1_right), each have shape (chi, bondRight)
        let mut right_0 = DMatrix::zeros(chi, bond_right);
        let mut right_1 = DMatrix::zeros(chi, bond_right);

        // U has shape (dL*2, chi), we do the row => alphaL, i split:
        for row in 0..(2 * bond_left) {
            let alpha_left = row / 2;
            let i = row % 2;
            for alpha_middle in 0..chi {
                let val = u_scaled[(row, alpha_middle)];
                if i == 0 {
                    left_0[(alpha_left, alpha_middle)] = val;
                } else {
                    left_1[(alpha_left, alpha_middle)] = val;
                }
            }
        }

        // V^T has shape (chi, 2 * bond_right), do:
        for row in 0..chi {
            for col in 0..(2 * bond_right) {
                let alpha_right = col / 2;
                let j = col % 2;
                let val = vt_scaled[(row, col)];
                if j == 0 {
                    right_0[(row, alpha_right)] = val;
                } else {
                    right_1[(row, alpha_right)] = val;
                }
            }
        }

        // Update MPS state
        self.tensors[left_site] = (left_0, left_1);
        self.tensors[right_site] = (right_0, right_1);

        // Update bond dims
        self.bond_dims[left_site] = (bond_left, chi);
        self.bond_dims[right_site] = (chi, bond_right);
    }

    fn apply_two_qubit_gate_nonadjacent<B: BasisIdx>(
        &mut self,
        config: &Config,
        matrix: &UnitaryMatrix,
        site1: usize,
        site2: usize,
    ) {
        // 1. Sort so we only handle q_low < q_high
        let (q_low, mut q_high) = if site1 < site2 {
            (site1, site2)
        } else {
            (site2, site1)
        };
        let orig_q_high = q_high;

        // 2. Move q_high left until it is directly next to q_low
        //    i.e. while q_high > q_low + 1, swap adjacent (q_high - 1, q_high)
        while q_high > q_low + 1 {
            self.apply_two_qubit_gate(
                config,
                &(Gate::<B>::new(GateDefn::Swap {
                    target1: q_high - 1,
                    target2: q_high,
                })
                .unitary()),
                q_high - 1,
                q_high,
            );
            q_high -= 1;
        }

        // Now q_high == q_low + 1 => adjacent

        // 3. Apply the actual two-qubit gate
        self.apply_two_qubit_gate(config, matrix, q_low, q_high);

        // 4. Move q_high back to its original position using SWAPs to the right
        while q_high < orig_q_high {
            self.apply_two_qubit_gate(
                config,
                &(Gate::<B>::new(GateDefn::Swap {
                    target1: q_high,
                    target2: q_high + 1,
                })
                .unitary()),
                q_high,
                q_high + 1,
            );
            q_high += 1;
        }
    }

    pub fn apply_gate<B: BasisIdx>(&mut self, config: &Config, gate: &Gate<B>) {
        match gate.defn {
            GateDefn::Hadamard(qindex)
            | GateDefn::S(qindex)
            | GateDefn::Sdg(qindex)
            | GateDefn::SqrtX(qindex)
            | GateDefn::SqrtXdg(qindex)
            | GateDefn::Tdg(qindex)
            | GateDefn::T(qindex)
            | GateDefn::X(qindex)
            | GateDefn::PauliY(qindex)
            | GateDefn::PauliZ(qindex)
            | GateDefn::Phase { target: qindex, .. }
            | GateDefn::RX { target: qindex, .. }
            | GateDefn::RY { target: qindex, .. }
            | GateDefn::RZ { target: qindex, .. }
            | GateDefn::U { target: qindex, .. } => {
                self.apply_single_qubit_gate(&gate.unitary(), qindex)
            }
            GateDefn::CZ { control, target }
            | GateDefn::CX { control, target }
            | GateDefn::CPhase {
                control, target, ..
            } => {
                let (left_site, right_site, mat) = if control < target {
                    (control, target, gate.unitary_rev())
                } else {
                    (target, control, gate.unitary())
                };

                if left_site + 1 == right_site {
                    self.apply_two_qubit_gate(config, &mat, left_site, right_site);
                } else {
                    self.apply_two_qubit_gate_nonadjacent::<B>(config, &mat, left_site, right_site);
                }
            }
            GateDefn::Swap { target1, target2 } => {
                let (left_site, right_site) = if target1 < target2 {
                    (target1, target2)
                } else {
                    (target2, target1)
                };
                if left_site + 1 == right_site {
                    self.apply_two_qubit_gate(config, &gate.unitary(), left_site, right_site);
                } else {
                    self.apply_two_qubit_gate_nonadjacent::<B>(
                        config,
                        &gate.unitary(),
                        left_site,
                        right_site,
                    );
                }
            }
            GateDefn::FSim {
                left: left_site,
                right: right_site,
                ..
            } => {
                if left_site + 1 == right_site {
                    self.apply_two_qubit_gate(config, &gate.unitary(), left_site, right_site);
                } else {
                    self.apply_two_qubit_gate_nonadjacent::<B>(
                        config,
                        &gate.unitary(),
                        left_site,
                        right_site,
                    );
                }
            }
            // We don't handle >= 3 qubit gates, they must have been decomposed already
            GateDefn::CSwap { .. } | GateDefn::CCX { .. } | GateDefn::Other { .. } => {
                log::warn!(
                    "Skipping gate {:?} as 3-qubit gates are not implemented by the MPS simulator.",
                    gate.defn
                );
            }
        };
    }
}
