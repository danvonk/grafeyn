use ndarray::*;
use ndarray_linalg::*;

use crate::{
    circuit::GateDefn,
    types::{BasisIdx, Complex},
};

// TODO: we should be able to switch to/from this and the other representations
#[derive(Debug)]
pub struct MPSState {
    pub tensors: Vec<ndarray::Array3<Complex>>,
    pub _bond_dims: Vec<(usize, usize)>,
    pub n_sites: usize,
}

impl MPSState {
    /// Create a new MPS (0,...,0) state
    pub fn singleton(num_qubits: usize) -> Self {
        let (tensors, bond_dims) = (0..num_qubits)
            .map(|_| {
                let mut site = Array::zeros((1, 2, 1));
                site[[0, 0, 0]] = Complex::new(1.0, 0.0);
                (site, (1, 1))
            })
            .unzip();

        Self {
            tensors,
            _bond_dims: bond_dims,
            n_sites: num_qubits,
        }
    }

    pub fn nonzeros<B: BasisIdx>(&self) -> Vec<(B, Complex)> {
        let n = self.n_sites;
        if n == 0 {
            // no qubits => empty state
            return vec![(B::zeros(), Complex::new(1.0, 0.0))];
        }

        // maintain a list of partial expansions: (bitstring_so_far, bond_vector)
        let mut partials: Vec<(B, Array1<Complex>)> = Vec::new();

        // start from site 0
        {
            let a_0 = &self.tensors[0]; // shape (1,2,d0_out)
            let (_, _, d_out) = a_0.dim();
            // BFS: for i=0..1
            for i in 0..2 {
                // build bond vec
                let mut bond_vec = Array1::zeros(d_out);
                for a_out in 0..d_out {
                    let c = a_0[[0, i, a_out]];
                    if c != Complex::default() {
                        bond_vec[a_out] = c;
                    }
                }
                // skip if entire bond_vec is zero
                if bond_vec.iter().all(|&x| x == Complex::default()) {
                    continue;
                }
                let bits = B::from_idx(i);
                partials.push((bits, bond_vec));
            }
        }

        // Go site by site
        for site_idx in 1..n {
            let tensor_k = &self.tensors[site_idx];
            let (d_in, _, d_out) = tensor_k.dim();

            let mut new_partials = Vec::with_capacity(partials.len() * 2);
            // for each partial expansion, spawn two new expansions for i=0,1
            for (old_bits, old_vec) in partials.into_iter() {
                for i in 0..2 {
                    let mut new_vec = Array1::zeros(d_out);
                    // contract old_vec with site_k[:, i, :]
                    for alpha_in in 0..d_in {
                        let coeff_in = old_vec[alpha_in];
                        if coeff_in != Complex::default() {
                            for alpha_out in 0..d_out {
                                new_vec[alpha_out] += coeff_in * tensor_k[[alpha_in, i, alpha_out]];
                            }
                        }
                    }
                    // skip if entire new_vec is zero
                    if new_vec.iter().all(|&x| x == Complex::default()) {
                        continue;
                    }
                    let new_bits = (old_bits.as_idx() << 1) | i;
                    new_partials.push((B::from_idx(new_bits), new_vec));
                }
            }
            partials = new_partials;
        }

        // Each partial bond vector is length=1 (assuming open boundary with bond_out=1).
        // So the amplitude is at index 0. If it's nonzero, store it.
        let mut results = Vec::new();
        for (bits, vec_ampl) in partials {
            let amp = vec_ampl[0];
            if amp != Complex::default() {
                results.push((bits, amp));
            }
        }

        results
    }

    pub fn _num_nonzeros<B: BasisIdx>(&self) -> usize {
        self.nonzeros::<B>().iter().count()
    }

    pub fn apply_single_qubit_gate(&mut self, matrix: &Array2<Complex>, site: usize) {
        println!("Apply gate to site {}", site);
        let a = &mut self.tensors[site];
        let (bond_in, _, bond_out) = a.dim();

        // 2. For each pair of bond indices (alpha_in, alpha_out),
        //    multiply [A(alpha_in, 0, alpha_out), A(alpha_in, 1, alpha_out)] by the 2x2 `matrix`.
        for alpha_in in 0..bond_in {
            for alpha_out in 0..bond_out {
                // Extract current amplitudes for |0> and |1>
                let v0 = a[[alpha_in, 0, alpha_out]];
                let v1 = a[[alpha_in, 1, alpha_out]];

                // matrix * (v0, v1)^T
                let w0 = matrix[[0, 0]] * v0 + matrix[[0, 1]] * v1;
                let w1 = matrix[[1, 0]] * v0 + matrix[[1, 1]] * v1;

                // Store updated amplitudes back in A
                a[[alpha_in, 0, alpha_out]] = w0;
                a[[alpha_in, 1, alpha_out]] = w1;
            }
        }
    }

    pub fn apply_two_qubit_gate(&mut self, matrix: &Array2<Complex>, site1: usize, site2: usize) {
        // 1. Sort site1, site2 so that site1 < site2
        let (left_site, right_site) = if site1 < site2 {
            (site1, site2)
        } else {
            (site2, site1)
        };
        // For simplicity, assume these are adjacent qubits
        assert_eq!(
            right_site,
            left_site + 1,
            "This apply_two_qubit_gate is only implemented for adjacent sites."
        );

        // 2. Retrieve MPS tensors for these sites
        let a_left = self.tensors[left_site].clone(); // shape = (dL, 2, dM)
        let a_right = self.tensors[right_site].clone(); // shape = (dM, 2, dR)

        let (dL, _, dM) = a_left.dim();
        let (dM2, _, dR) = a_right.dim();
        assert_eq!(
            dM, dM2,
            "Bond mismatch between site1's right bond and site2's left bond."
        );

        // 3. Combine the two site tensors into a rank-4 array combined(dL, 2, 2, dR)
        let mut combined = Array4::<Complex>::zeros((dL, 2, 2, dR));
        for alphaL in 0..dL {
            for alphaM in 0..dM {
                for alphaR in 0..dR {
                    for i in 0..2 {
                        for j in 0..2 {
                            combined[[alphaL, i, j, alphaR]] +=
                                a_left[[alphaL, i, alphaM]] * a_right[[alphaM, j, alphaR]];
                        }
                    }
                }
            }
        }

        // 4. Reshape the physical indices (i, j) from (2,2) -> 4 and multiply by matrix(4x4)
        //    combined_2d has shape (dL, 4, dR), we multiply along the axis of size=4
        let mut combined_2d = combined.into_shape((dL, 4, dR)).unwrap();

        // We'll multiply matrix(4x4) * combined_2d along the "4" axis.
        // One approach: loop over (dL, dR) and do a local dot product for each slice.
        // Or we can flatten the last two dims to (dL*dR, 4) if we prefer. We'll do the loop approach:

        for alphaL in 0..dL {
            for alphaR in 0..dR {
                let mut vec_in = Array::<Complex, _>::zeros((4,));
                // Copy out the 4 values
                for idx in 0..4 {
                    vec_in[idx] = combined_2d[[alphaL, idx, alphaR]];
                }

                // matrix(4x4) * vec_in(4x1) => out(4x1)
                let vec_out = matrix.dot(&vec_in);

                // Store back
                for idx in 0..4 {
                    combined_2d[[alphaL, idx, alphaR]] = vec_out[idx];
                }
            }
        }

        // 5. Reshape back to rank-4 => (dL, 2, 2, dR)
        let updated_4d = combined_2d.into_shape((dL, 2, 2, dR)).unwrap();

        // 6. Now we do an SVD to split updated_4d back into two rank-3 tensors.
        //    6a. Flatten to a 2D matrix for SVD:
        //        shape = (dL * 2, 2 * dR)
        let (dL2, _, _, dR2) = updated_4d.dim(); // should be (dL, 2, 2, dR)
        let reshaped_for_svd = updated_4d.into_shape((dL2 * 2, 2 * dR2)).unwrap(); // (bond_in*2, 2*bond_out)

        // 6b. Perform the SVD using ndarray-linalg: (U, S, V^H)
        //     SVDInto consumes the data and returns newly allocated U, S, V^H.
        //

        // let svd_res = reshaped_for_svd.svd_into(true, true).unwrap();
        let (u_opt, s, vt_opt) = reshaped_for_svd.svd_into(true, true).unwrap();

        let u = u_opt.expect("U was requested"); // shape (dL2*2, rank)
        let vt = vt_opt.expect("V^H was requested"); // shape (rank, 2*dR2)

        // 7. Potentially truncate the singular values if you have a max bond dimension
        let max_bond_dim = 64; // example
        let rank_full = s.len();
        let chi = std::cmp::min(rank_full, max_bond_dim);

        // Slice out the first `chi` columns/rows
        let s_trunc = s.slice(s![0..chi]).to_owned();
        let u_trunc = u.slice(s![.., 0..chi]).to_owned();
        let vt_trunc = vt.slice(s![0..chi, ..]).to_owned();

        // 8. Distribute singular values (sqrt) if you prefer balanced approach:
        let s_sqrt: Array1<Complex> = s_trunc.mapv(|x| Complex::new(x.sqrt(), 0.0));

        // Multiply columns of U by s_sqrt
        let mut u_svd = u_trunc.to_owned();
        for c in 0..chi {
            let scale = s_sqrt[c];
            // multiply column c by scale
            let mut col_mut = u_svd.column_mut(c);
            col_mut *= scale;
        }

        // Multiply rows of V by s_sqrt
        let mut v_svd = vt_trunc.to_owned();
        for r in 0..chi {
            let scale = s_sqrt[r];
            let mut row_mut = v_svd.row_mut(r);
            row_mut *= scale;
        }

        // 9. Reshape U => a_left_new with shape (dL2, 2, chi)
        //    U was shape (dL2*2, chi)
        let a_left_new = u_svd.into_shape((dL2, 2, chi)).unwrap();

        // 10. Reshape V => a_right_new with shape (chi, 2, dR2)
        //     V^H had shape (chi, 2*dR2)
        //     We'll interpret it as shape (chi, 2, dR2).
        let a_right_new = v_svd.into_shape((chi, 2, dR2)).unwrap();

        // 11. Store updated site tensors
        self.tensors[left_site] = a_left_new;
        self.tensors[right_site] = a_right_new;

        // 12. Update bond dimensions
        self._bond_dims[left_site] = (dL2, chi); // bond_in= dL2, new bond_out= chi
        self._bond_dims[right_site] = (chi, dR2); // new bond_in= chi, bond_out= dR2
    }

    pub fn apply_two_qubit_gate_nonadjacent(
        &mut self,
        matrix: &Array2<Complex>,
        site1: usize,
        site2: usize,
    ) {
        // 1. Sort so we only handle q_low < q_high
        let (mut q_low, mut q_high) = if site1 < site2 {
            (site1, site2)
        } else {
            (site2, site1)
        };
        let orig_q_high = q_high;

        // 2. Move q_high left until it is directly next to q_low
        //    i.e. while q_high > q_low + 1, swap adjacent (q_high - 1, q_high)
        while q_high > q_low + 1 {
            self.apply_two_qubit_gate(
                &(GateDefn::Swap {
                    target1: q_high - 1,
                    target2: q_high,
                }
                .gate_to_matrix()
                .unwrap()),
                q_high - 1,
                q_high,
            );
            q_high -= 1;
        }

        // Now q_high == q_low + 1 => adjacent

        // 3. Apply the actual two-qubit gate
        self.apply_two_qubit_gate(matrix, q_low, q_high);

        // 4. Move q_high back to its original position using SWAPs to the right
        while q_high < orig_q_high {
            self.apply_two_qubit_gate(
                &(GateDefn::Swap {
                    target1: q_high,
                    target2: q_high + 1,
                }
                .gate_to_matrix()
                .unwrap()),
                q_high,
                q_high + 1,
            );
            q_high += 1;
        }
    }
}
