use ndarray::*;

use crate::types::{BasisIdx, Complex};

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
}
