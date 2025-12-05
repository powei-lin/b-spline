use std::collections::HashMap;
use std::num::NonZero;
use std::sync::Arc;

use crate::common::{compute_base_coefficients, compute_blending_matrix};
use crate::traits::*;
use tiny_solver::loss_functions::HuberLoss;
use tiny_solver::manifold::{AutoDiffManifold, Manifold};
use tiny_solver::{self, GaussNewtonOptimizer, Optimizer, manifold::so3::SO3, na};

pub struct RvecManifold;
impl<T: na::RealField> AutoDiffManifold<T> for RvecManifold {
    fn plus(&self, x: na::DVectorView<T>, delta: na::DVectorView<T>) -> na::DVector<T> {
        (SO3::exp(x) * SO3::exp(delta)).log()
    }

    fn minus(&self, y: na::DVectorView<T>, x: na::DVectorView<T>) -> na::DVector<T> {
        let y_so3 = SO3::from_vec(y);
        let x_so3_inv = SO3::from_vec(x).inverse();
        (x_so3_inv * y_so3).log()
    }
}

impl Manifold for RvecManifold {
    fn tangent_size(&self) -> NonZero<usize> {
        NonZero::new(3).unwrap()
    }
}

#[derive(Debug, Clone)]
pub struct SO3Bspline<const N: usize> {
    timestamp_start_ns: u64,
    spacing_ns: u64,
    pub knots: Vec<[f64; 3]>,
    blending_matrix: na::DMatrix<f64>,
    first_derivative_bases: na::DVector<f64>,
}

pub fn so3_from_u_and_knots<T: na::RealField>(
    u: f64,
    knots: &[na::DVector<T>],
    blending_matrix: &na::DMatrix<f64>,
) -> SO3<T> {
    let n = knots.len();
    let uv = na::DVector::from_fn(n, |i, _| u.powi(i as i32));
    let kv = (blending_matrix * uv).cast::<T>();
    let mut r = knots[0].to_so3();
    for j in 1..n {
        let k = kv[j].clone();
        let r_j = knots[j].to_so3();
        let r_j_minus_inv = knots[j - 1].to_so3().inverse();
        let rj_vec = (r_j_minus_inv * r_j).log() * k;
        r = r * rj_vec.to_so3();
    }
    r
}

struct RotationCost {
    u: f64,
    rvec: na::Vector3<f64>,
    blending_matrix: na::DMatrix<f64>,
}
impl RotationCost {
    pub fn new(u: f64, rvec: &[f64; 3], blending_matrix: &na::DMatrix<f64>) -> Self {
        RotationCost {
            u,
            rvec: na::Vector3::new(rvec[0], rvec[1], rvec[2]),
            blending_matrix: blending_matrix.clone(),
        }
    }
}

impl<T: na::RealField> tiny_solver::factors::Factor<T> for RotationCost {
    fn residual_func(&self, params: &[na::DVector<T>]) -> na::DVector<T> {
        let r = so3_from_u_and_knots(self.u, params, &self.blending_matrix);
        let target = self.rvec.cast::<T>().to_dvec();
        let target = SO3::exp(target.as_view());
        (target.inverse() * r).log()
    }
}

impl<const N: usize> SO3Bspline<N> {
    fn fit(&mut self, timestamps_ns: &[u64], rvecs: &[[f64; 3]]) {
        let mut problem = tiny_solver::Problem::new();
        let mut initial_values = HashMap::new();
        for (&t_ns, rvec) in timestamps_ns.iter().zip(rvecs) {
            let (u, idx) = self.get_u_and_index(t_ns);
            let cost = RotationCost::new(u, rvec, &self.blending_matrix());
            let mut var_list = Vec::new();
            for i in idx..idx + N {
                let var_name = format!("r{}", i);
                if !initial_values.contains_key(&var_name) {
                    problem.set_variable_manifold(&var_name, Arc::new(RvecManifold));
                    let rvec = self.knots[i].to_dvec();
                    initial_values.insert(var_name.clone(), rvec);
                }
                var_list.push(var_name);
            }
            let var_list: Vec<_> = var_list.iter().map(|a| a.as_str()).collect();
            problem.add_residual_block(
                3,
                &var_list,
                Box::new(cost),
                Some(Box::new(HuberLoss::new(0.1))),
            );
        }
        let optimizer = GaussNewtonOptimizer::default();
        let result = optimizer.optimize(&problem, &initial_values, None).unwrap();
        self.knots = (0..self.knots.len())
            .map(|i| {
                let var_name = format!("r{}", i);
                let knot = result.get(&var_name).unwrap();
                [knot[0], knot[1], knot[2]]
            })
            .collect();
    }
    pub fn blending_matrix(&self) -> na::DMatrix<f64> {
        self.blending_matrix.clone()
    }

    pub fn first_derivative_bases(&self) -> na::DVector<f64> {
        self.first_derivative_bases.clone()
    }

    pub fn get_u_and_index(&self, timestamp_ns: u64) -> (f64, usize) {
        // println!("req: {}", timestamp_ns);
        let time_offset = timestamp_ns - self.timestamp_start_ns;
        // println!("offset {}", time_offset);
        // println!("mod {}", time_offset % self.spacing_ns);
        let u = (time_offset % self.spacing_ns) as f64 / self.spacing_ns as f64;
        // println!("u {}", u);
        let idx = time_offset / self.spacing_ns;
        (u, idx as usize)
    }

    pub fn get_rotation(&self, timestamp_ns: u64) -> SO3<f64> {
        let (u, i) = self.get_u_and_index(timestamp_ns);
        let knots: Vec<_> = (i..i + N).map(|idx| self.knots[idx].to_dvec()).collect();
        so3_from_u_and_knots(u, &knots, &self.blending_matrix())
    }

    pub fn get_velocity(&self, timestamp_ns: u64) -> na::Vector3<f64> {
        let (u, idx) = self.get_u_and_index(timestamp_ns);

        let coeff = &self.blending_matrix * Self::base_coeffs_with_time::<0>(u);
        let d_tn_s = self.spacing_ns as f64 / 1e9;
        let dcoeff = 1.0 / d_tn_s * &self.blending_matrix * Self::base_coeffs_with_time::<1>(u);
        let mut w = na::Vector3::<f64>::zeros();
        for j in 0..(N - 1) {
            let p0 = self.knots[idx + j].to_so3();
            let p1 = self.knots[idx + j + 1].to_so3();
            let r01 = p0.inverse() * p1;
            let delta = r01.log();
            let ww = ((-1.0 * &delta) * coeff[j + 1]).to_so3();
            w = &ww * w.as_view();
            w += delta * dcoeff[j + 1];
        }
        w
    }

    pub fn from_rotation_vectors(
        timestamps_ns: &[u64],
        rvecs: &[[f64; 3]],
        spacing_ns: u64,
    ) -> Self {
        if timestamps_ns.len() < N {
            panic!("timestams_ns should be larger than {}", N);
        } else if timestamps_ns.len() != rvecs.len() {
            panic!("timestameps_ns should have the same length as rvecs");
        }
        let mut prev_t = timestamps_ns[0];
        for &t in timestamps_ns {
            if t - prev_t > spacing_ns {
                panic!("spacing should be larger than all time stamp steps");
            }
            prev_t = t;
        }
        let num_of_knots =
            ((timestamps_ns.last().unwrap() - timestamps_ns[0]) / spacing_ns) as usize + N;
        let timestamp_start_ns = timestamps_ns[0];

        let mut current_idx = 0;
        let mut knots = Vec::new();
        for i in 0..num_of_knots {
            let knot_time_ns = timestamp_start_ns + i as u64 * spacing_ns;
            while current_idx < timestamps_ns.len() - 1 {
                if timestamps_ns[current_idx + 1] < knot_time_ns {
                    current_idx += 1;
                } else {
                    break;
                }
            }
            knots.push(rvecs[current_idx]);
        }

        let mut bspline = SO3Bspline {
            timestamp_start_ns,
            spacing_ns,
            knots,
            blending_matrix: compute_blending_matrix(N, true),
            first_derivative_bases: compute_base_coefficients(N).row(1).transpose(),
        };
        bspline.fit(timestamps_ns, rvecs);
        bspline
    }
    fn base_coeffs_with_time<const DERIVATIVE: usize>(u: f64) -> na::DVector<f64> {
        let mut res = na::DVector::zeros(N);
        let base_coefficients = compute_base_coefficients(N);
        if DERIVATIVE < N {
            res[DERIVATIVE] = base_coefficients[(DERIVATIVE, DERIVATIVE)];
            let mut ti = u;
            for j in (DERIVATIVE + 1)..N {
                res[j] = base_coefficients[(DERIVATIVE, j)] * ti;
                ti *= u;
            }
        }
        res
    }
}
