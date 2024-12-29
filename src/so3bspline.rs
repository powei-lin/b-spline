use crate::common::compute_blending_matrix_c;
use nalgebra as na;

pub struct SO3Bspline<const N: usize> {
    timestamp_start_ns: u64,
    spacing_ns: u64,
    pub knots: Vec<[f64; 3]>,
    blending_matrix: [[f64; N]; N],
}

// def rmat_from_u_knots(u: float, *knots) -> np.ndarray:
//     knots = list(knots)
//     n = len(knots)
//     uv = np.array([u**i for i in range(len(knots))]).reshape(-1, 1)
//     m = computeBlendingMatrix(n, True)
//     kv = m @ uv
//     r = Rotation.from_rotvec(knots[0]).as_matrix()
//     for j in range(1, len(knots)):
//         k = kv[j, 0]
//         r_j = Rotation.from_rotvec(knots[j]).as_matrix()
//         r_j_minus_inv = Rotation.from_rotvec(knots[j - 1]).as_matrix().T
//         rj_vec = k * Rotation.from_matrix(r_j_minus_inv @ r_j).as_rotvec()
//         r = r @ Rotation.from_rotvec(rj_vec).as_matrix()
//     # print(r)
//     return r

pub fn so3_from_u_and_knots<T: na::RealField>(u: T, knots: &[[T; 3]]) {
    let n = knots.len();
    let uv = na::DVector::from_fn(n, |i, _| u.clone().powi(i as i32));
}

impl<const N: usize> SO3Bspline<N> {
    // pub fn new() -> Self {
    //     SO3Bspline {
    //         knots: Vec::new(),
    //         blending_matrix: compute_blending_matrix_c(true),
    //     }
    // }
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
            println!("{} {} {}", t, prev_t, t - prev_t);
            if t - prev_t > spacing_ns {
                panic!("spacing should be larger than all time stamp steps");
            }
            prev_t = t;
        }
        let num_of_knots =
            ((timestamps_ns.last().unwrap() - timestamps_ns[0]) / spacing_ns) as usize + N - 1;
        let timestamp_start_ns = timestamps_ns[0];
        println!("{}", num_of_knots);

        let mut current_idx = 0;
        let mut knots = Vec::new();
        for i in 0..num_of_knots {
            let knot_time_ns = timestamp_start_ns + i as u64 * spacing_ns;
            if current_idx == timestamps_ns.len() - 1 {
                knots.push(rvecs[current_idx]);
            } else {
                if timestamps_ns[current_idx + 1] < knot_time_ns {
                    current_idx += 1;
                }
                knots.push(rvecs[current_idx]);
            }
        }

        SO3Bspline {
            timestamp_start_ns,
            spacing_ns,
            knots,
            blending_matrix: compute_blending_matrix_c(true),
        }
    }
}
