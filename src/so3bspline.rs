use crate::common::compute_blending_matrix_c;

pub struct SO3Bspline<const N: usize> {
    timestamp_start: u64,
    knots: Vec<[f64; 3]>,
    blending_matrix: [[f64; N]; N],
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
        println!("{}", num_of_knots);
        SO3Bspline {
            timestamp_start: timestamps_ns[0],
            knots: Vec::new(),
            blending_matrix: compute_blending_matrix_c(true),
        }
    }
}
