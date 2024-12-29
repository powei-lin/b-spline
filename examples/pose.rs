use std::collections::HashMap;

use b_spline::so3bspline::SO3Bspline;
use nalgebra as na;
use serde::{Deserialize, Serialize};

fn rtvec_to_transform(rvec: &[f64; 3], tvec: &[f64; 3]) -> na::Isometry3<f64> {
    na::Isometry3::new(
        na::Vector3::new(tvec[0], tvec[1], tvec[2]),
        na::Vector3::new(rvec[0], rvec[1], rvec[2]),
    )
}

pub fn na_isometry3_to_rerun_transform3d(transform: &na::Isometry3<f64>) -> rerun::Transform3D {
    let t = (
        transform.translation.x as f32,
        transform.translation.y as f32,
        transform.translation.z as f32,
    );
    let q_xyzw = (
        transform.rotation.quaternion().i as f32,
        transform.rotation.quaternion().j as f32,
        transform.rotation.quaternion().k as f32,
        transform.rotation.quaternion().w as f32,
    );
    rerun::Transform3D::from_translation_rotation(t, rerun::Quaternion::from_xyzw(q_xyzw.into()))
}

#[derive(Debug, Serialize, Deserialize)]
struct Pose {
    rvec: [f64; 3],
    tvec: [f64; 3],
}
fn main() {
    let recording = rerun::RecordingStreamBuilder::new("pose").spawn().unwrap();
    let contents = std::fs::read_to_string("examples/poses.json")
        .expect("Should have been able to read the file");
    let poses: HashMap<u64, Pose> = serde_json::from_str(&contents).unwrap();
    let mut rvecs_with_t: Vec<_> = poses.iter().map(|(k, v)| (*k, v.rvec)).collect();
    rvecs_with_t.sort_by(|a, b| a.0.cmp(&b.0));
    let (ts_ns, rvecs): (Vec<_>, Vec<_>) = rvecs_with_t.iter().map(|a| (a.0, a.1)).unzip();

    let spacing_ns = 250_000_000;

    let rotation_bspline = SO3Bspline::<5>::from_rotation_vectors(&ts_ns, &rvecs, spacing_ns);

    for p in poses {
        println!("{} {:?} {:?}", p.0, p.1.rvec, p.1.tvec);
        recording.set_time_nanos("stable", p.0 as i64);
        recording
            .log(
                "pose",
                &na_isometry3_to_rerun_transform3d(&rtvec_to_transform(&p.1.rvec, &p.1.tvec)),
            )
            .unwrap();
    }
}
