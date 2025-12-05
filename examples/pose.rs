use std::collections::HashMap;

use b_spline::so3bspline::SO3Bspline;
use serde::{Deserialize, Serialize};
use tiny_solver::na;

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
    env_logger::init();
    let contents = std::fs::read_to_string("examples/poses.json")
        .expect("Should have been able to read the file");
    let poses: HashMap<u64, Pose> = serde_json::from_str(&contents).unwrap();
    let mut rvecs_with_t: Vec<_> = poses.iter().map(|(k, v)| (*k, v.rvec)).collect();
    rvecs_with_t.sort_by(|a, b| a.0.cmp(&b.0));
    let mut pose_vec: Vec<_> = poses.iter().map(|(k, v)| (*k, v)).collect();
    pose_vec.sort_by(|a, b| a.0.cmp(&b.0));
    let (ts_ns, rvecs): (Vec<_>, Vec<_>) = rvecs_with_t.iter().map(|a| (a.0, a.1)).unzip();

    let spacing_ns = 250_000_000;

    let rotation_bspline = SO3Bspline::<5>::from_rotation_vectors(&ts_ns, &rvecs, spacing_ns);
    println!("knots: {}", rotation_bspline.knots.len());
    let recording = rerun::RecordingStreamBuilder::new("pose").spawn().unwrap();

    let start_t = pose_vec.first().unwrap().0;
    let end_t = pose_vec.last().unwrap().0;
    let delta_t = 10_000_000;
    let mut t = start_t;
    while t <= end_t {
        recording.set_time(
            "stable",
            rerun::TimeCell::from_timestamp_nanos_since_epoch(t as i64),
        );
        let r = rotation_bspline.get_rotation(t);
        let rvec = r.log();
        recording
            .log("/rvec/x/bpline", &rerun::Scalars::single(rvec[0]))
            .unwrap();
        recording
            .log("/rvec/y/bpline", &rerun::Scalars::single(rvec[1]))
            .unwrap();
        recording
            .log("/rvec/z/bpline", &rerun::Scalars::single(rvec[2]))
            .unwrap();
        let gyro = rotation_bspline.get_velocity(t);
        recording
            .log("/gyro/x", &rerun::Scalars::single(gyro[0]))
            .unwrap();
        recording
            .log("/gyro/y", &rerun::Scalars::single(gyro[1]))
            .unwrap();
        recording
            .log("/gyro/z", &rerun::Scalars::single(gyro[2]))
            .unwrap();
        t += delta_t;
    }

    for p in &pose_vec {
        recording.set_time(
            "stable",
            rerun::TimeCell::from_timestamp_nanos_since_epoch(p.0 as i64),
        );
        recording
            .log(
                "pose",
                &na_isometry3_to_rerun_transform3d(&rtvec_to_transform(&p.1.rvec, &p.1.tvec)),
            )
            .unwrap();
        let r = rotation_bspline.get_rotation(p.0);
        let rvec = r.log();
        let t_shift_for_vis = na::Isometry3::from_parts(
            na::Translation3::new(1.5, 0.0, 0.0),
            na::UnitQuaternion::identity(),
        );
        let t_bspline = rtvec_to_transform(&[rvec[0], rvec[1], rvec[2]], &p.1.tvec);
        recording
            .log(
                "pose_spline",
                &na_isometry3_to_rerun_transform3d(&(t_bspline * t_shift_for_vis)),
            )
            .unwrap();
        recording
            .log("/rvec/x/lidar", &rerun::Scalars::single(p.1.rvec[0]))
            .unwrap();
        recording
            .log("/rvec/y/lidar", &rerun::Scalars::single(p.1.rvec[1]))
            .unwrap();
        recording
            .log("/rvec/z/lidar", &rerun::Scalars::single(p.1.rvec[2]))
            .unwrap();
    }
}
