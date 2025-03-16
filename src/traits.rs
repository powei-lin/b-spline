use tiny_solver::manifold::so3::SO3;
use tiny_solver::na;

pub trait ToSO3<T: na::RealField> {
    fn to_rotation3(&self) -> na::Rotation3<T>;
    fn to_so3(&self) -> tiny_solver::manifold::so3::SO3<T>;
}

impl<T: na::RealField> ToSO3<T> for na::DVector<T> {
    fn to_rotation3(&self) -> na::Rotation3<T> {
        na::Rotation3::from_scaled_axis(na::Vector3::new(
            self[0].clone(),
            self[1].clone(),
            self[2].clone(),
        ))
    }
    fn to_so3(&self) -> tiny_solver::manifold::so3::SO3<T> {
        SO3::exp(self.as_view())
    }
}

impl<T: na::RealField> ToSO3<T> for [T; 3] {
    fn to_rotation3(&self) -> na::Rotation3<T> {
        na::Rotation3::from_scaled_axis(na::Vector3::new(
            self[0].clone(),
            self[1].clone(),
            self[2].clone(),
        ))
    }
    fn to_so3(&self) -> tiny_solver::manifold::so3::SO3<T> {
        SO3::exp(self.to_dvec().as_view())
    }
}
impl<T: na::RealField> ToSO3<T> for na::Vector3<T> {
    fn to_rotation3(&self) -> na::Rotation3<T> {
        na::Rotation3::from_scaled_axis(na::Vector3::new(
            self.x.clone(),
            self.y.clone(),
            self.z.clone(),
        ))
    }
    fn to_so3(&self) -> tiny_solver::manifold::so3::SO3<T> {
        SO3::exp(self.to_dvec().as_view())
    }
}

impl<T: na::RealField> ToDVec<T> for [T; 3] {
    fn to_dvec(&self) -> na::DVector<T> {
        na::DVector::<_>::from_vec_storage(na::VecStorage::new(
            na::Dyn(3usize),
            na::Const::<1>,
            vec![self[0].clone(), self[1].clone(), self[2].clone()],
        ))
    }
}

impl<T: na::RealField> ToDVec<T> for na::Vector3<T> {
    fn to_dvec(&self) -> na::DVector<T> {
        na::DVector::<_>::from_vec_storage(na::VecStorage::new(
            na::Dyn(3usize),
            na::Const::<1>,
            vec![self[0].clone(), self[1].clone(), self[2].clone()],
        ))
    }
}

pub trait ToDVec<T> {
    fn to_dvec(&self) -> na::DVector<T>;
}
