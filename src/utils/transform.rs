use crate::types::*;
use na::{Matrix3, Matrix4, Vector3};

#[allow(dead_code)]
pub fn rotate_about_start(el: &FrameElement, axis: Axis, angle: f64) -> FrameElement {
    let axis_angle = match axis {
        Axis::X => Vector3::x(),
        Axis::Y => Vector3::y(),
        Axis::Z => Vector3::z(),
    } * angle;

    let start = na::Point3::new(el.ends.0.x, el.ends.0.y, el.ends.0.z);
    let rotation = Matrix4::new_rotation_wrt_point(axis_angle, start);
    let end_point = Vector3::new(el.ends.1.x, el.ends.1.y, el.ends.1.z).to_homogeneous();
    let new_end_point = rotation * end_point;

    FrameElement {
        id: el.id.clone(),
        ends: (
            el.ends.0.clone(),
            crate::types::Point3 {
                x: new_end_point.x,
                y: new_end_point.y,
                z: new_end_point.z,
            },
        ),
    }
}

#[allow(dead_code)]
pub fn local_to_global(local: &Matrix3<f64>, transform: &Matrix3<f64>) -> Matrix3<f64> {
    transform * local
}

#[allow(dead_code)]
pub fn global_to_local(global: &Matrix3<f64>, transform: &Matrix3<f64>) -> Matrix3<f64> {
    transform.transpose() * global
}

/// Returns the 3D transformation matrix to go from a local to global coordinate system.
///
/// # Arguments
///
/// * `local_axes` - A 3x3 matrix with each of the local x, y, and z axes as column vectors
///
/// # Examples
/// TODO Add example
/// ```
///
/// ```
#[allow(dead_code)]
fn local_to_global_rotation_matrix(local_axes: &Matrix3<f64>) -> Matrix3<f64> {
    let local_x_axis = local_axes.column(0);
    let local_y_axis = local_axes.column(1);
    let local_z_axis = local_axes.column(2);

    Matrix3::new(
        local_x_axis.angle(&Vector3::x_axis()).cos(),
        local_x_axis.angle(&Vector3::y_axis()).cos(),
        local_x_axis.angle(&Vector3::z_axis()).cos(),
        local_y_axis.angle(&Vector3::x_axis()).cos(),
        local_y_axis.angle(&Vector3::y_axis()).cos(),
        local_y_axis.angle(&Vector3::z_axis()).cos(),
        local_z_axis.angle(&Vector3::x_axis()).cos(),
        local_z_axis.angle(&Vector3::y_axis()).cos(),
        local_z_axis.angle(&Vector3::z_axis()).cos(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn local_to_global_rotation_does_not_alter_global_axes() {
        let local_axes = Matrix3::<f64>::identity();
        let transformed = local_to_global_rotation_matrix(&local_axes);

        let expected = Matrix3::new(1., 0., 0., 0., 1., 0., 0., 0., 1.);

        assert!(transformed.relative_eq(&expected, 1e-6, 1e-6));
    }
}
