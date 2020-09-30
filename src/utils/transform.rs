use na::{Matrix, Matrix3, MatrixN, Vector3, U12};

// pub fn frame_element_stiff_matrix(element: )

/**
 * Returns the transformation matrix __Γ__ to go from world to local coordinates (force, displacement, stiffness)
 *
 * The 12x12 world-to-local transformation matrix takes the form of
 *
 * γ 0 0 0
 * 0 γ 0 0
 * 0 0 γ 0
 * 0 0 0 γ
 *
 * where __γ__ is the 3x3 world-to-local rotation transform matrix
 * and __0__ is the 3x3 zero matrix
 *
 * # Arguments
 *
 * `local_axes` - A 3x3 matrix with each of the member local x, y, and z axes as column vectors
 *
 */
pub fn world_to_local_transform(local_axes: &Matrix3<f64>) -> MatrixN<f64, U12> {
    let rotation: Matrix3<f64> = world_to_local_rotation(&local_axes);

    let mut transform = MatrixN::<f64, U12>::zeros();

    for i in 0..12 {
        if i < 3 {
            for j in 0..3 {
                transform[(i, j)] = rotation[(i, j)];
            }
        } else if i < 6 {
            for j in 3..6 {
                transform[(i, j)] = rotation[(i - 3, j - 3)];
            }
        } else if i < 9 {
            for j in 6..9 {
                transform[(i, j)] = rotation[(i - 6, j - 6)];
            }
        } else {
            for j in 9..12 {
                transform[(i, j)] = rotation[(i - 9, j - 9)];
            }
        }
    }

    transform
}

/**
 * Returns the rotation matrix __γ__ to go from world to local coordinate system
 * 
 * # Arguments
 * 
 * `local_axes` - A 3x3 matrix with each of the member local x, y, and z axes as column vectors.
 * They need not be basis vectors.
 * 
 */
#[rustfmt::skip]
fn world_to_local_rotation(local_axes: &Matrix3<f64>) -> Matrix3<f64> {

    // Ensure the  these are basis (unit) vectors
    let local_x = Matrix::normalize(&local_axes.column(0));
    let local_y = Matrix::normalize(&local_axes.column(1));
    let local_z = Matrix::normalize(&local_axes.column(2));

    let x = Vector3::<f64>::x_axis();
    let y = Vector3::<f64>::y_axis();
    let z = Vector3::<f64>::z_axis();

    Matrix3::new(local_x.dot(&x), local_x.dot(&y), local_x.dot(&z),
                 local_y.dot(&x), local_y.dot(&y), local_y.dot(&z),
                 local_z.dot(&x), local_z.dot(&y), local_z.dot(&z))
}
