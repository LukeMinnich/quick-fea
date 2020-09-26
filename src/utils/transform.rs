use na::{U12, Matrix, ArrayStorage, Matrix3, Vector3};

/**
 * 
 */
// pub fn frame_element_stiff_matrix(element: )

/**
 * Returns the transformation matrix __Γ__ to go from world to local vectors (force, displacement, stiffness)
 * 
 * The 12x12 world-to-local transformation matrix takes the form of
 * ```
 * γ 0 0 0
 * 0 γ 0 0
 * 0 0 γ 0
 * 0 0 0 γ
 * ```
 * where __γ__ is the 3x3 world-to-local rotation transform matrix
 * and __0__ is the 3x3 zero matrix
 * 
 * # Arguments
 * 
 * `local_axes` - A 3x3 matrix with each of the member local x, y, and z axes as column vectors
 * 
 */
pub fn world_to_local_force_transform(local_axes: &Matrix3<f64>) -> Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>> {

    let rotation: Matrix3<f64> = world_to_local_rotation(&local_axes);

    let mut transform = Matrix::<f64, U12, U12, ArrayStorage<f64, U12, U12>>::zeros();

    for i in 0..12 {
        if i < 3 {
            for j in 0..3 {
                transform[(i, j)] = rotation[(i, j)];
            }
        } else if i < 6 {
            for j in 3..6 {
                transform[(i, j)] = rotation[(i-3, j-3)];
            }
        } else if i < 9 {
            for j in 6..9 {
                transform[(i, j)] = rotation[(i-6, j-6)];
            }
        } else {
            for j in 9..12 {
                transform[(i, j)] = rotation[(i-9, j-9)];
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

#[cfg(test)]
mod tests {
    use super::*;

    mod matrix_structural_analysis_mcguire_2nd {
        use super::*;
        use na::{Matrix6};

        #[test]
        fn example_5_3() {
            let member_ab_local_stiffness = 200. * Matrix6::new(0.8, 0., 0., -0.8, 0., 0.,
                                                               0., 0.0048, 12., 0.,  -0.0048, 12.,
                                                               0., 12., 0.4e5, 0., -12., 0.2e5,
                                                               -0.8, 0., 0., 0.8, 0., 0.,
                                                               0., -0.0048, -12., 0., 0.0048, -12.,
                                                               0., 12., 0.2e5, 0., -12., 0.4e5);
            let member_ab_global_target = 200. * Matrix6::new(0.0048, 0., -12., -0.0048, 0., -12.,
                                                              0., 0.8, 0., 0., -0.8, 0.,
                                                              -12., 0., 0.4e5, 12., 0., 0.2e5,
                                                              -0.0048, 0., 122., 0.0048, 0., 12.,
                                                              0., -0.8, 0., 0., 0.8, 0.,
                                                              -12., 0., 0.2e5, 12., 0., 0.4e5);


            let member_ab_axes = Matrix3::new( 0., 1., 0.,
                                              -1., 0., 0.,
                                               0., 0., 1.);

            assert_abs_diff_eq!(1.0, 1.0, epsilon = 1e-6);
            assert_abs_diff_ne!(1.00000, 1.00001, epsilon = 1e-6);

            let stiffness_transform: Matrix6<f64> = local_to_world_stiffness_transform(&member_ab_axes)
        }
    }
}
