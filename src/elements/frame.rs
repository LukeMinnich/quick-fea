use crate::types::frame::{FrameElement, FrameEndRelease, FrameEndReleases};
use crate::utils::transform::world_to_local_transform;
use na::*;

pub fn transform_frame_local_to_world(
    frame: &FrameElement,
    m: &MatrixN<f64, U12>,
) -> MatrixN<f64, U12> {
    let transform = world_to_local_transform(&frame.geometry.local_axes);
    transform.transpose() * m
}
pub fn transform_frame_world_to_local(
    frame: &FrameElement,
    m: &MatrixN<f64, U12>,
) -> MatrixN<f64, U12> {
    let transform = world_to_local_transform(&frame.geometry.local_axes);
    transform * m
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub fn frame_element_stiffness_matrix(element: &FrameElement) -> MatrixN<f64, U12> {

    let section = &element.geometry.cross_section;
    let material = &element.material;

    let L = element.length_or_inf();
    let L2 = L * L;

    let (A, Iy, Iz, J) = (section.A, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = MatrixN::<f64, U12>::zeros();

    // Axial
    let index_map = Vector2::new(0, 6);
    let stiffness_matrix = E * A / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_2by2(&axial_action, &mut m);

    // Torsion
    let index_map = Vector2::new(3, 9);
    let stiffness_matrix = G * J / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_2by2(&axial_action, &mut m);

    // Bending About Z
    let index_map = Vector4::new(1, 5, 7, 11);
    let stiffness_matrix = E * Iz / L * Matrix4::new( 12. / L2,  6. / L, -12. / L2,  6. / L,
                                                       6. / L ,  4.    ,  -6. / L ,  2.    ,
                                                     -12. / L2, -6. / L,  12. / L2, -6. / L,
                                                       6. / L ,  2.    ,  -6. / L ,  4.    );
    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_4by4(&axial_action, &mut m);


    // Bending About Y
    let index_map = Vector4::new(2, 4, 8, 10);
    let stiffness_matrix = E * Iy / L * Matrix4::new( 12. / L2, -6. / L, -12. / L2,  -6. / L,
                                                      -6. / L ,  4.    ,   6. / L ,   2.    ,
                                                     -12. / L2,  6. / L,  12. / L2,   6. / L,
                                                      -6. / L ,  2.    ,   6. / L ,   4.    );
    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_4by4(&axial_action, &mut m);

    apply_end_releases(&mut m, element);    
    
    m
}

#[rustfmt::skip]
#[allow(non_snake_case)]
pub fn frame_element_with_shear_deformation_stiffness_matrix(element: &FrameElement) -> MatrixN<f64, U12> {

    let section = &element.geometry.cross_section;
    let material = &element.material;

    let L = element.length_or_inf();
    let L2 = L * L;

    let (A, Avy, Avz, Iy, Iz, J) = (section.A, section.Avy, section.Avz, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = MatrixN::<f64, U12>::zeros();

    // Axial
    let index_map = Vector2::new(0, 6);
    let stiffness_matrix = E * A / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_2by2(&axial_action, &mut m);

    // Torsion
    let index_map = Vector2::new(3, 9);
    let stiffness_matrix = G * J / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_2by2(&axial_action, &mut m);

    // Bending About Z
    let index_map = Vector4::new(7, 11, 1, 5);
    let eta = E * Iz / Avy / G;
    let stiffness_matrix = E * Iz / L / (L2 / 12. + eta) * Matrix4::new(     1.,       -L / 2.,     -1.,       -L / 2.,
                                                                        -L / 2., L2 / 3. + eta,  L / 2., L2 / 6. - eta,
                                                                            -1.,        L / 2.,      1.,        L / 2.,
                                                                        -L / 2., L2 / 6. - eta,  L / 2., L2 / 3. + eta);

    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_4by4(&axial_action, &mut m);


    // Bending About Y
    let index_map = Vector4::new(8, 10, 2, 4);
    let eta = E * Iy / Avz / G;
    let stiffness_matrix = E * Iy / L / (L2 / 12. + eta) * Matrix4::new(     1.,        L / 2.,     -1.,        L / 2.,
                                                                         L / 2., L2 / 3. + eta, -L / 2., L2 / 6. - eta,
                                                                            -1.,       -L / 2.,      1.,       -L / 2.,
                                                                         L / 2., L2 / 6. - eta, -L / 2., L2 / 3. + eta);

    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_4by4(&axial_action, &mut m);
    
    apply_end_releases(&mut m, element);
    
    m
}

/// Applies end releases if they occur
fn apply_end_releases(m: &mut MatrixN<f64, U12>, frame: &FrameElement) {
    let start: &FrameEndReleases = &(frame.start_releases);
    let end: &FrameEndReleases = &(frame.end_releases);

    if start.A == FrameEndRelease::Free || end.A == FrameEndRelease::Free {
        m[(0, 0)] = 0.;
        m[(6, 0)] = 0.;
        m[(0, 6)] = 0.;
        m[(6, 6)] = 0.;
    }

    if start.T == FrameEndRelease::Free || end.T == FrameEndRelease::Free {
        m[(3, 3)] = 0.;
        m[(3, 9)] = 0.;
        m[(9, 3)] = 0.;
        m[(9, 9)] = 0.;
    }

    // TODO Figure out this logic... there's moment and shear interaction.
    //      If it can't transmit shear through both ends, isn't that also effectively a release?
    if (start.Mz == FrameEndRelease::Free && end.Mz == FrameEndRelease::Free)
        || (start.Mz == FrameEndRelease::Free && start.Vy == FrameEndRelease::Free)
        || (end.Mz == FrameEndRelease::Free && end.Vy == FrameEndRelease::Free)
    {
        m[(1, 1)] = 0.;
        m[(1, 5)] = 0.;
        m[(1, 7)] = 0.;
        m[(1, 11)] = 0.;

        m[(5, 1)] = 0.;
        m[(5, 5)] = 0.;
        m[(5, 7)] = 0.;
        m[(5, 11)] = 0.;

        m[(7, 1)] = 0.;
        m[(7, 5)] = 0.;
        m[(7, 7)] = 0.;
        m[(7, 11)] = 0.;

        m[(11, 1)] = 0.;
        m[(11, 5)] = 0.;
        m[(11, 7)] = 0.;
        m[(11, 11)] = 0.;
    } else if start.Mz == FrameEndRelease::Free {
        unimplemented!();
    } else if end.Mz == FrameEndRelease::Free {
        unimplemented!();
    } else if start.Vy == FrameEndRelease::Free {
        unimplemented!();
    } else if end.Vy == FrameEndRelease::Free {
        unimplemented!();
    }

    // TODO Figure out this logic... there's moment and shear interaction.
    //      If it can't transmit shear through both ends, isn't that also effectively a release?
    if (start.My == FrameEndRelease::Free && end.My == FrameEndRelease::Free)
        || (start.My == FrameEndRelease::Free && start.Vz == FrameEndRelease::Free)
        || (end.My == FrameEndRelease::Free && end.Vz == FrameEndRelease::Free)
    {
        m[(2, 2)] = 0.;
        m[(2, 8)] = 0.;
        m[(8, 2)] = 0.;
        m[(8, 8)] = 0.;

        m[(2, 4)] = 0.;
        m[(2, 10)] = 0.;
        m[(8, 4)] = 0.;
        m[(8, 10)] = 0.;

        m[(4, 2)] = 0.;
        m[(4, 8)] = 0.;
        m[(10, 2)] = 0.;
        m[(10, 8)] = 0.;

        m[(4, 4)] = 0.;
        m[(4, 10)] = 0.;
        m[(10, 4)] = 0.;
        m[(10, 10)] = 0.;
    } else if start.My == FrameEndRelease::Free {
        unimplemented!();
    } else if end.My == FrameEndRelease::Free {
        unimplemented!();
    } else if start.Vz == FrameEndRelease::Free {
        unimplemented!();
    } else if end.Vz == FrameEndRelease::Free {
        unimplemented!();
    }
}

fn merge_single_action_into_complete_2by2(
    single: &SingleActionStiffnessMatrix2by2,
    complete: &mut MatrixN<f64, U12>,
) {
    let map = single.index_map;
    let source = single.stiffness_matrix;
    for i in 0..2 {
        for j in 0..2 {
            let target_i = map[i];
            let target_j = map[j];
            complete[(target_i, target_j)] = source[(i, j)];
        }
    }
}

fn merge_single_action_into_complete_4by4(
    single: &SingleActionStiffnessMatrix4by4,
    complete: &mut MatrixN<f64, U12>,
) {
    let map = single.index_map;
    let source = single.stiffness_matrix;
    for i in 0..4 {
        for j in 0..4 {
            let target_i = map[i];
            let target_j = map[j];
            complete[(target_i, target_j)] = source[(i, j)];
        }
    }
}

struct SingleActionStiffnessMatrix2by2 {
    index_map: Vector2<usize>,
    stiffness_matrix: Matrix2<f64>,
}

struct SingleActionStiffnessMatrix4by4 {
    index_map: Vector4<usize>,
    stiffness_matrix: Matrix4<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::frame::*;
    use crate::types::material::*;
    use crate::types::node::*;
    use crate::utils::ZERO_EPSILON;
    use crate::*;

    fn assert_matrix_symmetric_12by12(m: &MatrixN<f64, U12>) {
        for i in 0..12 {
            for j in 0..12 {
                assert_relative_eq!(m[(i, j)], m[(j, i)], max_relative = 1e-12);
            }
        }
    }

    fn assert_matrix_symmetric_18by18(m: &MatrixN<f64, U18>) {
        for i in 0..18 {
            for j in 0..18 {
                assert_relative_eq!(m[(i, j)], m[(j, i)], max_relative = 1e-12);
            }
        }
    }

    mod mcguire_matrix_structural_analysis_2nd_edition {
        use super::*;

        pub fn example_4_8_init() {
            let material = IsotropicMaterial::new(200., 0.3);
            add_node(Node {
                id: "a".to_string(),
                degrees_of_freedom: Vector6::from_iterator(0..6),
                coordinate: Point3::new(0., 0., 0.),
            });
            add_node(Node {
                id: "b".to_string(),
                degrees_of_freedom: Vector6::from_iterator(6..12),
                coordinate: Point3::new(8e3, 0., 0.),
            });
            add_node(Node {
                id: "c".to_string(),
                degrees_of_freedom: Vector6::from_iterator(12..18),
                coordinate: Point3::new(13e3, 0., 0.),
            });
            add_frame_element(FrameElement {
                id: "ab".to_string(),
                start_node_id: "a".to_string(),
                end_node_id: "b".to_string(),
                start_releases: FrameEndReleases::fully_fixed(),
                end_releases: FrameEndReleases::fully_fixed(),
                geometry: FrameGeometry {
                    cross_section: CrossSection {
                        A: 6e3,
                        Avy: 0.,
                        Avz: 0.,
                        J: 300e3,
                        Iy: 0.,
                        Iz: 200e6,
                    },
                    local_axes: Matrix3::identity(),
                },
                material: material.clone(),
            });
            add_frame_element(FrameElement {
                id: "bc".to_string(),
                start_node_id: "b".to_string(),
                end_node_id: "c".to_string(),
                start_releases: FrameEndReleases::fully_fixed(),
                end_releases: FrameEndReleases::fully_fixed(),
                geometry: FrameGeometry {
                    cross_section: CrossSection {
                        A: 4e3,
                        Avy: 0.,
                        Avz: 0.,
                        J: 100e3,
                        Iy: 0.,
                        Iz: 50e6,
                    },
                    local_axes: Matrix3::identity(),
                },
                material: material.clone(),
            });
        }

        #[test]
        pub fn example_4_8_part_1() {
            example_4_8_init();

            let member_ab = match get_frame_element_by_id("ab") {
                Some(x) => x,
                None => panic!(),
            };
            let member_bc = match get_frame_element_by_id("bc") {
                Some(x) => x,
                None => panic!(),
            };
            let local_ab = frame_element_stiffness_matrix(&member_ab) / 200.;
            let local_bc = frame_element_stiffness_matrix(&member_bc) / 200.;

            // Member ab

            // diagonals
            assert_relative_eq!(local_ab[(0, 0)], 0.750, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(1, 1)], 0.00469, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(3, 3)], 14.423, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(5, 5)], 1e5, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(6, 6)], 0.750, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(7, 7)], 0.00469, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(9, 9)], 14.423, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(11, 11)], 1e5, max_relative = 1e-3);

            // off-diagonal non-zeros
            assert_relative_eq!(local_ab[(0, 6)], -0.75, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(1, 5)], 18.75, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(1, 7)], -0.00469, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(1, 5)], 18.75, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(3, 9)], -14.423, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(5, 7)], -18.75, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(5, 11)], 0.5e5, max_relative = 1e-3);
            assert_relative_eq!(local_ab[(7, 11)], -18.75, max_relative = 1e-3);
            // row 0 zeros
            assert_abs_diff_eq!(local_ab[(0, 1)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(0, 3)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(0, 5)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(0, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(0, 9)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(0, 11)], 0., epsilon = ZERO_EPSILON);

            // row 1 zeros
            assert_abs_diff_eq!(local_ab[(1, 3)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(1, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(1, 9)], 0., epsilon = ZERO_EPSILON);

            // row 3 zeros
            assert_abs_diff_eq!(local_ab[(3, 5)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(3, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(3, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(3, 11)], 0., epsilon = ZERO_EPSILON);

            // row 5 zeros
            assert_abs_diff_eq!(local_ab[(5, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(5, 9)], 0., epsilon = ZERO_EPSILON);
            // row 6 zeros
            assert_abs_diff_eq!(local_ab[(6, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(6, 9)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_ab[(6, 11)], 0., epsilon = ZERO_EPSILON);

            // row 7 zeros
            assert_abs_diff_eq!(local_ab[(7, 9)], 0., epsilon = ZERO_EPSILON);

            // row 9 zeros
            assert_abs_diff_eq!(local_ab[(9, 11)], 0., epsilon = ZERO_EPSILON);

            assert_matrix_symmetric_12by12(&local_ab);

            // Member bc

            // diagonals
            assert_relative_eq!(local_bc[(0, 0)], 0.8, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(1, 1)], 0.0048, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(3, 3)], 7.692, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(5, 5)], 0.4e5, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(6, 6)], 0.8, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(7, 7)], 0.0048, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(9, 9)], 7.692, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(11, 11)], 0.4e5, max_relative = 1e-3);

            // off-diagonal non-zeros
            assert_relative_eq!(local_bc[(0, 6)], -0.8, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(1, 5)], 12., max_relative = 1e-3);
            assert_relative_eq!(local_bc[(1, 7)], -0.0048, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(1, 5)], 12., max_relative = 1e-3);
            assert_relative_eq!(local_bc[(3, 9)], -7.692, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(5, 7)], -12., max_relative = 1e-3);
            assert_relative_eq!(local_bc[(5, 11)], 0.2e5, max_relative = 1e-3);
            assert_relative_eq!(local_bc[(7, 11)], -12., max_relative = 1e-3);

            // row 0 zeros
            assert_abs_diff_eq!(local_bc[(0, 1)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(0, 3)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(0, 5)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(0, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(0, 9)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(0, 11)], 0., epsilon = ZERO_EPSILON);

            // row 1 zeros
            assert_abs_diff_eq!(local_bc[(1, 3)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(1, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(1, 9)], 0., epsilon = ZERO_EPSILON);

            // row 3 zeros
            assert_abs_diff_eq!(local_bc[(3, 5)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(3, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(3, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(3, 11)], 0., epsilon = ZERO_EPSILON);

            // row 5 zeros
            assert_abs_diff_eq!(local_bc[(5, 6)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(5, 9)], 0., epsilon = ZERO_EPSILON);
            // row 6 zeros
            assert_abs_diff_eq!(local_bc[(6, 7)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(6, 9)], 0., epsilon = ZERO_EPSILON);
            assert_abs_diff_eq!(local_bc[(6, 11)], 0., epsilon = ZERO_EPSILON);

            // row 7 zeros
            assert_abs_diff_eq!(local_bc[(7, 9)], 0., epsilon = ZERO_EPSILON);

            // row 9 zeros
            assert_abs_diff_eq!(local_bc[(9, 11)], 0., epsilon = ZERO_EPSILON);

            assert_matrix_symmetric_12by12(&local_bc);
        }

        #[test]
        pub fn example_4_8_part_2() {
            example_4_8_init();
            let member_ab = match get_frame_element_by_id("ab") {
                Some(x) => x,
                None => panic!(),
            };
            let member_bc = match get_frame_element_by_id("bc") {
                Some(x) => x,
                None => panic!(),
            };

            let local_ab = frame_element_stiffness_matrix(&member_ab) / 200.;
            let local_bc = frame_element_stiffness_matrix(&member_bc) / 200.;
            let world_ab = transform_frame_local_to_world(&member_ab, &local_ab);
            let world_bc = transform_frame_local_to_world(&member_bc, &local_bc);
            update_frame_element_stiffness(&member_ab, local_ab, world_ab);
            update_frame_element_stiffness(&member_bc, local_bc, world_bc);

            if let Ok(world) = analysis::assemble_world_stiffness_matrix() {
                assert_eq!(world.len(), 42); // The answer to life, the universe, and everything!

                let mut world_matrix = MatrixN::<f64, U18>::zeros();
                for (key, val) in world.iter() {
                    world_matrix[(key.0, key.1)] = *val;
                }

                // diagonals
                assert_relative_eq!(world_matrix[(0, 0)], 0.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(1, 1)], 0.00469, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(3, 3)], 14.423, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(5, 5)], 1e5, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(6, 6)], 1.55, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(7, 7)], 0.00949, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(9, 9)], 22.115, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(11, 11)], 1.4e5, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(12, 12)], 0.8, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(13, 13)], 0.0048, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(15, 15)], 7.692, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(17, 17)], 0.4e5, max_relative = 1e-3);
                // off-diagonal non-zeros
                assert_relative_eq!(world_matrix[(0, 6)], -0.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(1, 5)], 18.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(1, 7)], -0.00469, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(1, 11)], 18.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(3, 9)], -14.423, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(5, 7)], -18.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(5, 11)], 0.5e5, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(6, 12)], -0.8, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(7, 11)], -6.75, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(7, 13)], -0.0048, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(7, 17)], 12., max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(9, 15)], -7.692, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(11, 13)], -12., max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(11, 17)], 0.2e5, max_relative = 1e-3);
                assert_relative_eq!(world_matrix[(13, 17)], -12., max_relative = 1e-3);

                // row 0 zeros
                assert_abs_diff_eq!(world_matrix[(0, 1)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 3)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 5)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 7)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 9)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 11)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(0, 17)], 0., epsilon = ZERO_EPSILON);

                // row 1 zeros
                assert_abs_diff_eq!(world_matrix[(1, 3)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 6)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 9)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(1, 17)], 0., epsilon = ZERO_EPSILON);

                // row 3 zeros
                assert_abs_diff_eq!(world_matrix[(3, 5)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 6)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 7)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 11)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(3, 17)], 0., epsilon = ZERO_EPSILON);

                // row 5 zeros
                assert_abs_diff_eq!(world_matrix[(5, 6)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(5, 9)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(5, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(5, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(5, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(5, 17)], 0., epsilon = ZERO_EPSILON);

                // row 6 zeros
                assert_abs_diff_eq!(world_matrix[(6, 7)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(6, 9)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(6, 11)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(6, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(6, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(6, 17)], 0., epsilon = ZERO_EPSILON);

                // row 7 zeros
                assert_abs_diff_eq!(world_matrix[(7, 9)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(7, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(7, 15)], 0., epsilon = ZERO_EPSILON);

                // row 9 zeros
                assert_abs_diff_eq!(world_matrix[(9, 11)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(9, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(9, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(9, 17)], 0., epsilon = ZERO_EPSILON);

                // row 11 zeros
                assert_abs_diff_eq!(world_matrix[(11, 12)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(11, 15)], 0., epsilon = ZERO_EPSILON);

                // row 12 zeros
                assert_abs_diff_eq!(world_matrix[(12, 13)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(12, 15)], 0., epsilon = ZERO_EPSILON);
                assert_abs_diff_eq!(world_matrix[(12, 17)], 0., epsilon = ZERO_EPSILON);

                // row 13 zeros
                assert_abs_diff_eq!(world_matrix[(13, 15)], 0., epsilon = ZERO_EPSILON);

                // row 15 zeros
                assert_abs_diff_eq!(world_matrix[(15, 17)], 0., epsilon = ZERO_EPSILON);

                assert_matrix_symmetric_18by18(&world_matrix);
            }
        }
    }
}
