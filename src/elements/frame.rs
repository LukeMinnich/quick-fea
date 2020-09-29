use crate::types::frame::{FrameElement, FrameEndRelease, FrameEndReleases};
use na::{Matrix2, Matrix4, MatrixN, Vector2, Vector4, U12};

#[rustfmt::skip]
#[allow(non_snake_case)]
pub fn euler_bernoulli_element_stiffness_matrix(element: &FrameElement) -> MatrixN<f64, U12> {

    let section = &element.geometry.cross_section;
    let material = &element.material;

    let L = element.length_or_inf();
    let L2 = L * L;

    let (A, Iy, Iz, J) = (section.A, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = MatrixN::<f64, U12>::zeros();

    // Axial
    m[(0, 0)] =  A;
    m[(0, 6)] = -A;
    m[(6, 0)] = -A;
    m[(6, 6)] =  A;

    // Torsion
    m[(3, 3)] =  G * J / E;
    m[(3, 9)] = -G * J / E;
    m[(9, 3)] = -G * J / E;
    m[(9, 9)] =  G * J / E;

    // Bending About Z
    m[(1, 1)] =  12. * Iz / L2;
    m[(1, 7)] = -12. * Iz / L2;
    m[(7, 1)] = -12. * Iz / L2;
    m[(7, 7)] =  12. * Iz / L2;

    m[(1, 5)]  =  6. * Iz / L;
    m[(1, 11)] =  6. * Iz / L;
    m[(7, 5)]  = -6. * Iz / L;
    m[(7, 11)] = -6. * Iz / L;

    m[(5, 1)]  =  6. * Iz / L;
    m[(5, 7)]  = -6. * Iz / L;
    m[(11, 1)] =  6. * Iz / L;
    m[(11, 7)] = -6. * Iz / L;
    
    m[(5, 5)]   = 4. * Iz;
    m[(5, 11)]  = 2. * Iz;
    m[(11, 5)]  = 2. * Iz;
    m[(11, 11)] = 4. * Iz;

    // Bending About Y
    m[(2, 2)] =  12. * Iy / L2;
    m[(2, 8)] = -12. * Iy / L2;
    m[(8, 2)] = -12. * Iy / L2;
    m[(8, 8)] =  12. * Iy / L2;

    m[(2, 4)]  = -6. * Iy / L;
    m[(2, 10)] = -6. * Iy / L;
    m[(8, 4)]  =  6. * Iy / L;
    m[(8, 10)] =  6. * Iy / L;
    
    m[(4, 2)]  = -6. * Iy / L;
    m[(4, 8)]  =  6. * Iy / L;
    m[(10, 2)] = -6. * Iy / L;
    m[(10, 8)] =  6. * Iy / L;

    m[(4, 4)]   =  4. * Iy;
    m[(4, 10)]  =  2. * Iy;
    m[(10, 4)]  =  2. * Iy;
    m[(10, 10)] =  4. * Iy;

    apply_end_releases(&mut m, element);

    // Result must include matrix coefficients
    E / L * m
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
