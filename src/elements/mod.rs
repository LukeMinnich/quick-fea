use crate::types::*;
use na::{ArrayStorage, Matrix, Matrix2, Matrix4, Vector2, Vector4, U12};

#[rustfmt::skip]
pub fn euler_bernoulli_element_stiffness_matrix(element: &FrameElement) -> Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>>{

    let section = element.geometry.cross_section;
    let material = element.material;

    let L = element.geometry.length;
    let L2 = L * L;
    let L3 = L * L * L;

    let (A, Iy, Iz, J) = (section.A, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = Matrix::<f64, U12, U12, ArrayStorage<f64, U12, U12>>::zeros();

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

    // Result must include matrix coefficients
    E / L * m
}

#[rustfmt::skip]
pub fn frame_element_stiffness_matrix(
    element: &FrameElement
) -> Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>> {

    let section = element.geometry.cross_section;
    let material = element.material;

    let L = element.geometry.length;
    let L2 = L * L;

    let (A, Iy, Iz, J) = (section.A, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = Matrix::<f64, U12, U12, ArrayStorage<f64, U12, U12>>::zeros();

    // Axial
    let mut index_map = Vector2::new(0, 6);
    let stiffness_matrix = E * A / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U2(&axial_action, &mut m);

    // Torsion
    let mut index_map = Vector2::new(3, 9);
    let stiffness_matrix = G * J / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U2(&axial_action, &mut m);

    // Bending About Z
    let mut index_map = Vector4::new(1, 5, 7, 11);
    let stiffness_matrix = E * Iz / L * Matrix4::new( 12. / L2,  6. / L, -12. / L2,  6. / L,
                                                       6. / L ,  4.    ,  -6. / L ,  2.    ,
                                                     -12. / L2, -6. / L,  12. / L2, -6. / L,
                                                       6. / L ,  2.    ,  -6. / L ,  4.    );
    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U4(&axial_action, &mut m);


    // Bending About Y
    let mut index_map = Vector4::new(2, 4, 8, 10);
    let stiffness_matrix = E * Iy / L * Matrix4::new( 12. / L2, -6. / L, -12. / L2,  -6. / L,
                                                      -6. / L ,  4.    ,   6. / L ,   2.    ,
                                                     -12. / L2,  6. / L,  12. / L2,   6. / L,
                                                      -6. / L ,  2.    ,   6. / L ,   4.    );
    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U4(&axial_action, &mut m);
    
    m
}

pub fn frame_element_with_shear_deformation_stiffness_matrix(
    element: &FrameElement
) -> Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>> {

    let section = element.geometry.cross_section;
    let material = element.material;

    let L = element.geometry.length;
    let L2 = L * L;

    let (A, Avy, Avz, Iy, Iz, J) = (section.A, section.Avy, section.Avz, section.Iy, section.Iz, section.J);
    let (E, G) = (material.E, material.G);

    let mut m = Matrix::<f64, U12, U12, ArrayStorage<f64, U12, U12>>::zeros();

    // Axial
    let mut index_map = Vector2::new(0, 6);
    let stiffness_matrix = E * A / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U2(&axial_action, &mut m);

    // Torsion
    let mut index_map = Vector2::new(3, 9);
    let stiffness_matrix = G * J / L * Matrix2::new( 1., -1.,
                                                    -1.,  1.);
    let axial_action = SingleActionStiffnessMatrix2by2 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U2(&axial_action, &mut m);

    // Bending About Z
    let mut index_map = Vector4::new(7, 11, 1, 5);
    let eta = E * Iz / Avy / G;
    let stiffness_matrix = E * Iz / L / (L2 / 12. + eta) * Matrix4::new(     1.,       -L / 2.,     -1.,       -L / 2.,
                                                                        -L / 2., L2 / 3. + eta,  L / 2., L2 / 6. - eta,
                                                                            -1.,        L / 2.,      1.,        L / 2.,
                                                                        -L / 2., L2 / 6. - eta,  L / 2., L2 / 3. + eta);

    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U4(&axial_action, &mut m);


    // Bending About Y
    let mut index_map = Vector4::new(8, 10, 2, 4);
    let eta = E * Iy / Avz / G;
    let stiffness_matrix = E * Iy / L / (L2 / 12. + eta) * Matrix4::new(     1.,        L / 2.,     -1.,        L / 2.,
                                                                         L / 2., L2 / 3. + eta, -L / 2., L2 / 6. - eta,
                                                                            -1.,       -L / 2.,      1.,       -L / 2.,
                                                                         L / 2., L2 / 6. - eta, -L / 2., L2 / 3. + eta);

    let axial_action = SingleActionStiffnessMatrix4by4 { index_map, stiffness_matrix };
    merge_single_action_into_complete_U4(&axial_action, &mut m);
    
    m
}

fn merge_single_action_into_complete_U2(
    single: &SingleActionStiffnessMatrix2by2,
    complete: &mut Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>>,
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

fn merge_single_action_into_complete_U4(
    single: &SingleActionStiffnessMatrix4by4,
    complete: &mut Matrix<f64, U12, U12, ArrayStorage<f64, U12, U12>>,
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
