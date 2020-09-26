use na::{Matrix3, Point3};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Serialize, Deserialize, new)]
pub struct FrameElement {
    pub id: String,
    pub geometry: FrameElementGeometry,
    pub material: IsotropicMaterial,
}

#[derive(Serialize, Deserialize)]
pub struct FrameElementGeometry {
    pub ends: (Point3<f64>, Point3<f64>),
    pub local_axes: Matrix3<f64>,
    pub length: f64,
    pub cross_section: CrossSectionProperties,
}

// TODO change this function accept endpoints and a relative rotation about longitudinal axis
impl FrameElementGeometry {
    pub fn new(
        ends: (Point3<f64>, Point3<f64>),
        local_axes: Matrix3<f64>,
        cross_section: CrossSectionProperties,
    ) -> FrameElementGeometry {
        let length = na::distance(&ends.0, &ends.1);
        FrameElementGeometry {
            ends,
            local_axes,
            length,
            cross_section,
        }
    }
}

#[derive(Serialize, Deserialize, new)]
pub struct CrossSectionProperties {
    pub A: f64,
    pub Avy: f64,
    pub Avz: f64,
    pub J: f64,
    pub Iy: f64,
    pub Iz: f64,
}

#[derive(Serialize, Deserialize)]
pub struct IsotropicMaterial {
    pub E: f64,
    pub G: f64,
    pub nu: f64,
}

impl IsotropicMaterial {
    pub fn new(E: f64, nu: f64) -> IsotropicMaterial {
        let G = E / (2. * (1. + nu));
        IsotropicMaterial { E, G, nu }
    }
}
