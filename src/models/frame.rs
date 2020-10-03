use crate::models::material::*;
use crate::models::node::*;
use crate::*;
use na::{Matrix3, MatrixN, Point3, U12};

#[derive(Clone)]
pub struct FrameElement {
    pub id: String,
    pub start_node_id: String,
    pub end_node_id: String,
    pub start_releases: FrameEndReleases,
    pub end_releases: FrameEndReleases,
    pub geometry: FrameGeometry,
    pub material: IsotropicMaterial,
}

impl FrameElement {
    /// Returns the length if it can be determined from the start and end nodes
    /// or positive infinity otherwise.
    pub fn length_or_inf(&self) -> f64 {
        match self.length() {
            Some(x) => x,
            None => f64::INFINITY,
        }
    }
    /// Returns the length if it can be determined from the start and end nodes.
    pub fn length(&self) -> Option<f64> {
        let start: Point3<f64> = self.start_node()?.coordinate;
        let end: Point3<f64> = self.end_node()?.coordinate;
        Some(na::distance(&start, &end))
    }
    pub fn start_node(&self) -> Option<Node> {
        get_node_by_id(&self.start_node_id)
    }
    pub fn end_node(&self) -> Option<Node> {
        get_node_by_id(&self.end_node_id)
    }
}

#[derive(Clone)]
pub struct FrameGeometry {
    pub local_axes: Matrix3<f64>,
    pub cross_section: CrossSection,
}

// TODO change this function accept endpoints and a relative rotation about longitudinal axis
impl FrameGeometry {
    pub fn new(local_axes: Matrix3<f64>, cross_section: CrossSection) -> FrameGeometry {
        FrameGeometry {
            local_axes,
            cross_section,
        }
    }
}

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct FrameEndReleases {
    pub A: FrameEndRelease,
    pub Vy: FrameEndRelease,
    pub Vz: FrameEndRelease,
    pub T: FrameEndRelease,
    pub My: FrameEndRelease,
    pub Mz: FrameEndRelease,
}

#[allow(dead_code)]
impl FrameEndReleases {
    pub const fn fully_fixed() -> FrameEndReleases {
        FrameEndReleases {
            A: FrameEndRelease::Fixed,
            Vy: FrameEndRelease::Fixed,
            Vz: FrameEndRelease::Fixed,
            T: FrameEndRelease::Fixed,
            My: FrameEndRelease::Fixed,
            Mz: FrameEndRelease::Fixed,
        }
    }
    pub const fn pinned() -> FrameEndReleases {
        FrameEndReleases {
            A: FrameEndRelease::Fixed,
            Vy: FrameEndRelease::Fixed,
            Vz: FrameEndRelease::Fixed,
            T: FrameEndRelease::Fixed,
            My: FrameEndRelease::Free,
            Mz: FrameEndRelease::Free,
        }
    }
}

#[derive(PartialEq, Eq, Clone)]
pub enum FrameEndRelease {
    Fixed,
    Free,
}

#[derive(Clone)]
#[allow(non_snake_case)]
pub struct CrossSection {
    pub A: f64,
    pub Avy: f64,
    pub Avz: f64,
    pub J: f64,
    pub Iy: f64,
    pub Iz: f64,
}

pub struct FrameStiffness {
    pub local: MatrixN<f64, U12>,
    pub world: MatrixN<f64, U12>,
}
