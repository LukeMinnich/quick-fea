use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Serialize, Deserialize)]
pub enum Axis {
    X,
    Y,
    Z,
}

#[derive(Serialize, Deserialize)]
pub struct FrameElement {
    pub id: String,
    pub ends: (Point3, Point3),
}
