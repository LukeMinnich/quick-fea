use na::{Point3, Vector6};

#[derive(Clone)]
pub struct Node {
    pub id: String,
    pub coordinate: Point3<f64>,
    pub degrees_of_freedom: Vector6<usize>,
}
