extern crate nalgebra as na;
extern crate quick_fea_types as types;
extern crate serde;
extern crate wasm_bindgen;
#[macro_use]
extern crate approx;
#[macro_use]
extern crate lazy_static;

pub mod analysis;
pub mod elements;
pub mod models;
pub mod utils;

use crate::models::frame::{FrameElement, FrameStiffness};
use crate::models::node::Node;
use std::collections::HashMap;
use std::sync::RwLock;

lazy_static! {
    static ref ANALYSIS_DATA: RwLock<AnalysisData> = RwLock::new(AnalysisData::new());
    static ref ELEMENT_DATA: RwLock<ElementData> = RwLock::new(ElementData::new());
}

pub struct AnalysisData {
    pub applied_forces: Vec<f64>,
    pub frame_stiffnesses: HashMap<String, FrameStiffness>,
    pub world_stiffness: HashMap<(usize, usize), f64>,
}

impl AnalysisData {
    fn new() -> Self {
        AnalysisData {
            applied_forces: Vec::new(),
            frame_stiffnesses: HashMap::<String, FrameStiffness>::new(),
            world_stiffness: HashMap::<(usize, usize), f64>::new(),
        }
    }
}

pub struct ElementData {
    pub frames: HashMap<String, FrameElement>,
    pub nodes: HashMap<String, Node>,
}

impl ElementData {
    fn new() -> Self {
        ElementData {
            frames: HashMap::<String, FrameElement>::new(),
            nodes: HashMap::<String, Node>::new(),
        }
    }
}

pub fn add_node(node: Node) {
    ELEMENT_DATA
        .write()
        .unwrap()
        .nodes
        .insert(node.id.clone(), node);
}

// TODO is there any way to avoid the call to clone?
pub fn get_node_by_id(id: &str) -> Option<Node> {
    match ELEMENT_DATA.read() {
        Ok(data) => {
            let node: &Node = data.nodes.get(id)?;
            Some((*node).clone())
        }
        Err(_) => None,
    }
}

pub fn add_frame_element(frame: FrameElement) {
    ELEMENT_DATA
        .write()
        .unwrap()
        .frames
        .insert(frame.id.clone(), frame);
}

pub fn update_frame_element_stiffness(
    frame: &FrameElement,
    local: na::MatrixN<f64, na::U12>,
    world: na::MatrixN<f64, na::U12>,
) {
    ANALYSIS_DATA
        .write()
        .unwrap()
        .frame_stiffnesses
        .insert(frame.id.clone(), FrameStiffness { local, world });
}

// TODO is there any way to avoid the call to clone?
pub fn get_frame_element_by_id(id: &str) -> Option<FrameElement> {
    match ELEMENT_DATA.read() {
        Ok(data) => {
            let frame: &FrameElement = data.frames.get(id)?;
            Some((*frame).clone())
        }
        Err(_) => None,
    }
}
