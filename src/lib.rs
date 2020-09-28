extern crate nalgebra as na;
extern crate serde;
extern crate wasm_bindgen;
// #[macro_use]
// extern crate approx;
// #[macro_use]
// extern crate derive_new;
#[macro_use]
extern crate lazy_static;

pub mod analysis;
pub mod elements;
pub mod types;
pub mod utils;

use crate::types::frame::{FrameElement, FrameStiffness};
use crate::types::node::Node;
use std::collections::HashMap;

lazy_static! {
    static ref ANALYSIS_DATA: AnalysisData = AnalysisData::new();
    static ref ELEMENT_DATA: ElementData = ElementData::new();
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
