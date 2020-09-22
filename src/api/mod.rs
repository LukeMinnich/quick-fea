use crate::types::*;
use crate::utils::transform::*;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;

#[allow(dead_code)]
#[wasm_bindgen]
pub fn rotate_about_start_js(el: &JsValue, axis: &JsValue, angle: f64) -> JsValue {
    let el: FrameElement = el.into_serde().unwrap();
    let axis: Axis = axis.into_serde().unwrap();
    let rotated_el: FrameElement = rotate_about_start(&el, axis, angle);

    JsValue::from_serde(&rotated_el).unwrap()
}
