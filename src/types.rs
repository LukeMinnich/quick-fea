#[derive(Clone, PartialEq, ::prost::Message)]
pub struct IsotropicMaterial {
}
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct FrameElement {
    #[prost(string, tag="1")]
    pub id: std::string::String,
    #[prost(string, tag="2")]
    pub name: std::string::String,
    #[prost(string, tag="3")]
    pub start_node_id: std::string::String,
    #[prost(string, tag="4")]
    pub end_node_id: std::string::String,
    #[prost(message, optional, tag="5")]
    pub geometry: ::std::option::Option<FrameElementGeometry>,
    #[prost(message, optional, tag="6")]
    pub material: ::std::option::Option<IsotropicMaterial>,
    #[prost(message, optional, tag="7")]
    pub start_releases: ::std::option::Option<FrameElementEndRelease>,
    #[prost(message, optional, tag="8")]
    pub end_releases: ::std::option::Option<FrameElementEndRelease>,
}
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct FrameElementGeometry {
    #[prost(message, optional, tag="1")]
    pub cross_section: ::std::option::Option<FrameElementCrossSectionProperties>,
    #[prost(double, tag="2")]
    pub rotation: f64,
}
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct FrameElementCrossSectionProperties {
    #[prost(double, tag="1")]
    pub a: f64,
    #[prost(double, tag="2")]
    pub avy: f64,
    #[prost(double, tag="3")]
    pub avz: f64,
    #[prost(double, tag="4")]
    pub j: f64,
    #[prost(double, tag="5")]
    pub iy: f64,
    #[prost(double, tag="6")]
    pub iz: f64,
}
#[derive(Clone, PartialEq, ::prost::Message)]
pub struct FrameElementEndRelease {
    #[prost(enumeration="frame_element_end_release::Category", tag="1")]
    pub a: i32,
    #[prost(enumeration="frame_element_end_release::Category", tag="2")]
    pub vy: i32,
    #[prost(enumeration="frame_element_end_release::Category", tag="3")]
    pub vz: i32,
    #[prost(enumeration="frame_element_end_release::Category", tag="4")]
    pub t: i32,
    #[prost(enumeration="frame_element_end_release::Category", tag="5")]
    pub my: i32,
    #[prost(enumeration="frame_element_end_release::Category", tag="6")]
    pub mz: i32,
}
pub mod frame_element_end_release {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ::prost::Enumeration)]
    #[repr(i32)]
    pub enum Category {
        Fixed = 0,
        Free = 1,
    }
}
