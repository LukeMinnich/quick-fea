![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/LukeMinnich/quick-fea/build-master/master)

# quick-fea
A modern finite element library for structural analysis written entirely in Rust

## MVP Feature Set
* [ ] loads at nodes
* [x] 3D truss members
* [x] 3D frame members
* [ ] Fixed or free supports (in 6 dof)
* [ ] cross section input path definitions
* [x] cross section input section properties
* [ ] linear elastic static analysis
* [ ] AMV response of truss and frame elements
* [ ] Deflection response at nodes
* [ ] WebAssembly bindings
* [ ] Serde Serializable or Protobuf types

## Stretch Feature Set
* [ ] Loads on elements between nodes
* [ ] 3 point membrane triangular element
* [ ] 6 point membrane triangular element
* [ ] spring supports
* [ ] AMVD response of truss and frame elements
* [ ] Electron App build for Native Windows deployment
* [ ] Some cross section input path generators

## Wishlist Feature Set
* [ ] Undo / Redo
* [ ] Concurrent User Modeling Env
* [ ] Concurrent command line test runner
* [ ] Javascript Module support
* [ ] YAML definitions for specific FEA elements
* [ ] Save modeling changes
* [ ] Save analysis data
* [ ] Backwards compatibility with model versions and data versions
* [ ] Option to use old data version (based on previous build, e.g. v1.0.0 app)
* [ ] Transactional database 
* [ ] Custom sections
