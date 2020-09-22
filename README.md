# quick-fea
A modern finite element library for structural analysis written entirely in Rust

## MVP Feature Set
* loads at nodes
* 3D truss members
* 3D frame members
* Fixed or free supports (in 6 dof)
* cross section input path definitions
* cross section input section properties
* linear elastic static analysis
* AMV response of truss and frame elements
* Deflection response at nodes
* THREE.js viewer w/ setbacks at pinned ends
* UI to define member end points, pick cross section, set orientation, and releases
* Rust WASM bindings
* Use Hypar Elements defs for JS?
* Use Serde defs for JS?

## Stretch Feature Set
* Loads on elements between nodes
* 3 point membrane triangular element
* 6 point membrane triangular element
* spring supports
* AMVD response of truss and frame elements
* Electron App build for Native Windows deployment
* Some cross section input path generators

## Wishlist Feature Set
* Undo / Redo
* Concurrent User Modeling Env
* Concurrent command line test runner
* Javascript Module support
* YAML definitions for specific FEA elements
* Save modeling changes
* Save analysis data
* Backwards compatibility with model versions and data versions
* Option to use old data version (based on previous build, e.g. v1.0.0 app)
* Transactional database 
* Custom sections
