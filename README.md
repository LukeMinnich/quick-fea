# quick-fea

![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/LukeMinnich/quick-fea/build-master/master)

A modern finite element library for structural analysis written entirely in Rust.

### MVP Feature Set
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

### Stretch Feature Set
* [ ] Loads on elements between nodes
* [ ] 3 point membrane triangular element
* [ ] 6 point membrane triangular element
* [ ] spring supports
* [ ] AMVD response of truss and frame elements
* [ ] Electron App build for Native Windows deployment
* [ ] Some cross section input path generators

### Wishlist Feature Set
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


### 🛠️ Build with `wasm-pack build`

```
wasm-pack build
```

### 🔬 Test in Headless Browsers with `wasm-pack test`

```
wasm-pack test --headless --firefox
```

### 🎁 Publish to NPM with `wasm-pack publish`

```
wasm-pack publish
```


## 🔋 Batteries Included

* [`wasm-bindgen`](https://github.com/rustwasm/wasm-bindgen) for communicating
  between WebAssembly and JavaScript.
* [`console_error_panic_hook`](https://github.com/rustwasm/console_error_panic_hook)
  for logging panic messages to the developer console.
* [`wee_alloc`](https://github.com/rustwasm/wee_alloc), an allocator optimized
  for small code size.
