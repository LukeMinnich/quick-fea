use std::env;
use std::path::PathBuf;

fn main() {
    // Configure env variable for protoc compiled types
    let mut output_types_path: PathBuf = env::current_dir().unwrap();
    output_types_path.push("src");
    println!("cargo:rustc-env=OUT_TYPES={}", output_types_path.display());

    // Compile protoc types
    let mut prost_build = prost_build::Config::new();
    prost_build.out_dir("src");

    prost_build
        .compile_protos(
            &[
                "types/proto/utility.proto",
                "types/proto/material.proto",
                "types/proto/frame.proto",
            ],
            &["types/proto/"],
        )
        .unwrap();
}
