use std::{env, path::PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=heap.c");
    let os = env::var("CARGO_CFG_TARGET_OS").unwrap();
    let gnu = env::var("CARGO_CFG_TARGET_ENV").unwrap() == "gnu";
    let linkage = match os.as_str() {
        "macos" => "-dynamiclib",
        "linux" if gnu => "-shared",
        _ => return, // Native Criterion timings do not require the heap sampler.
    };
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("bitfront-heap.so");
    let status = cc::Build::new().get_compiler().to_command()
        .args([linkage, "-std=c11", "-fPIC", "-fno-builtin", "-pthread", "heap.c", "-o"])
        .arg(output).status().expect("could not compile the benchmark heap sampler");
    assert!(status.success(), "heap sampler compilation failed: {status}");
}
