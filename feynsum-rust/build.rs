use std::path::PathBuf;

fn main() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let futhark_backend = match option_env!("FUTHARK_BACKEND") {
        None => futhark_bindgen::Backend::C,
        Some("c") => futhark_bindgen::Backend::C,
        Some("multicore") => futhark_bindgen::Backend::Multicore,
        Some("cuda") => futhark_bindgen::Backend::CUDA,
        _ => panic!("Unknown FUTHARK_BACKEND"),
    };

    futhark_bindgen::build(
        futhark_backend,
        PathBuf::from(manifest_dir).join("futhark").join("lib.fut"),
        "futhark_lib.rs",
    )
}
