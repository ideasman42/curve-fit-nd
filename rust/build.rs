fn main() {
    let dims_max = std::env::var("DIMS_MAX").unwrap_or_else(|_| "32".to_string());
    println!("cargo:rerun-if-env-changed=DIMS_MAX");
    println!("cargo:rustc-env=DIMS_MAX={}", dims_max);
}
