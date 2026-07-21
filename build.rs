fn main() {
    // Only build the GMP shim when the gmp_cross_check feature is enabled.
    // This avoids requiring libgmp-dev for normal builds.
    #[cfg(feature = "gmp_cross_check")]
    {
        cc::Build::new()
            .file("tests/gmp_shim.c")
            .compile("gmp_shim");
        println!("cargo:rustc-link-lib=gmp");
    }
}
