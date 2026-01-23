use std::env;

fn main() {
    // Compile for neoverse n1
    env::set_var("CFLAGS", "-mcpu=neoverse-n1+sve");

    cc::Build::new()
        .file("src/sve/ntt_sve128.c")
        .include("src/sve/inc")
        .compile("ntt_sve128");

    println!("cargo::rerun-if-changed=src/sve/ntt_sve128.c");
    println!("cargo::rerun-if-changed=src/sve/inc/ntt_sve128.h");
}
