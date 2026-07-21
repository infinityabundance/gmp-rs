//! Benchmarks comparing gmp-rs vs raw GMP C calls.
//!
//! Run with:  cargo bench --features gmp_cross_check --bench gmp_comparison
//!
//! Requires:  libgmp-dev (or equivalent) installed on the system.
//!
//! These benchmarks quantify the overhead of gmp-rs's safe Rust API
//! relative to calling GMP's C API directly via FFI.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::ffi::CString;
use std::os::raw::c_char;

// ─────────────────────────────────────────────────────────────────────────
// GMP FFI bindings (same shim as cross-check tests)
// ─────────────────────────────────────────────────────────────────────────

extern "C" {
    fn gmp_add(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
    fn gmp_sub(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
    fn gmp_mul(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
}

const BUF_LEN: usize = 4096;

fn gmp_binop(
    op: unsafe extern "C" fn(*const c_char, *const c_char, *mut c_char, usize),
    a_str: &str,
    b_str: &str,
) {
    let a_c = CString::new(a_str).unwrap();
    let b_c = CString::new(b_str).unwrap();
    let mut buf = vec![0u8; BUF_LEN];
    unsafe {
        op(
            a_c.as_ptr(),
            b_c.as_ptr(),
            buf.as_mut_ptr() as *mut c_char,
            BUF_LEN,
        );
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Operands (from arithmetic.rs helpers)
// ─────────────────────────────────────────────────────────────────────────

use gmp_rs::Mpz;

fn op_single() -> Mpz {
    Mpz::from_u64(1_000_000_000_000)
}
fn op_two() -> Mpz {
    Mpz::from_u128((1u128 << 64) + 42)
}
fn op_four() -> Mpz {
    let two_pow_192 = Mpz::from_u64(1).try_mul_2exp(192).unwrap();
    two_pow_192.try_sub(&Mpz::from_u64(1)).unwrap()
}
fn op_eight() -> Mpz {
    let mut v = Mpz::from_u64(1);
    for _ in 0..7 {
        v = v.try_mul_2exp(64).unwrap();
        v = v.try_add(&Mpz::from_u64(u64::MAX)).unwrap();
    }
    v
}

fn mpz_to_string(mpz: &Mpz) -> String {
    let mut buf = [0u8; 4096];
    let len = mpz.write_decimal_buf(&mut buf);
    String::from_utf8_lossy(&buf[..len]).to_string()
}

// ─────────────────────────────────────────────────────────────────────────
// Comparison benchmarks
// ─────────────────────────────────────────────────────────────────────────

fn bench_add_gmp_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    let a_str = mpz_to_string(&a);
    let b_str = mpz_to_string(&b);
    c.bench_function("gmp/add_single_limb", |bench| {
        bench.iter(|| gmp_binop(gmp_add, black_box(&a_str), black_box(&b_str)))
    });
}

fn bench_add_ours_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    c.bench_function("ours/add_single_limb", |bench| {
        bench.iter(|| black_box(a.try_add(black_box(&b))))
    });
}

fn bench_add_gmp_eight(c: &mut Criterion) {
    let a = op_eight();
    let b = op_eight();
    let a_str = mpz_to_string(&a);
    let b_str = mpz_to_string(&b);
    c.bench_function("gmp/add_eight_limb", |bench| {
        bench.iter(|| gmp_binop(gmp_add, black_box(&a_str), black_box(&b_str)))
    });
}

fn bench_add_ours_eight(c: &mut Criterion) {
    let a = op_eight();
    let b = op_eight();
    c.bench_function("ours/add_eight_limb", |bench| {
        bench.iter(|| black_box(a.try_add(black_box(&b))))
    });
}

fn bench_mul_gmp_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    let a_str = mpz_to_string(&a);
    let b_str = mpz_to_string(&b);
    c.bench_function("gmp/mul_single_limb", |bench| {
        bench.iter(|| gmp_binop(gmp_mul, black_box(&a_str), black_box(&b_str)))
    });
}

fn bench_mul_ours_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    c.bench_function("ours/mul_single_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

fn bench_mul_gmp_eight(c: &mut Criterion) {
    let a = op_eight();
    let b = op_four();
    let a_str = mpz_to_string(&a);
    let b_str = mpz_to_string(&b);
    c.bench_function("gmp/mul_eight_x_four_limb", |bench| {
        bench.iter(|| gmp_binop(gmp_mul, black_box(&a_str), black_box(&b_str)))
    });
}

fn bench_mul_ours_eight(c: &mut Criterion) {
    let a = op_eight();
    let b = op_four();
    c.bench_function("ours/mul_eight_x_four_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

fn bench_sub_gmp_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    let a_str = mpz_to_string(&a);
    let b_str = mpz_to_string(&b);
    c.bench_function("gmp/sub_single_limb", |bench| {
        bench.iter(|| gmp_binop(gmp_sub, black_box(&a_str), black_box(&b_str)))
    });
}

fn bench_sub_ours_single(c: &mut Criterion) {
    let a = op_single();
    let b = op_single();
    c.bench_function("ours/sub_single_limb", |bench| {
        bench.iter(|| black_box(a.try_sub(black_box(&b))))
    });
}

criterion_group!(
    benches,
    // Single-limb add
    bench_add_gmp_single,
    bench_add_ours_single,
    // Eight-limb add
    bench_add_gmp_eight,
    bench_add_ours_eight,
    // Single-limb mul
    bench_mul_gmp_single,
    bench_mul_ours_single,
    // Multi-limb mul
    bench_mul_gmp_eight,
    bench_mul_ours_eight,
    // Single-limb sub
    bench_sub_gmp_single,
    bench_sub_ours_single,
);

criterion_main!(benches);
