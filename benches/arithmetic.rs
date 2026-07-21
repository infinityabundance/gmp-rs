//! Benchmarks for gmp-rs arithmetic operations.
//!
//! Run with:  cargo bench
//!
//! These benchmarks measure the throughput of core operations across
//! different operand sizes (1-limb, 2-limb, 4-limb, 8-limb).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gmp_rs::Mpz;

// ─────────────────────────────────────────────────────────────────────────
// Helper: construct operands of various sizes
// ─────────────────────────────────────────────────────────────────────────

fn single_limb_val(n: u64) -> Mpz {
    Mpz::from_u64(n)
}

fn two_limb_val() -> Mpz {
    // 2^64 + 42 (uses 2 limbs)
    Mpz::from_u128((1u128 << 64) + 42)
}

fn four_limb_val() -> Mpz {
    // Build a 4-limb value by shifting
    let mut v = Mpz::from_u64(1);
    for _ in 0..3 {
        v = v.try_mul_2exp(64).unwrap();
        v = v.try_add(&Mpz::from_u64(u64::MAX)).unwrap();
    }
    v
}

fn eight_limb_val() -> Mpz {
    // Build an 8-limb (512-bit) value
    let mut v = Mpz::from_u64(1);
    for _ in 0..7 {
        v = v.try_mul_2exp(64).unwrap();
        v = v.try_add(&Mpz::from_u64(u64::MAX)).unwrap();
    }
    v
}

// ─────────────────────────────────────────────────────────────────────────
// Addition benchmarks
// ─────────────────────────────────────────────────────────────────────────

fn bench_add_single_limb(c: &mut Criterion) {
    let a = single_limb_val(1_000_000_000_000);
    let b = single_limb_val(2_000_000_000_000);
    c.bench_function("add/single_limb", |bench| {
        bench.iter(|| black_box(a.try_add(black_box(&b))))
    });
}

fn bench_add_two_limb(c: &mut Criterion) {
    let a = two_limb_val();
    let b = two_limb_val();
    c.bench_function("add/two_limb", |bench| {
        bench.iter(|| black_box(a.try_add(black_box(&b))))
    });
}

fn bench_add_eight_limb(c: &mut Criterion) {
    let a = eight_limb_val();
    let b = eight_limb_val();
    c.bench_function("add/eight_limb", |bench| {
        bench.iter(|| black_box(a.try_add(black_box(&b))))
    });
}

// ─────────────────────────────────────────────────────────────────────────
// Multiplication benchmarks
// ─────────────────────────────────────────────────────────────────────────

fn bench_mul_single_limb(c: &mut Criterion) {
    let a = single_limb_val(1_000_000_000_000);
    let b = single_limb_val(2_000_000_000_000);
    c.bench_function("mul/single_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

fn bench_mul_two_limb(c: &mut Criterion) {
    let a = two_limb_val();
    let b = two_limb_val();
    c.bench_function("mul/two_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

fn bench_mul_four_limb(c: &mut Criterion) {
    let a = four_limb_val();
    let b = four_limb_val();
    c.bench_function("mul/four_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

fn bench_mul_eight_limb(c: &mut Criterion) {
    let a = eight_limb_val();
    let b = eight_limb_val();
    c.bench_function("mul/eight_limb", |bench| {
        bench.iter(|| black_box(a.try_mul(black_box(&b))))
    });
}

// ─────────────────────────────────────────────────────────────────────────
// Division benchmarks
// ─────────────────────────────────────────────────────────────────────────

fn bench_tdiv_qr_single_limb(c: &mut Criterion) {
    let a = single_limb_val(1_000_000_000_000_000);
    let b = single_limb_val(3_000_000);
    c.bench_function("tdiv_qr/single_limb_divisor", |bench| {
        bench.iter(|| black_box(a.tdiv_qr(black_box(&b))))
    });
}

fn bench_tdiv_qr_multi_limb(c: &mut Criterion) {
    let a = eight_limb_val();
    let b = four_limb_val();
    c.bench_function("tdiv_qr/multi_limb_divisor", |bench| {
        bench.iter(|| black_box(a.tdiv_qr(black_box(&b))))
    });
}

// ─────────────────────────────────────────────────────────────────────────
// Comparison benchmarks
// ─────────────────────────────────────────────────────────────────────────

fn bench_cmp_single_limb(c: &mut Criterion) {
    let a = single_limb_val(1_000_000);
    let b = single_limb_val(2_000_000);
    c.bench_function("cmp/single_limb", |bench| {
        bench.iter(|| black_box(a.cmp(black_box(&b))))
    });
}

fn bench_cmp_eight_limb(c: &mut Criterion) {
    let a = eight_limb_val();
    let b = eight_limb_val();
    c.bench_function("cmp/eight_limb", |bench| {
        bench.iter(|| black_box(a.cmp(black_box(&b))))
    });
}

// ─────────────────────────────────────────────────────────────────────────
// Group all benchmarks
// ─────────────────────────────────────────────────────────────────────────

criterion_group!(
    benches,
    // Addition
    bench_add_single_limb,
    bench_add_two_limb,
    bench_add_eight_limb,
    // Multiplication
    bench_mul_single_limb,
    bench_mul_two_limb,
    bench_mul_four_limb,
    bench_mul_eight_limb,
    // Division
    bench_tdiv_qr_single_limb,
    bench_tdiv_qr_multi_limb,
    // Comparison
    bench_cmp_single_limb,
    bench_cmp_eight_limb,
);

criterion_main!(benches);
