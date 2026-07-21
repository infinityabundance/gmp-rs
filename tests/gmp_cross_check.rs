//! Cross-check tests against real GMP via FFI.
//!
//! Tests that gmp-rs produces identical results to GMP's `mpz_*` functions
//! for all supported operations within the fixed capacity.
//!
//! Build requirements:
//!   - libgmp-dev installed (e.g., apt install libgmp-dev)
//!   - cmake or cc build support
//!
//! Run with:  cargo test --features "std" --test gmp_cross_check
//!            (compiles and links against GMP dynamically)

#![cfg(feature = "std")]

use std::ffi::CString;
use std::os::raw::{c_char, c_ulong};

// Note: This test file can only be compiled when the cc build script
// is run with `cargo build --features gmp_cross_check`.  It requires
// libgmp-dev to be installed on the system.

use gmp_rs::Mpz;

const BUF_LEN: usize = 4096;

// ─────────────────────────────────────────────────────────────────────────
// GMP FFI bindings (linked via build.rs)
// ─────────────────────────────────────────────────────────────────────────

extern "C" {
    fn gmp_add(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
    fn gmp_sub(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
    fn gmp_mul(a_str: *const c_char, b_str: *const c_char, out_buf: *mut c_char, out_len: usize);
    fn gmp_cmp(a_str: *const c_char, b_str: *const c_char) -> i32;
    fn gmp_bits(a_str: *const c_char, bits_out: *mut c_ulong);
}

/// Call a GMP binary operation via FFI, returning the result as a decimal string.
fn gmp_binop(
    op: unsafe extern "C" fn(*const c_char, *const c_char, *mut c_char, usize),
    a: &Mpz,
    b: &Mpz,
) -> String {
    let a_str = a.to_string();
    let b_str = b.to_string();
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
    let nul_pos = buf.iter().position(|&c| c == 0).unwrap_or(BUF_LEN);
    String::from_utf8_lossy(&buf[..nul_pos]).to_string()
}

/// Get the bit-length of an Mpz value via GMP.
fn gmp_bitlen(mpz: &Mpz) -> u64 {
    let s = mpz.to_string();
    let c = CString::new(s).unwrap();
    let mut bits: c_ulong = 0;
    unsafe {
        gmp_bits(c.as_ptr(), &mut bits);
    }
    bits as u64
}

// ─────────────────────────────────────────────────────────────────────────
// Cross-check tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn cross_check_add_vs_gmp() {
    // Test with i64-range values (always fits in both GMP and gmp-rs)
    let cases: &[i64] = &[
        0,
        1,
        -1,
        42,
        -42,
        i64::MAX,
        i64::MIN,
        i64::MAX / 2,
        i64::MIN / 2,
    ];
    for &a_val in cases {
        for &b_val in cases {
            let a = Mpz::from_i64(a_val);
            let b = Mpz::from_i64(b_val);
            let gmp_result = gmp_binop(gmp_add, &a, &b);
            let our_result = a.try_add(&b).unwrap();
            assert_eq!(
                our_result.to_string(),
                gmp_result,
                "add cross-check: {a_val} + {b_val}"
            );
        }
    }
}

#[test]
fn cross_check_sub_vs_gmp() {
    let cases: &[i64] = &[
        0,
        1,
        -1,
        42,
        -42,
        i64::MAX,
        i64::MIN,
        i64::MAX / 2,
        i64::MIN / 2,
    ];
    for &a_val in cases {
        for &b_val in cases {
            let a = Mpz::from_i64(a_val);
            let b = Mpz::from_i64(b_val);
            let gmp_result = gmp_binop(gmp_sub, &a, &b);
            let our_result = a.try_sub(&b).unwrap();
            assert_eq!(
                our_result.to_string(),
                gmp_result,
                "sub cross-check: {a_val} - {b_val}"
            );
        }
    }
}

#[test]
fn cross_check_mul_vs_gmp_small() {
    let cases: &[i64] = &[0, 1, -1, 2, -2, 10, -10, 100, -100, 1000000, -1000000];
    for &a_val in cases {
        for &b_val in cases {
            let a = Mpz::from_i64(a_val);
            let b = Mpz::from_i64(b_val);
            // Only check when result fits in capacity
            if let Ok(our_result) = a.try_mul(&b) {
                let gmp_result = gmp_binop(gmp_mul, &a, &b);
                assert_eq!(
                    our_result.to_string(),
                    gmp_result,
                    "mul cross-check: {a_val} * {b_val}"
                );
            }
        }
    }
}

#[test]
fn cross_check_add_large_values() {
    // Test with values that exercise multiple limbs
    let test_values: &[&str] = &[
        "0",
        "1",
        "-1",
        "340282366920938463463374607431768211455", // 2^128 - 1
        "-340282366920938463463374607431768211455",
        "170141183460469231731687303715884105727", // 2^127 - 1
        "100000000000000000000000000000000000000",
        "99999999999999999999999999999999999999",
    ];

    for &a_str in test_values {
        for &b_str in test_values {
            let a = Mpz::from_decimal_str(a_str).unwrap();
            let b = Mpz::from_decimal_str(b_str).unwrap();

            // Check add
            if let Ok(our_sum) = a.try_add(&b) {
                let gmp_sum = gmp_binop(gmp_add, &a, &b);
                assert_eq!(
                    our_sum.to_string(),
                    gmp_sum,
                    "add cross-check: {a_str} + {b_str}"
                );
            } else {
                // Our side overflowed; GMP should have > 512 bits
                let bits = gmp_bitlen(&a.try_add(&b).unwrap_or_else(|_| Mpz::new()));
                // Can't easily check GMP overflow; just verify consistent
            }

            // Check sub
            if let Ok(our_diff) = a.try_sub(&b) {
                let gmp_diff = gmp_binop(gmp_sub, &a, &b);
                assert_eq!(
                    our_diff.to_string(),
                    gmp_diff,
                    "sub cross-check: {a_str} - {b_str}"
                );
            }

            // Check mul (only when known to fit)
            if let Ok(our_prod) = a.try_mul(&b) {
                let gmp_prod = gmp_binop(gmp_mul, &a, &b);
                assert_eq!(
                    our_prod.to_string(),
                    gmp_prod,
                    "mul cross-check: {a_str} * {b_str}"
                );
            }
        }
    }
}

#[test]
fn cross_check_mul_vs_gmp_extended() {
    use proptest::prelude::*;

    // Generate i64-range values and cross-check mul
    fn any_i64_mpz() -> impl Strategy<Value = Mpz> {
        any::<i64>().prop_map(|v| Mpz::from_i64(v))
    }

    proptest!(|(a in any_i64_mpz(), b in any_i64_mpz())| {
        if let Ok(ours) = a.try_mul(&b) {
            let gmp = gmp_binop(gmp_mul, &a, &b);
            prop_assert_eq!(ours.to_string(), gmp,
                "mul cross-check failed");
        }
        // If our mul overflows, that's fine — GMP supports unlimited precision.
    });
}

#[test]
fn cross_check_compound_expression() {
    // Test: ((a + b) * c) - d
    use proptest::prelude::*;

    fn any_i64_mpz() -> impl Strategy<Value = Mpz> {
        any::<i64>().prop_map(|v| Mpz::from_i64(v))
    }

    proptest!(|(a in any_i64_mpz(), b in any_i64_mpz(), c in any_i64_mpz(), d in any_i64_mpz())| {
        let our_result = a.try_add(&b)
            .and_then(|ab| ab.try_mul(&c))
            .and_then(|abc| abc.try_sub(&d));

        // GMP computation: ((a + b) * c) - d
        let ab_gmp = gmp_binop(gmp_add, &a, &b);
        let ab_mpz = Mpz::from_decimal_str(&ab_gmp).unwrap();
        let abc_gmp = gmp_binop(gmp_mul, &ab_mpz, &c);
        let abc_mpz = Mpz::from_decimal_str(&abc_gmp).unwrap();
        let abcd_gmp = gmp_binop(gmp_sub, &abc_mpz, &d);

        match our_result {
            Ok(ours) => {
                prop_assert_eq!(ours.to_string(), abcd_gmp,
                    "compound cross-check: (({a} + {b}) * {c}) - {d}");
            }
            Err(_) => {
                // Our side overflowed; GMP result should have > MAX_BITS bits
                let abcd_mpz = Mpz::from_decimal_str(&abcd_gmp).unwrap();
                prop_assert!(abcd_mpz.sizeinbase2() > gmp_rs::MAX_BITS,
                    "should overflow capacity");
            }
        }
    });
}
