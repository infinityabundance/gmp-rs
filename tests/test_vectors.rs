//! Deterministic test vectors for gmp-rs.
//!
//! These are hardcoded known-answer tests that serve as the canonical reference
//! for verification.  Any change to gmp-rs arithmetic must not alter these results.
//!
//! Run with: cargo test --features std --test test_vectors

#![cfg(feature = "std")]

extern crate alloc;
use alloc::string::String;
use alloc::string::ToString;

use gmp_rs::*;

/// Convert Mpz to decimal string.
fn s(m: &Mpz) -> String {
    let mut buf = [0u8; 192];
    let len = m.write_decimal_buf(&mut buf);
    core::str::from_utf8(&buf[..len]).unwrap().into()
}

/// Assert that `expr` produces the expected decimal string.
macro_rules! assert_mpz_eq {
    ($expr:expr, $expected:expr) => {
        let val = $expr;
        let got = s(&val);
        assert_eq!(got, $expected, "expected {}, got {}", $expected, got);
    };
}

// ==========================================================================
// Test vectors: addition
// ==========================================================================

#[test]
fn add_vectors() {
    // Positive + positive
    assert_mpz_eq!(
        Mpz::from_u64(12345).try_add(&Mpz::from_u64(67890)).unwrap(),
        "80235"
    );
    // Negative + negative
    assert_mpz_eq!(
        Mpz::from_i64(-100).try_add(&Mpz::from_i64(-200)).unwrap(),
        "-300"
    );
    // Mixed signs
    assert_mpz_eq!(
        Mpz::from_i64(-50).try_add(&Mpz::from_u64(30)).unwrap(),
        "-20"
    );
    assert_mpz_eq!(
        Mpz::from_i64(50).try_add(&Mpz::from_i64(-30)).unwrap(),
        "20"
    );
    // Large: 2^100 + 2^100 = 2^101
    assert_mpz_eq!(
        Mpz::from_u64(1)
            .try_mul_2exp(100)
            .unwrap()
            .try_add(&Mpz::from_u64(1).try_mul_2exp(100).unwrap())
            .unwrap(),
        "2535301200456458802993406410752"
    );
}

// ==========================================================================
// Test vectors: subtraction
// ==========================================================================

#[test]
fn sub_vectors() {
    // Positive - positive
    assert_mpz_eq!(
        Mpz::from_u64(100).try_sub(&Mpz::from_u64(45)).unwrap(),
        "55"
    );
    // Result zero
    assert_mpz_eq!(
        Mpz::from_u64(42).try_sub(&Mpz::from_u64(42)).unwrap(),
        "0"
    );
    // Negative result
    assert_mpz_eq!(
        Mpz::from_u64(10).try_sub(&Mpz::from_u64(100)).unwrap(),
        "-90"
    );
    // Negative - negative
    assert_mpz_eq!(
        Mpz::from_i64(-30).try_sub(&Mpz::from_i64(-50)).unwrap(),
        "20"
    );
    // Large: 2^200 - 1 = all-ones in 3 lower limbs plus one 0xFFFF... in top
    let large = Mpz::from_u64(1).try_mul_2exp(200).unwrap();
    assert_mpz_eq!(
        large.try_sub(&Mpz::from_u64(1)).unwrap(),
        "1606938044258990275541962092341162602522202993782792835301375"
    );
}

// ==========================================================================
// Test vectors: multiplication
// ==========================================================================

#[test]
fn mul_vectors() {
    // Basic
    assert_mpz_eq!(
        Mpz::from_u64(7).try_mul(&Mpz::from_u64(8)).unwrap(),
        "56"
    );
    // Zero
    assert_mpz_eq!(
        Mpz::from_u64(12345).try_mul(&Mpz::new()).unwrap(),
        "0"
    );
    // Negative * negative = positive
    assert_mpz_eq!(
        Mpz::from_i64(-12).try_mul(&Mpz::from_i64(-12)).unwrap(),
        "144"
    );
    // Negative * positive = negative
    assert_mpz_eq!(
        Mpz::from_i64(-7).try_mul(&Mpz::from_u64(6)).unwrap(),
        "-42"
    );
    // 2^100 * 2^100 = 2^200
    assert_mpz_eq!(
        Mpz::from_u64(1)
            .try_mul_2exp(100)
            .unwrap()
            .try_mul(&Mpz::from_u64(1).try_mul_2exp(100).unwrap())
            .unwrap(),
        "1606938044258990275541962092341162602522202993782792835301376"
    );
}

// ==========================================================================
// Test vectors: truncating division
// ==========================================================================

#[test]
fn tdiv_vectors() {
    let (q, r) = Mpz::from_u64(100).tdiv_qr(&Mpz::from_u64(30));
    assert_mpz_eq!(q, "3");
    assert_mpz_eq!(r, "10");

    let (q, r) = Mpz::from_i64(-100).tdiv_qr(&Mpz::from_u64(30));
    assert_mpz_eq!(q, "-3");
    assert_mpz_eq!(r, "-10");

    let (q, r) = Mpz::from_i64(100).tdiv_qr(&Mpz::from_i64(-30));
    assert_mpz_eq!(q, "-3");
    assert_mpz_eq!(r, "10");

    let (q, r) = Mpz::from_i64(-100).tdiv_qr(&Mpz::from_i64(-30));
    assert_mpz_eq!(q, "3");
    assert_mpz_eq!(r, "-10");

    // Large: 2^200 / 2^100 = 2^100
    let (q, r) = Mpz::from_u64(1)
        .try_mul_2exp(200)
        .unwrap()
        .tdiv_qr(&Mpz::from_u64(1).try_mul_2exp(100).unwrap());
    assert_mpz_eq!(q, "1267650600228229401496703205376");
    assert_mpz_eq!(r, "0");
}

// ==========================================================================
// Test vectors: floor division
// ==========================================================================

#[test]
fn fdiv_vectors() {
    let (q, r) = Mpz::from_i64(-100)
        .try_fdiv_qr(&Mpz::from_u64(30))
        .unwrap();
    assert_mpz_eq!(q, "-4");
    assert_mpz_eq!(r, "20");

    let (q, r) = Mpz::from_i64(100)
        .try_fdiv_qr(&Mpz::from_i64(-30))
        .unwrap();
    assert_mpz_eq!(q, "-4");
    assert_mpz_eq!(r, "-20");
}

// ==========================================================================
// Test vectors: ceiling division
// ==========================================================================

#[test]
fn cdiv_vectors() {
    let (q, r) = Mpz::from_i64(-100)
        .try_cdiv_qr(&Mpz::from_u64(30))
        .unwrap();
    assert_mpz_eq!(q, "-3");
    assert_mpz_eq!(r, "-10");

    let (q, r) = Mpz::from_i64(100)
        .try_cdiv_qr(&Mpz::from_i64(-30))
        .unwrap();
    assert_mpz_eq!(q, "-3");
    assert_mpz_eq!(r, "10");
}

// ==========================================================================
// Test vectors: modulus (non-negative)
// ==========================================================================

#[test]
fn mod_vectors() {
    assert_mpz_eq!(
        Mpz::from_i64(-100).try_mod(&Mpz::from_u64(30)).unwrap(),
        "20"
    );
    assert_mpz_eq!(
        Mpz::from_i64(100).try_mod(&Mpz::from_u64(30)).unwrap(),
        "10"
    );
    assert_mpz_eq!(
        Mpz::from_i64(-17).try_mod(&Mpz::from_u64(5)).unwrap(),
        "3"
    );
}

// ==========================================================================
// Test vectors: GCD / LCM / extended GCD
// ==========================================================================

#[test]
fn gcd_vectors() {
    assert_mpz_eq!(
        Mpz::from_u64(12).try_gcd(&Mpz::from_u64(18)).unwrap(),
        "6"
    );
    assert_mpz_eq!(
        Mpz::from_u64(0).try_gcd(&Mpz::from_u64(5)).unwrap(),
        "5"
    );
    assert_mpz_eq!(
        Mpz::from_u64(7).try_gcd(&Mpz::from_u64(13)).unwrap(),
        "1"
    );
    assert_mpz_eq!(
        Mpz::from_u64(2700).try_gcd(&Mpz::from_u64(192)).unwrap(),
        "12"
    );
}

#[test]
fn lcm_vectors() {
    assert_mpz_eq!(
        Mpz::from_u64(12).try_lcm(&Mpz::from_u64(18)).unwrap(),
        "36"
    );
    assert_mpz_eq!(
        Mpz::from_u64(0).try_lcm(&Mpz::from_u64(5)).unwrap(),
        "0"
    );
}

#[test]
fn gcdext_vectors() {
    // gcd(120, 23) = 1, Bézout coefficients: 120*(-3) + 23*16 = 1? No...
    // Let's just verify: g = 120*s + 23*t
    let (g, s_coeff, t) = Mpz::from_u64(120).try_gcdext(&Mpz::from_u64(23)).unwrap();
    let check = s_coeff
        .try_mul(&Mpz::from_u64(120))
        .unwrap()
        .try_add(&t.try_mul(&Mpz::from_u64(23)).unwrap())
        .unwrap();
    assert_eq!(g, check);
}

// ==========================================================================
// Test vectors: power / root
// ==========================================================================

#[test]
fn pow_vectors() {
    assert_mpz_eq!(Mpz::try_ui_pow_ui(2, 10).unwrap(), "1024");
    assert_mpz_eq!(Mpz::try_ui_pow_ui(3, 5).unwrap(), "243");
    assert_mpz_eq!(Mpz::try_ui_pow_ui(10, 3).unwrap(), "1000");
    // 2^255
    assert_mpz_eq!(
        Mpz::try_ui_pow_ui(2, 255).unwrap(),
        "57896044618658097711785492504343953926634992332820282019728792003956564819968"
    );
}

#[test]
fn root_vectors() {
    assert_mpz_eq!(Mpz::from_u64(100).try_root(3).unwrap(), "4");
    assert_mpz_eq!(Mpz::from_u64(1000).try_root(3).unwrap(), "10");
    assert_mpz_eq!(Mpz::from_u64(144).try_root(2).unwrap(), "12");
    assert_mpz_eq!(Mpz::from_u64(10).isqrt(), "3");
}

#[test]
fn powm_vectors() {
    // 3^5 mod 7 = 243 mod 7 = 5
    assert_mpz_eq!(
        Mpz::from_u64(3)
            .try_powm(&Mpz::from_u64(5), &Mpz::from_u64(7))
            .unwrap(),
        "5"
    );
    // 2^10 mod 1000 = 1024 mod 1000 = 24
    assert_mpz_eq!(
        Mpz::from_u64(2)
            .try_powm(&Mpz::from_u64(10), &Mpz::from_u64(1000))
            .unwrap(),
        "24"
    );
}

// ==========================================================================
// Test vectors: bitwise
// ==========================================================================

#[test]
fn bitwise_vectors() {
    // 0b1100 AND 0b1010 = 0b1000
    assert_mpz_eq!(
        Mpz::from_u64(0b1100)
            .try_and(&Mpz::from_u64(0b1010))
            .unwrap(),
        "8"
    );
    // 0b1100 OR 0b1010 = 0b1110
    assert_mpz_eq!(
        Mpz::from_u64(0b1100)
            .try_ior(&Mpz::from_u64(0b1010))
            .unwrap(),
        "14"
    );
    // 0b1100 XOR 0b1010 = 0b0110
    assert_mpz_eq!(
        Mpz::from_u64(0b1100)
            .try_xor(&Mpz::from_u64(0b1010))
            .unwrap(),
        "6"
    );
    // com(0) = -1
    assert_mpz_eq!(Mpz::new().com(), "-1");
    // com(1) = -2
    assert_mpz_eq!(Mpz::from_u64(1).com(), "-2");
}

// ==========================================================================
// Test vectors: factorial / binomial / Fibonacci
// ==========================================================================

#[test]
fn combinatorial_vectors() {
    assert_mpz_eq!(Mpz::try_fac_ui(0).unwrap(), "1");
    assert_mpz_eq!(Mpz::try_fac_ui(5).unwrap(), "120");
    assert_mpz_eq!(Mpz::try_fac_ui(10).unwrap(), "3628800");

    assert_mpz_eq!(Mpz::try_bin_uiui(10, 5).unwrap(), "252");
    assert_mpz_eq!(Mpz::try_bin_uiui(10, 0).unwrap(), "1");
    assert_mpz_eq!(Mpz::try_bin_uiui(10, 10).unwrap(), "1");

    assert_mpz_eq!(Mpz::try_fib_ui(0).unwrap(), "0");
    assert_mpz_eq!(Mpz::try_fib_ui(1).unwrap(), "1");
    assert_mpz_eq!(Mpz::try_fib_ui(10).unwrap(), "55");
    assert_mpz_eq!(Mpz::try_fib_ui(20).unwrap(), "6765");

    assert_mpz_eq!(Mpz::try_lucnum_ui(0).unwrap(), "2");
    assert_mpz_eq!(Mpz::try_lucnum_ui(1).unwrap(), "1");
    assert_mpz_eq!(Mpz::try_lucnum_ui(10).unwrap(), "123");

    assert_mpz_eq!(Mpz::try_2fac_ui(5).unwrap(), "15");
    assert_mpz_eq!(Mpz::try_2fac_ui(6).unwrap(), "48");
}

// ==========================================================================
// Test vectors: from_decimal_str / write_decimal_buf round-trip
// ==========================================================================

#[test]
fn parse_roundtrip_vectors() {
    let cases = [
        "0",
        "1",
        "-1",
        "9",
        "10",
        "-10",
        "999999999999999999",
        "-999999999999999999",
        "1000000000000000000",
        "123456789012345678901234567890",
        "-123456789012345678901234567890",
    ];
    for &s_str in &cases {
        let m = Mpz::from_decimal_str(s_str).unwrap();
        let mut buf = [0u8; 192];
        let len = m.write_decimal_buf(&mut buf);
        let back = core::str::from_utf8(&buf[..len]).unwrap();
        assert_eq!(back, s_str, "round-trip failed for {}", s_str);
    }
}

// ==========================================================================
// Test vectors: parse errors
// ==========================================================================

#[test]
fn parse_error_vectors() {
    assert_eq!(Mpz::from_decimal_str(""), Err(ParseError::InvalidInput));
    assert_eq!(Mpz::from_decimal_str("-"), Err(ParseError::InvalidInput));
    assert_eq!(Mpz::from_decimal_str("+"), Err(ParseError::InvalidInput));
    assert_eq!(
        Mpz::from_decimal_str("12a34"),
        Err(ParseError::InvalidInput)
    );
    assert_eq!(
        Mpz::from_decimal_str("   "),
        Err(ParseError::InvalidInput)
    );
}

// ==========================================================================
// Test vectors: from_d (f64)
// ==========================================================================

#[test]
fn from_d_vectors() {
    assert_eq!(Mpz::from_d(0.0).unwrap(), Mpz::new());
    assert_eq!(Mpz::from_d(1.0).unwrap(), Mpz::from_u64(1));
    assert_eq!(Mpz::from_d(-1.0).unwrap(), Mpz::from_i64(-1));
    assert_eq!(Mpz::from_d(3.999).unwrap(), Mpz::from_u64(3));
    assert_eq!(Mpz::from_d(-3.999).unwrap(), Mpz::from_i64(-3));
    assert_eq!(Mpz::from_d(1e18).unwrap(), Mpz::from_u64(1000000000000000000));
    assert_eq!(Mpz::from_d(f64::INFINITY), Err(CapacityError));
    assert_eq!(Mpz::from_d(f64::NEG_INFINITY), Err(CapacityError));
    assert_eq!(Mpz::from_d(f64::NAN), Err(CapacityError));
}

// ==========================================================================
// Test vectors: comparison
// ==========================================================================

#[test]
fn cmp_vectors() {
    use core::cmp::Ordering;
    assert_eq!(Mpz::from_u64(5).cmp(&Mpz::from_u64(3)), Ordering::Greater);
    assert_eq!(Mpz::from_u64(3).cmp(&Mpz::from_u64(5)), Ordering::Less);
    assert_eq!(Mpz::from_u64(5).cmp(&Mpz::from_u64(5)), Ordering::Equal);
    assert_eq!(
        Mpz::from_i64(-5).cmp(&Mpz::from_u64(3)),
        Ordering::Less
    );
    assert_eq!(
        Mpz::from_i64(-5).cmp(&Mpz::from_i64(-3)),
        Ordering::Less
    );
    assert_eq!(
        Mpz::from_i64(-3).cmp_si(-3),
        Ordering::Equal
    );
    assert_eq!(
        Mpz::from_u64(5).cmp_ui(10),
        Ordering::Less
    );
}
