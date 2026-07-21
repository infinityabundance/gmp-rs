//! Boundary and stress tests — explicitly hammer the edges of the capacity.
//!
//! Run with:  cargo test --features std --test boundary_tests

#![cfg(feature = "std")]

use gmp_rs::{Mpz, MAX_BITS, MAX_DECIMAL_DIGITS};

// ─────────────────────────────────────────────────────────────────────────
// Capacity boundary: 2^511 is the highest power of 2 that fits (512 bits)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn two_pow_511_fits() {
    let v = Mpz::from_u64(1).try_mul_2exp(511).unwrap();
    assert_eq!(v.sizeinbase2(), 512);
}

#[test]
fn two_pow_512_overflows() {
    assert!(Mpz::from_u64(1).try_mul_2exp(512).is_err());
}

#[test]
fn max_mul_overflows() {
    // 2^255 * 2^255 = 2^510 which fits (511 bits)
    let v = Mpz::from_u64(1).try_mul_2exp(255).unwrap();
    assert!(v.try_mul(&v).is_ok());

    // 2^256 * 2^256 = 2^512 which overflows (needs 513 bits = 9 limbs)
    let c = Mpz::from_u64(1).try_mul_2exp(256).unwrap();
    assert!(c.try_mul(&c).is_err());
}

#[test]
fn max_add_near_boundary() {
    // 2^511 is 1 followed by 511 zeros. Adding any positive number stays in range.
    let v = Mpz::from_u64(1).try_mul_2exp(511).unwrap();
    let one = Mpz::from_u64(1);
    assert!(v.try_add(&one).is_ok()); // 2^511 + 1 fits (still 512 bits)
}

// ─────────────────────────────────────────────────────────────────────────
// Near-zero tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn zero_operations() {
    let z = Mpz::new();
    assert_eq!(z.sgn(), 0);
    assert_eq!(z.is_zero(), true);
    assert_eq!(z.to_string(), "0");

    assert_eq!(z.try_add(&z).unwrap(), z);
    assert_eq!(z.try_sub(&z).unwrap(), z);
    assert_eq!(z.try_mul(&z).unwrap(), z);

    let a = Mpz::from_i64(42);
    assert_eq!(z.try_mul(&a).unwrap(), z);
    assert_eq!(a.try_mul(&z).unwrap(), z);
    assert_eq!(a.try_add(&z).unwrap(), a);
    assert_eq!(z.try_add(&a).unwrap(), a);
    assert_eq!(a.try_sub(&z).unwrap(), a);
    assert_eq!(z.try_sub(&a).unwrap(), a.neg_to());
    assert_eq!(z.neg_to(), z);
    assert_eq!(z.neg_to().sgn(), 0);
}

#[test]
fn zero_string_form() {
    for s in &["0", "-0", "+0", "000"] {
        let m = Mpz::from_decimal_str(s).unwrap();
        assert_eq!(m.sgn(), 0, "from_decimal_str({s}) should be zero");
        assert_eq!(m.to_string(), "0", "from_decimal_str({s}) formats as 0");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Sign edge tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn sign_preservation_pos_neg() {
    let pos = Mpz::from_u64(42);
    let neg = Mpz::from_i64(-42);

    assert_eq!(
        neg.try_add(&Mpz::from_i64(-58)).unwrap(),
        Mpz::from_i64(100).neg_to(),
    );

    assert_eq!(
        neg.try_sub(&Mpz::from_i64(-58)).unwrap(),
        Mpz::from_i64(58).try_sub(&Mpz::from_i64(42)).unwrap(),
    );

    assert_eq!(
        neg.try_mul(&Mpz::from_i64(-58)).unwrap(),
        Mpz::from_i64(42).try_mul(&Mpz::from_i64(58)).unwrap(),
    );

    assert_eq!(
        neg.try_mul(&Mpz::from_i64(58)).unwrap(),
        Mpz::from_i64(42)
            .try_mul(&Mpz::from_i64(58))
            .unwrap()
            .neg_to(),
    );

    assert_eq!(
        pos.try_mul(&Mpz::from_i64(-58)).unwrap(),
        pos.try_mul(&Mpz::from_u64(58)).unwrap().neg_to(),
    );
}

#[test]
fn division_sign_edge_cases() {
    let (q, r) = Mpz::from_i64(-100).tdiv_qr(&Mpz::from_i64(30));
    assert_eq!(q.to_string(), "-3");
    assert_eq!(r.to_string(), "-10");

    let (q, r) = Mpz::from_i64(100).tdiv_qr(&Mpz::from_i64(-30));
    assert_eq!(q.to_string(), "-3");
    assert_eq!(r.to_string(), "10");

    let (q, r) = Mpz::from_i64(-100).tdiv_qr(&Mpz::from_i64(-30));
    assert_eq!(q.to_string(), "3");
    assert_eq!(r.to_string(), "-10");
}

// ─────────────────────────────────────────────────────────────────────────
// Stress: repeated operations
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn repeated_addition() {
    let mut sum = Mpz::new();
    for i in 1..=100 {
        sum = sum.try_add(&Mpz::from_u64(i)).unwrap();
    }
    assert_eq!(sum.to_string(), "5050");
}

#[test]
fn repeated_multiplication() {
    let mut prod = Mpz::from_u64(1);
    for _ in 0..10 {
        prod = prod.try_mul(&Mpz::from_u64(2)).unwrap();
    }
    assert_eq!(prod.to_string(), "1024");

    let mut fact = Mpz::from_u64(1);
    for i in 2..=10 {
        fact = fact.try_mul(&Mpz::from_u64(i)).unwrap();
    }
    assert_eq!(fact.to_string(), "3628800");
}

#[test]
fn oscillating_add_sub() {
    let mut val = Mpz::from_u64(1000);
    let inc = Mpz::from_u64(500);
    for _ in 0..50 {
        val = val.try_add(&inc).unwrap();
        val = val.try_sub(&inc).unwrap();
    }
    assert_eq!(val.to_string(), "1000");
}

// ─────────────────────────────────────────────────────────────────────────
// Capacity constants verification
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn constants_are_consistent() {
    assert_eq!(MAX_BITS, gmp_rs::MPZ_MAX_LIMBS * 64);
    assert_eq!(gmp_rs::LIMBS, gmp_rs::MPZ_MAX_LIMBS);
}

#[test]
fn max_decimal_digits_is_accurate() {
    // 10^154 should fit (it's ~511 bits)
    let s = "1".to_string() + &"0".repeat(MAX_DECIMAL_DIGITS);
    let m = Mpz::from_decimal_str(&s).unwrap();
    assert!(m.sizeinbase2() <= MAX_BITS);

    // 10^155 should overflow (~514 bits)
    let s2 = "1".to_string() + &"0".repeat(MAX_DECIMAL_DIGITS + 1);
    assert!(Mpz::from_decimal_str(&s2).is_err());
}
