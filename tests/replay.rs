//! Deterministic replay tests for known failure cases.
//!
//! When a property test or fuzz target finds a failure, the failure case
//! is added here as a concrete test so it never reappears.

use gmp_rs::Mpz;

fn check_add(a: &str, b: &str, expected: &str) {
    let a = Mpz::from_decimal_str(a).unwrap_or_else(|_| panic!("invalid a: {a}"));
    let b = Mpz::from_decimal_str(b).unwrap_or_else(|_| panic!("invalid b: {b}"));
    match a.try_add(&b) {
        Ok(r) => assert_eq!(r.to_string(), expected, "add({a}, {b})"),
        Err(_) => assert_eq!("CapacityError", expected, "add({a}, {b})"),
    }
}

fn check_sub(a: &str, b: &str, expected: &str) {
    let a = Mpz::from_decimal_str(a).unwrap();
    let b = Mpz::from_decimal_str(b).unwrap();
    match a.try_sub(&b) {
        Ok(r) => assert_eq!(r.to_string(), expected, "sub({a}, {b})"),
        Err(_) => assert_eq!("CapacityError", expected, "sub({a}, {b})"),
    }
}

fn check_mul(a: &str, b: &str, expected: &str) {
    let a = Mpz::from_decimal_str(a).unwrap();
    let b = Mpz::from_decimal_str(b).unwrap();
    match a.try_mul(&b) {
        Ok(r) => assert_eq!(r.to_string(), expected, "mul({a}, {b})"),
        Err(_) => assert_eq!("CapacityError", expected, "mul({a}, {b})"),
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Basic operations
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn replay_add_basic() {
    check_add("0", "0", "0");
    check_add("0", "1", "1");
    check_add("1", "0", "1");
    check_add("0", "-1", "-1");
    check_add("1", "1", "2");
    check_add("100", "200", "300");
    check_add("-5", "-3", "-8");
    check_add("-100", "50", "-50");
    check_add("50", "-100", "-50");
    check_add("-50", "100", "50");
}

#[test]
fn replay_sub_basic() {
    check_sub("5", "3", "2");
    check_sub("3", "5", "-2");
    check_sub("0", "5", "-5");
    check_sub("5", "0", "5");
    check_sub("-5", "-3", "-2");
    check_sub("-5", "3", "-8");
}

#[test]
fn replay_mul_basic() {
    check_mul("5", "3", "15");
    check_mul("-5", "3", "-15");
    check_mul("5", "-3", "-15");
    check_mul("-5", "-3", "15");
    check_mul("0", "5", "0");
    check_mul("5", "0", "0");
    check_mul("1000000000", "1000000000", "1000000000000000000");
}

#[test]
fn replay_add_large() {
    check_add("9999999999999999999", "1", "10000000000000000000");
    check_add(
        "12345678901234567890",
        "98765432109876543210",
        "111111111011111111100",
    );
}

#[test]
fn replay_mul_large() {
    check_mul(
        "1000000000000000000",
        "1000000000000000000",
        "1000000000000000000000000000000000000",
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Sign edge cases (programmatic, no giant strings)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn replay_sign_edge_cases() {
    // (-a) + (-b) == -(a + b)
    assert_eq!(
        Mpz::from_i64(-7).try_add(&Mpz::from_i64(-3)).unwrap(),
        Mpz::from_i64(7)
            .try_add(&Mpz::from_i64(3))
            .unwrap()
            .neg_to(),
    );

    // (-a) - (-b) == b - a
    assert_eq!(
        Mpz::from_i64(-7).try_sub(&Mpz::from_i64(-3)).unwrap(),
        Mpz::from_i64(3).try_sub(&Mpz::from_i64(7)).unwrap(),
    );

    // (-a) * (-b) == a * b
    assert_eq!(
        Mpz::from_i64(-7).try_mul(&Mpz::from_i64(-3)).unwrap(),
        Mpz::from_i64(7).try_mul(&Mpz::from_i64(3)).unwrap(),
    );

    // (-a) * b == -(a * b)
    assert_eq!(
        Mpz::from_i64(-7).try_mul(&Mpz::from_i64(3)).unwrap(),
        Mpz::from_i64(7)
            .try_mul(&Mpz::from_i64(3))
            .unwrap()
            .neg_to(),
    );
}

// ─────────────────────────────────────────────────────────────────────────
// Capacity boundary tests (using try_import)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn replay_capacity_boundary() {
    // Build max value: 2^511 * 2 - 1 = 2^512 - 1, but 2^511*2 = 2^512 overflows.
    // Instead, build via repeated shifts and adds up to the 512-bit boundary.

    // 2^511 should fit (512 bits)
    let two_pow_511 = Mpz::from_u64(1).try_mul_2exp(511).unwrap();
    assert_eq!(two_pow_511.sizeinbase2(), 512);

    // Build max = 2^512 - 1 by: (2^511 - 1)*2 + 1 = 2^512 - 1
    // 2^511 - 1 fits in 511 bits
    let max_minus_511 = two_pow_511.try_sub(&Mpz::from_u64(1)).unwrap();
    // 2^511 - 1 + 2^511 = 2^512 - 1 (but this overflows because 2^512 needs 9 limbs)
    // Actually let's just test the boundary directly:
    // A value with bit 511 set and all lower bits set = 2^512 - 1, which is max.
    // But 2^511 + (2^511 - 1) = 2^512 - 1... but 2^511 + 2^511 = 2^512 which overflows.
    // So we CAN'T construct the full max with shifts alone.

    // Instead, test: max-1 fits, and max-1 + 1 overflows if max-1 has 511 bits set
    // Actually, the safest test: a value known to exercise the boundary

    // 2^255 should easily fit (256 bits)
    let two_pow_255 = Mpz::from_u64(1).try_mul_2exp(255).unwrap();
    assert!(two_pow_255.sizeinbase2() <= 512);

    // 2^512 should overflow (requires 9 limbs)
    assert!(Mpz::from_u64(1).try_mul_2exp(512).is_err());

    // Multiplication of two 256-bit values: 2^256 * 2^256 = 2^512 which overflows
    let two_pow_256 = Mpz::from_u64(1).try_mul_2exp(256).unwrap();
    assert!(
        two_pow_256.try_mul(&two_pow_256).is_err(),
        "2^256 * 2^256 should overflow"
    );

    // But 2^255 * 2^255 = 2^510 which fits
    assert!(
        two_pow_255.try_mul(&two_pow_255).is_ok(),
        "2^255 * 2^255 should fit"
    );
}

#[test]
fn replay_negative_large_values() {
    // Build a large negative value and verify it works
    let large_positive = Mpz::from_u64(1).try_mul_2exp(510).unwrap(); // 2^510, fits
    let large_negative = large_positive.neg_to();

    assert_eq!(large_negative.sgn(), -1);
    // (-large) + large = 0
    assert_eq!(large_negative.try_add(&large_positive).unwrap(), Mpz::new());
}
