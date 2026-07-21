//! Capacity boundary tests — hammer the edges of the fixed-capacity model.
//!
//! Run with:  cargo test --test capacity_boundaries
//!
//! These tests prove:
//!   - No silent truncation or wraparound at the capacity boundary
//!   - Carry/borrow propagation is correct at every limb boundary
//!   - Decimal parsing uses the same overflow semantics as arithmetic
//!   - Zero normalization is correct
//!   - Sign handling is correct at the capacity extremes

use gmp_rs::{CapacityError, Mpz, LIMBS, MAX_BITS};

// ─────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────

/// Construct the maximum representable value: all limbs = u64::MAX, sign = +1.
fn max_mpz() -> Mpz {
    let limbs = [u64::MAX; LIMBS];
    Mpz::from_limbs_checked(1, &limbs).expect("max value should be valid")
}

/// Maximum value minus 1 in the least significant limb.
fn near_max_mpz() -> Mpz {
    let mut limbs = [u64::MAX; LIMBS];
    limbs[0] = u64::MAX - 1;
    Mpz::from_limbs_checked(1, &limbs).expect("near-max value should be valid")
}

/// A value with only the highest limb set (bit 511, value = 2^511).
fn high_limb_only() -> Mpz {
    let mut limbs = [0u64; LIMBS];
    limbs[LIMBS - 1] = 1u64 << 63;
    Mpz::from_limbs_checked(1, &limbs).expect("high-limb value should be valid")
}

/// Maximum negative value: all limbs = u64::MAX, sign = -1.
fn neg_max_mpz() -> Mpz {
    let limbs = [u64::MAX; LIMBS];
    Mpz::from_limbs_checked(-1, &limbs).expect("neg max value should be valid")
}

// ─────────────────────────────────────────────────────────────────────────
// 2. Addition overflow boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn add_overflow_boundary() {
    let max = max_mpz();
    let one = Mpz::from_i64(1);

    // max + 1 MUST overflow
    assert!(max.try_add(&one).is_err(), "max + 1 must overflow");

    // near_max + 1 MUST succeed
    let near = near_max_mpz();
    let res = near.try_add(&one).expect("near_max + 1 should fit");
    assert_eq!(
        res.to_string(),
        max.to_string(),
        "near_max + 1 should equal max"
    );
}

#[test]
fn add_max_with_negative_one() {
    let max = max_mpz();
    let neg_one = Mpz::from_i64(-1);
    let res = max.try_add(&neg_one).expect("max + (-1) should fit");
    assert_eq!(res, near_max_mpz(), "max + (-1) should equal near_max");
}

#[test]
fn add_two_large_values_barely_fit() {
    // 2^255 + 2^255 = 2^256 which fits
    let a = Mpz::from_u64(1).try_mul_2exp(255).unwrap();
    let b = Mpz::from_u64(1).try_mul_2exp(255).unwrap();
    let res = a.try_add(&b).expect("2^255 + 2^255 should fit");
    let expected = Mpz::from_u64(1).try_mul_2exp(256).unwrap();
    assert_eq!(res, expected);
}

// ─────────────────────────────────────────────────────────────────────────
// 3. Multiplication overflow boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn mul_overflow_boundary() {
    let max = max_mpz();
    let two = Mpz::from_i64(2);
    // max * 2 MUST overflow
    assert!(max.try_mul(&two).is_err(), "max * 2 must overflow");

    // near_max * 1 MUST fit
    let near = near_max_mpz();
    let one = Mpz::from_i64(1);
    let res = near.try_mul(&one).expect("near_max * 1 should fit");
    assert_eq!(res.to_string(), near.to_string());
}

#[test]
fn mul_high_limb_overflow() {
    // A value with only the highest limb set times 2 overflows
    let high = high_limb_only();
    let two = Mpz::from_i64(2);
    assert!(
        high.try_mul(&two).is_err(),
        "high-limb-only * 2 must overflow (carry into non-existent limb)"
    );

    // high * 1 fits
    let one = Mpz::from_i64(1);
    let res = high.try_mul(&one).unwrap();
    assert_eq!(res, high);
}

#[test]
fn mul_near_overflow_patterns() {
    let max = max_mpz();

    // max * 0 must be 0
    let zero = Mpz::new();
    assert_eq!(max.try_mul(&zero).unwrap(), zero);

    // max * 1 must be max
    let one = Mpz::from_i64(1);
    assert_eq!(max.try_mul(&one).unwrap(), max);

    // max * -1 must be -max
    let neg_one = Mpz::from_i64(-1);
    assert_eq!(max.try_mul(&neg_one).unwrap(), neg_max_mpz());
}

// ─────────────────────────────────────────────────────────────────────────
// 4. Decimal parsing overflow boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn parse_overflow_boundary() {
    let max = max_mpz();
    let max_str = max.to_string();

    // Parsing max MUST succeed
    let parsed = Mpz::from_decimal_str(&max_str).expect("max decimal must parse");
    assert_eq!(parsed.to_string(), max_str);

    // Construct max + 1 as a decimal string using manual addition
    let big_plus_one = {
        let mut digits: Vec<u8> = max_str.bytes().map(|b| b - b'0').collect();
        let mut carry = 1;
        for d in digits.iter_mut().rev() {
            let x = *d + carry;
            *d = x % 10;
            carry = x / 10;
        }
        if carry > 0 {
            digits.insert(0, carry);
        }
        digits
            .into_iter()
            .map(|d| (d + b'0') as char)
            .collect::<String>()
    };

    // Parsing max+1 MUST overflow
    assert!(
        Mpz::from_decimal_str(&big_plus_one).is_err(),
        "max+1 decimal string must overflow on parse"
    );
}

#[test]
fn parse_large_power_of_ten() {
    // 10^154 should fit (~511 bits)
    let s = "1".to_string() + &"0".repeat(154);
    assert!(Mpz::from_decimal_str(&s).is_ok(), "10^154 should parse");

    // 10^155 should overflow (~514 bits)
    let s2 = "1".to_string() + &"0".repeat(155);
    assert!(
        Mpz::from_decimal_str(&s2).is_err(),
        "10^155 should overflow on parse"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// 5. Sign-edge boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn sign_edge_boundaries() {
    let max = max_mpz();
    let neg_max = neg_max_mpz();

    // (-max) + max MUST be zero
    let zero = neg_max.try_add(&max).expect("(-max) + max should fit");
    assert_eq!(zero.to_string(), "0");
    assert_eq!(zero.sgn(), 0);

    // (-max) - max MUST overflow (result = -2*max, which exceeds capacity)
    assert!(neg_max.try_sub(&max).is_err(), "(-max) - max must overflow");
}

#[test]
fn sign_edge_add_neg_max() {
    let max = max_mpz();
    let neg_one = Mpz::from_i64(-1);

    // max + (-1) = max - 1, which fits
    let res = max.try_add(&neg_one).unwrap();
    assert_eq!(res, near_max_mpz());

    // (-max) + 1 = -(max - 1), which fits
    let neg_max = neg_max_mpz();
    let res = neg_max.try_add(&Mpz::from_i64(1)).unwrap();
    let expected = near_max_mpz().neg_to();
    assert_eq!(res, expected);
}

// ─────────────────────────────────────────────────────────────────────────
// 6. Zero-normalization boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn zero_normalization() {
    let zero = Mpz::from_i64(0);
    assert_eq!(zero.to_string(), "0");

    // Construct "negative zero" manually — must normalize to positive zero
    let limbs = [0u64; LIMBS]; // all zero limbs
    let neg_zero = Mpz::from_limbs_checked(-1, &limbs);
    // from_limbs_checked should reject this (sign=-1 with zero magnitude)
    assert!(neg_zero.is_none(), "negative zero should be rejected");

    // from_limbs_checked with sign=0 and zero limbs should succeed
    let actual_zero = Mpz::from_limbs_checked(0, &limbs);
    assert!(actual_zero.is_some());
    assert_eq!(actual_zero.unwrap().to_string(), "0");
}

#[test]
fn zero_is_always_canonical() {
    for s in &["0", "-0", "+0", "000"] {
        let m = Mpz::from_decimal_str(s).unwrap();
        assert_eq!(m.sgn(), 0);
        assert_eq!(m.to_string(), "0");
        // Operations on zero produce canonical zero
        assert_eq!(m.neg_to().sgn(), 0);
        assert_eq!(m.neg_to().to_string(), "0");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// 7. Limb-carry boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn full_limb_carry_propagation() {
    // Build a value where every limb is u64::MAX, then add 1.
    // This should overflow (since all 8 limbs carry).
    let max = max_mpz();
    assert!(
        max.try_add(&Mpz::from_i64(1)).is_err(),
        "all-limbs-carry add must overflow"
    );

    // Build a value with limb[0..6] = u64::MAX, limb[7] = 0,
    // then add 1.  The carry should propagate through all 7 limbs
    // and stop at limb[7] = 1.
    let mut limbs = [0u64; LIMBS];
    for i in 0..(LIMBS - 1) {
        limbs[i] = u64::MAX;
    }
    let val = Mpz::from_limbs_checked(1, &limbs).unwrap();
    let res = val.try_add(&Mpz::from_i64(1)).unwrap();
    // result should have limb[LIMBS-1] = 1, all others = 0
    let mut expected_limbs = [0u64; LIMBS];
    expected_limbs[LIMBS - 1] = 1;
    let expected = Mpz::from_limbs_checked(1, &expected_limbs).unwrap();
    assert_eq!(
        res, expected,
        "carry should propagate through 7 limbs to limb 7"
    );
}

#[test]
fn full_limb_borrow_propagation() {
    // Build a value with limb[0] = 0, all others = 0, then subtract 1.
    let val = Mpz::from_u64(0);
    let res = val.try_sub(&Mpz::from_i64(1)).unwrap();
    assert_eq!(res, Mpz::from_i64(-1));

    // Build 2^511 (only top bit set), subtract 1, get all lower bits set.
    let high = high_limb_only();
    let res = high.try_sub(&Mpz::from_i64(1)).unwrap();
    let mut expected_limbs = [u64::MAX; LIMBS];
    expected_limbs[LIMBS - 1] = (1u64 << 63) - 1;
    let expected = Mpz::from_limbs_checked(1, &expected_limbs).unwrap();
    assert_eq!(
        res, expected,
        "2^511 - 1 should have all 511 lower bits set"
    );
}

// ─────────────────────────────────────────────────────────────────────────
// 8. Subtraction overflow boundary tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn sub_overflow_boundary() {
    let max = max_mpz();
    let neg_max = neg_max_mpz();

    // max - (-1) = max + 1 which overflows
    let neg_one = Mpz::from_i64(-1);
    assert!(max.try_sub(&neg_one).is_err(), "max - (-1) must overflow");

    // (-max) - 1 = -(max + 1) which overflows
    assert!(
        neg_max.try_sub(&Mpz::from_i64(1)).is_err(),
        "(-max) - 1 must overflow"
    );

    // max - max = 0
    let zero = max.try_sub(&max).unwrap();
    assert_eq!(zero, Mpz::new());

    // (-max) - (-max) = 0
    let zero2 = neg_max.try_sub(&neg_max).unwrap();
    assert_eq!(zero2, Mpz::new());
}

// ─────────────────────────────────────────────────────────────────────────
// 9. Constant verification
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn constants_are_correct() {
    assert_eq!(MAX_BITS, LIMBS * 64);
    assert_eq!(LIMBS, 8);

    // Verify that from_limbs_checked validates correctly
    assert!(Mpz::from_limbs_checked(1, &[0; LIMBS]).is_none()); // zero mag with non-zero sign
    assert!(Mpz::from_limbs_checked(0, &[0; LIMBS]).is_some()); // proper zero
    assert!(Mpz::from_limbs_checked(1, &[42]).is_some()); // single limb
    assert!(Mpz::from_limbs_checked(-1, &[42]).is_some()); // negative single limb
    assert!(Mpz::from_limbs_checked(2, &[42]).is_none()); // invalid sign
    assert!(Mpz::from_limbs_checked(1, &[0; LIMBS + 1]).is_none()); // too many limbs
}
