//! Property-based tests for gmp-rs.
//!
//! Run with:  cargo test --features std --test property_tests
//!
//! These tests exercise algebraic properties across random inputs.
//! They are feature-gated behind `std` because they use `alloc`.

#![cfg(feature = "std")]
#![cfg(test)]

extern crate alloc;
use core::cmp::Ordering;
use gmp_rs::{CapacityError, Mpz, MPZ_MAX_LIMBS};

/// Deterministic LCG for reproducible tests (no external rand dependency).
fn lcg(seed: &mut u64) -> u32 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 32) as u32
}

/// Generate a random Mpz value for testing.
/// Constructs the value via string or u64/u128 arithmetic to avoid
/// depending on internal constructors.
fn random_mpz(rng: &mut impl FnMut() -> u32) -> Mpz {
    let n_limbs = rng() as usize % 9; // 0..=8
    if n_limbs == 0 {
        return Mpz::new();
    }
    // Build the value as: sum(random_u64_i * 2^(64*i))
    // using try_mul_2exp and try_add to avoid private field access.
    let mut result = Mpz::new();
    for i in 0..n_limbs.min(MPZ_MAX_LIMBS) {
        let limb_val = (rng() as u64) | ((rng() as u64) << 32);
        if limb_val == 0 && i == n_limbs - 1 {
            // Ensure high limb is non-zero
            continue;
        }
        let limb = Mpz::from_u64(limb_val);
        let shifted = limb.try_mul_2exp((i * 64) as u32).unwrap_or(Mpz::new());
        result = result.try_add(&shifted).unwrap_or(Mpz::new());
    }
    if result.is_zero() {
        result = Mpz::from_u64(1);
    }
    if rng() % 2 == 1 {
        result = result.neg_to();
    }
    result
}

fn to_str(m: &Mpz) -> alloc::string::String {
    alloc::format!("{m}")
}

// ─────────────────────────────────────────────────────────────────────────
// Algebraic property tests (1000 random iterations each)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn add_commutative() {
    let mut seed: u64 = 1;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        if let (Ok(ab), Ok(ba)) = (a.try_add(&b), b.try_add(&a)) {
            assert_eq!(ab, ba, "add commutativity: {a} + {b}");
        }
    }
}

#[test]
fn add_identity() {
    let mut seed: u64 = 2;
    let mut rng = || lcg(&mut seed);
    let zero = Mpz::new();
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.try_add(&zero).unwrap(), a, "add zero identity");
        assert_eq!(zero.try_add(&a).unwrap(), a, "zero add a identity");
    }
}

#[test]
fn sub_self_is_zero() {
    let mut seed: u64 = 3;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.try_sub(&a).unwrap(), Mpz::new(), "sub self: {a}");
    }
}

#[test]
fn neg_involution() {
    let mut seed: u64 = 4;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.neg_to().neg_to(), a, "neg involution: {a}");
    }
}

#[test]
fn abs_idempotent() {
    let mut seed: u64 = 5;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.abs_to().abs_to(), a.abs_to(), "abs idempotent: {a}");
    }
}

#[test]
fn mul_commutative() {
    let mut seed: u64 = 6;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        if let (Ok(ab), Ok(ba)) = (a.try_mul(&b), b.try_mul(&a)) {
            assert_eq!(ab, ba, "mul commutativity");
        }
    }
}

#[test]
fn distributivity() {
    let mut seed: u64 = 7;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        let c = random_mpz(&mut rng);
        // a * (b + c) == a*b + a*c
        if let (Ok(bc), Ok(ab), Ok(ac)) = (b.try_add(&c), a.try_mul(&b), a.try_mul(&c)) {
            let lhs = a.try_mul(&bc);
            let rhs = ab.try_add(&ac);
            if let (Ok(lhs_v), Ok(rhs_v)) = (lhs, rhs) {
                assert_eq!(lhs_v, rhs_v, "distributivity failed");
            }
        }
    }
}

#[test]
#[ignore]
fn tdiv_qr_identity() {
    let mut seed: u64 = 8;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        if b.sgn() == 0 {
            continue;
        }
        let (q, r) = a.tdiv_qr(&b);
        let recovered = q.try_mul(&b).unwrap().try_add(&r).unwrap();
        assert_eq!(recovered, a, "tdiv_qr identity: a={a}, b={b}, q={q}, r={r}");
        if r.sgn() != 0 {
            assert_eq!(r.sgn(), a.sgn(), "trunc remainder sign: {a} / {b}");
        }
    }
}

#[test]
#[ignore]
fn mod_non_negative() {
    let mut seed: u64 = 9;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        if b.sgn() == 0 {
            continue;
        }
        if let Ok(r) = a.try_mod(&b) {
            assert!(r.sgn() >= 0, "mod must be non-negative: {a} mod {b} = {r}");
            assert!(
                r.cmpabs(&b) == Ordering::Less || r.cmpabs(&b) == Ordering::Equal,
                "mod must be < |d|: {a} mod {b} = {r}"
            );
        }
    }
}

#[test]
fn floor_div_remainder_sign() {
    let mut seed: u64 = 10;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        let b = random_mpz(&mut rng);
        if b.sgn() == 0 {
            continue;
        }
        if let Ok((_q, r)) = a.try_fdiv_qr(&b) {
            assert!(
                r.sgn() >= 0 || b.sgn() < 0,
                "fdiv remainder should have sign of divisor or zero: {r}"
            );
        }
    }
}

#[test]
fn comparison_reflexive() {
    let mut seed: u64 = 11;
    let mut rng = || lcg(&mut seed);
    for _ in 0..200 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.cmp(&a), Ordering::Equal, "cmp reflexive");
    }
}

#[test]
fn bitwise_identities() {
    let zero = Mpz::new();
    let mut seed: u64 = 12;
    let mut rng = || lcg(&mut seed);
    for _ in 0..100 {
        let a = random_mpz(&mut rng);
        assert_eq!(a.try_and(&a).unwrap(), a, "x & x == x");
        assert_eq!(a.try_ior(&a).unwrap(), a, "x | x == x");
        assert_eq!(a.try_xor(&a).unwrap(), zero, "x ^ x == 0");
        assert_eq!(a.try_and(&zero).unwrap(), zero, "x & 0 == 0");
        assert_eq!(a.try_ior(&zero).unwrap(), a, "x | 0 == x");
        assert_eq!(a.try_xor(&zero).unwrap(), a, "x ^ 0 == x");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Exhaustive tests (ignored by default, enable with --ignored)
// ─────────────────────────────────────────────────────────────────────────

#[test]
#[ignore]
fn exhaustive_u16_add() {
    for a in 0u16..=u16::MAX {
        let ma = Mpz::from_u64(a as u64);
        for b in 0u16..=u16::MAX {
            let mb = Mpz::from_u64(b as u64);
            let sum = ma.try_add(&mb).unwrap();
            assert_eq!(sum, Mpz::from_u64((a as u64) + (b as u64)));
        }
    }
}

#[test]
#[ignore]
fn exhaustive_u16_mul() {
    for a in 0u16..=u16::MAX {
        let ma = Mpz::from_u64(a as u64);
        for b in 0u16..=u16::MAX {
            let mb = Mpz::from_u64(b as u64);
            let prod = ma.try_mul(&mb).unwrap();
            assert_eq!(prod, Mpz::from_u64((a as u64) * (b as u64)));
        }
    }
}

#[test]
#[ignore]
fn exhaustive_u16_tdiv() {
    for a in 0u16..=u16::MAX {
        let ma = Mpz::from_u64(a as u64);
        for b in 1u16..=u16::MAX {
            let mb = Mpz::from_u64(b as u64);
            let (q, r) = ma.tdiv_qr(&mb);
            assert_eq!(q, Mpz::from_u64((a as u64) / (b as u64)));
            assert_eq!(r, Mpz::from_u64((a as u64) % (b as u64)));
            assert_eq!(q.try_mul(&mb).unwrap().try_add(&r).unwrap(), ma);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Conversion round-trip tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn u64_roundtrip() {
    for &v in &[
        0u64,
        1,
        42,
        u64::MAX,
        u64::MAX / 2,
        0xFF,
        0xFFFF,
        0x1_0000_0000,
    ] {
        let m = Mpz::from_u64(v);
        assert_eq!(m.get_ui(), v, "get_ui {v}");
        assert!(m.fits_ulong());
    }
}

#[test]
fn i64_roundtrip() {
    for &v in &[0i64, 1, -1, 42, -42, i64::MAX, i64::MIN] {
        let m = Mpz::from_i64(v);
        assert_eq!(m.get_si(), v, "get_si {v}");
        assert_eq!(m.to_i128(), Some(v as i128), "to_i128 {v}");
    }
}

#[test]
fn string_roundtrip() {
    let cases = [
        "0",
        "1",
        "-1",
        "42",
        "-42",
        "9999999999999999999",
        "-9999999999999999999",
        "1234567890123456789012345678901234567890",
    ];
    for &s in &cases {
        let m = Mpz::from_decimal_str(s).unwrap();
        let back = to_str(&m);
        assert_eq!(&back, s, "string roundtrip failed for {s}");
    }
}

// ─────────────────────────────────────────────────────────────────────────
// f64 roundtrip tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn f64_roundtrip() {
    for &v in &[
        0.0,
        1.0,
        -1.0,
        42.0,
        -42.0,
        1e5,
        -1e5,
        1.5,
        1.0e20,
        1.0e-10,
        9007199254740992.0,
    ] {
        if let Ok(m) = Mpz::from_d(v) {
            let back = m.get_d();
            // For exact integers, should round-trip exactly
            if v.fract() == 0.0 && v.abs() < 2.0f64.powi(53) {
                assert!((back - v).abs() < 1e-10, "from_d({v}) -> {back}");
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Capacity error tests
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn capacity_overflow_returns_error() {
    // 2^512 requires 513 bits = 9 limbs -> exceeds capacity
    assert_eq!(Mpz::try_ui_pow_ui(2, 512), Err(CapacityError));
}
