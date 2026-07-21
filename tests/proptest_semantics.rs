//! Property-based tests for algebraic invariants using proptest.
//!
//! Run with:  cargo test --features "std" --test proptest_semantics
//!
//! Requires:  proptest = "1" (in dev-dependencies)

#![cfg(feature = "std")]

use gmp_rs::Mpz;
use proptest::prelude::*;

// ─────────────────────────────────────────────────────────────────────────
// proptest generators for Mpz
// ─────────────────────────────────────────────────────────────────────────

/// Generate an Mpz from a random i64 (always fits in capacity, single limb).
fn any_mpz_i64() -> impl Strategy<Value = Mpz> {
    any::<i64>().prop_map(|v| Mpz::from_i64(v))
}

/// Generate an Mpz from a random i128 (always fits, up to 2 limbs).
fn any_mpz_i128() -> impl Strategy<Value = Mpz> {
    any::<i128>().prop_map(|v| Mpz::from_i128(v))
}

/// Generate a small decimal string that fits within capacity.
fn any_mpz_str() -> impl Strategy<Value = String> {
    // Generate strings that definitely fit: up to ~40 digits
    prop::collection::vec(proptest::char::range('0', '9'), 1..20)
        .prop_map(|chars| {
            let s: String = chars.into_iter().collect();
            // Ensure no leading zeros (except for the value "0")
            if s.len() > 1 && s.starts_with('0') {
                "1".to_string() + &s[1..]
            } else {
                s
            }
        })
        .prop_flat_map(|digits| {
            prop_oneof![
                Just(digits.clone()),
                Just(format!("-{}", digits)),
                Just(format!("+{}", digits)),
            ]
        })
}

/// Generate a random i128-sized Mpz via decimal string roundtrip.
fn any_mpz_parseable() -> impl Strategy<Value = Mpz> {
    any_mpz_str().prop_filter_map("parseable decimal", |s| Mpz::from_decimal_str(&s).ok())
}

/// Generate a random Mpz from limb-level generation.
/// Uses try_import to construct from byte arrays.
fn any_mpz_limbs() -> impl Strategy<Value = Mpz> {
    prop::collection::vec(any::<u8>(), 0..65).prop_filter_map("construct via try_import", |bytes| {
        if bytes.is_empty() {
            return Some(Mpz::new());
        }
        // Use try_import with native-endian u64 chunks
        Mpz::try_import(
            bytes.len(),
            gmp_rs::Endian::Native,
            1,
            gmp_rs::Endian::Native,
            &bytes,
        )
        .ok()
    })
}

// ─────────────────────────────────────────────────────────────────────────
// Algebraic property tests
// ─────────────────────────────────────────────────────────────────────────

proptest! {
    // ── Add commutativity ──
    #[test]
    fn add_commutative_i64(a in any_mpz_i64(), b in any_mpz_i64()) {
        let ab = a.try_add(&b);
        let ba = b.try_add(&a);
        match (ab, ba) {
            (Ok(x), Ok(y)) => assert_eq!(x, y, "add commutativity failed"),
            (Err(_), Err(_)) => {}, // both overflow: acceptable
            _ => panic!("asymmetric overflow behaviour"),
        }
    }

    #[test]
    fn add_commutative_limb(a in any_mpz_limbs(), b in any_mpz_limbs()) {
        let ab = a.try_add(&b);
        let ba = b.try_add(&a);
        match (ab, ba) {
            (Ok(x), Ok(y)) => assert_eq!(x, y, "add commutativity (limb) failed"),
            (Err(_), Err(_)) => {},
            _ => panic!("asymmetric overflow behaviour"),
        }
    }

    // ── Add associativity ──
    #[test]
    fn add_associative_i64(a in any_mpz_i64(), b in any_mpz_i64(), c in any_mpz_i64()) {
        let left = a.try_add(&b).and_then(|ab| ab.try_add(&c));
        let right = b.try_add(&c).and_then(|bc| a.try_add(&bc));
        match (left, right) {
            (Ok(x), Ok(y)) => assert_eq!(x, y, "add associativity failed"),
            (Err(_), Err(_)) => {},
            _ => {}, // capacity differences acceptable
        }
    }

    #[test]
    fn add_associative_limb(a in any_mpz_limbs(), b in any_mpz_limbs(), c in any_mpz_limbs()) {
        let left = a.try_add(&b).and_then(|ab| ab.try_add(&c));
        let right = b.try_add(&c).and_then(|bc| a.try_add(&bc));
        match (left, right) {
            (Ok(x), Ok(y)) => assert_eq!(x, y),
            (Err(_), Err(_)) => {},
            _ => {},
        }
    }

    // ── Add identity ──
    #[test]
    fn add_identity(a in any_mpz_i64()) {
        let zero = Mpz::new();
        assert_eq!(a.try_add(&zero).unwrap(), a, "a + 0 == a");
        assert_eq!(zero.try_add(&a).unwrap(), a, "0 + a == a");
    }

    // ── Sub inverse ──
    #[test]
    fn sub_inverse(a in any_mpz_i64()) {
        assert_eq!(a.try_sub(&a).unwrap(), Mpz::new(), "a - a == 0");
        let neg = a.neg_to();
        assert_eq!(a.try_add(&neg).unwrap(), Mpz::new(), "a + (-a) == 0");
    }

    // ── Negation ──
    #[test]
    fn neg_involution(a in any_mpz_i64()) {
        assert_eq!(a.neg_to().neg_to(), a, "-(-a) == a");
    }

    // ── Mul commutativity ──
    #[test]
    fn mul_commutative_i64(a in any_mpz_i64(), b in any_mpz_i64()) {
        let ab = a.try_mul(&b);
        let ba = b.try_mul(&a);
        match (ab, ba) {
            (Ok(x), Ok(y)) => assert_eq!(x, y, "mul commutativity failed"),
            (Err(_), Err(_)) => {},
            _ => panic!("asymmetric overflow behaviour"),
        }
    }

    #[test]
    fn mul_commutative_limb(a in any_mpz_limbs(), b in any_mpz_limbs()) {
        let ab = a.try_mul(&b);
        let ba = b.try_mul(&a);
        match (ab, ba) {
            (Ok(x), Ok(y)) => assert_eq!(x, y),
            (Err(_), Err(_)) => {},
            _ => panic!("asymmetric overflow behaviour"),
        }
    }

    // ── Mul associativity ──
    #[test]
    fn mul_associative_i64(a in any_mpz_i64(), b in any_mpz_i64(), c in any_mpz_i64()) {
        let left = a.try_mul(&b).and_then(|ab| ab.try_mul(&c));
        let right = b.try_mul(&c).and_then(|bc| a.try_mul(&bc));
        match (left, right) {
            (Ok(x), Ok(y)) => assert_eq!(x, y),
            (Err(_), Err(_)) => {},
            _ => {},
        }
    }

    // ── Distributivity ──
    #[test]
    fn distributivity_i64(a in any_mpz_i64(), b in any_mpz_i64(), c in any_mpz_i64()) {
        let left = b.try_add(&c).and_then(|bc| a.try_mul(&bc));
        let right = a.try_mul(&b).and_then(|ab| {
            a.try_mul(&c).and_then(|ac| ab.try_add(&ac))
        });
        match (left, right) {
            (Ok(x), Ok(y)) => assert_eq!(x, y, "distributivity failed"),
            (Err(_), Err(_)) => {},
            _ => {}, // capacity differences
        }
    }

    // ── Ordering: if a < b then a + c < b + c ──
    #[test]
    fn ordering_add_monotonic(a in any_mpz_i64(), b in any_mpz_i64(), c in any_mpz_i64()) {
        if a.cmp(&b) == core::cmp::Ordering::Less {
            if let (Ok(ac), Ok(bc)) = (a.try_add(&c), b.try_add(&c)) {
                assert!(ac.cmp(&bc) == core::cmp::Ordering::Less
                        || ac.cmp(&bc) == core::cmp::Ordering::Equal,
                    "a < b ⇒ a + c ≤ b + c failed");
            }
        }
    }

    // ── Ordering: if c > 0 and a < b then a*c < b*c ──
    #[test]
    fn ordering_mul_monotonic(a in any_mpz_i64(), b in any_mpz_i64(), c in any_mpz_i64()) {
        if a.cmp(&b) == core::cmp::Ordering::Less && c.sgn() > 0 {
            if let (Ok(ac), Ok(bc)) = (a.try_mul(&c), b.try_mul(&c)) {
                assert_eq!(ac.cmp(&bc), core::cmp::Ordering::Less,
                    "a < b, c > 0 ⇒ a*c < b*c failed");
            }
        }
    }

    // ── Parse/format roundtrip ──
    #[test]
    fn parse_format_roundtrip_i128(n in any::<i128>()) {
        let s = n.to_string();
        if let Ok(m) = Mpz::from_decimal_str(&s) {
            // Roundtrip must produce the original string (canonical form)
            let mut buf = [0u8; 200];
            let len = m.write_decimal_buf(&mut buf);
            let formatted = core::str::from_utf8(&buf[..len]).unwrap();
            assert_eq!(formatted, s, "parse/format roundtrip: {s}");
        }
    }

    #[test]
    fn parse_format_roundtrip_str(s in any_mpz_str()) {
        let trimmed = s.trim_start_matches('+');
        if let Ok(m) = Mpz::from_decimal_str(&s) {
            let mut buf = [0u8; 200];
            let len = m.write_decimal_buf(&mut buf);
            let formatted = core::str::from_utf8(&buf[..len]).unwrap();
            // The canonical form should match (without leading +)
            let expected = if trimmed.starts_with('-') { trimmed } else { trimmed.trim_start_matches('+') };
            // Don't compare "0" with "-0" or "+0"
            if expected != "-0" && expected != "+0" {
                assert_eq!(formatted, expected, "parse/format roundtrip: '{s}' → '{formatted}'");
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Explicit edge-case tests (not generated)
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn zero_is_always_positive() {
    let z = Mpz::new();
    assert_eq!(z.sgn(), 0);
    assert_eq!(z.to_string(), "0");
    // Negating zero stays zero
    assert_eq!(z.neg_to().sgn(), 0);
    assert_eq!(z.neg_to().to_string(), "0");
}

#[test]
fn sign_after_operations() {
    // Positive * Positive = Positive
    assert!(Mpz::from_u64(3).try_mul(&Mpz::from_u64(4)).unwrap().sgn() > 0);
    // Positive * Negative = Negative
    assert!(Mpz::from_u64(3).try_mul(&Mpz::from_i64(-4)).unwrap().sgn() < 0);
    // Negative * Negative = Positive
    assert!(Mpz::from_i64(-3).try_mul(&Mpz::from_i64(-4)).unwrap().sgn() > 0);
}
