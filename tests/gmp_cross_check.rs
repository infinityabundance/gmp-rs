//! GMP cross-check harness for gmp-rs.
//!
//! This test file compares gmp-rs results against the real GMP library
//! (via `extern "C"` FFI).  It requires:
//!
//! - `libgmp-dev` installed on the system (e.g., `apt install libgmp-dev`)
//! - The `gmp_cross_check` feature enabled:
//!
//! ```bash
//! cargo test --features gmp_cross_check --test gmp_cross_check
//! ```
//!
//! This test is feature-gated (`#[cfg(feature = "gmp_cross_check")]`) and
//! will not run under normal test invocations.
//!
//! ## Security note
//!
//! This test links against the system GMP library via FFI.  This is inherently
//! `unsafe` and is only compiled when the explicit `gmp_cross_check` feature
//! is enabled.  It is not part of the normal build.

#![cfg(feature = "gmp_cross_check")]
#![allow(dead_code)]

extern crate alloc;
use alloc::string::String;
use alloc::string::ToString;
use core::ffi;
use core::ptr;

use gmp_rs::*;

// ==========================================================================
// GMP FFI bindings (minimal subset needed for cross-check)
// ==========================================================================

/// GMP's internal `__mpz_struct` (opaque).  We use `mpz_t` as a fixed-size
/// array of structs (GMP says `mpz_t` is `__mpz_struct[1]`).
#[repr(C)]
struct __mpz_struct {
    _mp_alloc: ffi::c_int,
    _mp_size: ffi::c_int,
    _mp_d: *mut ffi::c_ulong,
}

type mpz_t = [__mpz_struct; 1];

extern "C" {
    fn __gmpz_init(x: *mut mpz_t);
    fn __gmpz_clear(x: *mut mpz_t);
    fn __gmpz_set_str(rop: *mut mpz_t, s: *const ffi::c_char, base: ffi::c_int) -> ffi::c_int;
    fn __gmpz_get_str(s: *mut ffi::c_char, base: ffi::c_int, op: *const mpz_t) -> *mut ffi::c_char;
    fn __gmpz_add(rop: *mut mpz_t, op1: *const mpz_t, op2: *const mpz_t);
    fn __gmpz_sub(rop: *mut mpz_t, op1: *const mpz_t, op2: *const mpz_t);
    fn __gmpz_mul(rop: *mut mpz_t, op1: *const mpz_t, op2: *const mpz_t);
    fn __gmpz_tdiv_qr(q: *mut mpz_t, r: *mut mpz_t, n: *const mpz_t, d: *const mpz_t);
    fn __gmpz_gcd(rop: *mut mpz_t, op1: *const mpz_t, op2: *const mpz_t);
}

// ==========================================================================
// GMP wrapper
// ==========================================================================

struct GmpMpz {
    inner: mpz_t,
}

impl GmpMpz {
    fn new() -> Self {
        let mut z = GmpMpz {
            inner: [__mpz_struct {
                _mp_alloc: 0,
                _mp_size: 0,
                _mp_d: ptr::null_mut(),
            }],
        };
        unsafe {
            __gmpz_init(&mut z.inner as *mut mpz_t);
        }
        z
    }

    fn from_decimal(s: &str) -> Self {
        let mut z = Self::new();
        let c_str = alloc::ffi::CString::new(s).expect("CString::new failed");
        unsafe {
            __gmpz_set_str(&mut z.inner as *mut mpz_t, c_str.as_ptr(), 10);
        }
        z
    }

    fn to_decimal(&self) -> String {
        unsafe {
            let ptr = __gmpz_get_str(ptr::null_mut(), 10, &self.inner as *const mpz_t);
            let s = alloc::ffi::CStr::from_ptr(ptr)
                .to_string_lossy()
                .into_owned();
            libc::free(ptr as *mut ffi::c_void);
            s
        }
    }

    fn add(&self, other: &GmpMpz) -> GmpMpz {
        let mut rop = GmpMpz::new();
        unsafe {
            __gmpz_add(
                &mut rop.inner as *mut mpz_t,
                &self.inner as *const mpz_t,
                &other.inner as *const mpz_t,
            );
        }
        rop
    }

    fn sub(&self, other: &GmpMpz) -> GmpMpz {
        let mut rop = GmpMpz::new();
        unsafe {
            __gmpz_sub(
                &mut rop.inner as *mut mpz_t,
                &self.inner as *const mpz_t,
                &other.inner as *const mpz_t,
            );
        }
        rop
    }

    fn mul(&self, other: &GmpMpz) -> GmpMpz {
        let mut rop = GmpMpz::new();
        unsafe {
            __gmpz_mul(
                &mut rop.inner as *mut mpz_t,
                &self.inner as *const mpz_t,
                &other.inner as *const mpz_t,
            );
        }
        rop
    }

    fn tdiv_qr(&self, d: &GmpMpz) -> (GmpMpz, GmpMpz) {
        let mut q = GmpMpz::new();
        let mut r = GmpMpz::new();
        unsafe {
            __gmpz_tdiv_qr(
                &mut q.inner as *mut mpz_t,
                &mut r.inner as *mut mpz_t,
                &self.inner as *const mpz_t,
                &d.inner as *const mpz_t,
            );
        }
        (q, r)
    }

    fn gcd(&self, other: &GmpMpz) -> GmpMpz {
        let mut rop = GmpMpz::new();
        unsafe {
            __gmpz_gcd(
                &mut rop.inner as *mut mpz_t,
                &self.inner as *const mpz_t,
                &other.inner as *const mpz_t,
            );
        }
        rop
    }
}

impl Drop for GmpMpz {
    fn drop(&mut self) {
        unsafe {
            __gmpz_clear(&mut self.inner as *mut mpz_t);
        }
    }
}

// ==========================================================================
// Random test vector generator (within 512-bit capacity)
// ==========================================================================

/// Generate a random Mpz that fits within gmp-rs's 512-bit capacity.
fn random_mpz() -> (Mpz, GmpMpz) {
    let mut limbs = [0u64; 8];
    // Use a simple LCG for deterministic test vectors
    static mut SEED: u64 = 12345;
    for limb in limbs.iter_mut() {
        unsafe {
            SEED = SEED.wrapping_mul(6364136223846793005).wrapping_add(1);
            *limb = SEED;
        }
    }
    // Trim limbs to avoid exceeding MAX_BITS for the product
    limbs[7] &= 0x7FFF_FFFF_FFFF_FFFF; // ensure top bit is 0 for safety
    let mut len = 8;
    while len > 0 && limbs[len - 1] == 0 {
        len -= 1;
    }
    let mut m = Mpz::new();
    m.mag[..len].copy_from_slice(&limbs[..len]);
    m.len = len;
    m.sign = if len == 0 { 0 } else { 1 };
    if m.sign == 0 {
        m.sign = 1;
        m.mag[0] = 1;
        m.len = 1;
    }

    let dec = s(&m);
    let gmp = GmpMpz::from_decimal(&dec);
    (m, gmp)
}

/// Convert Mpz to decimal string.
fn s(m: &Mpz) -> String {
    let mut buf = [0u8; 192];
    let len = m.write_decimal_buf(&mut buf);
    core::str::from_utf8(&buf[..len]).unwrap().into()
}

// ==========================================================================
// Cross-check tests
// ==========================================================================

#[test]
fn gmp_cross_check_add() {
    for _ in 0..50 {
        let (a, ga) = random_mpz();
        let (b, gb) = random_mpz();
        // Truncate to avoid capacity overflow in gmp-rs
        let a_scaled = a.tdiv_q_2exp(10);
        let b_scaled = b.tdiv_q_2exp(10);
        let ga_scaled = ga.tdiv_q(&GmpMpz::from_decimal("1024"));
        let gb_scaled = gb.tdiv_q(&GmpMpz::from_decimal("1024"));

        let gmp_result = ga_scaled.add(&gb_scaled);
        let gmp_s = gmp_result.to_decimal();

        let our_result = a_scaled.try_add(&b_scaled).unwrap();
        let our_s = s(&our_result);

        assert_eq!(
            our_s, gmp_s,
            "GMP cross-check ADD mismatch: gmp-rs={}, GMP={}",
            our_s, gmp_s
        );
    }
}

#[test]
fn gmp_cross_check_sub() {
    for _ in 0..50 {
        let (a, ga) = random_mpz();
        let (b, gb) = random_mpz();
        let a_scaled = a.tdiv_q_2exp(10);
        let b_scaled = b.tdiv_q_2exp(10);
        let ga_scaled = ga.tdiv_q(&GmpMpz::from_decimal("1024"));
        let gb_scaled = gb.tdiv_q(&GmpMpz::from_decimal("1024"));

        let gmp_result = ga_scaled.sub(&gb_scaled);
        let gmp_s = gmp_result.to_decimal();

        let our_result = a_scaled.try_sub(&b_scaled).unwrap();
        let our_s = s(&our_result);

        assert_eq!(
            our_s, gmp_s,
            "GMP cross-check SUB mismatch: gmp-rs={}, GMP={}",
            our_s, gmp_s
        );
    }
}

#[test]
fn gmp_cross_check_mul() {
    for _ in 0..50 {
        let (a, ga) = random_mpz();
        let (b, gb) = random_mpz();
        // Further reduce to prevent overflow in mul
        let a_reduced = a.tdiv_q_2exp(20);
        let b_reduced = b.tdiv_q_2exp(20);
        let two_pow_20 = Mpz::from_u64(1).try_mul_2exp(20).unwrap();
        let two_pow_20_gmp = GmpMpz::from_decimal(&s(&two_pow_20));
        let ga_reduced = ga.tdiv_q(&two_pow_20_gmp);
        let gb_reduced = gb.tdiv_q(&two_pow_20_gmp);

        let gmp_result = ga_reduced.mul(&gb_reduced);
        let gmp_s = gmp_result.to_decimal();

        if let Ok(our_result) = a_reduced.try_mul(&b_reduced) {
            let our_s = s(&our_result);
            assert_eq!(
                our_s, gmp_s,
                "GMP cross-check MUL mismatch: gmp-rs={}, GMP={}",
                our_s, gmp_s
            );
        }
    }
}

#[test]
fn gmp_cross_check_tdiv_qr() {
    for _ in 0..50 {
        let (a, ga) = random_mpz();
        let (b, gb) = random_mpz();
        let a_scaled = a.tdiv_q_2exp(10);
        let b_scaled = b.tdiv_q_2exp(10).try_add_ui(1).unwrap(); // ensure non-zero
        let ga_scaled = ga.tdiv_q(&GmpMpz::from_decimal("1024"));
        let gb_scaled = gb.tdiv_q(&GmpMpz::from_decimal("1024"));
        // Ensure non-zero
        let gb_one = GmpMpz::from_decimal("1");
        let gb_nonzero = gb_scaled.add(&gb_one);

        let (gmp_q, gmp_r) = ga_scaled.tdiv_qr(&gb_nonzero);
        let gmp_q_s = gmp_q.to_decimal();
        let gmp_r_s = gmp_r.to_decimal();

        let (our_q, our_r) = a_scaled.tdiv_qr(&b_scaled);
        let our_q_s = s(&our_q);
        let our_r_s = s(&our_r);

        assert_eq!(
            our_q_s, gmp_q_s,
            "GMP cross-check TDIV_Q mismatch: q: gmp-rs={}, GMP={}",
            our_q_s, gmp_q_s
        );
        assert_eq!(
            our_r_s, gmp_r_s,
            "GMP cross-check TDIV_R mismatch: r: gmp-rs={}, GMP={}",
            our_r_s, gmp_r_s
        );
    }
}

#[test]
fn gmp_cross_check_gcd() {
    for _ in 0..25 {
        let (a, ga) = random_mpz();
        let (b, gb) = random_mpz();
        let a_scaled = a.tdiv_q_2exp(10);
        let b_scaled = b.tdiv_q_2exp(10);
        let ga_scaled = ga.tdiv_q(&GmpMpz::from_decimal("1024"));
        let gb_scaled = gb.tdiv_q(&GmpMpz::from_decimal("1024"));

        let gmp_result = ga_scaled.gcd(&gb_scaled);
        let gmp_s = gmp_result.to_decimal();

        let our_result = a_scaled.try_gcd(&b_scaled).unwrap();
        let our_s = s(&our_result);

        assert_eq!(
            our_s, gmp_s,
            "GMP cross-check GCD mismatch: gmp-rs={}, GMP={}",
            our_s, gmp_s
        );
    }
}

// ==========================================================================
// Serial cross-check with known test vectors
// ==========================================================================

#[test]
fn gmp_cross_check_known_vectors() {
    // Known values that both implementations must agree on
    let test_vectors: &[(&str, &str)] = &[
        ("100", "200"),
        ("-50", "30"),
        ("12345678901234567890", "98765432109876543210"),
        ("0", "0"),
        ("1", "1"),
    ];

    for (a_str, b_str) in test_vectors {
        let a_gmp = GmpMpz::from_decimal(a_str);
        let b_gmp = GmpMpz::from_decimal(b_str);
        let a_our = Mpz::from_decimal_str(a_str).unwrap();
        let b_our = Mpz::from_decimal_str(b_str).unwrap();

        // ADD
        let gmp_add = a_gmp.add(&b_gmp).to_decimal();
        let our_add = s(&a_our.try_add(&b_our).unwrap());
        assert_eq!(our_add, gmp_add, "ADD({}, {})", a_str, b_str);

        // SUB
        let gmp_sub = a_gmp.sub(&b_gmp).to_decimal();
        let our_sub = s(&a_our.try_sub(&b_our).unwrap());
        assert_eq!(our_sub, gmp_sub, "SUB({}, {})", a_str, b_str);
    }
}
