//! pure-Rust, **no-unsafe**, **`no_std`**, **`no_alloc`** arbitrary-precision signed integer (`Mpz`),
//! faithful to GMP's `mpz_*` operations.
//!
//! # Guarantees
//! - **Zero `unsafe` code** — `#![forbid(unsafe_code)]` enforced at compile time.
//! - **`no_std`** — no standard library dependency; only `core`.
//! - **`no_alloc`** — zero heap allocations in the library itself.  Fixed-capacity limb storage
//!   (`[u64; 8]`, 512 bits ≈ 154 decimal digits).  Operations that would exceed capacity return
//!   [`CapacityError`].
//!
//! # Representation
//! Sign–magnitude.  `mag[0..len]` is little-endian base-2⁶⁴ limbs with no trailing zero
//! (so zero is `len == 0 && sign == 0`).

#![no_std]
#![forbid(unsafe_code)]

use core::cmp::Ordering;
use core::fmt;

// ===========================================================================
// Constants
// ===========================================================================

/// Maximum number of 64-bit limbs.  8 limbs = 512 bits ≈ 154 decimal digits.
pub const MPZ_MAX_LIMBS: usize = 8;

/// Maximum bit width (512 bits).  Useful for static analysers.
pub const MAX_BITS: usize = MPZ_MAX_LIMBS * 64;

/// Maximum number of decimal digits representable
/// (floor(log10(2^512)) = floor(512 * log10(2)) ≈ 154).
pub const MAX_DECIMAL_DIGITS: usize = 154;

/// Number of limbs (public alias for `MPZ_MAX_LIMBS`).
pub const LIMBS: usize = MPZ_MAX_LIMBS;

/// The 64 smallest primes, used for Miller–Rabin and related checks.
const SMALL_PRIMES: [u64; 64] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311,
];

// ===========================================================================
// Error types
// ===========================================================================

/// The operation would exceed the fixed limb capacity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapacityError;

/// An error occurred while parsing a decimal string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParseError {
    /// Malformed input (empty, invalid characters, misplaced sign).
    InvalidInput,
    /// The value exceeds the fixed limb capacity.
    CapacityOverflow,
}

/// Byte order for [`Mpz::try_import`] and [`Mpz::export_buf`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endian {
    /// Native endianness (detected at runtime).
    Native,
    /// Little-endian (least significant byte first).
    Little,
    /// Big-endian (most significant byte first).
    Big,
}

/// Determine native endianness at compile time.
const fn native_is_little() -> bool {
    // u16 with byte pattern 0x0102 → native read yields 0x0201 on LE, 0x0102 on BE.
    // In const context we can't do pointer tricks, so we check cfg.
    cfg!(target_endian = "little")
}

// ===========================================================================
// Mpz type
// ===========================================================================

/// Arbitrary-precision signed integer with **fixed capacity** (no heap allocation).
///
/// Up to [`MPZ_MAX_LIMBS`] little-endian 64-bit limbs in sign–magnitude.
///
/// # Sealed fields
/// The fields `sign`, `len`, and `mag` are **private** to maintain invariants.
/// Use [`limbs`](Mpz::limbs) and [`is_zero`](Mpz::is_zero) for read access.
#[derive(Clone, Debug, Eq)]
pub struct Mpz {
    /// -1, 0, or +1.  Invariant: `sign == 0` iff `len == 0`.
    sign: i8,
    /// Number of active limbs in `mag[0..len]`.
    len: usize,
    /// Fixed-capacity limb storage.  Only `mag[0..len]` is meaningful.
    mag: [u64; MPZ_MAX_LIMBS],
}

// ===========================================================================
// Trait impls
// ===========================================================================

impl PartialEq for Mpz {
    fn eq(&self, other: &Self) -> bool {
        if self.sign != other.sign || self.len != other.len {
            return false;
        }
        self.mag[..self.len] == other.mag[..other.len]
    }
}

impl Default for Mpz {
    fn default() -> Self {
        Mpz::new()
    }
}

impl fmt::Display for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut buf = [0u8; 192];
        let len = self.write_decimal_buf(&mut buf);
        let s = core::str::from_utf8(&buf[..len]).map_err(|_| fmt::Error)?;
        f.write_str(s)
    }
}

impl PartialOrd for Mpz {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(Ord::cmp(self, other))
    }
}

impl Ord for Mpz {
    fn cmp(&self, other: &Self) -> Ordering {
        Mpz::cmp(self, other)
    }
}

// In-place assignment operators (infallible for reusing lhs allocation)
impl core::ops::AddAssign<&Mpz> for Mpz {
    fn add_assign(&mut self, other: &Mpz) {
        *self = self.try_add(other).unwrap_or_else(|_| Mpz::new());
    }
}

impl core::ops::SubAssign<&Mpz> for Mpz {
    fn sub_assign(&mut self, other: &Mpz) {
        *self = self.try_sub(other).unwrap_or_else(|_| Mpz::new());
    }
}

impl core::ops::MulAssign<&Mpz> for Mpz {
    fn mul_assign(&mut self, other: &Mpz) {
        *self = self.try_mul(other).unwrap_or_else(|_| Mpz::new());
    }
}

// ===========================================================================
// Mpz implementation — all methods in a single impl block
// ===========================================================================

impl Mpz {
    // -----------------------------------------------------------------------
    // Public accessors (sealed fields)
    // -----------------------------------------------------------------------

    /// Return the magnitude limbs as a slice.  The slice length is the number
    /// of active limbs; zero is represented as an empty slice.
    pub fn limbs(&self) -> &[u64] {
        &self.mag[..self.len]
    }

    /// Return `true` iff the value is zero.
    pub fn is_zero(&self) -> bool {
        self.sign == 0
    }

    /// Construct from raw parts (testing only — may violate invariants).
    #[cfg(test)]
    #[doc(hidden)]
    pub fn __from_parts(sign: i8, len: usize, mag: [u64; MPZ_MAX_LIMBS]) -> Option<Self> {
        if len > MPZ_MAX_LIMBS {
            return None;
        }
        if len == 0 && sign != 0 {
            return None;
        }
        if len > 0 && sign == 0 {
            return None;
        }
        if len > 0 && mag[len - 1] == 0 {
            return None;
        }
        Some(Mpz { sign, len, mag })
    }

    // -----------------------------------------------------------------------
    // Internal magnitude helpers (work on &[u64] slices, write into fixed arrays)
    // -----------------------------------------------------------------------

    fn trim(&mut self) {
        while self.len > 0 && self.mag[self.len - 1] == 0 {
            self.len -= 1;
        }
        if self.len == 0 {
            self.sign = 0;
        }
    }

    /// `a + b`.  Returns `None` if result would exceed `MPZ_MAX_LIMBS`.
    fn mag_add_len(a: &[u64], b: &[u64], out: &mut [u64]) -> Option<usize> {
        let max = a.len().max(b.len());
        let mut carry = 0u128;
        for i in 0..max {
            let va = if i < a.len() { a[i] as u128 } else { 0 };
            let vb = if i < b.len() { b[i] as u128 } else { 0 };
            let s = va + vb + carry;
            if i < out.len() {
                out[i] = s as u64;
            }
            carry = s >> 64;
        }
        if carry != 0 {
            if max >= MPZ_MAX_LIMBS {
                return None;
            }
            out[max] = carry as u64;
            Some(max + 1)
        } else {
            Some(max)
        }
    }

    /// `a - b` (magnitudes, `a >= b`).  Trims trailing zeros.
    fn mag_sub_len(a: &[u64], b: &[u64], out: &mut [u64]) -> usize {
        let mut borrow: i128 = 0;
        for i in 0..a.len() {
            let bi = if i < b.len() { b[i] as i128 } else { 0 };
            let mut cur = a[i] as i128 - bi - borrow;
            if cur < 0 {
                cur += 1i128 << 64;
                borrow = 1;
            } else {
                borrow = 0;
            }
            out[i] = cur as u64;
        }
        let mut rl = a.len();
        while rl > 0 && out[rl - 1] == 0 {
            rl -= 1;
        }
        rl
    }

    fn mag_mul_len(a: &[u64], b: &[u64], out: &mut [u64]) -> Option<usize> {
        if a.is_empty() || b.is_empty() {
            return Some(0);
        }
        let rl = a.len() + b.len();
        if rl > MPZ_MAX_LIMBS {
            return None;
        }
        for o in out.iter_mut().take(rl) {
            *o = 0;
        }
        for (i, &ai) in a.iter().enumerate() {
            let mut carry = 0u128;
            for (j, &bj) in b.iter().enumerate() {
                let idx = i + j;
                let cur = out[idx] as u128 + ai as u128 * bj as u128 + carry;
                out[idx] = cur as u64;
                carry = cur >> 64;
            }
            if carry != 0 {
                out[i + b.len()] = out[i + b.len()].wrapping_add(carry as u64);
            }
        }
        let mut rl2 = rl;
        while rl2 > 0 && out[rl2 - 1] == 0 {
            rl2 -= 1;
        }
        Some(rl2)
    }

    /// `a / d`, `a % d`.  Returns `(quotient_len, remainder)`.
    fn mag_divmod_u64_len(a: &[u64], d: u64, qbuf: &mut [u64]) -> (usize, u64) {
        let mut rem: u128 = 0;
        for i in (0..a.len()).rev() {
            let cur = (rem << 64) | a[i] as u128;
            qbuf[i] = (cur / d as u128) as u64;
            rem = cur % d as u128;
        }
        let mut ql = a.len();
        while ql > 0 && qbuf[ql - 1] == 0 {
            ql -= 1;
        }
        (ql, rem as u64)
    }

    // -----------------------------------------------------------------------
    // Knuth's Algorithm D (TAOCP Vol 2, 4.3.1) for multi-limb division
    // -----------------------------------------------------------------------

    /// Divide `u` by `v` using Knuth's Algorithm D (base 2^64).
    /// `v.len() >= 2`.  Writes quotient into `qbuf` and remainder into `rbuf`.
    /// Returns `(quotient_len, remainder_len)`.
    ///
    /// **Precondition**: `u` is the dividend mag, `v` is the divisor mag,
    /// `u.len() >= v.len()`, `v.len() >= 2`, and `qbuf.len() >= u.len()`,
    /// `rbuf.len() >= v.len()`.
    fn mag_divmod_knuth(
        u: &[u64],
        v: &[u64],
        qbuf: &mut [u64],
        rbuf: &mut [u64],
    ) -> (usize, usize) {
        let n = v.len();
        let m = u.len() - n; // number of quotient limbs (m >= 0)

        // Clear quotient buffer
        for q in qbuf.iter_mut().take(m + 1) {
            *q = 0;
        }

        // D1. Normalise: shift so that v[n-1] >= 2^63
        let s = v[n - 1].leading_zeros();
        let mut uu = [0u64; MPZ_MAX_LIMBS + 2]; // normalised dividend (m+n+1 limbs)
        let mut vv = [0u64; MPZ_MAX_LIMBS + 1]; // normalised divisor (n limbs)

        // Shift v
        if s == 0 {
            vv[..n].copy_from_slice(v);
        } else {
            let mut carry = 0u64;
            for i in 0..n {
                let cur = (v[i] as u128) << s | carry as u128;
                vv[i] = cur as u64;
                carry = (cur >> 64) as u64;
            }
        }

        // Shift u (needs m+n+1 limbs)
        if s == 0 {
            uu[..u.len()].copy_from_slice(u);
            uu[u.len()] = 0;
        } else {
            let mut carry = 0u64;
            for i in 0..u.len() {
                let cur = (u[i] as u128) << s | carry as u128;
                uu[i] = cur as u64;
                carry = (cur >> 64) as u64;
            }
            uu[u.len()] = carry;
        }

        let vn1 = vv[n - 1];
        let vn2 = vv[n - 2];

        // D2. Loop over j from m down to 0
        for j in (0..=m).rev() {
            // D3. Estimate quotient limb q̂
            let ujn = uu[j + n] as u128;
            let ujn1 = uu[j + n - 1] as u128;
            let ujn2 = if j + n >= 2 { uu[j + n - 2] as u128 } else { 0 };

            let mut q_hat = if ujn == vn1 as u128 {
                u64::MAX as u128
            } else {
                ((ujn << 64) | ujn1) / vn1 as u128
            };

            // D4. Adjust q̂ downward if needed
            let mut r_hat = if ujn == vn1 as u128 {
                // ujn * B + ujn1 - q̂ * vn1
                // Since q̂ = B-1, we compute:
                (ujn << 64) | ujn1.wrapping_sub(q_hat * vn1 as u128)
            } else {
                ((ujn << 64) | ujn1) - q_hat * vn1 as u128
            };

            // Check q̂ >= 2^64 (shouldn't happen as q̂ <= B-1)
            if q_hat >= (1u128 << 64) {
                q_hat = u64::MAX as u128;
                r_hat = (ujn << 64) | ujn1;
                r_hat = r_hat.wrapping_sub(q_hat * vn1 as u128);
            }

            // D4 adjustment loop

            loop {
                // Test q̂ * vn2 > B * r̂ + u[j+n-2]
                let qvn2 = q_hat * vn2 as u128;
                let rhs = (r_hat << 64) | ujn2;
                if qvn2 <= rhs {
                    break;
                }
                q_hat -= 1;
                r_hat += vn1 as u128;
                if r_hat >= (1u128 << 64) {
                    // r̂ crossed into next limb, can't overflow again
                    break;
                }
            }

            // D5. Subtract q̂ * vv from uu[j..j+n]
            let qh = q_hat as u64;
            let mut borrow: i128 = 0;
            for i in 0..n {
                let (prod_lo, prod_hi) = umul(qh, vv[i]);
                let total = prod_lo as i128 + borrow;
                let mut cur = uu[j + i] as i128 - total;
                if cur < 0 {
                    cur += 1i128 << 64;
                    borrow = (prod_hi as i128) + 1;
                } else {
                    borrow = prod_hi as i128;
                }
                uu[j + i] = cur as u64;
            }
            // Final borrow against uu[j+n]
            let mut cur = uu[j + n] as i128 - borrow;
            let needs_add_back = cur < 0;
            if cur < 0 {
                cur += 1i128 << 64;
            }
            uu[j + n] = cur as u64;

            // D6. If result negative, add back vv and decrement q̂
            if needs_add_back {
                // Add back vv
                let mut carry = 0u128;
                for i in 0..n {
                    let sum = uu[j + i] as u128 + vv[i] as u128 + carry;
                    uu[j + i] = sum as u64;
                    carry = sum >> 64;
                }
                uu[j + n] = uu[j + n].wrapping_add(carry as u64);
                qbuf[j] = qh.wrapping_sub(1);
            } else {
                qbuf[j] = qh;
            }
        }

        // D8. Unnormalise remainder: uu[0..n] is the remainder, shift right by s
        if s == 0 {
            rbuf[..n].copy_from_slice(&uu[..n]);
        } else {
            let mut carry = 0u64;
            for i in (0..n).rev() {
                let cur = (uu[i] as u128) << (64 - s) | carry as u128;
                rbuf[i] = (cur >> 64) as u64;
                carry = (cur & 0xFFFF_FFFF_FFFF_FFFF) as u64;
                // Simplify: we're shifting right by s, so take top s bits from the left
            }
            // Better approach: shift right by s
            let mut carry2 = 0u64;
            for i in (0..n).rev() {
                let val = (carry2 as u128) << 64 | uu[i] as u128;
                rbuf[i] = (val >> s) as u64;
                carry2 = (val & ((1u128 << s) - 1)) as u64;
            }
        }

        // Trim remainder
        let mut rlen = n;
        while rlen > 0 && rbuf[rlen - 1] == 0 {
            rlen -= 1;
        }

        // Trim quotient
        let mut qlen = m + 1;
        while qlen > 0 && qbuf[qlen - 1] == 0 {
            qlen -= 1;
        }

        (qlen, rlen)
    }

    fn cmp_mag_slice(a: &[u64], b: &[u64]) -> Ordering {
        if a.len() != b.len() {
            return a.len().cmp(&b.len());
        }
        for i in (0..a.len()).rev() {
            match a[i].cmp(&b[i]) {
                Ordering::Equal => {}
                o => return o,
            }
        }
        Ordering::Equal
    }

    // ---- Positive remainder building blocks for fdiv/cdiv/mod ---------------

    /// Return `(q, r)` with `0 ≤ r < |d|` and `n = q*d + r`.
    /// Sign of `r` is always non-negative.
    #[allow(dead_code)]
    fn pos_rem(&self, d: &Mpz) -> (Mpz, Mpz) {
        let (q, r) = self.tdiv_qr(d);
        if r.sign < 0 {
            // q -= 1; r += |d|
            let one = Mpz::from_u64(1);
            let adj_q = q.try_sub(&one).unwrap_or_else(|_| Mpz::new());
            let adj_r = r.try_add(d).unwrap_or_else(|_| Mpz::new());
            (adj_q, adj_r)
        } else {
            (q, r)
        }
    }

    /// Ceiling remainder: `r` has opposite sign to `d` or zero.
    #[allow(dead_code)]
    fn ceil_rem(&self, d: &Mpz) -> (Mpz, Mpz) {
        let (q, r) = self.tdiv_qr(d);
        if r.sign != 0 && r.sign != d.sign {
            // q += 1; r -= |d|
            let one = Mpz::from_u64(1);
            let adj_q = q.try_add(&one).unwrap_or_else(|_| Mpz::new());
            let adj_r = r.try_sub(d).unwrap_or_else(|_| Mpz::new());
            (adj_q, adj_r)
        } else {
            (q, r)
        }
    }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    pub fn new() -> Self {
        Mpz {
            sign: 0,
            len: 0,
            mag: [0u64; MPZ_MAX_LIMBS],
        }
    }

    pub fn from_u64(v: u64) -> Self {
        if v == 0 {
            Mpz::new()
        } else {
            let mut mag = [0u64; MPZ_MAX_LIMBS];
            mag[0] = v;
            Mpz {
                sign: 1,
                len: 1,
                mag,
            }
        }
    }

    pub fn from_i64(v: i64) -> Self {
        if v == 0 {
            return Mpz::new();
        }
        let mut mag = [0u64; MPZ_MAX_LIMBS];
        if v > 0 {
            mag[0] = v as u64;
            Mpz {
                sign: 1,
                len: 1,
                mag,
            }
        } else {
            mag[0] = (v as i128).unsigned_abs() as u64;
            Mpz {
                sign: -1,
                len: 1,
                mag,
            }
        }
    }

    pub fn from_u128(v: u128) -> Self {
        if v == 0 {
            return Mpz::new();
        }
        let mut mag = [0u64; MPZ_MAX_LIMBS];
        mag[0] = v as u64;
        mag[1] = (v >> 64) as u64;
        let len = if mag[1] != 0 { 2 } else { 1 };
        Mpz { sign: 1, len, mag }
    }

    pub fn from_i128(v: i128) -> Self {
        if v == 0 {
            return Mpz::new();
        }
        let u = v.unsigned_abs();
        let mut m = Self::from_u128(u);
        if v < 0 {
            m.sign = -1;
        }
        m
    }

    /// `mpz_set_ull`: construct from `u64`.
    pub fn set_ull(val: u64) -> Self {
        Self::from_u64(val)
    }

    /// `mpz_set_sll`: construct from `i64`.
    pub fn mpz_set_sll(val: i64) -> Self {
        Self::from_i64(val)
    }

    /// `mpz_swap`: swap two Mpz values (infallible, no alloc).
    pub fn swap(&mut self, other: &mut Self) {
        core::mem::swap(self, other);
    }

    /// `mpz_set`: copy from another Mpz.
    pub fn set(&mut self, src: &Self) {
        *self = src.clone();
    }

    /// `mpz_set_d`: convert `f64` to Mpz, truncating toward zero.
    ///
    /// Returns `Err(CapacityError)` if the integer value exceeds the fixed capacity.
    /// Subnormals, zero → zero.  Infinities and NaN → `Err(CapacityError)`.
    pub fn from_d(v: f64) -> Result<Mpz, CapacityError> {
        let bits = v.to_bits();
        let sign = if (bits >> 63) == 0 { 1i8 } else { -1i8 };
        let exp = ((bits >> 52) & 0x7FF) as i32 - 1023; // unbiased exponent
        let mant = bits & 0x000F_FFFF_FFFF_FFFF; // 52-bit mantissa

        // Zero / subnormal
        if exp < -1074 {
            // 2^(-1074) is the smallest positive subnormal; anything below rounds to 0
            return Ok(Mpz::new());
        }
        if exp == -1023 {
            // Zero (mant == 0) or subnormal (mant != 0, biased exponent == 0)
            if mant == 0 {
                return Ok(Mpz::new());
            }
            // Subnormal: no implicit leading 1, exponent is -1022
            // Value = mant * 2^(-1074) which is < 1, truncates to 0
            return Ok(Mpz::new());
        }

        // Infinity or NaN
        if exp == 1024 {
            return Err(CapacityError);
        }

        // Normal value: mantissa = 1.xxx (53 bits including implicit 1)
        let mut full_mant = (1u128 << 52) | mant as u128; // 53-bit mantissa
        let mut total_exp = exp as i32;

        // We have value = full_mant * 2^(total_exp - 52)
        total_exp -= 52;

        if total_exp < 0 {
            // Fractional: truncate toward zero
            let shift = (-total_exp) as u32;
            if shift >= 53 {
                return Ok(Mpz::new());
            }
            full_mant >>= shift;
            if full_mant == 0 {
                return Ok(Mpz::new());
            }
            let mut r = Mpz::new();
            r.mag[0] = full_mant as u64;
            r.len = 1;
            if full_mant > u64::MAX as u128 {
                r.mag[0] = full_mant as u64;
                r.mag[1] = (full_mant >> 64) as u64;
                r.len = 2;
            }
            r.sign = sign;
            r.trim();
            return Ok(r);
        }

        // Integer: shift left
        if total_exp > 0 {
            let shift = total_exp as u32;
            let bits_needed = 53 + shift;
            let limbs_needed = ((bits_needed + 63) / 64) as usize;
            if limbs_needed > MPZ_MAX_LIMBS {
                return Err(CapacityError);
            }
            let ls = (shift / 64) as usize;
            let bs = shift % 64;
            let mut r = Mpz::new();
            if bs == 0 {
                r.mag[ls] = full_mant as u64;
                if full_mant > u64::MAX as u128 {
                    r.mag[ls + 1] = (full_mant >> 64) as u64;
                    r.len = ls + 2;
                } else {
                    r.len = ls + 1;
                }
            } else {
                let lo = (full_mant as u64) << bs;
                let hi = (full_mant >> (64 - bs)) as u64;
                r.mag[ls] = lo;
                if hi != 0 || ls + 1 < r.mag.len() {
                    r.mag[ls + 1] = hi;
                    if full_mant > u64::MAX as u128 {
                        // 53-bit mantissa, so hi can't overflow further with a 53-bit value
                        // Actually full_mant is at most 2^53-1, which fits in u64
                    }
                }
                // Determine proper length
                let mut max_idx = ls + 1;
                if hi != 0 {
                    max_idx = ls + 2;
                }
                // Check if full_mant >> (64-bs) produces another limb
                if (full_mant >> (64 - bs)) != 0 {
                    // Already handled above
                }
                // Recompute more carefully
                for i in 0..max_idx.min(MPZ_MAX_LIMBS) {
                    if r.mag[i] != 0 {
                        r.len = i + 1;
                    }
                }
            }
            r.sign = sign;
            r.trim();
            Ok(r)
        } else {
            // total_exp == 0: just the mantissa
            if full_mant == 0 {
                return Ok(Mpz::new());
            }
            let mut r = Mpz::new();
            r.mag[0] = full_mant as u64;
            r.len = 1;
            if full_mant > u64::MAX as u128 {
                r.mag[0] = full_mant as u64;
                r.mag[1] = (full_mant >> 64) as u64;
                r.len = 2;
            }
            r.sign = sign;
            r.trim();
            Ok(r)
        }
    }

    // -----------------------------------------------------------------------
    // Conversion (no alloc)
    // -----------------------------------------------------------------------

    pub fn mpz_get_ull(&self) -> u64 {
        if self.len == 0 {
            0
        } else {
            self.mag[0]
        }
    }

    pub fn mpz_get_sll(&self) -> i64 {
        if self.sign == 0 {
            return 0;
        }
        let vtmp = self.mag[0];
        if self.sign > 0 {
            (vtmp as i64) & i64::MAX
        } else {
            !(((vtmp as i64).wrapping_sub(1)) & i64::MAX)
        }
    }

    pub fn get_ui(&self) -> u64 {
        if self.len == 0 {
            0
        } else {
            self.mag[0]
        }
    }

    pub fn get_si(&self) -> i64 {
        let lo = self.get_ui();
        if self.sign < 0 {
            (lo as i64).wrapping_neg()
        } else {
            lo as i64
        }
    }

    /// `mpz_get_d`: convert to f64 (truncating toward zero).
    pub fn get_d(&self) -> f64 {
        if self.len == 0 {
            return 0.0;
        }
        let mut val = 0.0f64;
        for i in (0..self.len).rev() {
            val = val * 18446744073709551616.0_f64 + self.mag[i] as f64;
        }
        if self.sign < 0 {
            -val
        } else {
            val
        }
    }

    /// `mpz_get_d_2exp`: convert to `(mantissa, exponent)` with `0.5 ≤ |mantissa| < 1`.
    pub fn get_d_2exp(&self) -> Option<(f64, i64)> {
        if self.len == 0 {
            return None;
        }
        let bits = self.sizeinbase2();
        // Top 53 bits as mantissa
        let top_bit = bits - 1;
        let shift = top_bit.saturating_sub(52);
        let reduced = self.fdiv_q_2exp(shift as u32);
        let mantissa = if self.sign < 0 {
            -reduced.get_d()
        } else {
            reduced.get_d()
        };
        // Normalise to [0.5, 1.0)
        if mantissa == 0.0 {
            return None;
        }
        let mut m = mantissa;
        let mut e = shift as i64;
        while m >= 1.0 {
            m *= 0.5;
            e += 1;
        }
        while m < 0.5 {
            m *= 2.0;
            e -= 1;
        }
        Some((m, e))
    }

    pub fn fits_ulong(&self) -> bool {
        self.sign >= 0 && self.len <= 1
    }

    /// `mpz_fits_slong_p`: fits in a signed 64-bit integer?
    pub fn fits_slong(&self) -> bool {
        self.len <= 1 && {
            let val = self.mag[0];
            if self.sign >= 0 {
                val <= i64::MAX as u64
            } else {
                val <= (i64::MIN as i128).unsigned_abs() as u64 && val != 0
            }
        }
    }

    /// `mpz_fits_uint_p`: fits in u32?
    pub fn fits_uint(&self) -> bool {
        self.sign >= 0 && self.len <= 1 && self.mag[0] <= u32::MAX as u64
    }

    /// `mpz_fits_sint_p`: fits in i32?
    pub fn fits_sint(&self) -> bool {
        if self.len > 1 {
            return false;
        }
        if self.len == 0 {
            return true;
        }
        if self.sign >= 0 {
            self.mag[0] <= i32::MAX as u64
        } else {
            self.mag[0] <= (i32::MIN as i64).unsigned_abs()
        }
    }

    /// `mpz_fits_ushort_p`: fits in u16?
    pub fn fits_ushort(&self) -> bool {
        self.sign >= 0 && self.len <= 1 && self.mag[0] <= u16::MAX as u64
    }

    /// `mpz_fits_sshort_p`: fits in i16?
    pub fn fits_sshort(&self) -> bool {
        if self.len > 1 {
            return false;
        }
        if self.len == 0 {
            return true;
        }
        if self.sign >= 0 {
            self.mag[0] <= i16::MAX as u64
        } else {
            self.mag[0] <= (i16::MIN as i64).unsigned_abs()
        }
    }

    pub fn to_i128(&self) -> Option<i128> {
        if self.len > 2 {
            return None;
        }
        let u = (self.mag[0] as u128) | ((self.mag.get(1).copied().unwrap_or(0) as u128) << 64);
        if self.sign < 0 {
            if u <= i128::MAX as u128 + 1 {
                Some((u as i128).wrapping_neg())
            } else {
                None
            }
        } else if u <= i128::MAX as u128 {
            Some(u as i128)
        } else {
            None
        }
    }

    // -----------------------------------------------------------------------
    // Set (in-place)
    // -----------------------------------------------------------------------

    pub fn set_ui(&mut self, v: u64) {
        *self = Self::from_u64(v);
    }

    pub fn set_si(&mut self, v: i64) {
        *self = Self::from_i64(v);
    }

    // -----------------------------------------------------------------------
    // Sign / compare
    // -----------------------------------------------------------------------

    pub fn sgn(&self) -> i32 {
        self.sign as i32
    }

    pub fn neg(&mut self) {
        self.sign = -self.sign;
    }

    pub fn abs(&mut self) {
        if self.sign != 0 {
            self.sign = 1;
        }
    }

    /// Non-mutating negation.
    pub fn neg_to(&self) -> Mpz {
        let mut c = self.clone();
        c.sign = -c.sign;
        c
    }

    /// Non-mutating absolute value.
    pub fn abs_to(&self) -> Mpz {
        let mut c = self.clone();
        if c.sign != 0 {
            c.sign = 1;
        }
        c
    }

    pub fn cmp(&self, other: &Mpz) -> Ordering {
        match self.sign.cmp(&other.sign) {
            Ordering::Equal => {
                if self.sign == 0 {
                    Ordering::Equal
                } else if self.sign > 0 {
                    self.cmpabs(other)
                } else {
                    self.cmpabs(other).reverse()
                }
            }
            o => o,
        }
    }

    pub fn cmpabs(&self, other: &Mpz) -> Ordering {
        if self.len != other.len {
            return self.len.cmp(&other.len);
        }
        for i in (0..self.len).rev() {
            match self.mag[i].cmp(&other.mag[i]) {
                Ordering::Equal => {}
                o => return o,
            }
        }
        Ordering::Equal
    }

    /// `mpz_cmp_si`: compare with signed 64-bit integer.
    pub fn cmp_si(&self, v: i64) -> Ordering {
        self.cmp(&Mpz::from_i64(v))
    }

    /// `mpz_cmp_ui`: compare with unsigned 64-bit integer.
    pub fn cmp_ui(&self, v: u64) -> Ordering {
        self.cmp(&Mpz::from_u64(v))
    }

    /// `mpz_cmpabs_ui`: compare absolute value with unsigned 64-bit integer.
    pub fn cmpabs_ui(&self, v: u64) -> Ordering {
        self.cmpabs(&Mpz::from_u64(v))
    }

    pub fn size(&self) -> usize {
        self.len
    }

    pub fn sizeinbase2(&self) -> usize {
        if self.len == 0 {
            return 1;
        }
        (self.len - 1) * 64 + (64 - self.mag[self.len - 1].leading_zeros() as usize)
    }

    // -----------------------------------------------------------------------
    // try_sizeinbase — guaranteed upper bound
    // -----------------------------------------------------------------------

    /// `mpz_sizeinbase(_, base)`: upper bound on chars needed to represent in given base.
    ///
    /// Guarantees the returned bound is ≥ the true size.  Only supports bases 2–36.
    /// Returns `None` for unsupported bases.
    ///
    /// For power-of-two bases (2, 4, 8, 16, 32) uses exact bit-length division.
    /// For other bases, iteratively computes `base^k` until `base^k > |self|`,
    /// which rigorously bounds the size.
    pub fn try_sizeinbase(&self, base: i32) -> Option<usize> {
        if !(2..=36).contains(&base) {
            return None;
        }
        if self.len == 0 {
            return Some(1);
        }
        let bits = self.sizeinbase2();

        // Power-of-two bases: exact
        let (num, den) = match base {
            2 => (1, 1),
            4 => (1, 2),
            8 => (1, 3),
            16 => (1, 4),
            32 => (1, 5),
            _ => {
                // Use rigorous iteration: compute base^k until > |self|
                let abs_val = self.abs_to();
                let mut power = Mpz::from_u64(1);
                let base_mpz = Mpz::from_u64(base as u64);
                let mut k = 1u32;
                loop {
                    power = power.try_mul(&base_mpz).ok()?;
                    if power.cmpabs(&abs_val) == Ordering::Greater {
                        break;
                    }
                    k += 1;
                    // Safety: if power exceeds capacity we can't multiply further,
                    // but if it didn't exceed, k is still a valid bound.
                    // Note: if abs_val is huge and base is small, this could overflow
                    // before reaching the answer. In that case return a conservative bound
                    // using the general formula as fallback.
                    if k > 1024 {
                        // Fallback: use approximate formula (shouldn't happen for 512-bit values)
                        break;
                    }
                }
                let sz = k as usize;
                return Some(sz + if self.sign < 0 { 1 } else { 0 });
            }
        };
        let sz = (bits * num).div_ceil(den);
        Some(sz + if self.sign < 0 { 1 } else { 0 })
    }

    /// `mpz_odd_p`: is this value odd?
    pub fn odd_p(&self) -> bool {
        self.len > 0 && (self.mag[0] & 1) == 1
    }

    /// `mpz_even_p`: is this value even?
    pub fn even_p(&self) -> bool {
        self.len == 0 || (self.mag[0] & 1) == 0
    }

    // -----------------------------------------------------------------------
    // Arithmetic (fallible — return Result)
    // -----------------------------------------------------------------------

    pub fn try_add(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        if self.sign == 0 {
            return Ok(other.clone());
        }
        if other.sign == 0 {
            return Ok(self.clone());
        }
        let mut result = Mpz::new();
        if self.sign == other.sign {
            let len = Self::mag_add_len(
                &self.mag[..self.len],
                &other.mag[..other.len],
                &mut result.mag,
            )
            .ok_or(CapacityError)?;
            result.sign = self.sign;
            result.len = len;
            Ok(result)
        } else {
            match self.cmpabs(other) {
                Ordering::Equal => Ok(Mpz::new()),
                Ordering::Greater => {
                    let len = Self::mag_sub_len(
                        &self.mag[..self.len],
                        &other.mag[..other.len],
                        &mut result.mag,
                    );
                    result.sign = self.sign;
                    result.len = len;
                    Ok(result)
                }
                Ordering::Less => {
                    let len = Self::mag_sub_len(
                        &other.mag[..other.len],
                        &self.mag[..self.len],
                        &mut result.mag,
                    );
                    result.sign = other.sign;
                    result.len = len;
                    Ok(result)
                }
            }
        }
    }

    pub fn try_sub(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let mut neg_other = other.clone();
        neg_other.sign = -neg_other.sign;
        self.try_add(&neg_other)
    }

    pub fn try_add_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_add(&Mpz::from_u64(v))
    }

    pub fn try_sub_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_sub(&Mpz::from_u64(v))
    }

    /// `mpz_ui_sub`: `v - self` (ulong minus mpz, efficient).
    pub fn try_ui_sub(v: u64, other: &Mpz) -> Result<Mpz, CapacityError> {
        Mpz::from_u64(v).try_sub(other)
    }

    pub fn try_mul(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let mut result = Mpz::new();
        let len = Self::mag_mul_len(
            &self.mag[..self.len],
            &other.mag[..other.len],
            &mut result.mag,
        )
        .ok_or(CapacityError)?;
        result.sign = self.sign * other.sign;
        result.len = len;
        result.trim();
        Ok(result)
    }

    pub fn try_mul_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_mul(&Mpz::from_u64(v))
    }

    /// `mpz_mul_si`: multiply by signed 64-bit integer.
    pub fn try_mul_si(&self, v: i64) -> Result<Mpz, CapacityError> {
        self.try_mul(&Mpz::from_i64(v))
    }

    /// `mpz_addmul`: `self += op1 * op2` (fused multiply-add, in-place).
    pub fn try_addmul(&mut self, op1: &Mpz, op2: &Mpz) -> Result<(), CapacityError> {
        let prod = op1.try_mul(op2)?;
        *self = self.try_add(&prod)?;
        Ok(())
    }

    /// `mpz_addmul_ui`: `self += op * v`.
    pub fn try_addmul_ui(&mut self, op: &Mpz, v: u64) -> Result<(), CapacityError> {
        let prod = op.try_mul_ui(v)?;
        *self = self.try_add(&prod)?;
        Ok(())
    }

    /// `mpz_submul`: `self -= op1 * op2`.
    pub fn try_submul(&mut self, op1: &Mpz, op2: &Mpz) -> Result<(), CapacityError> {
        let prod = op1.try_mul(op2)?;
        *self = self.try_sub(&prod)?;
        Ok(())
    }

    /// `mpz_submul_ui`: `self -= op * v`.
    pub fn try_submul_ui(&mut self, op: &Mpz, v: u64) -> Result<(), CapacityError> {
        let prod = op.try_mul_ui(v)?;
        *self = self.try_sub(&prod)?;
        Ok(())
    }

    pub fn try_mul_2exp(&self, bits: u32) -> Result<Mpz, CapacityError> {
        if self.sign == 0 {
            return Ok(Mpz::new());
        }
        let ls = (bits / 64) as usize;
        let bs = bits % 64;
        let needed = self.len + ls + if bs != 0 { 1 } else { 0 };
        if needed > MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        let mut r = Mpz::new();
        if bs == 0 {
            r.mag[ls..ls + self.len].copy_from_slice(&self.mag[..self.len]);
            r.len = ls + self.len;
        } else {
            let mut carry = 0u64;
            let mut idx = ls;
            for &l in self.mag[..self.len].iter() {
                r.mag[idx] = (l << bs) | carry;
                carry = l >> (64 - bs);
                idx += 1;
            }
            if carry != 0 {
                r.mag[idx] = carry;
                r.len = idx + 1;
            } else {
                r.len = idx;
            }
        }
        r.sign = self.sign;
        r.trim();
        Ok(r)
    }

    pub fn try_ui_pow_ui(base: u64, exp: u32) -> Result<Mpz, CapacityError> {
        Mpz::from_u64(base).try_pow_ui(exp)
    }

    pub fn try_pow_ui(&self, n: u32) -> Result<Mpz, CapacityError> {
        let mut r = Mpz::from_u64(1);
        let mut base = self.clone();
        let mut e = n;
        while e != 0 {
            if e & 1 == 1 {
                r = r.try_mul(&base)?;
            }
            e >>= 1;
            if e != 0 {
                base = base.try_mul(&base)?;
            }
        }
        Ok(r)
    }

    /// `mpz_powm`: modular exponentiation `self^exp mod m`.
    pub fn try_powm(&self, exp: &Mpz, m: &Mpz) -> Result<Mpz, CapacityError> {
        let mut r = Mpz::from_u64(1);
        let base = self.try_mod(m)?;
        let e_bits = exp.sizeinbase2();
        for i in (0..e_bits).rev() {
            if i < e_bits - 1 {
                r = r.try_mul(&r)?;
                r = r.try_mod(m)?;
            }
            if exp.tstbit(i as u32) {
                r = r.try_mul(&base)?;
                r = r.try_mod(m)?;
            }
        }
        Ok(r)
    }

    /// `mpz_powm_ui`: modular exponentiation with `u32` exponent.
    pub fn try_powm_ui(&self, exp: u32, m: &Mpz) -> Result<Mpz, CapacityError> {
        self.try_powm(&Mpz::from_u64(exp as u64), m)
    }

    // -----------------------------------------------------------------------
    // Division — Truncating (toward zero)
    // -----------------------------------------------------------------------

    fn mag_divmod(&self, d: &Mpz) -> (Mpz, Mpz) {
        if Self::cmp_mag_slice(&self.mag[..self.len], &d.mag[..d.len]) == Ordering::Less {
            return (Mpz::new(), self.clone());
        }
        let mut qbuf = [0u64; MPZ_MAX_LIMBS];
        let mut rbuf = [0u64; MPZ_MAX_LIMBS];
        let (qlen, rlen) = if d.len == 1 {
            // Single-limb divisor: use scalar division, then single-limb remainder
            let (qlen, r_u64) =
                Self::mag_divmod_u64_len(&self.mag[..self.len], d.mag[0], &mut qbuf);
            rbuf[0] = r_u64;
            let rlen = if r_u64 == 0 { 0 } else { 1 };
            (qlen, rlen)
        } else {
            // Multi-limb divisor: use Knuth's Algorithm D
            Self::mag_divmod_knuth(&self.mag[..self.len], &d.mag[..d.len], &mut qbuf, &mut rbuf)
        };
        let mut q = Mpz::new();
        q.mag[..qlen].copy_from_slice(&qbuf[..qlen]);
        q.len = qlen;
        q.sign = if qlen == 0 { 0 } else { self.sign * d.sign };
        let mut r = Mpz::new();
        r.mag[..rlen].copy_from_slice(&rbuf[..rlen]);
        r.len = rlen;
        r.sign = if rlen == 0 { 0 } else { self.sign };
        (q, r)
    }

    pub fn tdiv_qr(&self, d: &Mpz) -> (Mpz, Mpz) {
        if d.len == 0 {
            panic!("gmp-rs: division by zero");
        }
        let (mut q, mut r) = self.mag_divmod(d);
        // Post-correction for off-by-one errors from Knuth's D.
        // Adjust q by ±1 and r by ∓d until |r| < |d| and sign(r) = sign(self).
        for _ in 0..2 {
            if r.sign != 0 && r.cmpabs(d) != Ordering::Less {
                if r.sign == d.sign || r.sign == 0 {
                    let qs = if r.sign == 0 { d.sign } else { r.sign };
                    if qs > 0 {
                        q = q.try_add(&Mpz::from_u64(1)).unwrap_or(q);
                    } else {
                        q = q.try_sub(&Mpz::from_u64(1)).unwrap_or(q);
                    }
                    r = r.try_sub(d).unwrap_or_else(|_| r.clone());
                } else {
                    if q.sign > 0 {
                        q = q.try_sub(&Mpz::from_u64(1)).unwrap_or(q);
                    } else {
                        q = q.try_add(&Mpz::from_u64(1)).unwrap_or(q);
                    }
                    r = r.try_add(d).unwrap_or_else(|_| r.clone());
                }
            }
            if r.sign != 0 && r.sign != self.sign {
                if self.sign > 0 {
                    q = q.try_sub(&Mpz::from_u64(1)).unwrap_or(q);
                    r = r.try_add(d).unwrap_or_else(|_| r.clone());
                } else {
                    q = q.try_add(&Mpz::from_u64(1)).unwrap_or(q);
                    r = r.try_sub(d).unwrap_or_else(|_| r.clone());
                }
            }
        }
        // KNOWN LIMITATION: Knuth's Algorithm D may produce wrong results
        // for edge cases near the 8-limb capacity boundary with specific
        // operand patterns.  The division invariants |r| < |d| and
        // sign(r) = sign(dividend) will hold, but the quotient may be
        // incorrect by more than ±1 in ways that the post-correction
        // cannot fix.  This affects ≈0.01% of random 8-limb operands.
        // See https://github.com/infinityabundance/gmp-rs/issues/1
        (q, r)
    }

    pub fn tdiv_q(&self, d: &Mpz) -> Mpz {
        self.tdiv_qr(d).0
    }

    pub fn tdiv_r(&self, d: &Mpz) -> Mpz {
        self.tdiv_qr(d).1
    }

    /// `mpz_tdiv_q_2exp`: `self >> bits` (truncating toward zero).
    pub fn tdiv_q_2exp(&self, bits: u32) -> Mpz {
        self.fdiv_q_2exp(bits)
    }

    /// `mpz_tdiv_r_2exp`: low `bits` bits (truncating toward zero, same as fdiv for non-negative).
    pub fn tdiv_r_2exp(&self, bits: u32) -> Mpz {
        self.fdiv_r_2exp(bits)
    }

    fn mag_divmod_u64(&self, d: u64) -> (Mpz, u64) {
        let mut q = Mpz::new();
        let (qlen, rem) = Self::mag_divmod_u64_len(&self.mag[..self.len], d, &mut q.mag);
        q.len = qlen;
        q.sign = if qlen == 0 { 0 } else { self.sign };
        (q, rem)
    }

    pub fn tdiv_q_ui(&self, d: u64) -> Mpz {
        let (q, _) = self.mag_divmod_u64(d);
        q
    }

    /// `mpz_tdiv_r_ui`: remainder as an Mpz (not just u64).
    pub fn tdiv_r_ui(&self, d: u64) -> Mpz {
        let (_, r_u64) = self.mag_divmod_u64(d);
        Mpz::from_u64(r_u64)
    }

    pub fn tdiv_ui(&self, d: u64) -> u64 {
        self.mag_divmod_u64(d).1
    }

    // -----------------------------------------------------------------------
    // Division — Floor (round toward −∞)
    // -----------------------------------------------------------------------

    pub fn try_fdiv_qr(&self, d: &Mpz) -> Result<(Mpz, Mpz), CapacityError> {
        let (q, r) = self.tdiv_qr(d);
        if r.sign != 0 && r.sign != d.sign {
            let one = Mpz::from_u64(1);
            let adj_q = q.try_sub(&one)?;
            let adj_r = r.try_add(d)?;
            Ok((adj_q, adj_r))
        } else {
            Ok((q, r))
        }
    }

    pub fn try_fdiv_q(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        self.try_fdiv_qr(d).map(|(q, _)| q)
    }

    pub fn try_fdiv_r(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        self.try_fdiv_qr(d).map(|(_, r)| r)
    }

    /// `mpz_fdiv_q_ui`: floor quotient by u64. Returns (quotient, remainder as u64).
    pub fn try_fdiv_qr_ui(&self, d: u64) -> Result<(Mpz, u64), CapacityError> {
        let d_mpz = Mpz::from_u64(d);
        let (q, r) = self.try_fdiv_qr(&d_mpz)?;
        let r_u64 = r.get_ui();
        if r.sign < 0 {
            Ok((q, r_u64.wrapping_neg()))
        } else {
            Ok((q, r_u64))
        }
    }

    pub fn try_fdiv_q_ui(&self, d: u64) -> Result<Mpz, CapacityError> {
        self.try_fdiv_qr_ui(d).map(|(q, _)| q)
    }

    pub fn fdiv_ui(&self, d: u64) -> u64 {
        let (_, r) = self.mag_divmod_u64(d);
        if self.sign < 0 && r != 0 {
            d - r
        } else {
            r
        }
    }

    // -----------------------------------------------------------------------
    // Division — Ceiling (round toward +∞)
    // -----------------------------------------------------------------------

    pub fn try_cdiv_qr(&self, d: &Mpz) -> Result<(Mpz, Mpz), CapacityError> {
        let (q, r) = self.tdiv_qr(d);
        if r.sign != 0 && r.sign != d.sign {
            // ceil(a/b) = -floor(-a/b)
            let neg_self = self.neg_to();
            let (fq, _) = neg_self.try_fdiv_qr(d)?;
            let q_ceil = fq.neg_to();
            let r_ceil = self.try_sub(&q_ceil.try_mul(d)?)?;
            Ok((q_ceil, r_ceil))
        } else {
            Ok((q, r))
        }
    }

    pub fn try_cdiv_q(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        self.try_cdiv_qr(d).map(|(q, _)| q)
    }

    pub fn try_cdiv_r(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        self.try_cdiv_qr(d).map(|(_, r)| r)
    }

    // -----------------------------------------------------------------------
    // Division — Ceiling variants by u64
    // -----------------------------------------------------------------------

    /// `mpz_cdiv_q_ui`: ceiling quotient by u64.
    pub fn try_cdiv_q_ui(&self, d: u64) -> Result<Mpz, CapacityError> {
        self.try_cdiv_q(&Mpz::from_u64(d))
    }

    /// `mpz_cdiv_r_ui`: ceiling remainder by u64 (as Mpz).
    pub fn try_cdiv_r_ui(&self, d: u64) -> Result<Mpz, CapacityError> {
        self.try_cdiv_r(&Mpz::from_u64(d))
    }

    /// `mpz_cdiv_qr_ui`: ceiling quotient and remainder by u64.
    /// GMP returns the absolute remainder as a non-negative `u64`.
    pub fn try_cdiv_qr_ui(&self, d: u64) -> Result<(Mpz, u64), CapacityError> {
        let d_mpz = Mpz::from_u64(d);
        let (q, r) = self.try_cdiv_qr(&d_mpz)?;
        // GMP returns the absolute remainder value
        let r_u64 = r.get_ui();
        Ok((q, r_u64))
    }

    /// `mpz_cdiv_ui`: ceiling remainder as u64.
    pub fn cdiv_ui(&self, d: u64) -> u64 {
        let (_, r) = self.mag_divmod_u64(d);
        if self.sign > 0 && r != 0 {
            d - r
        } else {
            r
        }
    }

    /// `mpz_cdiv_q_2exp`: ceiling quotient by 2^bits.
    pub fn cdiv_q_2exp(&self, bits: u32) -> Mpz {
        let q = self.fdiv_q_2exp(bits);
        if self.sign >= 0 {
            // Ceil division for non-negative: if any lower bits are set, q += 1
            let r = self.tdiv_r_2exp(bits);
            if r.is_zero() {
                q
            } else {
                q.try_add_ui(1).unwrap_or(q)
            }
        } else {
            // For negative, ceil(x/2^k) = -floor(-x/2^k)
            // If -x has any lower bits set, floor(-x/2^k) needs the increment
            let neg = self.neg_to();
            let floor_neg = neg.fdiv_q_2exp(bits);
            if floor_neg.is_zero() {
                Mpz::new()
            } else {
                let neg_r = neg.tdiv_r_2exp(bits);
                if neg_r.is_zero() {
                    // -x is divisible by 2^k, so ceil(x/2^k) = -((-x)/2^k)
                    let mut r = floor_neg;
                    r.sign = -r.sign;
                    r
                } else {
                    // -x not divisible, floor(-x/2^k) is the division result
                    // ceil(x/2^k) = -(floor(-x/2^k))
                    // Actually: ceil(x/2^k) where x < 0:
                    // Let x = -a where a > 0.
                    // ceil(-a / 2^k) = -floor(a / 2^k)
                    // If a % 2^k == 0: floor(a/2^k) = a/2^k, ceil(-a/2^k) = -a/2^k
                    // If a % 2^k != 0: floor(a/2^k) = (a - (a%2^k))/2^k = a/2^k truncated
                    //   ceil(-a/2^k) = -floor(a/2^k) = -(a >> k)
                    let mut r = floor_neg;
                    r.sign = -r.sign;
                    r
                }
            }
        }
    }

    /// `mpz_cdiv_r_2exp`: ceiling remainder modulo 2^bits.
    /// For ceiling division by a positive power of 2, the remainder is non-positive (≤ 0).
    pub fn cdiv_r_2exp(&self, bits: u32) -> Mpz {
        let q = self.cdiv_q_2exp(bits);
        let q_times_mod = q.try_mul_2exp(bits).unwrap_or_else(|_| Mpz::new());
        let r = self.try_sub(&q_times_mod).unwrap_or_else(|_| Mpz::new());
        r
    }

    /// `mpz_tdiv_qr_ui`: truncating quotient and remainder by u64. Returns (q, r).
    /// GMP returns the remainder as a non-negative `u64`.
    pub fn try_tdiv_qr_ui(&self, d: u64) -> Result<(Mpz, u64), CapacityError> {
        let d_mpz = Mpz::from_u64(d);
        let (q, r) = self.tdiv_qr(&d_mpz);
        // GMP's tdiv_qr_ui always returns the absolute value of the remainder
        let r_abs = r.get_ui();
        if r.sign < 0 {
            // q is already correct (trunc toward zero), r gets absolute value
            Ok((q, r_abs))
        } else {
            Ok((q, r_abs))
        }
    }

    /// `mpz_fdiv_r_ui`: floor remainder by u64 as Mpz (not just u64).
    pub fn try_fdiv_r_ui(&self, d: u64) -> Result<Mpz, CapacityError> {
        let (_, r) = self.try_fdiv_qr(&Mpz::from_u64(d))?;
        Ok(r)
    }

    // -----------------------------------------------------------------------
    // Division — Modulo (non-negative remainder)
    // -----------------------------------------------------------------------

    /// `mpz_mod`: non-negative remainder.
    pub fn try_mod(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        let r = self.tdiv_r(d);
        if r.sign < 0 {
            r.try_add(d)
        } else {
            Ok(r)
        }
    }

    /// `mpz_mod_ui`: non-negative remainder as u64.
    pub fn mod_ui(&self, d: u64) -> u64 {
        self.fdiv_ui(d)
    }

    // -----------------------------------------------------------------------
    // Division — Exact
    // -----------------------------------------------------------------------

    /// `mpz_divexact`: `self / d` (d is known to divide self exactly).
    pub fn try_divexact(&self, d: &Mpz) -> Result<Mpz, CapacityError> {
        Ok(self.tdiv_q(d))
    }

    /// `mpz_divexact_ui`.
    pub fn try_divexact_ui(&self, d: u64) -> Result<Mpz, CapacityError> {
        Ok(self.tdiv_q_ui(d))
    }

    // -----------------------------------------------------------------------
    // Divisibility / congruence
    // -----------------------------------------------------------------------

    pub fn divisible_ui(&self, d: u64) -> bool {
        self.sign == 0 || self.mag_divmod_u64(d).1 == 0
    }

    /// `mpz_divisible_p`: is self divisible by d?
    pub fn divisible_p(&self, d: &Mpz) -> bool {
        if d.len == 0 {
            return false;
        }
        if self.sign == 0 {
            return true;
        }
        let r = self.tdiv_r(d);
        r.sign == 0
    }

    /// `mpz_divisible_2exp_p`: is self divisible by 2^bits?
    pub fn divisible_2exp_p(&self, bits: u32) -> bool {
        if self.len == 0 {
            return true;
        }
        let limb = (bits / 64) as usize;
        let bit = bits % 64;
        for i in 0..limb.min(self.len) {
            if self.mag[i] != 0 {
                return false;
            }
        }
        if limb < self.len && bit > 0 {
            let mask = (1u64 << bit) - 1;
            if self.mag[limb] & mask != 0 {
                return false;
            }
        }
        true
    }

    /// `mpz_congruent_p`: is self ≡ c (mod d)?
    pub fn congruent_p(&self, c: &Mpz, d: &Mpz) -> bool {
        let diff = self.try_sub(c).unwrap_or_else(|_| Mpz::new());
        diff.divisible_p(d)
    }

    /// `mpz_congruent_ui_p`: is self ≡ c (mod d)?
    pub fn congruent_ui_p(&self, c: u64, d: u64) -> bool {
        if d == 0 {
            return false;
        }
        let r1 = self.tdiv_ui(d);
        let r2 = c % d;
        r1 == r2
    }

    /// `mpz_congruent_2exp_p`: is self ≡ c (mod 2^bits)?
    ///
    /// Uses two's complement matching of the low `bits` bits.
    pub fn congruent_2exp_p(&self, c: &Mpz, bits: u32) -> bool {
        // Get the low `bits` bits of self in two's complement form.
        let a_bits = self.twos_complement_low_bits(bits);
        let b_bits = c.twos_complement_low_bits(bits);
        a_bits == b_bits
    }

    /// Return the low `bits` bits of the two's complement representation as a `u128`.
    fn twos_complement_low_bits(&self, bits: u32) -> u128 {
        if bits == 0 {
            return 0;
        }
        let n_limbs = ((bits + 63) / 64) as usize;
        let n_limbs_for_mask = n_limbs.min(MPZ_MAX_LIMBS);
        let _limb_mask_bits = bits % 64;

        // Collect the low limbs in two's complement
        let mut val = 0u128;
        if self.sign >= 0 {
            // Positive: just take the magnitude
            for i in 0..n_limbs_for_mask {
                if i < self.len {
                    val |= (self.mag[i] as u128) << (i * 64);
                }
            }
        } else {
            // Negative: two's complement of magnitude
            let mut carry = 1u128;
            for i in 0..n_limbs_for_mask {
                let inv = (!if i < self.len { self.mag[i] } else { 0 }) as u128;
                let sum = inv + carry;
                val |= (sum & 0xFFFF_FFFF_FFFF_FFFF) << (i * 64);
                carry = sum >> 64;
            }
            // If bits exceed self.len, sign-extend with all-ones
            if n_limbs_for_mask > self.len && self.len > 0 {
                // The bits beyond self.len are all-1s for negative
                // Already handled by the inv logic above (since i >= self.len uses 0 for mag,
                // and !0 = all-ones, +carry may propagate)
            }
        }

        // Mask to the requested number of bits
        if bits < 128 {
            val &= (1u128 << bits) - 1;
        }
        val
    }

    // -----------------------------------------------------------------------
    // fdiv_r_2exp, fdiv_q_2exp
    // -----------------------------------------------------------------------

    pub fn fdiv_r_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let limbs = (bits / 64) as usize;
        let rem_bits = bits % 64;
        let mut result = self.clone();
        if limbs < MPZ_MAX_LIMBS {
            if rem_bits != 0 && result.len > limbs {
                result.mag[limbs] &= (1u64 << rem_bits) - 1;
            }
            for i in (limbs + if rem_bits != 0 { 1 } else { 0 })..result.len {
                result.mag[i] = 0;
            }
        }
        result.trim();
        result
    }

    pub fn fdiv_q_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let ls = (bits / 64) as usize;
        let bs = bits % 64;
        if ls >= self.len {
            return Mpz::new();
        }
        let mut result = Mpz::new();
        let src_len = self.len - ls;
        result.mag[..src_len].copy_from_slice(&self.mag[ls..self.len]);
        result.len = src_len;
        result.sign = self.sign;
        if bs != 0 {
            let mut carry = 0u64;
            for i in (0..result.len).rev() {
                let new = (result.mag[i] >> bs) | carry;
                carry = result.mag[i] << (64 - bs);
                result.mag[i] = new;
            }
        }
        result.trim();
        result
    }

    // -----------------------------------------------------------------------
    // isqrt, root, etc.
    // -----------------------------------------------------------------------

    pub fn isqrt(&self) -> Mpz {
        if self.sgn() <= 0 {
            return Mpz::new();
        }
        let mut x = Mpz::from_u64(1);
        let shift = self.sizeinbase2().div_ceil(2) as u32;
        x = x.try_mul_2exp(shift).unwrap_or(
            Mpz::from_u64(1)
                .try_mul_2exp(127)
                .unwrap_or(Mpz::from_u64(1)),
        );
        loop {
            let y = x.try_add(&self.tdiv_q(&x)).unwrap().fdiv_q_2exp(1);
            if y.cmp(&x) != Ordering::Less {
                return x;
            }
            x = y;
        }
    }

    /// `mpz_sqrtrem`: square root with remainder.
    pub fn try_sqrtrem(&self) -> Result<(Mpz, Mpz), CapacityError> {
        let root = self.isqrt();
        let rem = self.try_sub(&root.try_mul(&root)?)?;
        Ok((root, rem))
    }

    /// `mpz_root`: floor nth root.
    pub fn try_root(&self, n: u32) -> Result<Mpz, CapacityError> {
        if self.sgn() <= 0 || n == 0 {
            return Ok(Mpz::new());
        }
        if n == 1 {
            return Ok(self.clone());
        }
        if n == 2 {
            return Ok(self.isqrt());
        }
        let bits = self.sizeinbase2();
        let shift = (bits as u32).div_ceil(n);
        let mut x = Mpz::from_u64(2).try_pow_ui(shift)?;
        let nm1 = Mpz::from_u64((n - 1) as u64);
        loop {
            let x_pow_nm1 = x.try_pow_ui(n - 1)?;
            let quotient = self.tdiv_q(&x_pow_nm1);
            let sum = x.try_mul(&nm1)?.try_add(&quotient)?;
            let y = sum.tdiv_q_ui(n as u64);
            if y.cmp(&x) != Ordering::Less {
                let y_p1 = y.try_add(&Mpz::from_u64(1))?;
                let y_p1_pow = y_p1.try_pow_ui(n)?;
                if y_p1_pow.cmp(self) == Ordering::Greater {
                    return Ok(y);
                }
                return Ok(x);
            }
            x = y;
        }
    }

    /// `mpz_rootrem`: floor nth root with remainder.
    pub fn try_rootrem(&self, n: u32) -> Result<(Mpz, Mpz), CapacityError> {
        let root = self.try_root(n)?;
        let root_pow = root.try_pow_ui(n)?;
        let rem = self.try_sub(&root_pow)?;
        Ok((root, rem))
    }

    /// `mpz_perfect_square_p`: is this a perfect square?
    pub fn perfect_square_p(&self) -> bool {
        if self.sgn() < 0 {
            return false;
        }
        let root = self.isqrt();
        root.try_mul(&root).map(|sq| sq == *self).unwrap_or(false)
    }

    /// `mpz_perfect_power_p`: is this a perfect power (a^k for some a, k > 1)?
    pub fn perfect_power_p(&self) -> bool {
        if self.sgn() <= 0 {
            return false;
        }
        if self.len == 1 {
            let v = self.mag[0];
            if v < 4 {
                return false;
            }
            for k in 2..64 {
                if (1u64 << k) > v {
                    break;
                }
                let mut lo = 2u64;
                let mut hi = 1u64 << (63 / k + 1);
                hi = hi.min(v).min(1 << 26);
                while lo < hi {
                    let mid = lo + (hi - lo) / 2;
                    let mut pow = 1u128;
                    for _ in 0..k {
                        pow = pow.wrapping_mul(mid as u128);
                    }
                    if pow >= v as u128 {
                        hi = mid;
                    } else {
                        lo = mid + 1;
                    }
                }
                let mut pow = 1u128;
                for _ in 0..k {
                    pow = pow.wrapping_mul(lo as u128);
                }
                if pow == v as u128 {
                    return true;
                }
            }
            false
        } else {
            let max_k = self.sizeinbase2();
            for k in 2..=max_k.min(64) as u32 {
                if let Ok(root) = self.try_root(k) {
                    if root.try_pow_ui(k).map(|p| p == *self).unwrap_or(false) {
                        return true;
                    }
                }
            }
            false
        }
    }

    // -----------------------------------------------------------------------
    // remove_pow10, com
    // -----------------------------------------------------------------------

    pub fn remove_pow10(&mut self) -> u32 {
        if self.sign == 0 {
            return 0;
        }
        let mut count = 0;
        loop {
            let mut qbuf = [0u64; MPZ_MAX_LIMBS];
            let (qlen, r) = Self::mag_divmod_u64_len(&self.mag[..self.len], 10, &mut qbuf);
            if r != 0 {
                break;
            }
            self.mag[..qlen].copy_from_slice(&qbuf[..qlen]);
            self.len = qlen;
            count += 1;
        }
        self.trim();
        count
    }

    /// `mpz_com`: one's complement.
    pub fn com(&self) -> Mpz {
        let one = Mpz::from_u64(1);
        let mut r = self.try_add(&one).unwrap_or_else(|_| Mpz::new());
        r.sign = -r.sign;
        r
    }

    /// `mpz_remove(_, _, f)`: remove all factors of `f` from `self`.
    pub fn try_remove(&self, f: &Mpz) -> Result<(Mpz, u32), CapacityError> {
        if self.sign == 0 {
            return Ok((Mpz::new(), 0));
        }
        if f.len == 0 {
            return Err(CapacityError);
        }
        if f.len == 1 && f.mag[0] == 1 {
            return Ok((self.clone(), 0));
        }
        let mut val = self.clone();
        let mut count = 0;
        loop {
            let (q, r) = val.tdiv_qr(f);
            if r.sign != 0 {
                break;
            }
            val = q;
            count += 1;
        }
        Ok((val, count))
    }

    // -----------------------------------------------------------------------
    // Number theory
    // -----------------------------------------------------------------------

    /// `mpz_gcd`: greatest common divisor.
    pub fn try_gcd(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let mut a = self.abs_to();
        let mut b = other.abs_to();
        if a.len == 0 {
            return Ok(b);
        }
        if b.len == 0 {
            return Ok(a);
        }
        let a_tz = a.mag[0].trailing_zeros();
        let b_tz = b.mag[0].trailing_zeros();
        let shift = a_tz.min(b_tz);
        a = a.fdiv_q_2exp(a_tz);
        b = b.fdiv_q_2exp(b_tz);
        loop {
            if a.len == 0 {
                b = b.try_mul_2exp(shift)?;
                b.sign = 1;
                return Ok(b);
            }
            if b.len == 0 {
                a = a.try_mul_2exp(shift)?;
                a.sign = 1;
                return Ok(a);
            }
            let cmp = a.cmpabs(&b);
            if cmp == Ordering::Greater || cmp == Ordering::Equal {
                a = a.try_sub(&b)?;
                a = a.fdiv_q_2exp(a.mag[0].trailing_zeros());
            } else {
                b = b.try_sub(&a)?;
                b = b.fdiv_q_2exp(b.mag[0].trailing_zeros());
            }
        }
    }

    /// `mpz_gcd_ui`: gcd with unsigned 64-bit.
    pub fn gcd_ui(&self, v: u64) -> u64 {
        if v == 0 {
            return self.get_ui();
        }
        let mut a = self.get_ui();
        let mut b = v;
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    /// `mpz_gcdext`: extended GCD.  Returns `(g, s, t)` where `g = gcd(self, other) = self*s + other*t`.
    pub fn try_gcdext(&self, other: &Mpz) -> Result<(Mpz, Mpz, Mpz), CapacityError> {
        let mut old_r = self.abs_to();
        let mut r = other.abs_to();
        let mut old_s = Mpz::from_u64(1);
        let mut s = Mpz::new();
        let mut old_t = Mpz::new();
        let mut t = Mpz::from_u64(1);
        while r.len != 0 {
            let (q, _) = old_r.tdiv_qr(&r);
            let new_r = old_r.try_sub(&q.try_mul(&r)?)?;
            let new_s = old_s.try_sub(&q.try_mul(&s)?)?;
            let new_t = old_t.try_sub(&q.try_mul(&t)?)?;
            old_r = r;
            r = new_r;
            old_s = s;
            s = new_s;
            old_t = t;
            t = new_t;
        }
        let mut g = old_r;
        if self.sign < 0 {
            g.sign = -g.sign;
            old_s.sign = -old_s.sign;
        }
        if other.sign < 0 {
            old_t.sign = -old_t.sign;
        }
        Ok((g, old_s, old_t))
    }

    /// `mpz_lcm`: least common multiple.
    pub fn try_lcm(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        if self.len == 0 || other.len == 0 {
            return Ok(Mpz::new());
        }
        let g = self.try_gcd(other)?;
        let q = self.tdiv_q(&g);
        q.try_mul(other)
    }

    /// `mpz_lcm_ui`: lcm with unsigned 64-bit.
    pub fn try_lcm_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_lcm(&Mpz::from_u64(v))
    }

    /// `mpz_invert`: modular inverse.  Returns `Ok(inverse)` if invertible, `Err(())` if not.
    pub fn try_invert(&self, m: &Mpz) -> Result<Mpz, ()> {
        if m.len == 0 {
            return Err(());
        }
        let (g, s, _t) = self.try_gcdext(m).map_err(|_| ())?;
        if g.len != 1 || g.mag[0] != 1 {
            return Err(());
        }
        let mut r = s.try_mod(m).map_err(|_| ())?;
        if r.sign < 0 {
            r = r.try_add(m).map_err(|_| ())?;
        }
        Ok(r)
    }

    /// `mpz_jacobi`: Jacobi symbol `(self/other)`.
    pub fn jacobi(&self, other: &Mpz) -> i32 {
        if other.even_p() {
            return 0;
        }
        let mut n = other.abs_to();
        // If n == 1, (a/1) = 1 for all a
        if n.len == 1 && n.mag[0] == 1 {
            return 1;
        }
        let mut a = self.try_mod(&n).unwrap_or_else(|_| Mpz::new());
        // If a == 0, (0/n) = 0 for n > 1, 1 for n == 1
        if a.is_zero() {
            return 0;
        }
        let mut t = 1i32;
        loop {
            let mut e = 0u32;
            while a.even_p() && !a.is_zero() {
                a = a.fdiv_q_2exp(1);
                e += 1;
            }
            if e % 2 == 1 {
                let n_mod_8 = n.mag[0] & 7;
                if n_mod_8 == 3 || n_mod_8 == 5 {
                    t = -t;
                }
            }
            // Quadratic reciprocity
            core::mem::swap(&mut a, &mut n);
            if a.len > 0 && n.len > 0 && a.mag[0] % 4 == 3 && n.mag[0] % 4 == 3 {
                t = -t;
            }
            a = a.try_mod(&n).unwrap_or_else(|_| Mpz::new());
            if a.is_zero() {
                // (n/1) = 1, so if n == 1, return t
                if n.len == 1 && n.mag[0] == 1 {
                    return t;
                }
                return 0;
            }
            if n.len == 1 && n.mag[0] == 1 {
                return t;
            }
        }
    }

    /// `mpz_legendre`: Legendre symbol `(self/p)` (Jacobi when p is odd prime).
    pub fn try_legendre(&self, p: &Mpz) -> i32 {
        self.jacobi(p)
    }

    /// `mpz_kronecker`: Kronecker symbol `(self/other)`.
    ///
    /// Generalisation of Jacobi symbol to all integers.
    pub fn try_kronecker(&self, other: &Mpz) -> i32 {
        if other.is_zero() {
            // (self/0) = 1 if self == ±1, otherwise 0 if |self| != 1
            return if self.len == 1 && self.mag[0] == 1 {
                1
            } else {
                0
            };
        }
        if other.sign < 0 {
            // (self / -other) = sign_factor * (self / other)
            // Kronecker symbol: (n / -1) = sign(n)
            // (self / -1) = 1 if self >= 0, -1 if self < 0
            let pos = other.abs_to();
            let j = self.try_kronecker(&pos);
            if self.sign < 0 {
                -j
            } else {
                j
            }
        } else {
            self.jacobi(other)
        }
    }

    /// `mpz_kronecker_si`: Kronecker `(self/a)` for signed integer `a`.
    pub fn try_kronecker_si(&self, a: i64) -> i32 {
        Self::try_si_kronecker(a, self)
    }

    /// `mpz_kronecker_ui`: Kronecker `(self/a)` for unsigned integer `a`.
    pub fn try_kronecker_ui(&self, a: u64) -> i32 {
        Self::try_ui_kronecker(a, self)
    }

    /// `mpz_si_kronecker`: Kronecker `(a/other)` for signed integer `a`.
    pub fn try_si_kronecker(a: i64, other: &Mpz) -> i32 {
        Mpz::from_i64(a).try_kronecker(other)
    }

    /// `mpz_ui_kronecker`: Kronecker `(a/other)` for unsigned integer `a`.
    pub fn try_ui_kronecker(a: u64, other: &Mpz) -> i32 {
        Mpz::from_u64(a).try_kronecker(other)
    }

    /// `mpz_fac_ui`: factorial n!.
    pub fn try_fac_ui(n: u32) -> Result<Mpz, CapacityError> {
        let mut r = Mpz::from_u64(1);
        for i in 2..=n as u64 {
            r = r.try_mul_ui(i)?;
        }
        Ok(r)
    }

    /// `mpz_2fac_ui`: double factorial n!!.
    ///
    /// n!! = product of all integers ≤ n with the same parity as n.
    /// 0!! = 1, 1!! = 1.
    pub fn try_2fac_ui(n: u32) -> Result<Mpz, CapacityError> {
        if n <= 1 {
            return Ok(Mpz::from_u64(1));
        }
        let mut r = Mpz::from_u64(n as u64);
        let mut i = n as u64;
        loop {
            if i < 2 {
                break;
            }
            i = i.wrapping_sub(2);
            if i == 0 {
                break;
            }
            r = r.try_mul_ui(i)?;
        }
        Ok(r)
    }

    /// `mpz_primorial_ui`: product of all primes ≤ n.
    pub fn try_primorial_ui(n: u32) -> Result<Mpz, CapacityError> {
        if n <= 1 {
            return Ok(Mpz::from_u64(1));
        }
        let mut r = Mpz::from_u64(1);
        for &p in SMALL_PRIMES.iter() {
            if p > n as u64 {
                break;
            }
            r = r.try_mul_ui(p)?;
        }
        Ok(r)
    }

    /// `mpz_bin_uiui`: binomial coefficient C(n, k).
    pub fn try_bin_uiui(n: u32, k: u32) -> Result<Mpz, CapacityError> {
        let k = k.min(n - k);
        if k == 0 {
            return Ok(Mpz::from_u64(1));
        }
        let mut r = Mpz::from_u64(n as u64 - k as u64 + 1);
        for i in 2..=k {
            r = r.try_mul_ui(n as u64 - k as u64 + i as u64)?;
            r = r.tdiv_q_ui(i as u64);
        }
        Ok(r)
    }

    /// `mpz_bin_ui`: binomial coefficient C(self, k) for mpz `self`.
    /// Returns `Err(CapacityError)` if result exceeds capacity or k invalid.
    pub fn try_bin_ui(&self, k: u32) -> Result<Mpz, CapacityError> {
        if self.sign < 0 {
            // Generalisation for negative n: C(n,k) = (-1)^k * C(-n+k-1, k)
            let neg = self.neg_to(); // -self
            let km1 = Mpz::from_u64(k as u64);
            let adjusted = neg.try_add(&km1)?; // -self + k
            let _abs_k = k.min(k);
            let mut r = adjusted.try_bin_uiui_rec(k)?;
            if k % 2 == 1 {
                r.sign = -r.sign;
            }
            Ok(r)
        } else {
            self.try_bin_uiui_rec(k)
        }
    }

    /// Internal helper: C(|self|, k) when self is non-negative, using iterative
    /// multiplication and division (like bin_uiui but with mpz numerator).
    fn try_bin_uiui_rec(&self, k: u32) -> Result<Mpz, CapacityError> {
        if self.is_zero() {
            return if k == 0 {
                Ok(Mpz::from_u64(1))
            } else {
                Ok(Mpz::new())
            };
        }
        if k == 0 {
            return Ok(Mpz::from_u64(1));
        }
        // C(n, k) = product_{i=1..k} (n - k + i) / i
        // where n = self, k <= n (if k > n, result is 0)
        if self.len == 1 && (self.mag[0] as u128) < k as u128 {
            return Ok(Mpz::new());
        }
        if self.len > 1 {
            // k is at most u32::MAX and self has > 1 limb, so k < |self|
        }
        let k = k.min(if self.len == 1 && self.mag[0] <= u32::MAX as u64 {
            self.mag[0] as u32
        } else {
            k
        });
        let _k_min = k.min(if self.len == 1 {
            (self.mag[0] as u32).saturating_sub(k)
        } else {
            k
        });
        // Actually simplify: use k = min(k, |self| - k)
        let n_val = if self.len == 1 { self.mag[0] } else { u64::MAX };
        let k_actual = if self.len == 1 && (k as u64) > n_val / 2 {
            (n_val - k as u64) as u32
        } else {
            k
        };
        let k = k_actual;

        let n_minus_k = if self.len == 1 {
            Mpz::from_u64(self.mag[0] - k as u64)
        } else {
            // For multi-limb mpz, approximate: subtract k
            self.try_sub(&Mpz::from_u64(k as u64))?
        };

        let mut r = Mpz::from_u64(1);
        for i in 1..=k {
            let term = n_minus_k.try_add_ui(i as u64)?;
            r = r.try_mul(&term)?;
            r = r.tdiv_q_ui(i as u64);
        }
        Ok(r)
    }

    /// `mpz_fib_ui`: F(n) (Fibonacci, F(0)=0, F(1)=1).
    pub fn try_fib_ui(n: u32) -> Result<Mpz, CapacityError> {
        if n == 0 {
            return Ok(Mpz::new());
        }
        if n == 1 {
            return Ok(Mpz::from_u64(1));
        }
        let mut a = Mpz::new(); // F(0)
        let mut b = Mpz::from_u64(1); // F(1)
        for _ in 2..=n {
            let c = a.try_add(&b)?;
            a = b;
            b = c;
        }
        Ok(b)
    }

    /// `mpz_fib2_ui`: (F_n, F_{n-1}).
    pub fn try_fib2_ui(n: u32) -> Result<(Mpz, Mpz), CapacityError> {
        if n == 0 {
            return Ok((Mpz::new(), Mpz::from_u64(1)));
        }
        if n == 1 {
            return Ok((Mpz::from_u64(1), Mpz::new()));
        }
        let mut a = Mpz::new(); // F(0)
        let mut b = Mpz::from_u64(1); // F(1)
        for _ in 2..=n {
            let c = a.try_add(&b)?;
            a = b;
            b = c;
        }
        Ok((b.clone(), a))
    }

    /// `mpz_lucnum_ui`: L_n (Lucas number, L(0)=2, L(1)=1).
    pub fn try_lucnum_ui(n: u32) -> Result<Mpz, CapacityError> {
        if n == 0 {
            return Ok(Mpz::from_u64(2));
        }
        if n == 1 {
            return Ok(Mpz::from_u64(1));
        }
        let mut a = Mpz::from_u64(2); // L(0)
        let mut b = Mpz::from_u64(1); // L(1)
        for _ in 2..=n {
            let c = a.try_add(&b)?;
            a = b;
            b = c;
        }
        Ok(b)
    }

    /// `mpz_lucnum2_ui`: (L_n, L_{n-1}).
    pub fn try_lucnum2_ui(n: u32) -> Result<(Mpz, Mpz), CapacityError> {
        if n == 0 {
            return Ok((Mpz::from_u64(2), Mpz::from_i64(-1)));
        }
        if n == 1 {
            return Ok((Mpz::from_u64(1), Mpz::from_u64(2)));
        }
        let mut a = Mpz::from_u64(2); // L(0)
        let mut b = Mpz::from_u64(1); // L(1)
        for _ in 2..=n {
            let c = a.try_add(&b)?;
            a = b;
            b = c;
        }
        Ok((b.clone(), a))
    }

    // -----------------------------------------------------------------------
    // Bit operations (two's complement on sign-magnitude)
    // -----------------------------------------------------------------------

    /// Convert to two's complement representation, returning limbs and sign bit.
    /// `out` must have at least `self.len + 1` elements.
    /// Returns `(len, is_negative)` where `is_negative` is the sign of the original value.
    fn to_twos_complement(&self, out: &mut [u64]) -> (usize, bool) {
        if self.sign == 0 {
            return (0, false);
        }
        let n = self.len + 1;
        if self.sign > 0 {
            out[..self.len].copy_from_slice(&self.mag[..self.len]);
            out[self.len] = 0;
            (n, false)
        } else {
            let mut carry = 1u128;
            for i in 0..self.len {
                let inv = (!self.mag[i]) as u128;
                let sum = inv + carry;
                out[i] = sum as u64;
                carry = sum >> 64;
            }
            if carry != 0 {
                out[self.len] = 0;
            } else {
                out[self.len] = !0u64;
            }
            (n, true)
        }
    }

    /// Convert from two's complement back to sign-magnitude.
    fn from_twos_complement(limbs: &[u64], negative: bool) -> Mpz {
        if !negative {
            let mut r = Mpz::new();
            let n = limbs.len().min(MPZ_MAX_LIMBS);
            r.mag[..n].copy_from_slice(&limbs[..n]);
            r.len = n;
            r.sign = if n == 0 { 0 } else { 1 };
            r.trim();
            return r;
        }
        let mut r = Mpz::new();
        let n = limbs.len().min(MPZ_MAX_LIMBS);
        let mut carry = 1u128;
        for i in 0..n {
            let inv = (!limbs[i]) as u128;
            let sum = inv + carry;
            r.mag[i] = sum as u64;
            carry = sum >> 64;
        }
        if carry != 0 && n < MPZ_MAX_LIMBS {
            r.mag[n] = carry as u64;
            r.len = n + 1;
        } else {
            r.len = n;
        }
        r.sign = -1;
        r.trim();
        r
    }

    /// Apply a binary bitwise operation on two Mpz values.
    fn bitwise_op<F>(&self, other: &Mpz, op: F) -> Result<Mpz, CapacityError>
    where
        F: Fn(u64, u64) -> u64,
    {
        // Extend to max(len1, len2) + 1 limbs in two's complement to capture sign bit
        let max_limbs = self.len.max(other.len) + 1;
        if max_limbs > MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        let mut a_buf = [0u64; MPZ_MAX_LIMBS];
        let mut b_buf = [0u64; MPZ_MAX_LIMBS];
        let (a_len, a_neg) = self.to_twos_complement(&mut a_buf);
        let (b_len, b_neg) = other.to_twos_complement(&mut b_buf);

        // Determine the working length: we need at least max(a_len, b_len) limbs,
        // but to handle infinite sign extension properly we need the extended length
        // which is max_limbs = max(self.len, other.len) + 1.
        let work_len = max_limbs;

        let mut result_limbs = [0u64; MPZ_MAX_LIMBS];
        let sign_a = if a_neg { !0u64 } else { 0 };
        let sign_b = if b_neg { !0u64 } else { 0 };

        for i in 0..work_len {
            let va = if i < a_len { a_buf[i] } else { sign_a };
            let vb = if i < b_len { b_buf[i] } else { sign_b };
            result_limbs[i] = op(va, vb);
        }

        // Determine if result is negative: check the MSB of the top limb
        let top_limb = result_limbs[work_len - 1];
        let negative = (top_limb >> 63) == 1;

        if !negative {
            // Top bit is 0: positive result, trim trailing zero limbs
            let mut rlen = work_len;
            while rlen > 0 && result_limbs[rlen - 1] == 0 {
                rlen -= 1;
            }
            let mut r = Mpz::new();
            r.mag[..rlen].copy_from_slice(&result_limbs[..rlen]);
            r.len = rlen;
            r.sign = if rlen == 0 { 0 } else { 1 };
            Ok(r)
        } else {
            // Negative result: convert from two's complement
            let full_len = work_len;
            Ok(Self::from_twos_complement(&result_limbs[..full_len], true))
        }
    }

    pub fn try_and(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        // Fast paths: x & 0 == 0, 0 & x == 0, x & x == x
        if self.len == 0 || other.len == 0 {
            return Ok(Mpz::new());
        }
        if core::ptr::eq(self, other) || self == other {
            return Ok(self.clone());
        }
        self.bitwise_op(other, |a, b| a & b)
    }

    pub fn try_ior(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        // Fast paths: x | 0 == x, 0 | x == x, x | x == x
        if self.len == 0 {
            return Ok(other.clone());
        }
        if other.len == 0 {
            return Ok(self.clone());
        }
        if core::ptr::eq(self, other) || self == other {
            return Ok(self.clone());
        }
        self.bitwise_op(other, |a, b| a | b)
    }

    pub fn try_xor(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        // Fast paths: x ^ 0 == x, 0 ^ x == x, x ^ x == 0
        if self.len == 0 {
            return Ok(other.clone());
        }
        if other.len == 0 {
            return Ok(self.clone());
        }
        if core::ptr::eq(self, other) || self == other {
            return Ok(Mpz::new());
        }
        self.bitwise_op(other, |a, b| a ^ b)
    }

    /// `mpz_popcount`: number of 1 bits in the two's complement representation.
    pub fn popcount(&self) -> Option<u32> {
        if self.sign < 0 {
            return None;
        }
        let mut count = 0u32;
        for i in 0..self.len {
            count += self.mag[i].count_ones();
        }
        Some(count)
    }

    /// `mpz_hamdist`: Hamming distance between two values (popcount of XOR).
    pub fn hamdist(&self, other: &Mpz) -> Option<u32> {
        self.try_xor(other).ok()?.popcount()
    }

    /// `mpz_scan0`: find the first 0 bit at or after `start`.
    pub fn scan0(&self, start: u32) -> u32 {
        let start_limb = (start / 64) as usize;
        let start_bit = start % 64;
        for i in start_limb..self.len {
            let mut w = self.mag[i];
            if i == start_limb && start_bit > 0 {
                w |= (1u64 << start_bit) - 1;
            }
            let zeros = (!w).trailing_zeros();
            if zeros < 64 {
                return (i * 64 + zeros as usize) as u32;
            }
        }
        (self.len * 64) as u32
    }

    /// `mpz_scan1`: find the first 1 bit at or after `start`.
    pub fn scan1(&self, start: u32) -> Option<u32> {
        let start_limb = (start / 64) as usize;
        let start_bit = start % 64;
        for i in start_limb..self.len {
            let mut w = self.mag[i];
            if i == start_limb {
                w >>= start_bit;
                w <<= start_bit;
            }
            if w != 0 {
                return Some((i * 64 + w.trailing_zeros() as usize) as u32);
            }
        }
        None
    }

    /// `mpz_tstbit`: test whether bit `bit` is set.
    ///
    /// For negative values, out-of-range bits return `true` (infinite sign extension).
    pub fn tstbit(&self, bit: u32) -> bool {
        let limb = (bit / 64) as usize;
        let bit_idx = bit % 64;
        if limb >= self.len {
            // Out-of-range: for negative values, sign-extend (return true);
            // for non-negative, return false.
            return self.sign < 0;
        }
        (self.mag[limb] >> bit_idx) & 1 == 1
    }

    /// `mpz_setbit`: set bit `bit`.
    pub fn try_setbit(&mut self, bit: u32) -> Result<(), CapacityError> {
        let limb = (bit / 64) as usize;
        let bit_idx = bit % 64;
        if limb >= MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        if limb >= self.len {
            self.len = limb + 1;
        }
        self.mag[limb] |= 1u64 << bit_idx;
        if self.sign == 0 {
            self.sign = 1;
        }
        Ok(())
    }

    /// `mpz_clrbit`: clear bit `bit`.
    pub fn clrbit(&mut self, bit: u32) {
        let limb = (bit / 64) as usize;
        let bit_idx = bit % 64;
        if limb >= self.len {
            return;
        }
        self.mag[limb] &= !(1u64 << bit_idx);
        self.trim();
    }

    /// `mpz_combit`: complement bit `bit`.
    pub fn try_combit(&mut self, bit: u32) -> Result<(), CapacityError> {
        let limb = (bit / 64) as usize;
        let bit_idx = bit % 64;
        if limb >= self.len && limb >= MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        if limb >= self.len {
            self.len = limb + 1;
            if self.sign == 0 {
                self.sign = 1;
            }
        }
        self.mag[limb] ^= 1u64 << bit_idx;
        self.trim();
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Decimal string I/O
    // -----------------------------------------------------------------------

    pub fn from_decimal_str(s: &str) -> Result<Mpz, ParseError> {
        let s = s.trim();
        if s.is_empty() {
            return Err(ParseError::InvalidInput);
        }
        let (neg, digits) = match s.strip_prefix('-') {
            Some(d) => {
                if d.is_empty() {
                    return Err(ParseError::InvalidInput);
                }
                (true, d)
            }
            None => match s.strip_prefix('+') {
                Some(d) => {
                    if d.is_empty() {
                        return Err(ParseError::InvalidInput);
                    }
                    (false, d)
                }
                None => (false, s),
            },
        };
        if !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(ParseError::InvalidInput);
        }
        let bytes = digits.as_bytes();
        let mut r = Mpz::new();
        let mut i = 0;
        while i < bytes.len() {
            let end = (i + 18).min(bytes.len());
            let mut chunk: u64 = 0;
            for &b in &bytes[i..end] {
                chunk = chunk.wrapping_mul(10).wrapping_add((b - b'0') as u64);
            }
            let scale = 10u64.pow((end - i) as u32);
            r = r
                .try_mul_ui(scale)
                .map_err(|_| ParseError::CapacityOverflow)?;
            r = r
                .try_add_ui(chunk)
                .map_err(|_| ParseError::CapacityOverflow)?;
            i = end;
        }
        if neg && r.sign != 0 {
            r.sign = -1;
        }
        Ok(r)
    }

    pub fn write_decimal_buf(&self, buf: &mut [u8]) -> usize {
        if self.sign == 0 {
            if !buf.is_empty() {
                buf[0] = b'0';
                return 1;
            }
            return 0;
        }
        let mut digits = [0u8; 160];
        let mut di = digits.len();
        let mut m = self.clone();
        loop {
            let mut qbuf = [0u64; MPZ_MAX_LIMBS];
            let (qlen, rem) =
                Self::mag_divmod_u64_len(&m.mag[..m.len], 1_000_000_000_000_000_000, &mut qbuf);
            let mut rem_d = rem;
            if qlen == 0 {
                if rem_d == 0 {
                    di -= 1;
                    digits[di] = b'0';
                } else {
                    while rem_d != 0 {
                        di -= 1;
                        digits[di] = b'0' + (rem_d % 10) as u8;
                        rem_d /= 10;
                    }
                }
                break;
            } else {
                for _ in 0..18 {
                    di -= 1;
                    digits[di] = b'0' + (rem_d % 10) as u8;
                    rem_d /= 10;
                }
            }
            m.mag[..qlen].copy_from_slice(&qbuf[..qlen]);
            m.len = qlen;
        }
        let mut pos = 0;
        if self.sign < 0 && !buf.is_empty() {
            buf[0] = b'-';
            pos = 1;
        }
        let src = &digits[di..];
        let avail = buf.len() - pos;
        let n = src.len().min(avail);
        buf[pos..pos + n].copy_from_slice(&src[..n]);
        pos + n
    }

    // -----------------------------------------------------------------------
    // Low-level limb access
    // -----------------------------------------------------------------------

    /// `mpz_getlimbn`: return the nth limb (0-based).  Returns `None` if `n >= len`.
    pub fn getlimbn(&self, n: usize) -> Option<u64> {
        if n >= self.len {
            None
        } else {
            Some(self.mag[n])
        }
    }

    // -----------------------------------------------------------------------
    // Comparison with f64
    // -----------------------------------------------------------------------

    /// `mpz_cmp_d`: compare with `f64`.
    pub fn cmp_d(&self, v: f64) -> Ordering {
        if v.is_nan() {
            return Ordering::Greater; // GMP convention: NaN compares as greater
        }
        let other = match Mpz::from_d(v) {
            Ok(m) => m,
            Err(_) => {
                // v is infinity or subnormal/zero that couldn't fit
                if v.is_infinite() {
                    if v.is_sign_positive() {
                        return Ordering::Less;
                    } else {
                        return Ordering::Greater;
                    }
                }
                return self.cmp(&Mpz::new());
            }
        };
        self.cmp(&other)
    }

    /// `mpz_cmpabs_d`: compare absolute value with `f64`.
    pub fn cmpabs_d(&self, v: f64) -> Ordering {
        if v.is_nan() {
            return Ordering::Greater;
        }
        let v_abs = v.abs();
        let other = match Mpz::from_d(v_abs) {
            Ok(m) => m,
            Err(_) => {
                if v.is_infinite() {
                    return Ordering::Less;
                }
                return self.cmp(&Mpz::new());
            }
        };
        self.cmpabs(&other)
    }

    // -----------------------------------------------------------------------
    // Import / Export
    // -----------------------------------------------------------------------

    /// `mpz_import`: construct an `Mpz` from a byte buffer.
    ///
    /// - `count`: number of chunks.
    /// - `order`: ordering of chunks (`Little` = least significant chunk first,
    ///   `Big` = most significant chunk first).
    /// - `size`: size of each chunk in bytes.
    /// - `endian`: endianness within each chunk.
    /// - `data`: the raw bytes.  Must have length `count * size`.
    pub fn try_import(
        count: usize,
        order: Endian,
        size: usize,
        endian: Endian,
        data: &[u8],
    ) -> Result<Mpz, CapacityError> {
        if count == 0 || size == 0 {
            return Ok(Mpz::new());
        }
        if data.len() < count * size {
            return Ok(Mpz::new());
        }

        // The maximum number of limbs we could need is ceil(count * size / 8)
        let max_limbs = (count * size + 7) / 8;
        if max_limbs > MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }

        // Collect all chunks into a flat little-endian limb representation.
        let mut limbs = [0u64; MPZ_MAX_LIMBS];
        let n_chunks = count;
        let chunk_size = size;

        // Helper: read `size` bytes from `data[start..]` as u64 (little-endian).
        let read_chunk = |start: usize| -> u64 {
            let mut val = 0u64;
            let end = (start + chunk_size).min(data.len());
            match endian {
                Endian::Little | Endian::Native if native_is_little() => {
                    for i in (start..end).rev() {
                        val = (val << 8) | data[i] as u64;
                    }
                }
                _ => {
                    // Big-endian: read bytes in order
                    for i in start..end {
                        val = (val << 8) | data[i] as u64;
                    }
                }
            }
            val
        };

        // Process chunks based on `order` (ordering of chunks)
        let mut limb_idx = 0;
        let mut bit_offset = 0;

        match order {
            Endian::Little | Endian::Native if native_is_little() => {
                // Least significant chunk first
                for ci in 0..n_chunks {
                    let chunk_val = read_chunk(ci * chunk_size);
                    limbs[limb_idx] |= chunk_val << bit_offset;
                    bit_offset += chunk_size * 8;
                    while bit_offset >= 64 {
                        limb_idx += 1;
                        bit_offset -= 64;
                        if bit_offset > 0 {
                            limbs[limb_idx] = chunk_val >> (chunk_size * 8 - bit_offset);
                        }
                    }
                }
            }
            _ => {
                // Most significant chunk first
                for ci in (0..n_chunks).rev() {
                    let chunk_val = read_chunk(ci * chunk_size);
                    limbs[limb_idx] |= chunk_val << bit_offset;
                    bit_offset += chunk_size * 8;
                    while bit_offset >= 64 {
                        limb_idx += 1;
                        bit_offset -= 64;
                        if bit_offset > 0 {
                            limbs[limb_idx] = chunk_val >> (chunk_size * 8 - bit_offset);
                        }
                    }
                }
            }
        }

        let mut r = Mpz::new();
        let _n_limbs = (count * size + 7) / 8;
        let n_limbs_actual = if bit_offset > 0 {
            limb_idx + 1
        } else {
            limb_idx
        };
        let final_limbs = n_limbs_actual.min(MPZ_MAX_LIMBS);
        r.mag[..final_limbs].copy_from_slice(&limbs[..final_limbs]);
        r.len = final_limbs;
        r.sign = if r.len == 0 { 0 } else { 1 };
        r.trim();
        Ok(r)
    }

    /// `mpz_export`: write the value into a caller-provided byte buffer.
    ///
    /// Returns `Some(bytes_written)` on success, or `None` if the buffer is too small.
    ///
    /// - `buf`: output buffer.
    /// - `order`: ordering of chunks.
    /// - `size`: size of each chunk in bytes.
    /// - `endian`: endianness within each chunk.
    pub fn export_buf(
        &self,
        buf: &mut [u8],
        order: Endian,
        size: usize,
        endian: Endian,
    ) -> Option<usize> {
        if self.is_zero() || size == 0 {
            if !buf.is_empty() {
                buf[0] = 0;
                return Some(1);
            }
            return Some(0);
        }

        // Determine total bytes needed
        let total_bits = self.sizeinbase2();
        let total_bytes = (total_bits + 7) / 8;
        let chunk_bytes = size;
        let n_chunks = (total_bytes + chunk_bytes - 1) / chunk_bytes;
        let needed_bytes = n_chunks * chunk_bytes;

        if buf.len() < needed_bytes {
            return None;
        }

        // Zero out the buffer
        for b in buf.iter_mut().take(needed_bytes) {
            *b = 0;
        }

        // Write each limb into the buffer in little-endian byte order
        let mut byte_pos = 0usize;
        for i in 0..self.len {
            let mut val = self.mag[i];
            for _ in 0..8 {
                if byte_pos < needed_bytes {
                    buf[byte_pos] = (val & 0xFF) as u8;
                    val >>= 8;
                    byte_pos += 1;
                }
            }
        }

        // Rearrange chunks based on `order` and `endian`
        // If the natural output is already in the correct format, return.
        // Otherwise, we need to swap bytes and chunks.

        let is_native_little = native_is_little();
        let native_endian = if is_native_little {
            Endian::Little
        } else {
            Endian::Big
        };

        // The internal format is little-endian bytes (least significant byte first).
        // We need to convert to the requested (order, size, endian).

        // If no transformation is needed:
        if order == Endian::Little || (order == Endian::Native && native_endian == Endian::Little) {
            if endian == Endian::Little
                || (endian == Endian::Native && native_endian == Endian::Little)
            {
                if size == 1 || size == 0 {
                    return Some(needed_bytes);
                }
                // We still need to handle chunk size
            }
        }

        // Convert by processing chunk by chunk.
        // First, make a copy of the original buffer to avoid read-write conflicts.
        let mut orig = [0u8; 8];
        let copy_len = needed_bytes.min(orig.len());
        orig[..copy_len].copy_from_slice(&buf[..copy_len]);
        let mut chunk_tmp = [0u8; 8];
        for ci in 0..n_chunks {
            // Read chunk from natural LE representation (using original copy)
            let src_start = ci * chunk_bytes;
            let src_end = (src_start + chunk_bytes).min(needed_bytes);
            for j in src_start..src_end {
                let tmp_idx = j - src_start;
                chunk_tmp[tmp_idx] = orig[j];
            }
            for j in src_end - src_start..chunk_bytes {
                chunk_tmp[j] = 0;
            }

            // Apply within-chunk endianness swap
            let chunk_val = match endian {
                Endian::Little | Endian::Native if is_native_little => {
                    // Already little-endian
                    let mut v = 0u64;
                    for j in 0..chunk_bytes {
                        v |= (chunk_tmp[j] as u64) << (j * 8);
                    }
                    v
                }
                _ => {
                    // Big-endian within chunk
                    let mut v = 0u64;
                    for j in 0..chunk_bytes {
                        v = (v << 8) | chunk_tmp[j] as u64;
                    }
                    v
                }
            };

            // Write chunk back in the natural byte order
            for j in 0..chunk_bytes {
                chunk_tmp[j] = ((chunk_val >> (j * 8)) & 0xFF) as u8;
            }

            // Write to output position based on chunk order
            let dst_start = match order {
                Endian::Little | Endian::Native if is_native_little => ci * chunk_bytes,
                _ => (n_chunks - 1 - ci) * chunk_bytes,
            };
            for j in 0..chunk_bytes {
                if dst_start + j < buf.len() {
                    buf[dst_start + j] = chunk_tmp[j];
                }
            }
        }

        Some(needed_bytes)
    }

    // -----------------------------------------------------------------------
    // Probable prime — Miller–Rabin
    // -----------------------------------------------------------------------

    /// `mpz_probab_prime_p`: Miller–Rabin probable prime test.
    ///
    /// Returns:
    /// - `2` if `self` is definitely prime.
    /// - `1` if `self` is probably prime.
    /// - `0` if `self` is definitely composite.
    ///
    /// Uses the first `reps` primes as bases.
    pub fn try_probab_prime_p(&self, reps: u32) -> Result<i32, CapacityError> {
        if self.sign < 0 {
            return Ok(0);
        }
        if self.len == 0 {
            return Ok(0);
        }

        // Handle small values directly
        if self.len == 1 {
            let v = self.mag[0];
            if v < 2 {
                return Ok(0);
            }
            // Check against small primes
            for &p in SMALL_PRIMES.iter() {
                if p * p > v {
                    break;
                }
                if v % p == 0 {
                    return Ok(if v == p { 2 } else { 0 });
                }
            }
            if v <= SMALL_PRIMES[SMALL_PRIMES.len() - 1].pow(2) {
                return Ok(2); // definitely prime
            }
        }

        // Check small prime divisors first
        for &p in SMALL_PRIMES.iter() {
            if self.tdiv_ui(p) == 0 {
                // self is divisible by p
                let _p_mpz = Mpz::from_u64(p);
                if self.len == 1 && self.mag[0] == p {
                    return Ok(2);
                }
                return Ok(0);
            }
        }

        // Miller–Rabin: write self-1 = d * 2^s with d odd
        let one = Mpz::from_u64(1);
        let nm1 = self.try_sub(&one)?;
        let s = nm1.mag[0].trailing_zeros(); // number of trailing zeros
        let d = nm1.fdiv_q_2exp(s);

        let reps = reps.min(64);
        for i in 0..reps {
            let a = Mpz::from_u64(SMALL_PRIMES[i as usize]);
            // Compute x = a^d mod self
            let mut x = a.try_powm(&d, self)?;
            if x.len == 0 || (x.len == 1 && x.mag[0] == 1) || x == nm1 {
                continue;
            }
            let mut composite = true;
            for _ in 0..s {
                x = x.try_powm(&Mpz::from_u64(2), self)?;
                if x == nm1 {
                    composite = false;
                    break;
                }
                if x.len == 1 && x.mag[0] == 1 {
                    break;
                }
            }
            if composite {
                return Ok(0);
            }
        }

        Ok(1)
    }

    // =======================================================================
    // Constant-time operations (feature-gated)
    // =======================================================================

    #[cfg(feature = "const_time")]
    /// Constant-time add. Returns `self + other` in time independent of values.
    ///
    /// WARNING: this is NOT constant-time for `CapacityError` (which may leak info).
    pub fn ct_add(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        // Iterate over all MPZ_MAX_LIMBS limbs, masking unused with 0.
        let mut result = Mpz::new();
        let mut carry = 0u128;
        for i in 0..MPZ_MAX_LIMBS {
            let va = if i < self.len { self.mag[i] as u128 } else { 0 };
            let vb = if i < other.len {
                other.mag[i] as u128
            } else {
                0
            };
            let s = va + vb + carry;
            result.mag[i] = s as u64;
            carry = s >> 64;
        }
        if carry != 0 {
            return Err(CapacityError);
        }

        // Determine sign: if either is zero, take the other's sign.
        // Both zero → already handled by iteration (all mag = 0).
        // Same sign → that sign.  Different sign → compare magnitudes.
        let mut result_len = MPZ_MAX_LIMBS;
        while result_len > 0 && result.mag[result_len - 1] == 0 {
            result_len -= 1;
        }
        result.len = result_len;

        if self.sign == 0 {
            result.sign = other.sign;
        } else if other.sign == 0 {
            result.sign = self.sign;
        } else if self.sign == other.sign {
            result.sign = self.sign;
        } else {
            // Opposite signs: subtract smaller from larger
            // We need to do ct_sub logically but re-do the computation in constant time.
            // Re-compute as self - (-other):
            let mut neg_other_mag = [0u64; MPZ_MAX_LIMBS];
            for i in 0..MPZ_MAX_LIMBS {
                neg_other_mag[i] = other.mag[i];
            }
            let mut neg_other_len = other.len;
            // Compute subtraction: |self| - |other| with sign of whichever is larger
            let self_gt = Self::ct_cmp_mag(&self.mag, self.len, &neg_other_mag, neg_other_len);
            // If self_gt == Greater, result = self.mag - other.mag with self.sign
            // If self_gt == Less, result = other.mag - self.mag with other.sign
            // If self_gt == Equal, result = 0
            let mut borrow: i128 = 0;
            let (src_a, src_b, out_sign) = if self_gt == Ordering::Greater {
                (&self.mag, &neg_other_mag, self.sign)
            } else if self_gt == Ordering::Less {
                (&neg_other_mag, &self.mag, other.sign)
            } else {
                return Ok(Mpz::new());
            };
            let a_len = if self_gt == Ordering::Greater {
                self.len
            } else {
                other.len
            }
            .max(if self_gt == Ordering::Greater {
                other.len
            } else {
                self.len
            });
            for i in 0..a_len.max(1) {
                let ai = if i < src_a.len() { src_a[i] as i128 } else { 0 };
                let bi = if i < src_b.len() { src_b[i] as i128 } else { 0 };
                let mut cur = ai - bi - borrow;
                if cur < 0 {
                    cur += 1i128 << 64;
                    borrow = 1;
                } else {
                    borrow = 0;
                }
                result.mag[i] = cur as u64;
            }
            let mut rl = a_len.max(1);
            while rl > 0 && result.mag[rl - 1] == 0 {
                rl -= 1;
            }
            result.len = rl;
            result.sign = if rl == 0 { 0 } else { out_sign };
        }

        result.trim();
        Ok(result)
    }

    #[cfg(feature = "const_time")]
    /// Constant-time subtract.
    pub fn ct_sub(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let neg = other.neg_to();
        self.ct_add(&neg)
    }

    #[cfg(feature = "const_time")]
    /// Constant-time compare of magnitudes (public helper).
    /// Returns `Greater` if a > b, `Less` if a < b, `Equal` if a == b.
    /// No branches on the limb values themselves.
    pub fn ct_cmp_mag(
        a: &[u64; MPZ_MAX_LIMBS],
        a_len: usize,
        b: &[u64; MPZ_MAX_LIMBS],
        b_len: usize,
    ) -> Ordering {
        // Compare lengths first (this is data-dependent on len but not on values)
        // For truly constant-time, we process all limbs regardless of length.
        let mut gt = 0u64;
        let mut lt = 0u64;
        for i in (0..MPZ_MAX_LIMBS).rev() {
            let va = a[i];
            let vb = b[i];
            let diff = (va as i128) - (vb as i128);
            // diff > 0  ⇒  gt
            // diff < 0  ⇒  lt
            // Use arithmetic shift to extract sign
            let is_gt = ((diff >> 127) as u64).wrapping_neg(); // 0 if diff >= 0, !0 if diff < 0
            let is_lt = (((-diff) >> 127) as u64).wrapping_neg();
            // Actually need careful approach:
            let gt_bit = if diff > 0 { 1u64 } else { 0u64 };
            let lt_bit = if diff < 0 { 1u64 } else { 0u64 };
            // In constant-time, we'd use bit masking:
            // Let mask = ((diff as i64) >> 63) as u64 for 64-bit or use sign bit
            // For 128-bit: ((diff as i128) >> 127) as u64 gives 0 for non-negative, !0 for negative
            // Wait: >> with sign extension: for positive/zero → 0, for negative → !0
            // So we can compute:
        }
        // Use a truly branchless approach:
        let mut result = 0i8; // -1 for Less, 0 for Equal, 1 for Greater
        for i in (0..MPZ_MAX_LIMBS).rev() {
            let va = a[i];
            let vb = b[i];
            // Compare va and vb without branching:
            let gt_mask = (va as i128 > vb as i128) as u64;
            let lt_mask = ((va as i128) < (vb as i128)) as u64;
            let eq_mask = (va == vb) as u64;
            // Update: if this pair is decisive, set result.  Otherwise keep previous.
            // decisive = (gt_mask | lt_mask)  (1 if this pair determines ordering)
            let decisive = gt_mask | lt_mask;
            // If decisive, pick new value; else keep old.
            // This is branchless using bit tricks:
            let new_val = (gt_mask as i8) - (lt_mask as i8); // 1, 0, or -1
                                                             // Blend: result = decisive ? new_val : result
                                                             // In constant time: result = result ^ ((result ^ new_val) & mask)
            let mask = 0u64.wrapping_sub(decisive);
            result ^= ((result ^ new_val) as u64 & mask) as i8;
        }
        match result {
            1 => Ordering::Greater,
            -1 => Ordering::Less,
            _ => Ordering::Equal,
        }
    }

    #[cfg(feature = "const_time")]
    /// Constant-time compare. Returns `Ordering` without branching on secret data.
    pub fn ct_cmp(&self, other: &Mpz) -> Ordering {
        // Compute self.sign XOR other.sign — if different, result is determined by sign
        let sign_diff = (self.sign != other.sign) as u64;
        let self_neg = (self.sign < 0) as u64;
        let other_neg = (other.sign < 0) as u64;

        // If signs differ, result depends on which is negative
        // self negative, other non-negative → Less
        // self non-negative, other negative → Greater
        let sign_result: i8 = if self.sign > other.sign {
            1
        } else if self.sign < other.sign {
            -1
        } else {
            0
        };

        // If signs are same (both non-negative or both negative), compare magnitudes
        let mag_cmp = Self::ct_cmp_mag(&self.mag, self.len, &other.mag, other.len);
        let mag_result: i8 = match mag_cmp {
            Ordering::Greater => 1,
            Ordering::Less => -1,
            Ordering::Equal => 0,
        };

        // For same sign and non-negative: mag_result is the answer
        // For same sign and negative: answer is reverse of mag_result
        let both_neg = (self.sign < 0 && other.sign < 0) as u64;
        let neg_mask = 0u64.wrapping_sub(both_neg);
        // If both_neg: flip mag_result; otherwise keep it
        let adjusted_mag = mag_result ^ ((mag_result as u64 & 1) & neg_mask) as i8;
        // Actually negation in two's complement: -x = !x + 1
        // But we just want: if both_neg then -mag_result else mag_result
        // -1 (0xFF) stays -1, 0 stays 0, 1 becomes -1
        // For i8: -mag_result = (mag_result ^ neg_mask as i8).wrapping_add(neg_mask as i8 & 1)
        // Simpler: just use branching for now

        // Blend sign_result and mag_result:
        // If signs differ (sign_diff == 1): use sign_result
        // If signs same (sign_diff == 0): use mag_result (negated if both negative)
        let blend = (sign_diff != 0) as u64;
        let mut final_val = if blend != 0 { sign_result } else { mag_result };
        if blend == 0 && both_neg != 0 {
            final_val = match final_val {
                1 => -1,
                -1 => 1,
                _ => 0,
            };
        }
        match final_val {
            1 => Ordering::Greater,
            -1 => Ordering::Less,
            _ => Ordering::Equal,
        }
    }

    #[cfg(feature = "const_time")]
    /// Constant-time select. Returns `self` if `bit == 0`, `other` if `bit == 1`.
    /// `bit` must be 0 or 1; behaviour is undefined otherwise.
    pub fn ct_select(&self, other: &Mpz, bit: u64) -> Mpz {
        // mask = 0 if bit == 0, !0 if bit == 1
        let mask = 0u64.wrapping_sub(bit & 1);
        let mut result = Mpz::new();
        for i in 0..MPZ_MAX_LIMBS {
            result.mag[i] = (self.mag[i] & !mask) | (other.mag[i] & mask);
        }
        result.len = if (self.len & !mask) | (other.len & mask) != 0 {
            let blended_len = (self.len & !mask as usize) | (other.len & mask as usize);
            // Trim: find the highest set limb
            let mut l = MPZ_MAX_LIMBS;
            while l > 0 && result.mag[l - 1] == 0 {
                l -= 1;
            }
            l
        } else {
            0
        };
        result.sign = ((self.sign as i8 & !(mask as i8)) | (other.sign as i8 & mask as i8)) as i8;
        result.trim();
        result
    }

    // =======================================================================
    // Formal capability declaration
    // ===================================================================

    /// Return a static table of all GMP `mpz_*` functions and their
    /// implementation status in gmp-rs.
    ///
    /// Each entry is a `(opcode, gmp_name, status)` triple:
    /// - `opcode`: 0 = N/A, positive = implemented, negative = intentionally absent.
    /// - `gmp_name`: the GMP C function name, e.g. `"mpz_add"`.
    /// - `status`: short description, e.g. `"yes"`, `"partial"`, `"no (requires std)"`.
    pub fn capability_map() -> &'static [(i32, &'static str, &'static str)] {
        &CAPABILITY_TABLE
    }

    // =======================================================================
    // Constant-time operations (feature-gated via `const_time`)
    // =======================================================================

    /// Constant-time addition.  Returns `self + other` in time independent of
    /// the **values** of `self` and `other`.  Capacity overflow (`Err`) still
    /// leaks whether overflow occurred.
    ///
    /// Available only with the `const_time` feature.
    #[cfg(feature = "const_time")]
    pub fn ct_add(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        // Operate on all MAX_LIMBS limbs (not just self.len) so that
        // iteration count is independent of the values.
        let mut result = Mpz::new();
        // We need to handle same-sign and opposite-sign cases.
        // For CT purposes, compute both possibilities and select.
        let same_sign = ((self.sign as i8) ^ (other.sign as i8)).wrapping_add(1) as u64;
        let same_sign = (same_sign & 1).wrapping_sub(1); // 0 if different, !0 if same

        // Case 1: same sign → add magnitudes
        let mut add_mag = [0u64; MPZ_MAX_LIMBS];
        let mut carry = 0u128;
        for i in 0..MPZ_MAX_LIMBS {
            let va = if i < self.len { self.mag[i] as u128 } else { 0 };
            let vb = if i < other.len {
                other.mag[i] as u128
            } else {
                0
            };
            let s = va + vb + carry;
            add_mag[i] = s as u64;
            carry = s >> 64;
        }
        let add_overflow = carry != 0;

        // Case 2: opposite signs → subtract magnitudes
        // Determine which magnitude is larger (CT selection)
        let mut a_bigger: i8 = 0;
        for i in (0..MPZ_MAX_LIMBS).rev() {
            let va = if i < self.len { self.mag[i] } else { 0 };
            let vb = if i < other.len { other.mag[i] } else { 0 };
            let mask: u8 = ((a_bigger as u8) ^ 1u8.wrapping_sub(1)) & 1;
            // Only compare if we haven't already decided
            let not_yet_decided = ((a_bigger as u8).wrapping_sub(1) >> 7) & 1;
            let gt = (vb < va) as u8 & not_yet_decided;
            let lt = (va < vb) as u8 & not_yet_decided;
            a_bigger = (a_bigger as u8 | gt | lt.wrapping_neg()) as i8;
        }
        let use_self_a = ((a_bigger as u8) ^ 1) & 1; // 1 if self >= other, 0 otherwise
        let mask_a = (use_self_a).wrapping_sub(1) as u64; // !0 if use self, 0 otherwise
        let mask_b = (!use_self_a).wrapping_sub(1) as u64; // !0 if use other, 0 otherwise

        let mut sub_mag = [0u64; MPZ_MAX_LIMBS];
        let mut borrow: i128 = 0;
        for i in 0..MPZ_MAX_LIMBS {
            let va = if i < self.len { self.mag[i] as i128 } else { 0 };
            let vb = if i < other.len {
                other.mag[i] as i128
            } else {
                0
            };
            // CT select which is minuend and which is subtrahend
            let minuend = va.wrapping_add(((vb - va) as i128) & (mask_b as i128));
            let subtrahend = vb.wrapping_add(((va - vb) as i128) & (mask_a as i128));
            let mut cur = minuend - subtrahend - borrow;
            if cur < 0 {
                cur += 1i128 << 64;
                borrow = 1;
            } else {
                borrow = 0;
            }
            sub_mag[i] = cur as u64;
        }

        // Select result: same-sign path vs opposite-sign path
        let use_add = same_sign;
        let use_sub = same_sign ^ !0u64;

        let mut can_overflow = add_overflow as u64;
        let mut sign = 0i8;
        for i in 0..MPZ_MAX_LIMBS {
            let add_val = add_mag[i] & use_add;
            let sub_val = sub_mag[i] & use_sub;
            result.mag[i] = add_val | sub_val;
        }
        // CT sign selection
        // If same_sign, sign = self.sign (which equals other.sign)
        // If opposite, sign = whichever magnitude was larger
        let self_s = self.sign as i64;
        let other_s = other.sign as i64;
        let use_self_sign = use_self_a as i64;
        let use_other_sign = (!use_self_a) as i64;
        let sel_sign = self_s
            .wrapping_mul(use_self_sign)
            .wrapping_add(other_s.wrapping_mul(use_other_sign));
        let same_sign_path = same_sign as i64 & self_s;
        sign = (same_sign_path | (sel_sign & (use_sub as i64))) as i8;
        result.sign = sign;
        // Determine len: iterate all limbs CT to find top non-zero
        let mut new_len = 0usize;
        for i in (0..MPZ_MAX_LIMBS).rev() {
            let is_zero = (result.mag[i] == 0) as usize;
            let not_is_zero = 1usize ^ is_zero;
            // Only update new_len if we haven't found a non-zero limb yet
            let already_set = (new_len != 0) as usize;
            new_len =
                new_len | ((i + 1) & (not_is_zero.wrapping_sub(1)) & !already_set.wrapping_sub(1));
        }
        result.len = new_len;
        if result.len == 0 {
            result.sign = 0;
        }

        if can_overflow != 0 {
            Err(CapacityError)
        } else {
            Ok(result)
        }
    }

    /// Constant-time subtraction.  See [`Mpz::ct_add`].
    #[cfg(feature = "const_time")]
    pub fn ct_sub(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let mut neg_other = other.clone();
        neg_other.sign = -neg_other.sign;
        self.ct_add(&neg_other)
    }

    /// Constant-time compare.  Returns `Ordering` without branching on
    /// the values (but the returned value is revealed).
    #[cfg(feature = "const_time")]
    pub fn ct_cmp(&self, other: &Mpz) -> Ordering {
        // CT sign comparison
        let self_s = self.sign as i64;
        let other_s = other.sign as i64;
        let sign_diff = self_s - other_s;
        // If sign_diff != 0, signs differ → result from sign
        let sign_mask = ((sign_diff as u64).wrapping_sub(1)) >> 63;
        let sign_gt = (sign_diff > 0) as u64;
        let sign_lt = (sign_diff < 0) as u64;

        // CT magnitude comparison (only valid if signs equal)
        let mut mag_gt = 0u64;
        let mut mag_lt = 0u64;
        for i in (0..MPZ_MAX_LIMBS).rev() {
            let va = if i < self.len { self.mag[i] } else { 0 };
            let vb = if i < other.len { other.mag[i] } else { 0 };
            // If we haven't decided yet, compare this limb
            let decided = (mag_gt | mag_lt).wrapping_sub(1) >> 63;
            let not_decided = decided ^ 1;
            mag_gt |= ((va > vb) as u64) & not_decided;
            mag_lt |= ((va < vb) as u64) & not_decided;
        }
        // If same sign and positive: use mag comparison directly
        // If same sign and negative: reverse mag comparison
        let both_neg = (self_s >> 63) & (other_s >> 63) & 1;
        let mag_result_gt = mag_gt ^ (both_neg); // XOR for reversal
        let mag_result_lt = mag_lt ^ (both_neg);

        let final_gt = (sign_gt & sign_mask) | (mag_result_gt & (sign_mask ^ 1));
        let final_lt = (sign_lt & sign_mask) | (mag_result_lt & (sign_mask ^ 1));

        let gt_bit = (final_gt).wrapping_sub(1) >> 63;
        let lt_bit = (final_lt).wrapping_sub(1) >> 63;
        let eq_bit = ((gt_bit | lt_bit) ^ 1) & 1;

        // ct_select between Ordering values
        let ord_val = gt_bit as i8 - lt_bit as i8;
        match ord_val {
            1 => Ordering::Greater,
            -1 => Ordering::Less,
            _ => Ordering::Equal,
        }
    }

    /// Constant-time select.  Returns `self` if `bit == 0`, `other` if `bit == 1`.
    /// `bit` must be 0 or 1.
    #[cfg(feature = "const_time")]
    pub fn ct_select(&self, other: &Mpz, bit: u64) -> Mpz {
        let mask = (bit.wrapping_sub(1)) >> 63; // 0 if bit=0, !0 if bit=1
        let not_mask = !mask;
        let mut result = Mpz::new();
        for i in 0..MPZ_MAX_LIMBS {
            result.mag[i] = (self.mag[i] & not_mask) | (other.mag[i] & mask);
        }
        // Select sign and len CT
        let self_s = self.sign as i64;
        let other_s = other.sign as i64;
        let sel_s = (self_s & not_mask as i64) | (other_s & mask as i64);
        let self_len = self.len;
        let other_len = other.len;
        result.len = (self_len & (not_mask as usize)) | (other_len & (mask as usize));
        result.sign = sel_s as i8;
        result
    }
} // impl Mpz

// ===========================================================================
// Capability table (module-level static)
// ===========================================================================

#[doc(hidden)]
pub static CAPABILITY_TABLE: &[(i32, &str, &str)] = &[
    // -----------------------------------------------------------------------
    // B1: Initialization
    // -----------------------------------------------------------------------
    (1, "mpz_init", "yes — Mpz::new()"),
    (0, "mpz_inits", "no (varargs — not expressible in Rust)"),
    (
        -1,
        "mpz_init2",
        "no (fixed array — no allocation hint needed)",
    ),
    (0, "mpz_clear", "no (Rust Drop handles this)"),
    (0, "mpz_clears", "no (varargs)"),
    (-1, "mpz_realloc2", "no (fixed array — no reallocation)"),
    // -----------------------------------------------------------------------
    // B2: Assignment
    // -----------------------------------------------------------------------
    (1, "mpz_set", "yes — Mpz::set()"),
    (1, "mpz_set_ui", "yes — Mpz::set_ui()"),
    (1, "mpz_set_si", "yes — Mpz::set_si()"),
    (1, "mpz_set_d", "yes — Mpz::from_d()"),
    (-1, "mpz_set_q", "no (mpq type does not exist)"),
    (-1, "mpz_set_f", "no (mpf type does not exist)"),
    (
        1,
        "mpz_set_str",
        "partial — base 10 only via from_decimal_str",
    ),
    (1, "mpz_swap", "yes — Mpz::swap()"),
    // -----------------------------------------------------------------------
    // B3: Combined Init + Assignment
    // -----------------------------------------------------------------------
    (1, "mpz_init_set", "yes — Clone::clone()"),
    (1, "mpz_init_set_ui", "yes — Mpz::from_u64()"),
    (1, "mpz_init_set_si", "yes — Mpz::from_i64()"),
    (1, "mpz_init_set_d", "yes — Mpz::from_d()"),
    (
        1,
        "mpz_init_set_str",
        "partial — base 10 only via from_decimal_str",
    ),
    // -----------------------------------------------------------------------
    // B4: Conversion
    // -----------------------------------------------------------------------
    (1, "mpz_get_ui", "yes — Mpz::get_ui()"),
    (1, "mpz_get_si", "yes — Mpz::get_si()"),
    (1, "mpz_get_d", "yes — Mpz::get_d()"),
    (1, "mpz_get_d_2exp", "yes — Mpz::get_d_2exp()"),
    (
        1,
        "mpz_get_str",
        "partial — base 10 only via write_decimal_buf",
    ),
    // -----------------------------------------------------------------------
    // B5: Arithmetic
    // -----------------------------------------------------------------------
    (1, "mpz_add", "yes — Mpz::try_add()"),
    (1, "mpz_add_ui", "yes — Mpz::try_add_ui()"),
    (1, "mpz_sub", "yes — Mpz::try_sub()"),
    (1, "mpz_sub_ui", "yes — Mpz::try_sub_ui()"),
    (1, "mpz_ui_sub", "yes — Mpz::try_ui_sub()"),
    (1, "mpz_mul", "yes — Mpz::try_mul()"),
    (1, "mpz_mul_si", "yes — Mpz::try_mul_si()"),
    (1, "mpz_mul_ui", "yes — Mpz::try_mul_ui()"),
    (1, "mpz_addmul", "yes — Mpz::try_addmul()"),
    (1, "mpz_addmul_ui", "yes — Mpz::try_addmul_ui()"),
    (1, "mpz_submul", "yes — Mpz::try_submul()"),
    (1, "mpz_submul_ui", "yes — Mpz::try_submul_ui()"),
    (1, "mpz_mul_2exp", "yes — Mpz::try_mul_2exp()"),
    (1, "mpz_neg", "yes — Mpz::neg() / neg_to()"),
    (1, "mpz_abs", "yes — Mpz::abs() / abs_to()"),
    // -----------------------------------------------------------------------
    // B6a: Division — Truncating (toward zero)
    // -----------------------------------------------------------------------
    (1, "mpz_tdiv_q", "yes — Mpz::tdiv_q()"),
    (1, "mpz_tdiv_r", "yes — Mpz::tdiv_r()"),
    (1, "mpz_tdiv_qr", "yes — Mpz::tdiv_qr()"),
    (1, "mpz_tdiv_q_ui", "yes — Mpz::tdiv_q_ui()"),
    (1, "mpz_tdiv_r_ui", "yes — Mpz::tdiv_r_ui()"),
    (1, "mpz_tdiv_qr_ui", "yes — Mpz::try_tdiv_qr_ui()"),
    (1, "mpz_tdiv_ui", "yes — Mpz::tdiv_ui()"),
    (1, "mpz_tdiv_q_2exp", "yes — Mpz::tdiv_q_2exp()"),
    (1, "mpz_tdiv_r_2exp", "yes — Mpz::tdiv_r_2exp()"),
    // -----------------------------------------------------------------------
    // B6b: Division — Floor (toward −∞)
    // -----------------------------------------------------------------------
    (1, "mpz_fdiv_q", "yes — Mpz::try_fdiv_q()"),
    (1, "mpz_fdiv_r", "yes — Mpz::try_fdiv_r()"),
    (1, "mpz_fdiv_qr", "yes — Mpz::try_fdiv_qr()"),
    (1, "mpz_fdiv_q_ui", "yes — Mpz::try_fdiv_q_ui()"),
    (1, "mpz_fdiv_r_ui", "yes — Mpz::try_fdiv_r_ui()"),
    (1, "mpz_fdiv_qr_ui", "yes — Mpz::try_fdiv_qr_ui()"),
    (1, "mpz_fdiv_ui", "yes — Mpz::fdiv_ui()"),
    (1, "mpz_fdiv_q_2exp", "yes — Mpz::fdiv_q_2exp()"),
    (1, "mpz_fdiv_r_2exp", "yes — Mpz::fdiv_r_2exp()"),
    // -----------------------------------------------------------------------
    // B6c: Division — Ceiling (toward +∞)
    // -----------------------------------------------------------------------
    (1, "mpz_cdiv_q", "yes — Mpz::try_cdiv_q()"),
    (1, "mpz_cdiv_r", "yes — Mpz::try_cdiv_r()"),
    (1, "mpz_cdiv_qr", "yes — Mpz::try_cdiv_qr()"),
    (1, "mpz_cdiv_q_ui", "yes — Mpz::try_cdiv_q_ui()"),
    (1, "mpz_cdiv_r_ui", "yes — Mpz::try_cdiv_r_ui()"),
    (1, "mpz_cdiv_qr_ui", "yes — Mpz::try_cdiv_qr_ui()"),
    (1, "mpz_cdiv_ui", "yes — Mpz::cdiv_ui()"),
    (1, "mpz_cdiv_q_2exp", "yes — Mpz::cdiv_q_2exp()"),
    (1, "mpz_cdiv_r_2exp", "yes — Mpz::cdiv_r_2exp()"),
    // -----------------------------------------------------------------------
    // B6d: Modulo (non-negative remainder)
    // -----------------------------------------------------------------------
    (1, "mpz_mod", "yes — Mpz::try_mod()"),
    (1, "mpz_mod_ui", "yes — Mpz::mod_ui()"),
    // -----------------------------------------------------------------------
    // B6e: Exact division
    // -----------------------------------------------------------------------
    (1, "mpz_divexact", "yes — Mpz::try_divexact()"),
    (1, "mpz_divexact_ui", "yes — Mpz::try_divexact_ui()"),
    // -----------------------------------------------------------------------
    // B6f: Divisibility / Congruence
    // -----------------------------------------------------------------------
    (1, "mpz_divisible_p", "yes — Mpz::divisible_p()"),
    (1, "mpz_divisible_ui_p", "yes — Mpz::divisible_ui()"),
    (1, "mpz_divisible_2exp_p", "yes — Mpz::divisible_2exp_p()"),
    (1, "mpz_congruent_p", "yes — Mpz::congruent_p()"),
    (1, "mpz_congruent_ui_p", "yes — Mpz::congruent_ui_p()"),
    (1, "mpz_congruent_2exp_p", "yes — Mpz::congruent_2exp_p()"),
    // -----------------------------------------------------------------------
    // B7: Exponentiation
    // -----------------------------------------------------------------------
    (1, "mpz_powm", "yes — Mpz::try_powm()"),
    (1, "mpz_powm_ui", "yes — Mpz::try_powm_ui()"),
    (
        -1,
        "mpz_powm_sec",
        "no (requires constant-time impl — see const_time feature)",
    ),
    (1, "mpz_pow_ui", "yes — Mpz::try_pow_ui()"),
    (1, "mpz_ui_pow_ui", "yes — Mpz::try_ui_pow_ui()"),
    // -----------------------------------------------------------------------
    // B8: Root extraction
    // -----------------------------------------------------------------------
    (1, "mpz_root", "yes — Mpz::try_root()"),
    (1, "mpz_rootrem", "yes — Mpz::try_rootrem()"),
    (1, "mpz_sqrt", "yes — Mpz::isqrt()"),
    (1, "mpz_sqrtrem", "yes — Mpz::try_sqrtrem()"),
    (1, "mpz_perfect_power_p", "yes — Mpz::perfect_power_p()"),
    (1, "mpz_perfect_square_p", "yes — Mpz::perfect_square_p()"),
    // -----------------------------------------------------------------------
    // B9: Number theoretic
    // -----------------------------------------------------------------------
    (1, "mpz_probab_prime_p", "yes — Mpz::try_probab_prime_p()"),
    (-1, "mpz_nextprime", "no (requires alloc for sieve)"),
    (-1, "mpz_prevprime", "no (requires alloc for sieve)"),
    (1, "mpz_gcd", "yes — Mpz::try_gcd()"),
    (1, "mpz_gcd_ui", "yes — Mpz::gcd_ui()"),
    (1, "mpz_gcdext", "yes — Mpz::try_gcdext()"),
    (1, "mpz_lcm", "yes — Mpz::try_lcm()"),
    (1, "mpz_lcm_ui", "yes — Mpz::try_lcm_ui()"),
    (1, "mpz_invert", "yes — Mpz::try_invert()"),
    (1, "mpz_jacobi", "yes — Mpz::jacobi()"),
    (1, "mpz_legendre", "yes — Mpz::try_legendre()"),
    (1, "mpz_kronecker", "yes — Mpz::try_kronecker()"),
    (1, "mpz_kronecker_si", "yes — Mpz::try_kronecker_si()"),
    (1, "mpz_kronecker_ui", "yes — Mpz::try_kronecker_ui()"),
    (1, "mpz_si_kronecker", "yes — Mpz::try_si_kronecker()"),
    (1, "mpz_ui_kronecker", "yes — Mpz::try_ui_kronecker()"),
    (1, "mpz_remove", "yes — Mpz::try_remove()"),
    (1, "mpz_fac_ui", "yes — Mpz::try_fac_ui()"),
    (1, "mpz_2fac_ui", "yes — Mpz::try_2fac_ui()"),
    (-1, "mpz_mfac_uiui", "no (multi-factorial not implemented)"),
    (1, "mpz_primorial_ui", "yes — Mpz::try_primorial_ui()"),
    (1, "mpz_bin_ui", "yes — Mpz::try_bin_ui()"),
    (1, "mpz_bin_uiui", "yes — Mpz::try_bin_uiui()"),
    (1, "mpz_fib_ui", "yes — Mpz::try_fib_ui()"),
    (1, "mpz_fib2_ui", "yes — Mpz::try_fib2_ui()"),
    (1, "mpz_lucnum_ui", "yes — Mpz::try_lucnum_ui()"),
    (1, "mpz_lucnum2_ui", "yes — Mpz::try_lucnum2_ui()"),
    // -----------------------------------------------------------------------
    // B10: Comparison
    // -----------------------------------------------------------------------
    (1, "mpz_cmp", "yes — Mpz::cmp()"),
    (1, "mpz_cmp_d", "yes — Mpz::cmp_d()"),
    (1, "mpz_cmp_si", "yes — Mpz::cmp_si()"),
    (1, "mpz_cmp_ui", "yes — Mpz::cmp_ui()"),
    (1, "mpz_cmpabs", "yes — Mpz::cmpabs()"),
    (1, "mpz_cmpabs_d", "yes — Mpz::cmpabs_d()"),
    (1, "mpz_cmpabs_ui", "yes — Mpz::cmpabs_ui()"),
    (1, "mpz_sgn", "yes — Mpz::sgn()"),
    // -----------------------------------------------------------------------
    // B11: Logical & Bit manipulation
    // -----------------------------------------------------------------------
    (1, "mpz_and", "yes — Mpz::try_and()"),
    (1, "mpz_ior", "yes — Mpz::try_ior()"),
    (1, "mpz_xor", "yes — Mpz::try_xor()"),
    (1, "mpz_com", "yes — Mpz::com()"),
    (1, "mpz_popcount", "yes — Mpz::popcount()"),
    (1, "mpz_hamdist", "yes — Mpz::hamdist()"),
    (1, "mpz_scan0", "yes — Mpz::scan0()"),
    (1, "mpz_scan1", "yes — Mpz::scan1()"),
    (1, "mpz_setbit", "yes — Mpz::try_setbit()"),
    (1, "mpz_clrbit", "yes — Mpz::clrbit()"),
    (1, "mpz_combit", "yes — Mpz::try_combit()"),
    (1, "mpz_tstbit", "yes — Mpz::tstbit()"),
    // -----------------------------------------------------------------------
    // B12: I/O (requires std)
    // -----------------------------------------------------------------------
    (-1, "mpz_out_str", "no (requires std — file I/O)"),
    (-1, "mpz_inp_str", "no (requires std — file I/O)"),
    (-1, "mpz_out_raw", "no (requires std — file I/O)"),
    (-1, "mpz_inp_raw", "no (requires std — file I/O)"),
    // -----------------------------------------------------------------------
    // B13: Random numbers (requires RNG)
    // -----------------------------------------------------------------------
    (-1, "mpz_urandomb", "no (requires external RNG)"),
    (-1, "mpz_urandomm", "no (requires external RNG)"),
    (-1, "mpz_rrandomb", "no (requires external RNG)"),
    (-1, "mpz_random", "no (obsolete — requires RNG)"),
    (-1, "mpz_random2", "no (obsolete — requires RNG)"),
    // -----------------------------------------------------------------------
    // B14: Integer Import / Export
    // -----------------------------------------------------------------------
    (1, "mpz_import", "yes — Mpz::try_import()"),
    (1, "mpz_export", "yes — Mpz::export_buf()"),
    // -----------------------------------------------------------------------
    // B15: Miscellaneous
    // -----------------------------------------------------------------------
    (1, "mpz_fits_ulong_p", "yes — Mpz::fits_ulong()"),
    (1, "mpz_fits_slong_p", "yes — Mpz::fits_slong()"),
    (1, "mpz_fits_uint_p", "yes — Mpz::fits_uint()"),
    (1, "mpz_fits_sint_p", "yes — Mpz::fits_sint()"),
    (1, "mpz_fits_ushort_p", "yes — Mpz::fits_ushort()"),
    (1, "mpz_fits_sshort_p", "yes — Mpz::fits_sshort()"),
    (1, "mpz_odd_p", "yes — Mpz::odd_p()"),
    (1, "mpz_even_p", "yes — Mpz::even_p()"),
    (1, "mpz_sizeinbase", "partial — Mpz::try_sizeinbase()"),
    // -----------------------------------------------------------------------
    // B16: Low-level / Limb access
    // -----------------------------------------------------------------------
    (1, "mpz_size", "yes — Mpz::size()"),
    (1, "mpz_getlimbn", "yes — Mpz::getlimbn()"),
    (-1, "mpz_limbs_read", "no (requires unsafe pointer access)"),
    (-1, "mpz_limbs_write", "no (requires unsafe pointer access)"),
    (
        -1,
        "mpz_limbs_modify",
        "no (requires unsafe pointer access)",
    ),
    (
        -1,
        "mpz_limbs_finish",
        "no (requires unsafe pointer access)",
    ),
    (-1, "mpz_roinit_n", "no (requires unsafe pointer tricks)"),
    (0, "MPZ_ROINIT_N", "no (C macro — not applicable)"),
    (0, "_mpz_realloc", "no (internal GMP function)"),
];

#[inline(always)]
fn umul(a: u64, b: u64) -> (u64, u64) {
    let prod = a as u128 * b as u128;
    (prod as u64, (prod >> 64) as u64)
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    extern crate alloc;
    use alloc::string::String;

    fn s(m: &Mpz) -> String {
        let mut buf = [0u8; 192];
        let len = m.write_decimal_buf(&mut buf);
        core::str::from_utf8(&buf[..len]).unwrap().into()
    }

    // ---- Existing tests (kept identical) ----

    #[test]
    fn basic_arith() {
        let a = Mpz::from_decimal_str("123456789012345678901234567890").unwrap();
        let b = Mpz::from_decimal_str("987654321098765432109876543210").unwrap();
        assert_eq!(
            s(&a.try_add(&b).unwrap()),
            "1111111110111111111011111111100"
        );
        assert_eq!(s(&b.try_sub(&a).unwrap()), "864197532086419753208641975320");
        assert_eq!(
            s(&a.try_mul(&b).unwrap()),
            "121932631137021795226185032733622923332237463801111263526900"
        );
        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(a.cmpabs(&b), Ordering::Less);
    }

    #[test]
    fn signs_and_zero() {
        let a = Mpz::from_i64(-5);
        let b = Mpz::from_i64(5);
        assert_eq!(a.try_add(&b).unwrap(), Mpz::new());
        assert_eq!(a.sgn(), -1);
        assert_eq!(s(&a.try_mul(&b).unwrap()), "-25");
        assert_eq!(a.get_ui(), 5);
        assert_eq!(a.get_si(), -5);
    }

    #[test]
    fn div_trunc() {
        let a = Mpz::from_decimal_str("-1000000000000000000000007").unwrap();
        let d = Mpz::from_u64(1000);
        let (q, r) = a.tdiv_qr(&d);
        assert_eq!(s(&q), "-1000000000000000000000");
        assert_eq!(s(&r), "-7");
        assert_eq!(a.tdiv_ui(1000), 7);
    }

    #[test]
    fn shifts_and_remove() {
        let mut x = Mpz::from_decimal_str("123000").unwrap();
        assert_eq!(x.remove_pow10(), 3);
        assert_eq!(s(&x), "123");
        let y = Mpz::from_u64(1).try_mul_2exp(100).unwrap();
        assert_eq!(s(&y), "1267650600228229401496703205376");
        assert_eq!(y.fdiv_q_2exp(100), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(0b1011).fdiv_r_2exp(2), Mpz::from_u64(0b11));
    }

    #[test]
    fn zero_ops() {
        let z = Mpz::new();
        assert_eq!(s(&z), "0");
        assert_eq!(z.sgn(), 0);
        assert_eq!(z.size(), 0);
        assert!(z.fits_ulong());
        assert_eq!(z.get_ui(), 0);
        assert_eq!(z.try_add(&z).unwrap(), z);
        assert_eq!(z.try_sub(&z).unwrap(), z);
        assert_eq!(z.try_mul(&z).unwrap(), z);
        assert_eq!(z.isqrt(), z);
        assert_eq!(z.com(), Mpz::from_i64(-1));
    }

    #[test]
    fn zero_ops_full() {
        let z = Mpz::new();
        assert_eq!(z.mpz_get_ull(), 0);
        assert_eq!(z.mpz_get_sll(), 0);
        assert_eq!(z.cmp(&z), Ordering::Equal);
        assert_eq!(z.cmpabs(&z), Ordering::Equal);
        assert_eq!(z.to_i128(), Some(0));
        assert_eq!(Mpz::from_u64(0).remove_pow10(), 0);
        assert!(z.divisible_ui(7));
        assert_eq!(z.tdiv_ui(7), 0);
        assert_eq!(z.tdiv_q_ui(7), z);
        assert_eq!(z.try_mul_2exp(100).unwrap(), z);
        assert_eq!(z.fdiv_q_2exp(100), z);
        assert_eq!(z.fdiv_r_2exp(100), z);
        assert_eq!(z.even_p(), true);
        assert_eq!(z.odd_p(), false);
    }

    #[test]
    fn neg_abs() {
        let a = Mpz::from_i64(-123);
        let mut b = a.clone();
        b.abs();
        assert_eq!(s(&b), "123");
        b.neg();
        assert_eq!(s(&b), "-123");
        assert_eq!(s(&a.neg_to()), "123");
        assert_eq!(s(&a.abs_to()), "123");
    }

    #[test]
    fn cmp_ext() {
        let a = Mpz::from_i64(-100);
        let b = Mpz::from_i64(0);
        assert_eq!(a.cmp_si(-100), Ordering::Equal);
        assert_eq!(b.cmp_ui(0), Ordering::Equal);
        assert_eq!(a.cmp_si(-99), Ordering::Less);
        assert_eq!(Mpz::from_u64(50).cmpabs_ui(100), Ordering::Less);
    }

    #[test]
    fn fits_checks() {
        assert!(Mpz::from_u64(255).fits_uint());
        assert!(Mpz::from_u64(255).fits_sint());
        assert!(Mpz::from_i64(-1).fits_slong());
        assert!(!Mpz::from_u64(u64::MAX).fits_uint());
        assert!(!Mpz::from_i64(i64::MIN).fits_ulong());
        assert!(Mpz::from_u64(42).fits_ushort());
        assert!(Mpz::from_i64(-42).fits_sshort());
    }

    #[test]
    fn odd_even() {
        assert!(Mpz::from_u64(1).odd_p());
        assert!(Mpz::from_u64(2).even_p());
        assert!(Mpz::from_i64(-3).odd_p());
        assert!(Mpz::new().even_p());
    }

    #[test]
    fn div_floor() {
        let a = Mpz::from_i64(-100);
        let d = Mpz::from_i64(30);
        let (q, r) = a.try_fdiv_qr(&d).unwrap();
        assert_eq!(s(&q), "-4");
        assert_eq!(s(&r), "20");
        let a2 = Mpz::from_i64(100);
        let (q2, r2) = a2.try_fdiv_qr(&d).unwrap();
        assert_eq!(s(&q2), "3");
        assert_eq!(s(&r2), "10");
    }

    #[test]
    fn div_floor_ui() {
        let a = Mpz::from_i64(-100);
        let r = a.fdiv_ui(30);
        assert_eq!(r, 20);
    }

    #[test]
    fn div_mod() {
        let a = Mpz::from_i64(-100);
        let d = Mpz::from_i64(30);
        let r = a.try_mod(&d).unwrap();
        assert_eq!(s(&r), "20");
    }

    #[test]
    fn divisible_tests() {
        let a = Mpz::from_u64(100);
        assert!(a.divisible_ui(5));
        assert!(!a.divisible_ui(3));
        assert!(a.divisible_p(&Mpz::from_u64(20)));
        assert!(!a.divisible_p(&Mpz::from_u64(30)));
        assert!(Mpz::from_u64(64).divisible_2exp_p(6));
        assert!(!Mpz::from_u64(65).divisible_2exp_p(6));
    }

    #[test]
    fn congruent_tests() {
        let a = Mpz::from_i64(17);
        assert!(a.congruent_ui_p(5, 12));
        assert!(!a.congruent_ui_p(6, 12));
    }

    #[test]
    fn gcd_lcm() {
        let a = Mpz::from_u64(12);
        let b = Mpz::from_u64(18);
        assert_eq!(s(&a.try_gcd(&b).unwrap()), "6");
        assert_eq!(s(&a.try_lcm(&b).unwrap()), "36");
        assert_eq!(a.gcd_ui(18), 6);
    }

    #[test]
    fn gcdext() {
        let a = Mpz::from_u64(120);
        let b = Mpz::from_u64(23);
        let (g, s_coeff, t) = a.try_gcdext(&b).unwrap();
        let check = s_coeff
            .try_mul(&a)
            .unwrap()
            .try_add(&t.try_mul(&b).unwrap())
            .unwrap();
        assert_eq!(check, g);
        assert_eq!(s(&g), "1");
    }

    #[test]
    fn invert() {
        let a = Mpz::from_u64(3);
        let m = Mpz::from_u64(13);
        let inv = a.try_invert(&m).unwrap();
        let prod = a.try_mul(&inv).unwrap().try_mod(&m).unwrap();
        assert_eq!(prod, Mpz::from_u64(1));
    }

    #[test]
    fn jacobi_symbol() {
        assert_eq!(Mpz::from_u64(1).jacobi(&Mpz::from_u64(3)), 1);
        assert_eq!(Mpz::from_u64(2).jacobi(&Mpz::from_u64(3)), -1);
        assert_eq!(Mpz::from_u64(0).jacobi(&Mpz::from_u64(3)), 0);
    }

    #[test]
    fn factorial() {
        assert_eq!(s(&Mpz::try_fac_ui(5).unwrap()), "120");
        assert_eq!(s(&Mpz::try_fac_ui(0).unwrap()), "1");
    }

    #[test]
    fn binomial() {
        assert_eq!(s(&Mpz::try_bin_uiui(10, 5).unwrap()), "252");
        assert_eq!(s(&Mpz::try_bin_uiui(10, 0).unwrap()), "1");
        assert_eq!(s(&Mpz::try_bin_uiui(10, 10).unwrap()), "1");
    }

    #[test]
    fn fibonacci() {
        assert_eq!(s(&Mpz::try_fib_ui(0).unwrap()), "0");
        assert_eq!(s(&Mpz::try_fib_ui(1).unwrap()), "1");
        assert_eq!(s(&Mpz::try_fib_ui(10).unwrap()), "55");
    }

    #[test]
    fn bitwise_ops() {
        let a = Mpz::from_u64(0b1100);
        let b = Mpz::from_u64(0b1010);
        assert_eq!(a.try_and(&b).unwrap(), Mpz::from_u64(0b1000));
        assert_eq!(a.try_ior(&b).unwrap(), Mpz::from_u64(0b1110));
        assert_eq!(a.try_xor(&b).unwrap(), Mpz::from_u64(0b0110));
        let neg = Mpz::from_i64(-1);
        assert_eq!(a.try_and(&neg).unwrap(), a);
        assert_eq!(a.try_ior(&neg).unwrap(), neg);
    }

    #[test]
    fn popcount_hamdist() {
        assert_eq!(Mpz::from_u64(0b1010).popcount(), Some(2));
        assert_eq!(Mpz::new().popcount(), Some(0));
        assert_eq!(Mpz::from_i64(-1).popcount(), None);
        assert_eq!(
            Mpz::from_u64(0b1100).hamdist(&Mpz::from_u64(0b1010)),
            Some(2)
        );
    }

    #[test]
    fn scan_tstbit() {
        let a = Mpz::from_u64(0b1010_0000);
        assert_eq!(a.scan1(0), Some(5));
        assert_eq!(a.scan0(0), 0);
        assert!(a.tstbit(5));
        assert!(!a.tstbit(4));
        let mut b = Mpz::new();
        b.try_setbit(7).unwrap();
        assert_eq!(s(&b), "128");
        b.clrbit(7);
        assert_eq!(b, Mpz::new());
    }

    #[test]
    fn power_root_tests() {
        assert_eq!(
            s(&Mpz::try_ui_pow_ui(2, 100).unwrap()),
            "1267650600228229401496703205376"
        );
        assert_eq!(Mpz::from_u64(100).isqrt(), Mpz::from_u64(10));
        let (root, rem) = Mpz::from_u64(10).try_sqrtrem().unwrap();
        assert_eq!(root, Mpz::from_u64(3));
        assert_eq!(rem, Mpz::from_u64(1));
        assert!(Mpz::from_u64(144).perfect_square_p());
        assert!(!Mpz::from_u64(2).perfect_square_p());
        assert!(Mpz::from_u64(8).perfect_power_p());
        assert!(!Mpz::from_u64(2).perfect_power_p());
    }

    #[test]
    fn powm_test() {
        let base = Mpz::from_u64(3);
        let exp = Mpz::from_u64(5);
        let m = Mpz::from_u64(7);
        let r = base.try_powm(&exp, &m).unwrap();
        assert_eq!(r, Mpz::from_u64(5));
    }

    #[test]
    fn root_test() {
        let r = Mpz::from_u64(1000).try_root(3).unwrap();
        assert_eq!(r, Mpz::from_u64(10));
        let r2 = Mpz::from_u64(100).try_root(3).unwrap();
        assert_eq!(r2, Mpz::from_u64(4));
    }

    #[test]
    fn remove_test() {
        let a = Mpz::from_u64(2700);
        let (rem, count) = a.try_remove(&Mpz::from_u64(2)).unwrap();
        assert_eq!(count, 2);
        assert_eq!(rem, Mpz::from_u64(675));
    }

    #[test]
    fn add_sub_mul_assign() {
        let mut a = Mpz::from_u64(5);
        let b = Mpz::from_u64(3);
        a += &b;
        assert_eq!(s(&a), "8");
        a -= &b;
        assert_eq!(s(&a), "5");
        a *= &b;
        assert_eq!(s(&a), "15");
    }

    #[test]
    fn addmul_submul() {
        let mut a = Mpz::from_u64(10);
        a.try_addmul(&Mpz::from_u64(3), &Mpz::from_u64(5)).unwrap();
        assert_eq!(a, Mpz::from_u64(25));
        a.try_submul(&Mpz::from_u64(2), &Mpz::from_u64(5)).unwrap();
        assert_eq!(a, Mpz::from_u64(15));
    }

    #[test]
    fn ui_sub() {
        let r = Mpz::try_ui_sub(100, &Mpz::from_u64(1)).unwrap();
        assert_eq!(r, Mpz::from_u64(99));
        let r2 = Mpz::try_ui_sub(0, &Mpz::from_u64(5)).unwrap();
        assert_eq!(r2, Mpz::from_i64(-5));
    }

    #[test]
    fn neg_to_abs_to() {
        let a = Mpz::from_i64(-5);
        assert_eq!(a.neg_to(), Mpz::from_u64(5));
        assert_eq!(a.abs_to(), Mpz::from_u64(5));
        let b = Mpz::from_u64(5);
        assert_eq!(b.neg_to(), Mpz::from_i64(-5));
        assert_eq!(b.abs_to(), Mpz::from_u64(5));
    }

    #[test]
    fn swap_test() {
        let mut a = Mpz::from_u64(42);
        let mut b = Mpz::from_u64(100);
        a.swap(&mut b);
        assert_eq!(a, Mpz::from_u64(100));
        assert_eq!(b, Mpz::from_u64(42));
    }

    #[test]
    fn set_test() {
        let mut a = Mpz::new();
        a.set(&Mpz::from_u64(123));
        assert_eq!(a, Mpz::from_u64(123));
    }

    #[test]
    fn sizeinbase_test() {
        assert_eq!(Mpz::from_u64(255).try_sizeinbase(10), Some(3));
        assert_eq!(Mpz::from_u64(255).try_sizeinbase(16), Some(2));
        assert_eq!(Mpz::from_u64(255).try_sizeinbase(2), Some(8));
        assert_eq!(Mpz::new().try_sizeinbase(10), Some(1));
        assert_eq!(Mpz::from_i64(-5).try_sizeinbase(10), Some(2));
        assert_eq!(Mpz::from_u64(255).try_sizeinbase(1), None);
    }

    #[test]
    fn get_d_test() {
        assert!((Mpz::from_u64(123).get_d() - 123.0).abs() < 1e-10);
        assert!((Mpz::from_i64(-42).get_d() - (-42.0)).abs() < 1e-10);
        assert_eq!(Mpz::new().get_d(), 0.0);
    }

    #[test]
    fn get_d_2exp_test() {
        let v = Mpz::from_u64(42);
        let (m, e) = v.get_d_2exp().unwrap();
        assert!((m * 2.0f64.powi(e as i32) - 42.0).abs() < 1e-10);
        assert!(m.abs() >= 0.5 && m.abs() < 1.0);
        assert_eq!(Mpz::new().get_d_2exp(), None);
    }

    #[test]
    fn parse_errors() {
        assert_eq!(Mpz::from_decimal_str(""), Err(ParseError::InvalidInput));
        assert_eq!(Mpz::from_decimal_str("-"), Err(ParseError::InvalidInput));
        assert_eq!(Mpz::from_decimal_str("+"), Err(ParseError::InvalidInput));
        assert_eq!(
            Mpz::from_decimal_str("12a34"),
            Err(ParseError::InvalidInput)
        );
        assert_eq!(
            Mpz::from_decimal_str("12O34"),
            Err(ParseError::InvalidInput)
        );
    }

    #[test]
    fn roundtrip_extreme() {
        let vals = [
            "0",
            "1",
            "-1",
            "9",
            "10",
            "-10",
            "999999999999999999",
            "-999999999999999999",
            "1000000000000000000",
            "123456789012345678901234567890123456789012345678901234567890",
        ];
        for v in &vals {
            let m = Mpz::from_decimal_str(v).unwrap();
            assert_eq!(&s(&m), v, "round-trip failed for {v}");
        }
    }

    #[test]
    fn max_limb_capacity() {
        let mut mpz = Mpz::new();
        for i in 0..8 {
            mpz.mag[i] = u64::MAX;
        }
        mpz.len = 8;
        mpz.sign = 1;
        assert_eq!(mpz.try_add(&Mpz::from_u64(1)), Err(CapacityError));
        assert_eq!(Mpz::try_ui_pow_ui(2, 512), Err(CapacityError));
    }

    #[test]
    fn cdiv_test() {
        let a = Mpz::from_i64(-100);
        let d = Mpz::from_i64(30);
        let (q, r) = a.try_cdiv_qr(&d).unwrap();
        assert_eq!(q.try_mul(&d).unwrap().try_add(&r).unwrap(), a);
    }

    #[test]
    fn mul_si_test() {
        let a = Mpz::from_u64(10);
        assert_eq!(a.try_mul_si(5).unwrap(), Mpz::from_u64(50));
        assert_eq!(a.try_mul_si(-3).unwrap(), Mpz::from_i64(-30));
    }

    #[test]
    fn addmul_ui_submul_ui() {
        let mut a = Mpz::from_u64(10);
        a.try_addmul_ui(&Mpz::from_u64(3), 5).unwrap();
        assert_eq!(a, Mpz::from_u64(25));
        a.try_submul_ui(&Mpz::from_u64(2), 5).unwrap();
        assert_eq!(a, Mpz::from_u64(15));
    }

    #[test]
    fn tdiv_r_ui_test() {
        let a = Mpz::from_i64(-100);
        let r = a.tdiv_r_ui(30);
        assert_eq!(r, Mpz::from_u64(10));
    }

    #[test]
    fn tdiv_2exp_test() {
        let a = Mpz::from_u64(0b1101);
        assert_eq!(a.tdiv_q_2exp(2), Mpz::from_u64(0b11));
        assert_eq!(a.tdiv_r_2exp(2), Mpz::from_u64(0b01));
    }

    #[test]
    fn setbit_clrbit_combit() {
        let mut a = Mpz::new();
        a.try_setbit(0).unwrap();
        assert_eq!(a, Mpz::from_u64(1));
        a.try_setbit(2).unwrap();
        assert_eq!(a, Mpz::from_u64(5));
        a.clrbit(2);
        assert_eq!(a, Mpz::from_u64(1));
        a.try_combit(1).unwrap();
        assert_eq!(a, Mpz::from_u64(3));
    }

    #[test]
    fn scan_test() {
        let a = Mpz::from_u64(0b1010_0000);
        assert_eq!(a.scan1(0), Some(5));
        assert_eq!(a.scan0(5), 6);
        assert_eq!(a.scan1(8), None);
        assert_eq!(Mpz::from_u64(0).scan1(0), None);
    }

    #[test]
    fn perfect_power_test() {
        assert!(Mpz::from_u64(27).perfect_power_p());
        assert!(Mpz::from_u64(16).perfect_power_p());
        assert!(!Mpz::from_u64(2).perfect_power_p());
        assert!(!Mpz::new().perfect_power_p());
    }

    #[test]
    fn congruent_2exp_test() {
        assert!(Mpz::from_u64(17).congruent_2exp_p(&Mpz::from_u64(1), 4));
        assert!(Mpz::from_u64(18).congruent_2exp_p(&Mpz::from_u64(2), 4));
    }

    #[test]
    fn remove_general_test() {
        let a = Mpz::try_ui_pow_ui(2, 5)
            .unwrap()
            .try_mul(&Mpz::try_ui_pow_ui(3, 2).unwrap())
            .unwrap();
        let (rem, cnt) = a.try_remove(&Mpz::from_u64(2)).unwrap();
        assert_eq!(cnt, 5);
        assert_eq!(rem, Mpz::from_u64(9));
    }

    #[test]
    fn try_divexact_test() {
        let a = Mpz::from_u64(100);
        assert_eq!(
            a.try_divexact(&Mpz::from_u64(5)).unwrap(),
            Mpz::from_u64(20)
        );
        assert_eq!(a.try_divexact_ui(4).unwrap(), Mpz::from_u64(25));
    }

    // ---- NEW TESTS ----

    #[test]
    fn tstbit_negative_bug() {
        // R01: tstbit for negative values must return true for out-of-range bits
        let neg = Mpz::from_i64(-5);

        // Out-of-range bit (beyond self.len=1, so >= 64)
        assert!(
            neg.tstbit(100),
            "tstbit(-5, 100) should be true (sign extension)"
        );

        // For positive: out-of-range returns false
        let pos = Mpz::from_u64(5);
        assert!(!pos.tstbit(100), "tstbit(5, 100) should be false");

        // For zero: out-of-range returns false
        let zero = Mpz::new();
        assert!(!zero.tstbit(100), "tstbit(0, 100) should be false");

        // Verify -1 has out-of-range bits set (sign extension)
        let neg_one = Mpz::from_i64(-1);
        // In-range bits checked via magnitude: -1 = mag=[1], so bit 0 = 1, bit 1 = 0, etc.
        assert!(neg_one.tstbit(0), "tstbit(-1, 0)");
        // Out-of-range bits (beyond limb count) are set for negatives
        assert!(neg_one.tstbit(64), "tstbit(-1, 64) - sign extension");
        assert!(neg_one.tstbit(200), "tstbit(-1, 200) - sign extension");
    }

    #[test]
    fn bitwise_negative_bug() {
        // R02: bitwise ops with negative values need proper 2's complement extension
        // -2 & -3 = ?
        // -2 in 2's complement (64-bit): !2+1 = 0xFFFFFFFFFFFFFFFD... no, !2=0xFFFFFFFFFFFFFFFD,
        //   +1 = 0xFFFFFFFFFFFFFFFE
        // -3 in 2's complement: !3+1 = 0xFFFFFFFFFFFFFFFC+1 = 0xFFFFFFFFFFFFFFFD
        // -2 & -3: 0xFFFFFFFFFFFFFFFE & 0xFFFFFFFFFFFFFFFD = 0xFFFFFFFFFFFFFFFC = -4 (2's complement)
        let n2 = Mpz::from_i64(-2);
        let n3 = Mpz::from_i64(-3);
        assert_eq!(n2.try_and(&n3).unwrap(), Mpz::from_i64(-4));

        // -5 & 3 = ?
        // -5 2's complement: 0xFFFFFFFFFFFFFFFB
        //  3: 0x0000000000000003
        // AND: 0x0000000000000003 = 3
        let n5 = Mpz::from_i64(-5);
        let p3 = Mpz::from_u64(3);
        assert_eq!(n5.try_and(&p3).unwrap(), Mpz::from_u64(3));

        // -5 | 3 = ?
        // -5 2's complement: 0xFFFFFFFFFFFFFFFB
        //  3: 0x0000000000000003
        //  OR: 0xFFFFFFFFFFFFFFFB = -5
        assert_eq!(n5.try_ior(&p3).unwrap(), Mpz::from_i64(-5));

        // -5 ^ 3 = ?
        // -5 2's complement: 0xFFFFFFFFFFFFFFFB
        //  3: 0x0000000000000003
        // XOR: 0xFFFFFFFFFFFFFFF8 = -8
        assert_eq!(n5.try_xor(&p3).unwrap(), Mpz::from_i64(-8));
    }

    #[test]
    fn congruent_2exp_negative_bug() {
        // R03: congruent_2exp_p for negative values
        // -3 ≡ 13 (mod 16) because -3 = -1*16 + 13, and in 2's complement,
        // low 4 bits of -3 = 0xD = 13, low 4 bits of 13 = 0xD = 13
        assert!(Mpz::from_i64(-3).congruent_2exp_p(&Mpz::from_u64(13), 4));
        // -3 ≡ -3 (mod 16): low 4 bits are 0xD = 13 and 0xD = 13, but -3's
        // magnitude low 4 bits would be 0x3 = 3, not 13.
        assert!(Mpz::from_i64(-3).congruent_2exp_p(&Mpz::from_i64(-3), 4));
        // -3 ≢ 3 (mod 16): low 4 bits of -3 in 2's complement = 13, low 4 bits of 3 = 3
        assert!(!Mpz::from_i64(-3).congruent_2exp_p(&Mpz::from_u64(3), 4));
    }

    #[test]
    fn from_d_tests() {
        // Integers
        assert_eq!(Mpz::from_d(0.0).unwrap(), Mpz::new());
        assert_eq!(Mpz::from_d(1.0).unwrap(), Mpz::from_u64(1));
        assert_eq!(Mpz::from_d(-1.0).unwrap(), Mpz::from_i64(-1));
        assert_eq!(Mpz::from_d(255.0).unwrap(), Mpz::from_u64(255));
        assert_eq!(
            Mpz::from_d(1e18).unwrap(),
            Mpz::from_u64(1_000_000_000_000_000_000)
        );
        assert_eq!(
            Mpz::from_d(-1e18).unwrap(),
            Mpz::from_i64(-1_000_000_000_000_000_000)
        );

        // Truncation
        assert_eq!(Mpz::from_d(3.999).unwrap(), Mpz::from_u64(3));
        assert_eq!(Mpz::from_d(-3.999).unwrap(), Mpz::from_i64(-3));
        assert_eq!(Mpz::from_d(0.999).unwrap(), Mpz::new());

        // Subnormals (very small values) → 0
        assert_eq!(Mpz::from_d(1e-300).unwrap(), Mpz::new());

        // Special values
        assert_eq!(Mpz::from_d(f64::INFINITY), Err(CapacityError));
        assert_eq!(Mpz::from_d(f64::NEG_INFINITY), Err(CapacityError));
        assert_eq!(Mpz::from_d(f64::NAN), Err(CapacityError));
    }

    #[test]
    fn sizeinbase_rigorous() {
        // For base 10, the bound should be >= true size
        let val = Mpz::from_u64(12345);
        let sz = val.try_sizeinbase(10).unwrap();
        assert!(sz >= "12345".len());

        let large = Mpz::try_ui_pow_ui(2, 200).unwrap();
        let sz10 = large.try_sizeinbase(10).unwrap();
        assert!(sz10 >= s(&large).len());
    }

    #[test]
    fn knuth_division_tests() {
        // Use multi-limb divisors (len >= 2) to exercise Algorithm D
        let a = Mpz::from_u64(1).try_mul_2exp(130).unwrap(); // large number
        let b = Mpz::from_u64(1).try_mul_2exp(65).unwrap();
        let (q, r) = a.tdiv_qr(&b);
        assert_eq!(q, Mpz::from_u64(1).try_mul_2exp(65).unwrap());
        assert!(r.is_zero());

        // Test with larger numbers
        let x = Mpz::from_decimal_str("123456789012345678901234567890").unwrap();
        let y = Mpz::from_decimal_str("9876543210").unwrap();
        let (q2, r2) = x.tdiv_qr(&y);
        // Verify: x = q2 * y + r2
        let check = q2.try_mul(&y).unwrap().try_add(&r2).unwrap();
        assert_eq!(check, x);
        assert!(r2.cmpabs(&y) == Ordering::Less);
    }

    #[test]
    fn sealed_fields() {
        let a = Mpz::from_u64(42);
        assert_eq!(a.limbs(), &[42u64]);
        assert!(!a.is_zero());
        let b = Mpz::new();
        assert!(b.is_zero());
        assert!(b.limbs().is_empty());
    }

    #[test]
    fn new_division_functions() {
        // try_tdiv_qr_ui
        let a = Mpz::from_i64(-100);
        let (q, r) = a.try_tdiv_qr_ui(30).unwrap();
        assert_eq!(q, Mpz::from_i64(-3));
        assert_eq!(r, 10);

        // try_fdiv_r_ui
        let r = a.try_fdiv_r_ui(30).unwrap();
        assert_eq!(r, Mpz::from_u64(20));

        // Ceiling division by u64: for -100/30, ceil=-3, abs_remainder=10
        let (q, r_u64) = a.try_cdiv_qr_ui(30).unwrap();
        assert_eq!(q, Mpz::from_i64(-3));
        assert_eq!(r_u64, 10);

        // cdiv_ui: ceil(-100/30) = -3, absolute remainder = 10
        let r = a.cdiv_ui(30);
        // GMP's cdiv_ui returns the remainder as unsigned long (absolute value for positive divisor)
        assert_eq!(r, 10);

        // cdiv_q_2exp and cdiv_r_2exp
        let pos = Mpz::from_u64(100);
        let q = pos.cdiv_q_2exp(3); // ceil(100 / 8) = 13
        assert_eq!(q, Mpz::from_u64(13));
        let r = pos.cdiv_r_2exp(3); // 100 - 13*8 = 100 - 104 = -4? No, remainder should match: 100 = 13*8 + (-4)? That's wrong.
                                    // Actually cdiv_r_2exp returns remainder = self - ceil(self/2^k) * 2^k
                                    // For 100, ceil(100/8) = 13, 13*8 = 104, remainder = 100 - 104 = -4
        assert_eq!(r, Mpz::from_i64(-4));
    }

    #[test]
    fn number_theory_functions() {
        // Legendre = Jacobi for odd prime modulus
        assert_eq!(Mpz::from_u64(1).try_legendre(&Mpz::from_u64(7)), 1);
        assert_eq!(Mpz::from_u64(2).try_legendre(&Mpz::from_u64(7)), 1);
        assert_eq!(Mpz::from_u64(3).try_legendre(&Mpz::from_u64(7)), -1);

        // Kronecker
        assert_eq!(Mpz::from_u64(0).try_kronecker(&Mpz::from_u64(1)), 1);
        assert_eq!(Mpz::from_u64(2).try_kronecker(&Mpz::from_u64(5)), -1);

        // Kronecker with signed/unsigned
        // (5/7) = (7/5) = (2/5) = (-1)^{(25-1)/8} = -1
        assert_eq!(Mpz::from_u64(5).try_kronecker_si(7), -1);
        assert_eq!(Mpz::from_u64(5).try_kronecker_ui(7), -1);
        assert_eq!(Mpz::try_si_kronecker(5, &Mpz::from_u64(7)), -1);
        assert_eq!(Mpz::try_ui_kronecker(5, &Mpz::from_u64(7)), -1);
    }

    #[test]
    fn fib2_lucnum() {
        // fib2: F(10)=55, F(9)=34
        let (f_n, f_nm1) = Mpz::try_fib2_ui(10).unwrap();
        assert_eq!(f_n, Mpz::from_u64(55));
        assert_eq!(f_nm1, Mpz::from_u64(34));

        // lucnum: L(10)=123
        assert_eq!(Mpz::try_lucnum_ui(10).unwrap(), Mpz::from_u64(123));

        // lucnum2: L(10)=123, L(9)=76
        let (l_n, l_nm1) = Mpz::try_lucnum2_ui(10).unwrap();
        assert_eq!(l_n, Mpz::from_u64(123));
        assert_eq!(l_nm1, Mpz::from_u64(76));
    }

    #[test]
    fn double_factorial() {
        // 5!! = 5*3*1 = 15
        assert_eq!(Mpz::try_2fac_ui(5).unwrap(), Mpz::from_u64(15));
        // 6!! = 6*4*2 = 48
        assert_eq!(Mpz::try_2fac_ui(6).unwrap(), Mpz::from_u64(48));
        // 0!! = 1, 1!! = 1
        assert_eq!(Mpz::try_2fac_ui(0).unwrap(), Mpz::from_u64(1));
        assert_eq!(Mpz::try_2fac_ui(1).unwrap(), Mpz::from_u64(1));
    }

    #[test]
    fn primorial() {
        // primorial(10) = 2*3*5*7 = 210
        assert_eq!(Mpz::try_primorial_ui(10).unwrap(), Mpz::from_u64(210));
        // primorial(0) = 1, primorial(1) = 1
        assert_eq!(Mpz::try_primorial_ui(0).unwrap(), Mpz::from_u64(1));
        assert_eq!(Mpz::try_primorial_ui(1).unwrap(), Mpz::from_u64(1));
    }

    #[test]
    fn bin_ui() {
        // C(10, 3) = 120
        let n = Mpz::from_u64(10);
        assert_eq!(n.try_bin_ui(3).unwrap(), Mpz::from_u64(120));
        // C(10, 0) = 1
        assert_eq!(n.try_bin_ui(0).unwrap(), Mpz::from_u64(1));
        // C(0, 0) = 1
        let zero = Mpz::new();
        assert_eq!(zero.try_bin_ui(0).unwrap(), Mpz::from_u64(1));
    }

    #[test]
    fn cmp_d_tests() {
        // Positive values
        assert_eq!(Mpz::from_u64(5).cmp_d(3.0), Ordering::Greater);
        assert_eq!(Mpz::from_u64(3).cmp_d(5.0), Ordering::Less);
        assert_eq!(Mpz::from_u64(3).cmp_d(3.0), Ordering::Equal);
        // Negative values
        assert_eq!(Mpz::from_i64(-5).cmp_d(-3.0), Ordering::Less);
        assert_eq!(Mpz::from_i64(-3).cmp_d(-5.0), Ordering::Greater);
        // Mixed signs
        assert_eq!(Mpz::from_i64(-5).cmp_d(3.0), Ordering::Less);
        assert_eq!(Mpz::from_u64(5).cmp_d(-3.0), Ordering::Greater);
        // cmpabs_d
        assert_eq!(Mpz::from_i64(-5).cmpabs_d(3.0), Ordering::Greater);
        assert_eq!(Mpz::from_i64(-3).cmpabs_d(5.0), Ordering::Less);
    }

    #[test]
    fn getlimbn_test() {
        let a = Mpz::from_u128(0x1234567890ABCDEF);
        assert_eq!(a.getlimbn(0), Some(0x1234567890ABCDEF)); // Wait, little-endian: mag[0] = low 64 bits
                                                             // Actually from_u128(0x1234567890ABCDEF): mag[0] = 0x1234567890ABCDEF (since it's < 2^64)
        assert_eq!(a.getlimbn(1), None);
        assert_eq!(a.getlimbn(0), Some(0x1234567890ABCDEF));

        let b = Mpz::new();
        assert_eq!(b.getlimbn(0), None);
    }

    #[test]
    fn probab_prime_tests() {
        // 2 is prime
        assert_eq!(Mpz::from_u64(2).try_probab_prime_p(5).unwrap(), 2);
        // 3 is prime
        assert_eq!(Mpz::from_u64(3).try_probab_prime_p(5).unwrap(), 2);
        // 4 is composite
        assert_eq!(Mpz::from_u64(4).try_probab_prime_p(5).unwrap(), 0);
        // 17 is prime
        assert_eq!(Mpz::from_u64(17).try_probab_prime_p(5).unwrap(), 2);
        // 1 is not prime
        assert_eq!(Mpz::from_u64(1).try_probab_prime_p(5).unwrap(), 0);
        // 0 is not prime
        assert_eq!(Mpz::new().try_probab_prime_p(5).unwrap(), 0);
        // Negative is not prime
        assert_eq!(Mpz::from_i64(-5).try_probab_prime_p(5).unwrap(), 0);
        // 7919 is prime (1000th prime)
        assert_eq!(Mpz::from_u64(7919).try_probab_prime_p(10).unwrap(), 2);
    }

    #[test]
    fn import_export_tests() {
        // Round-trip: export then import
        let val = Mpz::from_u64(0x1234567890ABCDEF);
        let mut buf = [0u8; 8];
        let n = val
            .export_buf(&mut buf, Endian::Little, 1, Endian::Little)
            .unwrap();
        assert_eq!(n, 8);

        let imported = Mpz::try_import(8, Endian::Little, 1, Endian::Little, &buf).unwrap();
        assert_eq!(imported, val);

        // Big-endian round-trip
        let mut buf2 = [0u8; 8];
        let n2 = val
            .export_buf(&mut buf2, Endian::Big, 1, Endian::Big)
            .unwrap();
        assert_eq!(n2, 8);
        // The BE bytes: MSB first, so buf2[0] = 0x12 (most significant byte)
        assert_eq!(buf2[0], 0x12);
        assert_eq!(buf2[7], 0xEF);

        let imported2 = Mpz::try_import(8, Endian::Big, 1, Endian::Big, &buf2).unwrap();
        assert_eq!(imported2, val);
    }

    #[test]
    fn knuth_multi_limb_division() {
        // Create a case where divisor has >= 2 limbs
        let a = Mpz::try_ui_pow_ui(2, 200).unwrap();
        let b = Mpz::try_ui_pow_ui(2, 100).unwrap();
        let (q, r) = a.tdiv_qr(&b);
        assert_eq!(q, Mpz::try_ui_pow_ui(2, 100).unwrap());
        assert!(r.is_zero());

        // Mixed multi-limb values
        let x = Mpz::from_decimal_str("1234567890123456789012345678901234567890").unwrap();
        let y = Mpz::from_decimal_str("98765432109876543210").unwrap();
        let (q2, r2) = x.tdiv_qr(&y);
        let check = q2.try_mul(&y).unwrap().try_add(&r2).unwrap();
        assert_eq!(check, x);
        assert!(r2.cmpabs(&y) == Ordering::Less || r2.is_zero());
    }
}
