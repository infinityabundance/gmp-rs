//! A pure-Rust, **no-unsafe**, **`no_std`**, **`no_alloc`** arbitrary-precision signed integer (`Mpz`),
//! faithful to the GMP `mpz_*` operations.
//!
//! # No-alloc guarantee
//!
//! This crate performs **zero heap allocations**. The internal limb storage is a fixed-capacity
//! array of `MPZ_MAX_LIMBS` (8) 64-bit words, giving a maximum precision of 512 bits
//! (≈ 154 decimal digits). Operations that would exceed this capacity return
//! [`CapacityError`](enum.CapacityError.html).
//!
//! # Representation
//!
//! Sign–magnitude. `mag` is little-endian base-2⁶⁴ limbs, stored in `mag[0..len]` with no trailing
//! zero limb (so the zero value is `len == 0` and `sign == 0`).

#![no_std]
#![forbid(unsafe_code)]

use core::cmp::Ordering;
use core::fmt;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Maximum number of 64-bit limbs in the internal representation.
///
/// 8 limbs = 512 bits ≈ 154 decimal digits. This is a compile-time constant;
/// operations that would exceed this capacity return [`CapacityError`].
///
/// This is sufficient for the GnuCOBOL decimal use case (max 38 digits ≈ 2 limbs)
/// with generous headroom.
pub const MPZ_MAX_LIMBS: usize = 8;

// ---------------------------------------------------------------------------
// Error types
// ---------------------------------------------------------------------------

/// The operation would exceed the fixed limb capacity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapacityError;

/// An error occurred while parsing a decimal string.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParseError {
    /// The input string is malformed (empty, invalid characters, misplaced sign).
    InvalidInput,
    /// The value parsed correctly but exceeds the fixed limb capacity.
    CapacityOverflow,
}

// ---------------------------------------------------------------------------
// Mpz type
// ---------------------------------------------------------------------------

/// An arbitrary-precision signed integer with **fixed capacity** (no heap allocation).
///
/// Internally a sign–magnitude representation with up to [`MPZ_MAX_LIMBS`] little-endian
/// 64-bit limbs.
#[derive(Clone, Debug, Eq)]
pub struct Mpz {
    /// -1, 0, or +1. Invariant: `sign == 0` iff `len == 0`.
    sign: i8,
    /// Number of active limbs in `mag[0..len]`. No trailing zero limb.
    len: usize,
    /// Fixed-capacity limb storage. Only `mag[0..len]` is meaningful.
    mag: [u64; MPZ_MAX_LIMBS],
}

// ---------------------------------------------------------------------------
// PartialEq
// ---------------------------------------------------------------------------

impl PartialEq for Mpz {
    fn eq(&self, other: &Self) -> bool {
        if self.sign != other.sign || self.len != other.len {
            return false;
        }
        self.mag[..self.len] == other.mag[..other.len]
    }
}

// ---------------------------------------------------------------------------
// Default
// ---------------------------------------------------------------------------

impl Default for Mpz {
    fn default() -> Self {
        Mpz::new()
    }
}

// ---------------------------------------------------------------------------
// Display (uses a fixed stack buffer, never allocates)
// ---------------------------------------------------------------------------

impl fmt::Display for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Max decimal digits for 512-bit value: ceil(512 * log10(2)) ≈ 155
        // + sign byte = 156. Use 192 for safety.
        let mut buf = [0u8; 192];
        let len = self.write_decimal_buf(&mut buf);
        let s = core::str::from_utf8(&buf[..len]).map_err(|_| fmt::Error)?;
        f.write_str(s)
    }
}

// ---------------------------------------------------------------------------
// PartialOrd / Ord
// ---------------------------------------------------------------------------

impl PartialOrd for Mpz {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(core::cmp::Ord::cmp(self, other))
    }
}

impl Ord for Mpz {
    fn cmp(&self, other: &Self) -> Ordering {
        // Re-use the public cmp method
        Mpz::cmp(self, other)
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Methods
// ---------------------------------------------------------------------------

impl Mpz {
    // ---- Construction ----

    /// Zero value.
    pub fn new() -> Self {
        Mpz {
            sign: 0,
            len: 0,
            mag: [0u64; MPZ_MAX_LIMBS],
        }
    }

    /// Construct from an unsigned 64-bit value (always fits).
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

    /// Construct from a signed 64-bit value (always fits).
    pub fn from_i64(v: i64) -> Self {
        if v == 0 {
            return Mpz::new();
        }
        if v > 0 {
            let mut mag = [0u64; MPZ_MAX_LIMBS];
            mag[0] = v as u64;
            Mpz {
                sign: 1,
                len: 1,
                mag,
            }
        } else {
            let mut mag = [0u64; MPZ_MAX_LIMBS];
            mag[0] = (v as i128).unsigned_abs() as u64;
            Mpz {
                sign: -1,
                len: 1,
                mag,
            }
        }
    }

    /// Construct from an unsigned 128-bit value (always fits — max 2 limbs).
    pub fn from_u128(v: u128) -> Self {
        if v == 0 {
            return Mpz::new();
        }
        let mut mag = [0u64; MPZ_MAX_LIMBS];
        let lo = v as u64;
        let hi = (v >> 64) as u64;
        mag[0] = lo;
        if hi != 0 {
            mag[1] = hi;
            Mpz {
                sign: 1,
                len: 2,
                mag,
            }
        } else {
            Mpz {
                sign: 1,
                len: 1,
                mag,
            }
        }
    }

    /// Construct from a signed 128-bit value (always fits — max 2 limbs).
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

    /// `mpz_set_ull`: construct from an unsigned 64-bit value.
    pub fn set_ull(val: u64) -> Self {
        Self::from_u64(val)
    }

    /// `mpz_set_sll`: construct from a signed 64-bit value.
    pub fn mpz_set_sll(val: i64) -> Self {
        Self::from_i64(val)
    }

    // ---- Conversion (no alloc) ----

    /// `mpz_get_ull`: the low 64-bit limb of the magnitude.
    pub fn mpz_get_ull(&self) -> u64 {
        if self.len == 0 {
            0
        } else {
            self.mag[0]
        }
    }

    /// `mpz_get_sll`: reconstruct a signed 64-bit host integer from the low
    /// limb and the sign, mirroring the C bit-for-bit.
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

    /// `mpz_get_ui`: the low 64 bits of the absolute value.
    pub fn get_ui(&self) -> u64 {
        if self.len == 0 {
            0
        } else {
            self.mag[0]
        }
    }

    /// `mpz_get_si`: low bits with sign.
    pub fn get_si(&self) -> i64 {
        let lo = self.get_ui();
        if self.sign < 0 {
            (lo as i64).wrapping_neg()
        } else {
            lo as i64
        }
    }

    /// `mpz_fits_ulong_p`: non-negative and a single limb.
    pub fn fits_ulong(&self) -> bool {
        self.sign >= 0 && self.len <= 1
    }

    /// The value as `i128` if it fits (≤ 2 limbs), else `None`.
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

    // ---- Set (in-place) ----

    /// `mpz_set_ui`.
    pub fn set_ui(&mut self, v: u64) {
        *self = Self::from_u64(v);
    }

    /// `mpz_set_si`.
    pub fn set_si(&mut self, v: i64) {
        *self = Self::from_i64(v);
    }

    // ---- Sign / compare ----

    /// `mpz_sgn`.
    pub fn sgn(&self) -> i32 {
        self.sign as i32
    }

    /// `mpz_neg` (in place).
    pub fn neg(&mut self) {
        self.sign = -self.sign;
    }

    /// `mpz_abs` (in place).
    pub fn abs(&mut self) {
        if self.sign != 0 {
            self.sign = 1;
        }
    }

    /// `mpz_cmp`: signed compare.
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

    /// `mpz_cmpabs`: compare absolute values.
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

    /// `mpz_size`: number of limbs.
    pub fn size(&self) -> usize {
        self.len
    }

    /// `mpz_sizeinbase(_, 2)`: number of significant bits (1 for zero).
    pub fn sizeinbase2(&self) -> usize {
        if self.len == 0 {
            return 1;
        }
        (self.len - 1) * 64 + (64 - self.mag[self.len - 1].leading_zeros() as usize)
    }

    // ---- Internal magnitude ops (all operate on &[u64] slices, write into fixed arrays) ----

    fn trim(&mut self) {
        while self.len > 0 && self.mag[self.len - 1] == 0 {
            self.len -= 1;
        }
        if self.len == 0 {
            self.sign = 0;
        }
    }

    /// `a + b` (magnitudes). Returns `None` if the result would exceed `MPZ_MAX_LIMBS` limbs.
    fn mag_add_len(a: &[u64], b: &[u64], out: &mut [u64]) -> Option<usize> {
        let max = a.len().max(b.len());
        if max >= MPZ_MAX_LIMBS {
            // Check if there will be a carry into an extra limb
            if a.len() == MPZ_MAX_LIMBS || b.len() == MPZ_MAX_LIMBS {
                // Need to detect carry out of the top limb
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
                    return None; // would overflow MPZ_MAX_LIMBS
                }
                return Some(max);
            }
        }
        let mut carry = 0u128;
        for i in 0..max {
            let va = if i < a.len() { a[i] as u128 } else { 0 };
            let vb = if i < b.len() { b[i] as u128 } else { 0 };
            let s = va + vb + carry;
            out[i] = s as u64;
            carry = s >> 64;
        }
        let mut result_len = max;
        if carry != 0 {
            if max >= MPZ_MAX_LIMBS {
                return None;
            }
            out[max] = carry as u64;
            result_len = max + 1;
        }
        Some(result_len)
    }

    /// `a - b` (magnitudes, `a >= b` assumed).
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
        // Trim trailing zeros
        let mut result_len = a.len();
        while result_len > 0 && out[result_len - 1] == 0 {
            result_len -= 1;
        }
        result_len
    }

    /// `a * b` (magnitudes). Returns `None` if the result would exceed `MPZ_MAX_LIMBS` limbs.
    fn mag_mul_len(a: &[u64], b: &[u64], out: &mut [u64]) -> Option<usize> {
        if a.is_empty() || b.is_empty() {
            return Some(0);
        }
        let result_limbs = a.len() + b.len();
        if result_limbs > MPZ_MAX_LIMBS {
            return None;
        }
        // Zero out the output
        for o in out.iter_mut().take(result_limbs) {
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
                let idx = i + b.len();
                out[idx] = out[idx].wrapping_add(carry as u64);
            }
        }
        // Trim trailing zeros
        let mut result_len = result_limbs;
        while result_len > 0 && out[result_len - 1] == 0 {
            result_len -= 1;
        }
        Some(result_len)
    }

    /// Divide magnitude by a single u64. Returns `(quotient_len, remainder)`.
    /// `qbuf` must have at least `a.len()` elements.
    fn mag_divmod_u64_len(a: &[u64], d: u64, qbuf: &mut [u64]) -> (usize, u64) {
        let mut rem: u128 = 0;
        for i in (0..a.len()).rev() {
            let cur = (rem << 64) | a[i] as u128;
            qbuf[i] = (cur / d as u128) as u64;
            rem = cur % d as u128;
        }
        // Trim quotient
        let mut qlen = a.len();
        while qlen > 0 && qbuf[qlen - 1] == 0 {
            qlen -= 1;
        }
        (qlen, rem as u64)
    }

    // ---- Arithmetic (fallible — return Result) ----

    /// `mpz_add`: `self + other`.
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
            let cmp = self.cmpabs(other);
            match cmp {
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

    /// `mpz_sub`: `self - other`.
    pub fn try_sub(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        let neg = Mpz {
            sign: -other.sign,
            ..other.clone()
        };
        self.try_add(&neg)
    }

    /// `mpz_add_ui`: `self + v`.
    pub fn try_add_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_add(&Mpz::from_u64(v))
    }

    /// `mpz_sub_ui`: `self - v`.
    pub fn try_sub_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        self.try_sub(&Mpz::from_u64(v))
    }

    /// `mpz_mul`: `self * other`.
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

    /// `mpz_mul_ui`: `self * v`.
    pub fn try_mul_ui(&self, v: u64) -> Result<Mpz, CapacityError> {
        let v_mpz = Mpz::from_u64(v);
        // Avoid clone: construct a reference to a temporary
        self.try_mul(&v_mpz)
    }

    /// `mpz_mul_2exp`: `self << bits`.
    pub fn try_mul_2exp(&self, bits: u32) -> Result<Mpz, CapacityError> {
        if self.sign == 0 {
            return Ok(Mpz::new());
        }
        let limb_shift = (bits / 64) as usize;
        let bit_shift = bits % 64;
        let needed = self.len + limb_shift + if bit_shift != 0 { 1 } else { 0 };
        if needed > MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        let mut result = Mpz::new();
        for i in 0..limb_shift {
            result.mag[i] = 0;
        }
        if bit_shift == 0 {
            result.mag[limb_shift..limb_shift + self.len].copy_from_slice(&self.mag[..self.len]);
            result.len = limb_shift + self.len;
        } else {
            let mut carry = 0u64;
            let mut idx = limb_shift;
            for &l in self.mag[..self.len].iter() {
                result.mag[idx] = (l << bit_shift) | carry;
                carry = l >> (64 - bit_shift);
                idx += 1;
            }
            if carry != 0 {
                result.mag[idx] = carry;
                result.len = idx + 1;
            } else {
                result.len = idx;
            }
        }
        result.sign = self.sign;
        result.trim();
        Ok(result)
    }

    /// `mpz_ui_pow_ui`: `base^exp` using exponentiation by squaring.
    pub fn try_ui_pow_ui(base: u64, exp: u32) -> Result<Mpz, CapacityError> {
        let b = Mpz::from_u64(base);
        b.try_pow_ui(exp)
    }

    /// `mpz_pow_ui`: `self^n` by exponentiation by squaring.
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

    // ---- Division ----

    /// Divide magnitude by a single u64, returning `(quotient, remainder)`.
    fn mag_divmod_u64(&self, d: u64) -> (Mpz, u64) {
        let mut q = Mpz::new();
        let (qlen, rem) = Self::mag_divmod_u64_len(&self.mag[..self.len], d, &mut q.mag);
        q.len = qlen;
        q.sign = if qlen == 0 { 0 } else { self.sign };
        (q, rem)
    }

    /// `mpz_tdiv_q_ui`: truncated quotient by a u64.
    pub fn tdiv_q_ui(&self, d: u64) -> Mpz {
        let (q, _) = self.mag_divmod_u64(d);
        q
    }

    /// `mpz_tdiv_ui`: the absolute remainder mod `d`.
    pub fn tdiv_ui(&self, d: u64) -> u64 {
        self.mag_divmod_u64(d).1
    }

    /// `mpz_divisible_ui_p`.
    pub fn divisible_ui(&self, d: u64) -> bool {
        self.sign == 0 || self.mag_divmod_u64(d).1 == 0
    }

    /// Schoolbook long division of magnitudes (binary).
    /// Returns `(quotient_mag_len, remainder_mag_len)`.
    /// `qbuf` must have `a.len()` elements.
    fn mag_divmod_len(a: &[u64], d: &[u64], qbuf: &mut [u64], rbuf: &mut [u64]) -> (usize, usize) {
        let nbits = a.len() * 64;
        // Zero out qbuf
        for o in qbuf.iter_mut().take(a.len()) {
            *o = 0;
        }
        let mut rem = Mpz::new();
        // dm = d
        for i in (0..nbits).rev() {
            // rem = (rem << 1) | bit_i(a)
            rem = rem.try_mul_2exp(1).unwrap(); // can't overflow: 1 extra bit at a time
            let bit = (a[i / 64] >> (i % 64)) & 1;
            if bit != 0 {
                // rem.mag[0] |= 1;
                if rem.len == 0 {
                    rem.mag[0] = 1;
                    rem.len = 1;
                    rem.sign = 1;
                } else {
                    rem.mag[0] |= 1;
                }
            }
            // Compare rem.mag with d
            if Self::cmp_mag_slice(&rem.mag[..rem.len], d) != Ordering::Less {
                // rem = rem - d: need to copy because mag_sub_len needs &mut and & immut
                let mut tmp = [0u64; MPZ_MAX_LIMBS];
                let new_len = Self::mag_sub_len(&rem.mag[..rem.len], d, &mut tmp);
                rem.mag[..new_len].copy_from_slice(&tmp[..new_len]);
                rem.len = new_len;
                if rem.len == 0 {
                    rem.sign = 0;
                }
                qbuf[i / 64] |= 1u64 << (i % 64);
            }
        }
        // Trim quotient
        let mut qlen = a.len();
        while qlen > 0 && qbuf[qlen - 1] == 0 {
            qlen -= 1;
        }
        // Copy remainder
        rbuf[..rem.len].copy_from_slice(&rem.mag[..rem.len]);
        (qlen, rem.len)
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

    /// Full truncated division `self / d`. Returns `(quotient, remainder)`.
    pub fn tdiv_qr(&self, d: &Mpz) -> (Mpz, Mpz) {
        debug_assert!(d.sign != 0, "division by zero");
        if Self::cmp_mag_slice(&self.mag[..self.len], &d.mag[..d.len]) == Ordering::Less {
            return (Mpz::new(), self.clone());
        }
        let mut qbuf = [0u64; MPZ_MAX_LIMBS];
        let mut rbuf = [0u64; MPZ_MAX_LIMBS];
        let (qlen, rlen) =
            Self::mag_divmod_len(&self.mag[..self.len], &d.mag[..d.len], &mut qbuf, &mut rbuf);
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

    /// `mpz_tdiv_q`.
    pub fn tdiv_q(&self, d: &Mpz) -> Mpz {
        self.tdiv_qr(d).0
    }

    /// `mpz_tdiv_r`.
    pub fn tdiv_r(&self, d: &Mpz) -> Mpz {
        self.tdiv_qr(d).1
    }

    /// `mpz_fdiv_r_2exp`: the low `bits` bits.
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
            // Zero out higher limbs
            for i in (limbs + if rem_bits != 0 { 1 } else { 0 })..result.len {
                result.mag[i] = 0;
            }
        }
        result.trim();
        result
    }

    /// `mpz_fdiv_q_2exp`: `self >> bits`.
    pub fn fdiv_q_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let limb_shift = (bits / 64) as usize;
        let bit_shift = bits % 64;
        if limb_shift >= self.len {
            return Mpz::new();
        }
        let mut result = Mpz::new();
        let src_len = self.len - limb_shift;
        result.mag[..src_len].copy_from_slice(&self.mag[limb_shift..self.len]);
        result.len = src_len;
        result.sign = self.sign;
        if bit_shift != 0 {
            let mut carry = 0u64;
            for i in (0..result.len).rev() {
                let new = (result.mag[i] >> bit_shift) | carry;
                carry = result.mag[i] << (64 - bit_shift);
                result.mag[i] = new;
            }
        }
        result.trim();
        result
    }

    /// `mpz_sqrt`: floor of integer square root (Newton iteration).
    pub fn isqrt(&self) -> Mpz {
        if self.sgn() <= 0 {
            return Mpz::new();
        }
        // Initial over-estimate: 2^ceil(bits/2)
        let mut x = Mpz::from_u64(1);
        let shift = self.sizeinbase2().div_ceil(2) as u32;
        // try_mul_2exp can't fail for a single-limb value shifted by <= 256
        x = x
            .try_mul_2exp(shift)
            .unwrap_or(Mpz::from_u64(1).try_mul_2exp(255).unwrap());
        loop {
            let y = x.try_add(&self.tdiv_q(&x)).unwrap().fdiv_q_2exp(1);
            if y.cmp(&x) != Ordering::Less {
                return x;
            }
            x = y;
        }
    }

    /// `mpz_remove(_, _, 10)`: divide out all factors of ten.
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

    /// `mpz_com`: one's complement (`-self - 1`).
    pub fn com(&self) -> Mpz {
        let one = Mpz::from_u64(1);
        self.try_add(&one)
            .map(|mut m| {
                m.sign = -m.sign;
                m
            })
            .unwrap_or_else(|_| Mpz::new())
    }

    // ---- Decimal string I/O (no alloc) ----

    /// Parse a decimal string into an `Mpz`.
    ///
    /// Returns `ParseError::InvalidInput` for malformed strings and
    /// `ParseError::CapacityOverflow` if the value exceeds [`MPZ_MAX_LIMBS`] limbs.
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
        // Validate that all characters are ASCII digits
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

    /// Write the decimal representation into a `&mut [u8]` buffer.
    /// Returns the number of bytes written (the string is NOT null-terminated).
    ///
    /// The buffer should be at least 160 bytes for the maximum 512-bit value.
    /// If the buffer is too small, the value is truncated on the high end
    /// (buffer overflow is not possible; the method will simply not write
    /// the full number, but this is not recommended).
    pub fn write_decimal_buf(&self, buf: &mut [u8]) -> usize {
        if self.sign == 0 {
            if !buf.is_empty() {
                buf[0] = b'0';
                return 1;
            }
            return 0;
        }
        // Collect digits in a temporary stack buffer (right-to-left)
        let mut digits = [0u8; 160];
        let mut di = digits.len();
        // Divide by 10 repeatedly using 18-digit chunks
        let mut m = self.clone();
        loop {
            let mut qbuf = [0u64; MPZ_MAX_LIMBS];
            let (qlen, rem) =
                Self::mag_divmod_u64_len(&m.mag[..m.len], 1_000_000_000_000_000_000, &mut qbuf);
            // Write rem (up to 18 zero-padded digits)
            let mut rem_d = rem;
            if qlen == 0 {
                // Last chunk — no zero-padding
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
                // Zero-padded to 18 digits
                for _ in 0..18 {
                    di -= 1;
                    digits[di] = b'0' + (rem_d % 10) as u8;
                    rem_d /= 10;
                }
            }
            // Update m
            m.mag[..qlen].copy_from_slice(&qbuf[..qlen]);
            m.len = qlen;
        }
        // Write sign
        let mut pos = 0;
        if self.sign < 0 {
            if buf.is_empty() {
                return 0;
            }
            buf[0] = b'-';
            pos = 1;
        }
        // Now write from digits[di..] to buf
        let src = &digits[di..];
        let avail = buf.len() - pos;
        let n = src.len().min(avail);
        buf[pos..pos + n].copy_from_slice(&src[..n]);
        pos + n
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: format an Mpz as a string using the stack buffer.
    fn s(m: &Mpz) -> alloc::string::String {
        let mut buf = [0u8; 192];
        let len = m.write_decimal_buf(&mut buf);
        core::str::from_utf8(&buf[..len]).unwrap().into()
    }

    // The test module needs alloc for String comparison in tests.
    // But the library itself is no_alloc.
    extern crate alloc;

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
    fn division_truncating() {
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
        assert_eq!(
            Mpz::try_ui_pow_ui(10, 20).unwrap(),
            Mpz::from_decimal_str("100000000000000000000").unwrap()
        );
    }

    #[test]
    fn sizeinbase_and_str_roundtrip() {
        assert_eq!(Mpz::from_u64(255).sizeinbase2(), 8);
        assert_eq!(Mpz::new().sizeinbase2(), 1);
        for v in ["0", "-1", "42", "-1000000000000000000000000000000001"] {
            assert_eq!(s(&Mpz::from_decimal_str(v).unwrap()), v);
        }
    }

    #[test]
    fn sll_ull_roundtrip() {
        for &v in &[0i64, 1, -1, 42, -42, i64::MAX, i64::MIN + 1] {
            assert_eq!(Mpz::mpz_set_sll(v).mpz_get_sll(), v, "sll round-trip {v}");
        }
        assert_eq!(Mpz::mpz_set_sll(0).mpz_get_ull(), 0);
        assert_eq!(Mpz::mpz_set_sll(123456789).mpz_get_ull(), 123456789);
        assert_eq!(Mpz::mpz_set_sll(-1000).mpz_get_ull(), 1000);
    }

    #[test]
    fn zero_operations() {
        let mut z = Mpz::new();
        assert_eq!(s(&z), "0");
        assert_eq!(z.sgn(), 0);
        assert_eq!(z.size(), 0);
        assert!(z.fits_ulong());
        assert_eq!(z.get_ui(), 0);
        assert_eq!(z.get_si(), 0);
        assert_eq!(z.mpz_get_ull(), 0);
        assert_eq!(z.mpz_get_sll(), 0);
        assert_eq!(z.cmp(&z), Ordering::Equal);
        assert_eq!(z.cmpabs(&z), Ordering::Equal);
        assert_eq!(z.to_i128(), Some(0));
        assert_eq!(z.remove_pow10(), 0);
        assert_eq!(z.try_add(&z).unwrap(), z);
        assert_eq!(z.try_sub(&z).unwrap(), z);
        assert_eq!(z.try_mul(&z).unwrap(), z);
        assert!(z.divisible_ui(7));
        assert_eq!(z.tdiv_ui(7), 0);
        assert_eq!(z.tdiv_q_ui(7), z);
        assert_eq!(z.try_mul_2exp(100).unwrap(), z);
        assert_eq!(z.fdiv_q_2exp(100), z);
        assert_eq!(z.fdiv_r_2exp(100), z);
        assert_eq!(z.isqrt(), z);
        assert_eq!(z.com(), Mpz::from_i64(-1));
    }

    #[test]
    fn neg_abs() {
        let a = Mpz::from_i64(-123);
        let mut b = a.clone();
        b.abs();
        assert_eq!(s(&b), "123");
        b.neg();
        assert_eq!(s(&b), "-123");
        b.neg();
        assert_eq!(s(&b), "123");
        let mut c = Mpz::from_u64(42);
        c.abs();
        assert_eq!(s(&c), "42");
    }

    #[test]
    fn partial_ord_and_ord() {
        let a = Mpz::from_i64(-100);
        let b = Mpz::from_i64(0);
        let c = Mpz::from_i64(50);
        let d = Mpz::from_i64(50);
        assert!(a < b);
        assert!(b < c);
        assert!(c > a);
        assert!(c >= d);
        assert!(c <= d);
        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(b.cmp(&a), Ordering::Greater);
        assert_eq!(c.cmp(&d), Ordering::Equal);
    }

    #[test]
    fn conversions() {
        for v in &[0u64, 1, 255, u64::MAX / 2, u64::MAX] {
            let m = Mpz::from_u64(*v);
            assert_eq!(m.get_ui(), *v);
            if *v <= i64::MAX as u64 {
                assert_eq!(m.get_si(), *v as i64);
            }
            assert!(m.fits_ulong());
        }
        for v in &[0i64, 1, -1, i64::MAX, i64::MIN, 1234567890123456789] {
            let m = Mpz::from_i64(*v);
            assert_eq!(m.to_i128(), Some(*v as i128));
        }
        let big = 0xdeadbeef_cafebabe_12345678_9abcdef0u128;
        let m = Mpz::from_u128(big);
        assert_eq!(m.len, 2);
        assert_eq!(m.to_i128(), None);
        assert_eq!(Mpz::from_i128(i128::MAX).to_i128(), Some(i128::MAX));
        assert_eq!(Mpz::from_i128(i128::MIN).to_i128(), Some(i128::MIN));
    }

    #[test]
    fn display_format() {
        let cases = [
            ("0", "0"),
            ("-1", "-1"),
            ("12345", "12345"),
            ("-99999999999999999999", "-99999999999999999999"),
        ];
        for (src, expected) in &cases {
            let m = Mpz::from_decimal_str(src).unwrap();
            assert_eq!(alloc::format!("{m}"), *expected);
        }
    }

    #[test]
    fn large_multiplication() {
        let a = Mpz::from_decimal_str("123456789012345678901234567890").unwrap();
        let b = Mpz::from_decimal_str("987654321098765432109876543210").unwrap();
        assert_eq!(
            s(&a.try_mul(&b).unwrap()),
            "121932631137021795226185032733622923332237463801111263526900"
        );
        let ten_30 = Mpz::try_ui_pow_ui(10, 30).unwrap();
        let sq = ten_30.try_mul(&ten_30).unwrap();
        let expected = alloc::format!(
            "1{}",
            alloc::vec!['0'; 60]
                .iter()
                .collect::<alloc::string::String>()
        );
        assert_eq!(s(&sq), expected);
    }

    #[test]
    fn power_tests() {
        let p = Mpz::try_ui_pow_ui(2, 100).unwrap();
        assert_eq!(s(&p), "1267650600228229401496703205376");
        let p = Mpz::try_ui_pow_ui(3, 10).unwrap();
        assert_eq!(p, Mpz::from_u64(59049));
        let base = Mpz::from_u64(5);
        let p = base.try_pow_ui(15).unwrap();
        assert_eq!(s(&p), "30517578125");
        assert_eq!(Mpz::try_ui_pow_ui(0, 0).unwrap(), Mpz::from_u64(1));
    }

    #[test]
    fn division_edge_cases() {
        let a = Mpz::from_decimal_str("-12345678901234567890").unwrap();
        let one = Mpz::from_u64(1);
        let (q, r) = a.tdiv_qr(&one);
        assert_eq!(q, a);
        assert_eq!(s(&r), "0");
        let (q, r) = a.tdiv_qr(&a);
        assert_eq!(s(&q), "1");
        assert_eq!(s(&r), "0");
        let (q, r) = a.tdiv_qr(&Mpz::from_u64(2));
        assert_eq!(s(&q), "-6172839450617283945");
        assert_eq!(s(&r), "0");
        assert_eq!(s(&Mpz::from_u64(100).tdiv_q_ui(3)), "33");
        assert_eq!(Mpz::from_u64(100).tdiv_ui(3), 1);
        assert!(!Mpz::from_u64(100).divisible_ui(3));
        assert!(Mpz::from_u64(100).divisible_ui(5));
        assert!(Mpz::new().divisible_ui(42));
    }

    #[test]
    fn isqrt_tests() {
        assert_eq!(Mpz::from_u64(0).isqrt(), Mpz::from_u64(0));
        assert_eq!(Mpz::from_u64(1).isqrt(), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(4).isqrt(), Mpz::from_u64(2));
        assert_eq!(Mpz::from_u64(9).isqrt(), Mpz::from_u64(3));
        assert_eq!(Mpz::from_u64(144).isqrt(), Mpz::from_u64(12));
        assert_eq!(Mpz::from_u64(2).isqrt(), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(8).isqrt(), Mpz::from_u64(2));
        assert_eq!(Mpz::from_u64(15).isqrt(), Mpz::from_u64(3));
        assert_eq!(Mpz::from_u64(99).isqrt(), Mpz::from_u64(9));
        assert_eq!(Mpz::from_i64(-1).isqrt(), Mpz::new());
        let big = Mpz::from_decimal_str("15241578750190521027815549000000000000000").unwrap();
        let r = big.isqrt();
        let r_plus_1 = r.try_add(&Mpz::from_u64(1)).unwrap();
        assert!(r.try_mul(&r).unwrap().cmp(&big) != Ordering::Greater);
        assert!(r_plus_1.try_mul(&r_plus_1).unwrap().cmp(&big) == Ordering::Greater);
    }

    #[test]
    fn com_tests() {
        assert_eq!(Mpz::new().com(), Mpz::from_i64(-1));
        assert_eq!(Mpz::from_u64(5).com(), Mpz::from_i64(-6));
        assert_eq!(Mpz::from_i64(-5).com(), Mpz::from_i64(4));
    }

    #[test]
    fn parse_edge_cases() {
        assert_eq!(s(&Mpz::from_decimal_str("+42").unwrap()), "42");
        assert_eq!(s(&Mpz::from_decimal_str("  42  ").unwrap()), "42");
        assert_eq!(s(&Mpz::from_decimal_str("1").unwrap()), "1");
        assert_eq!(
            s(&Mpz::from_decimal_str("123456789012345678").unwrap()),
            "123456789012345678"
        );
        assert_eq!(
            s(&Mpz::from_decimal_str("1234567890123456789").unwrap()),
            "1234567890123456789"
        );
        let s36 = "123456789012345678901234567890123456";
        assert_eq!(s(&Mpz::from_decimal_str(s36).unwrap()), s36);
        // Parse errors
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
    fn arithmetic_with_carry_and_borrow() {
        let a = Mpz::from_u128(u64::MAX as u128);
        let b = Mpz::from_u128(1);
        let c = a.try_add(&b).unwrap();
        assert_eq!(c, Mpz::from_u128((u64::MAX as u128) + 1));
        assert_eq!(c.len, 2);
        let big = Mpz::from_u128(1u128 << 64);
        let one = Mpz::from_u64(1);
        let d = big.try_sub(&one).unwrap();
        assert_eq!(d, Mpz::from_u128(u64::MAX as u128));
        assert_eq!(d.len, 1);
        let e = Mpz::from_u128(u128::MAX);
        let f = Mpz::from_u128(u128::MAX);
        let g = e.try_add(&f).unwrap();
        assert_eq!(g.len, 3);
        assert_eq!(s(&g), "680564733841876926926749214863536422910");
    }

    #[test]
    fn shift_operations() {
        assert_eq!(
            Mpz::from_u64(42).try_mul_2exp(0).unwrap(),
            Mpz::from_u64(42)
        );
        assert_eq!(
            Mpz::from_u64(1).try_mul_2exp(64).unwrap(),
            Mpz::from_u128(1u128 << 64)
        );
        assert_eq!(
            Mpz::from_u64(1).try_mul_2exp(63).unwrap(),
            Mpz::from_u64(1u64 << 63)
        );
        let x = Mpz::from_u128((1u128 << 100) - 1);
        assert_eq!(x.fdiv_q_2exp(99), Mpz::from_u64(1));
        assert_eq!(x.fdiv_r_2exp(99), Mpz::from_u128((1u128 << 99) - 1));
        assert_eq!(Mpz::from_u64(1).fdiv_q_2exp(128), Mpz::new());
    }

    #[test]
    fn from_i64_negative_boundaries() {
        assert_eq!(Mpz::from_i64(i64::MIN).to_i128(), Some(i64::MIN as i128));
        assert_eq!(Mpz::from_i64(i64::MIN).get_si(), i64::MIN);
        assert_eq!(Mpz::from_i64(i64::MIN).sgn(), -1);
        assert!(!Mpz::from_i64(i64::MIN).fits_ulong());
    }

    #[test]
    fn cmp_and_cmpabs() {
        let a = Mpz::from_i64(-50);
        let b = Mpz::from_i64(30);
        let c = Mpz::from_i64(-30);
        assert_eq!(a.cmpabs(&b), Ordering::Greater);
        assert_eq!(b.cmpabs(&c), Ordering::Equal);
        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(c.cmp(&b), Ordering::Less);
        assert_eq!(a.cmp(&c), Ordering::Less);
    }

    #[test]
    fn set_ull_and_mpz_set_sll_consistency() {
        for v in &[0u64, 1, 42, u64::MAX] {
            assert_eq!(Mpz::set_ull(*v), Mpz::from_u64(*v));
        }
        for v in &[0i64, 1, -1, 42, -42, i64::MAX, i64::MIN + 1] {
            assert_eq!(Mpz::mpz_set_sll(*v), Mpz::from_i64(*v));
        }
    }

    #[test]
    fn sizeinbase2() {
        assert_eq!(Mpz::from_u64(1u64 << 63).sizeinbase2(), 64);
        assert_eq!(Mpz::from_u64(u64::MAX).sizeinbase2(), 64);
        let two_64 = Mpz::from_u128(1u128 << 64);
        assert_eq!(two_64.sizeinbase2(), 65);
        assert_eq!(Mpz::new().sizeinbase2(), 1);
    }

    #[test]
    fn remove_pow10_edge_cases() {
        let mut x = Mpz::from_decimal_str("12345").unwrap();
        assert_eq!(x.remove_pow10(), 0);
        assert_eq!(s(&x), "12345");
        let mut y = Mpz::from_u64(0);
        assert_eq!(y.remove_pow10(), 0);
        let mut z = Mpz::from_i64(-3000);
        assert_eq!(z.remove_pow10(), 3);
        assert_eq!(s(&z), "-3");
    }

    #[test]
    fn to_i128_limits() {
        assert_eq!(Mpz::from_i128(i128::MAX).to_i128(), Some(i128::MAX));
        assert_eq!(Mpz::from_i128(i128::MIN).to_i128(), Some(i128::MIN));
        let big3 = Mpz::from_decimal_str("340282366920938463463374607431768211456").unwrap(); // 2^128
        assert_eq!(big3.to_i128(), None);
        let neg_big3 = Mpz::from_decimal_str("-340282366920938463463374607431768211456").unwrap();
        assert_eq!(neg_big3.to_i128(), None);
    }

    #[test]
    fn capacity_error_on_overflow() {
        // Mpz::try_ui_pow_ui(2, 512) needs 513 bits = 9 limbs, which exceeds MPZ_MAX_LIMBS=8
        assert_eq!(Mpz::try_ui_pow_ui(2, 512), Err(CapacityError));
        // try_mul_2exp on a value that would exceed the limit
        let big = Mpz::from_u128(u128::MAX);
        assert_eq!(big.try_mul_2exp(400), Err(CapacityError));
    }

    #[test]
    fn decimal_string_roundtrip_extreme() {
        let vals = [
            "0",
            "1",
            "-1",
            "9",
            "10",
            "-10",
            "999999999999999999",  // 18 digits
            "-999999999999999999", // 18 digits
            "1000000000000000000", // 19 digits
            "123456789012345678901234567890123456789012345678901234567890", // 60 digits
        ];
        for v in &vals {
            let m = Mpz::from_decimal_str(v).unwrap();
            let back = s(&m);
            assert_eq!(&back, v, "round-trip failed for {v}");
        }
    }

    #[test]
    fn max_limb_capacity() {
        // Create a value that uses all 8 limbs: 2^(8*64) - 1
        let mut mpz = Mpz::new();
        for i in 0..8 {
            mpz.mag[i] = u64::MAX;
        }
        mpz.len = 8;
        mpz.sign = 1;
        assert_eq!(mpz.len, 8);
        // Adding 1 would need a 9th limb and should error
        let one = Mpz::from_u64(1);
        assert_eq!(mpz.try_add(&one), Err(CapacityError));
    }
}
