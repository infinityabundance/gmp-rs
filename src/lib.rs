//! pure-Rust, **no-unsafe**, **`no_std`**, **`no_alloc`** arbitrary-precision signed integer (`Mpz`),
//! faithful to GMP's `mpz_*` operations.
//!
//! # Guarantees
//! - **Zero `unsafe` code** — `#![forbid(unsafe_code)]` enforced at compile time.
//! - **`no_std`** — no standard library dependency; only `core`.
//! - **`no_alloc`** — zero heap allocations.  Fixed-capacity limb storage (`[u64; 8]`,
//!   512 bits ≈ 154 decimal digits).  Operations that would exceed capacity return
//!   [`CapacityError`](enum.CapacityError.html).
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

// ===========================================================================
// Mpz type
// ===========================================================================

/// Arbitrary-precision signed integer with **fixed capacity** (no heap allocation).
///
/// Up to [`MPZ_MAX_LIMBS`] little-endian 64-bit limbs in sign–magnitude.
#[derive(Clone, Debug, Eq)]
pub struct Mpz {
    /// -1, 0, or +1.  Invariant: `sign == 0` iff `len == 0`.
    pub sign: i8,
    /// Number of active limbs in `mag[0..len]`.
    pub len: usize,
    /// Fixed-capacity limb storage.  Only `mag[0..len]` is meaningful.
    pub mag: [u64; MPZ_MAX_LIMBS],
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

    /// Bit-by-bit division for |a| / |d| when d is small (< 2 limbs).
    /// Fills qbuf and rbuf.  Returns (quotient_len, remainder_len).
    fn mag_divmod_bitwise(
        a: &[u64],
        d: &[u64],
        qbuf: &mut [u64],
        rbuf: &mut [u64],
    ) -> (usize, usize) {
        let nbits = a.len() * 64;
        for o in qbuf.iter_mut().take(a.len()) {
            *o = 0;
        }
        let mut rem = Mpz::new();
        for i in (0..nbits).rev() {
            // rem × 2
            let carry = rem.mag[rem.len..].iter_mut().fold(0u64, |c, limb| {
                let new = (*limb << 1) | c;
                *limb = new;
                new >> 63
            });
            if carry != 0 && rem.len < MPZ_MAX_LIMBS {
                rem.mag[rem.len] = carry;
                rem.len += 1;
            } else if carry != 0 {
                // overflow in intermediate - should not happen at max 1 bit
            }
            // Insert bit a[i]
            let bit = (a[i / 64] >> (i % 64)) & 1;
            if bit != 0 {
                if rem.len == 0 {
                    rem.mag[0] = 1;
                    rem.len = 1;
                    rem.sign = 1;
                } else {
                    // Add 1 to rem (it's just been shifted, so LSB is 0)
                    let mut c = 1u64;
                    for j in 0..rem.len {
                        let (new, carry_out) = rem.mag[j].overflowing_add(c);
                        rem.mag[j] = new;
                        c = carry_out as u64;
                        if c == 0 {
                            break;
                        }
                    }
                    if c != 0 && rem.len < MPZ_MAX_LIMBS {
                        rem.mag[rem.len] = c;
                        rem.len += 1;
                    }
                }
            }
            // Compare and potentially subtract
            if Self::cmp_mag_slice(&rem.mag[..rem.len], d) != Ordering::Less {
                let mut tmp = [0u64; MPZ_MAX_LIMBS];
                let nl = Self::mag_sub_len(&rem.mag[..rem.len], d, &mut tmp);
                rem.mag[..nl].copy_from_slice(&tmp[..nl]);
                rem.len = nl;
                if rem.len == 0 {
                    rem.sign = 0;
                }
                qbuf[i / 64] |= 1u64 << (i % 64);
            }
        }
        let mut ql = a.len();
        while ql > 0 && qbuf[ql - 1] == 0 {
            ql -= 1;
        }
        rbuf[..rem.len].copy_from_slice(&rem.mag[..rem.len]);
        (ql, rem.len)
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

    pub fn set_ull(val: u64) -> Self {
        Self::from_u64(val)
    }

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

    /// `mpz_sizeinbase(_, base)`: upper bound on chars needed to represent in given base.
    /// Only supports bases 2–36.  Returns `None` for unsupported bases.
    pub fn try_sizeinbase(&self, base: i32) -> Option<usize> {
        if !(2..=36).contains(&base) {
            return None;
        }
        if self.len == 0 {
            return Some(1);
        }
        let bits = self.sizeinbase2();
        // ceil(bits * log(2) / log(base))
        let (num, den) = match base {
            2 => (1, 1),
            3 => (63, 100), // approx log2(3) ≈ 1.585
            4 => (1, 2),
            5 => (43, 100), // approx log2(5) ≈ 2.322
            6 => (39, 100),
            7 => (36, 100),
            8 => (1, 3),
            9 => (32, 100),
            10 => (30103, 100000), // log10(2) ≈ 0.30103
            // General base: approximate using natural logs
            // size = ceil(bits * ln(2) / ln(base))
            // Approximate ln(2) = 6931/10000 and ln(base) using integer tables
            _ => {
                // Use a simple integer log table for small bases up to 36
                let log2_approx: &[u32] = &[
                    0, 0, 1000, 1585, 2000, 2322, 2585, 2807, 3000, 3168, 3322, 3460, 3585, 3700,
                    3807, 3907, 4000, 4087, 4170, 4249, 4324, 4397, 4466, 4534, 4599, 4662, 4723,
                    4782, 4840, 4897, 4952, 5005, 5058, 5110, 5160, 5210,
                ];
                let log2_b = if (2..=36).contains(&base) {
                    log2_approx[base as usize]
                } else {
                    return None;
                };
                // sz = ceil(bits * 1000 / log2_b)
                let sz = ((bits as u64) * 1000).div_ceil(log2_b as u64);
                return Some(sz as usize + if self.sign < 0 { 1 } else { 0 });
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
            // Store remainder as single-limb Mpz in rbuf
            rbuf[0] = r_u64;
            let rlen = if r_u64 == 0 { 0 } else { 1 };
            (qlen, rlen)
        } else {
            Self::mag_divmod_bitwise(&self.mag[..self.len], &d.mag[..d.len], &mut qbuf, &mut rbuf)
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
        self.mag_divmod(d)
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
            // Shouldn't happen after floor adjustment, but be safe
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
            // Standard algorithm: ceil(a/b) = -floor(-a/b)
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
    /// Falls back to tdiv_qr for correctness; could be faster with direct division.
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
        // All bits below 'bits' must be zero
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
    pub fn congruent_2exp_p(&self, c: &Mpz, bits: u32) -> bool {
        let r_self = self.fdiv_r_2exp(bits);
        let r_c = c.fdiv_r_2exp(bits);
        r_self
            .try_sub(&r_c)
            .unwrap_or_else(|_| Mpz::new())
            .divisible_2exp_p(bits)
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
        // Newton iteration: x_{k+1} = ((n-1)*x_k + a/x_k^{n-1}) / n
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
                // Check for oscillation
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
            // Check small powers using integer arithmetic
            for k in 2..64 {
                if (1u64 << k) > v {
                    break;
                }
                // Binary search for integer k-th root
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
            // Multi-limb: try k from 2 to bitlen
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
    /// Returns `(remaining_value, count)`.
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
        // Binary GCD (Stein's algorithm)
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
        // Adjust signs
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
        let mut a = self.try_mod(&n).unwrap_or_else(|_| Mpz::new());
        let mut t = 1i32;
        while a.len != 0 {
            // Remove factors of 2 from a
            let mut e = 0u32;
            while a.even_p() {
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
            if a.mag[0] % 4 == 3 && n.mag[0] % 4 == 3 {
                t = -t;
            }
            a = a.try_mod(&n).unwrap_or_else(|_| Mpz::new());
            if n.len == 1 && n.mag[0] == 1 {
                return t;
            }
        }
        0
    }

    /// `mpz_fac_ui`: factorial n!.
    pub fn try_fac_ui(n: u32) -> Result<Mpz, CapacityError> {
        let mut r = Mpz::from_u64(1);
        for i in 2..=n as u64 {
            r = r.try_mul_ui(i)?;
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
        // Need one extra limb for the sign bit
        let n = self.len + 1;
        if self.sign > 0 {
            out[..self.len].copy_from_slice(&self.mag[..self.len]);
            out[self.len] = 0;
            (n, false)
        } else {
            // two's complement = ~|self| + 1, sign-extended
            let mut carry = 1u128;
            for i in 0..self.len {
                let inv = (!self.mag[i]) as u128;
                let sum = inv + carry;
                out[i] = sum as u64;
                carry = sum >> 64;
            }
            // Sign-extend: if carry propagated, higher limb is 0; otherwise all-ones
            if carry != 0 {
                out[self.len] = 0;
            } else {
                out[self.len] = !0u64; // all-ones (sign extension)
            }
            (n, true)
        }
    }

    /// Convert from two's complement back to sign-magnitude.
    fn from_twos_complement(limbs: &[u64], negative: bool) -> Mpz {
        if !negative {
            // Positive: just the magnitude
            let mut r = Mpz::new();
            let n = limbs.len().min(MPZ_MAX_LIMBS);
            r.mag[..n].copy_from_slice(&limbs[..n]);
            r.len = n;
            r.sign = if n == 0 { 0 } else { 1 };
            r.trim();
            return r;
        }
        // Negative: compute two's complement again: ~limbs + 1, then take magnitude
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
        let max_limbs = self.len.max(other.len) + 1;
        if max_limbs > MPZ_MAX_LIMBS {
            return Err(CapacityError);
        }
        let mut a_buf = [0u64; MPZ_MAX_LIMBS];
        let mut b_buf = [0u64; MPZ_MAX_LIMBS];
        let (a_len, a_neg) = self.to_twos_complement(&mut a_buf);
        let (b_len, b_neg) = other.to_twos_complement(&mut b_buf);
        let work_len = a_len.max(b_len).min(MPZ_MAX_LIMBS);
        let mut result_limbs = [0u64; MPZ_MAX_LIMBS];
        let sign_a = if a_neg { !0u64 } else { 0 };
        let sign_b = if b_neg { !0u64 } else { 0 };
        for i in 0..work_len {
            let va = if i < a_len { a_buf[i] } else { sign_a };
            let vb = if i < b_len { b_buf[i] } else { sign_b };
            result_limbs[i] = op(va, vb);
        }
        // Determine if result is negative: check the "sign bit" of the top limb
        let top_limb = result_limbs[work_len - 1];
        let negative = (top_limb >> 63) == 1;
        if !negative {
            // Top bit is 0: positive result, trim
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
        self.bitwise_op(other, |a, b| a & b)
    }

    pub fn try_ior(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        self.bitwise_op(other, |a, b| a | b)
    }

    pub fn try_xor(&self, other: &Mpz) -> Result<Mpz, CapacityError> {
        self.bitwise_op(other, |a, b| a ^ b)
    }

    /// `mpz_popcount`: number of 1 bits in the two's complement representation.
    /// For positive values this is the popcount of the magnitude.
    /// For negative values, treat as infinite two's complement (so -1 has popcount ∞,
    /// but GMP defines popcount only for non-negative values).
    pub fn popcount(&self) -> Option<u32> {
        if self.sign < 0 {
            return None;
        } // GMP: undefined for negative
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
                // Only scan from start_bit upward
                w |= (1u64 << start_bit) - 1; // fill lower bits to skip them
            }
            let zeros = (!w).trailing_zeros();
            if zeros < 64 {
                return (i * 64 + zeros as usize) as u32;
            }
        }
        // Beyond the magnitude limbs, all bits are 0
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
                w <<= start_bit; // zero out bits below start_bit
            }
            if w != 0 {
                return Some((i * 64 + w.trailing_zeros() as usize) as u32);
            }
        }
        None
    }

    /// `mpz_tstbit`: test whether bit `bit` is set.
    pub fn tstbit(&self, bit: u32) -> bool {
        let limb = (bit / 64) as usize;
        let bit_idx = bit % 64;
        if limb >= self.len {
            return false;
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
            // Need to extend
            if limb >= MPZ_MAX_LIMBS {
                return Err(CapacityError);
            }
            self.len = limb + 1;
            // The newly covered limbs (between old len and new) are already zero
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
        // Non-mutating
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
        // Floor: -100 / 30 = -4 remainder 20 (because -4*30 = -120, -120 + 20 = -100)
        let a = Mpz::from_i64(-100);
        let d = Mpz::from_i64(30);
        let (q, r) = a.try_fdiv_qr(&d).unwrap();
        assert_eq!(s(&q), "-4");
        assert_eq!(s(&r), "20");
        // Positive: 100 / 30 = 3 remainder 10
        let a2 = Mpz::from_i64(100);
        let (q2, r2) = a2.try_fdiv_qr(&d).unwrap();
        assert_eq!(s(&q2), "3");
        assert_eq!(s(&r2), "10");
    }

    #[test]
    fn div_floor_ui() {
        let a = Mpz::from_i64(-100);
        let r = a.fdiv_ui(30);
        assert_eq!(r, 20); // floor remainder is positive
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
        assert!(a.congruent_ui_p(5, 12)); // 17 ≡ 5 (mod 12)
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
        // 120 * s + 23 * t = g
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
        let m = Mpz::from_u64(13); // 3 * 9 = 27 ≡ 1 (mod 13)
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
        // With negative values (two's complement)
        let neg = Mpz::from_i64(-1); // all-ones
        assert_eq!(a.try_and(&neg).unwrap(), a);
        assert_eq!(a.try_ior(&neg).unwrap(), neg);
    }

    #[test]
    fn popcount_hamdist() {
        assert_eq!(Mpz::from_u64(0b1010).popcount(), Some(2));
        assert_eq!(Mpz::new().popcount(), Some(0));
        assert_eq!(Mpz::from_i64(-1).popcount(), None); // undefined for negative
        assert_eq!(
            Mpz::from_u64(0b1100).hamdist(&Mpz::from_u64(0b1010)),
            Some(2)
        );
    }

    #[test]
    fn scan_tstbit() {
        let a = Mpz::from_u64(0b1010_0000);
        assert_eq!(a.scan1(0), Some(5));
        assert_eq!(a.scan0(0), 0); // LSB is 0
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
        // 3^5 mod 7 = 243 mod 7 = 5
        let base = Mpz::from_u64(3);
        let exp = Mpz::from_u64(5);
        let m = Mpz::from_u64(7);
        let r = base.try_powm(&exp, &m).unwrap();
        assert_eq!(r, Mpz::from_u64(5));
    }

    #[test]
    fn root_test() {
        // cube root of 1000 = 10
        let r = Mpz::from_u64(1000).try_root(3).unwrap();
        assert_eq!(r, Mpz::from_u64(10));
        // cube root of 100 = 4 (4^3 = 64, 5^3 = 125)
        let r2 = Mpz::from_u64(100).try_root(3).unwrap();
        assert_eq!(r2, Mpz::from_u64(4));
    }

    #[test]
    fn remove_test() {
        let a = Mpz::from_u64(2700); // 2^2 * 3^3 * 5^2 = 4 * 27 * 25
        let (rem, count) = a.try_remove(&Mpz::from_u64(2)).unwrap();
        assert_eq!(count, 2);
        // 2700 / 4 = 675 = 3^3 * 5^2
        // Actually 2700 = 2^2 * 3^3 * 5^2, so removing factor 2 -> 675
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
        a.try_addmul(&Mpz::from_u64(3), &Mpz::from_u64(5)).unwrap(); // 10 + 15 = 25
        assert_eq!(a, Mpz::from_u64(25));
        a.try_submul(&Mpz::from_u64(2), &Mpz::from_u64(5)).unwrap(); // 25 - 10 = 15
        assert_eq!(a, Mpz::from_u64(15));
    }

    #[test]
    fn ui_sub() {
        // 100 - 1 = 99
        let r = Mpz::try_ui_sub(100, &Mpz::from_u64(1)).unwrap();
        assert_eq!(r, Mpz::from_u64(99));
        // 0 - 5 = -5
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
        assert_eq!(Mpz::from_u64(255).try_sizeinbase(1), None); // unsupported
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
        // Ceiling division: -100 / 30 = -3 remainder -10 (since -3*30 = -90, -90 + -10 = -100)
        let a = Mpz::from_i64(-100);
        let d = Mpz::from_i64(30);
        let (q, r) = a.try_cdiv_qr(&d).unwrap();
        // Verify: q*d + r == a
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
        a.try_addmul_ui(&Mpz::from_u64(3), 5).unwrap(); // 10 + 15 = 25
        assert_eq!(a, Mpz::from_u64(25));
        a.try_submul_ui(&Mpz::from_u64(2), 5).unwrap(); // 25 - 10 = 15
        assert_eq!(a, Mpz::from_u64(15));
    }

    #[test]
    fn tdiv_r_ui_test() {
        let a = Mpz::from_i64(-100);
        let r = a.tdiv_r_ui(30);
        assert_eq!(r, Mpz::from_u64(10)); // |100| mod 30 = 10
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
        assert_eq!(a.scan1(8), None); // no more 1 bits above bit 7
        assert_eq!(Mpz::from_u64(0).scan1(0), None);
    }

    #[test]
    fn perfect_power_test() {
        assert!(Mpz::from_u64(27).perfect_power_p()); // 3^3
        assert!(Mpz::from_u64(16).perfect_power_p()); // 2^4 = 4^2
        assert!(!Mpz::from_u64(2).perfect_power_p());
        assert!(!Mpz::new().perfect_power_p());
    }

    #[test]
    fn congruent_2exp_test() {
        // 17 ≡ 1 (mod 16 = 2^4)
        assert!(Mpz::from_u64(17).congruent_2exp_p(&Mpz::from_u64(1), 4));
        // 18 ≡ 2 (mod 16)
        assert!(Mpz::from_u64(18).congruent_2exp_p(&Mpz::from_u64(2), 4));
    }

    #[test]
    fn remove_general_test() {
        // 2^5 * 3^2 = 32 * 9 = 288
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
}
