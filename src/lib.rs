//! A pure-Rust, no-unsafe, no_std arbitrary-precision signed integer (`Mpz`), faithful to the GMP
//! `mpz_*` operations.
//!
//! This crate implements a subset of GMP's `mpz` surface — enough to serve as the integer
//! foundation for higher-level decimal / numeric ports — with **zero `unsafe` code** and **no
//! dependency on the standard library** (only the `alloc` crate is needed for `Vec<u64>` limbs).
//!
//! # Representation
//!
//! Sign–magnitude. `mag` is little-endian base-2⁶⁴ limbs with no trailing zero limb (so the zero
//! value is `mag.is_empty()` and `sign == 0`).
//!
//! # Capabilities
//!
//! - Construction from `u64`, `i64`, `u128`, `i128`, and decimal strings
//! - Output as decimal string or to `i128` (when it fits)
//! - Signed comparison, `add`, `sub`, `mul`, truncated division
//! - `pow_ui`, `ui_pow_ui`, `mul_2exp`, `fdiv_q/r_2exp`
//! - `isqrt` (integer floor square root)
//! - `remove_pow10` (factor out trailing decimal zeros)
//! - `com` (one's complement)

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

use alloc::string::String;
use alloc::string::ToString;
use alloc::vec::Vec;
use core::cmp::Ordering;
use core::fmt;

/// An arbitrary-precision signed integer (`mpz_t`).
#[derive(Clone, Debug, Eq)]
pub struct Mpz {
    /// -1, 0, or +1. Invariant: `sign == 0` iff `mag` is empty.
    sign: i8,
    /// magnitude, little-endian 64-bit limbs, no trailing zeros.
    mag: Vec<u64>,
}

// ---- Display ----

impl fmt::Display for Mpz {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_decimal_string())
    }
}

// ---- PartialEq ----

impl PartialEq for Mpz {
    fn eq(&self, other: &Self) -> bool {
        self.sign == other.sign && self.mag == other.mag
    }
}

// ---- Default ----

impl Default for Mpz {
    fn default() -> Self {
        Mpz::new()
    }
}

// ---- Operator overloads ----

impl core::ops::Add for &Mpz {
    type Output = Mpz;
    fn add(self, rhs: Self) -> Mpz {
        Mpz::add(self, rhs)
    }
}
impl core::ops::Add<&Mpz> for Mpz {
    type Output = Mpz;
    fn add(self, rhs: &Mpz) -> Mpz {
        Mpz::add(&self, rhs)
    }
}
impl core::ops::Add<Mpz> for &Mpz {
    type Output = Mpz;
    fn add(self, rhs: Mpz) -> Mpz {
        self.add(&rhs)
    }
}
impl core::ops::Add for Mpz {
    type Output = Mpz;
    fn add(self, rhs: Mpz) -> Mpz {
        self.add(&rhs)
    }
}

impl core::ops::Sub for &Mpz {
    type Output = Mpz;
    fn sub(self, rhs: Self) -> Mpz {
        Mpz::sub(self, rhs)
    }
}
impl core::ops::Sub<&Mpz> for Mpz {
    type Output = Mpz;
    fn sub(self, rhs: &Mpz) -> Mpz {
        Mpz::sub(&self, rhs)
    }
}
impl core::ops::Sub<Mpz> for &Mpz {
    type Output = Mpz;
    fn sub(self, rhs: Mpz) -> Mpz {
        self.sub(&rhs)
    }
}
impl core::ops::Sub for Mpz {
    type Output = Mpz;
    fn sub(self, rhs: Mpz) -> Mpz {
        self.sub(&rhs)
    }
}

impl core::ops::Neg for &Mpz {
    type Output = Mpz;
    fn neg(self) -> Mpz {
        let mut c = self.clone();
        Mpz::neg(&mut c);
        c
    }
}
impl core::ops::Neg for Mpz {
    type Output = Mpz;
    fn neg(self) -> Mpz {
        let mut c = self;
        Mpz::neg(&mut c);
        c
    }
}

// ---- PartialOrd / Ord ----

impl PartialOrd for Mpz {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(core::cmp::Ord::cmp(self, other))
    }
}
impl Ord for Mpz {
    fn cmp(&self, other: &Self) -> Ordering {
        Mpz::cmp(self, other)
    }
}

// ---- From conversions ----

impl From<u64> for Mpz {
    fn from(v: u64) -> Self {
        Mpz::from_u64(v)
    }
}
impl From<i64> for Mpz {
    fn from(v: i64) -> Self {
        Mpz::from_i64(v)
    }
}
impl From<u128> for Mpz {
    fn from(v: u128) -> Self {
        Mpz::from_u128(v)
    }
}
impl From<i128> for Mpz {
    fn from(v: i128) -> Self {
        Mpz::from_i128(v)
    }
}

// ---- Methods ----

impl Mpz {
    /// `mpz_init` / a freshly-`mpz_init2`-d value: zero.
    pub fn new() -> Self {
        Mpz {
            sign: 0,
            mag: Vec::new(),
        }
    }

    /// `mpz_set_ull (dest, val)` (numeric.c:182, COB_EXPERIMENTAL): set to an unsigned 64-bit host
    /// integer. GnuCOBOL writes `_mp_d[0] = val & GMP_NUMB_MASK` and `_mp_size = (val != 0)` (a single
    /// limb where `GMP_LIMB_BITS >= 64`, which holds here); the magnitude is one 64-bit limb of `val`.
    pub fn set_ull(val: u64) -> Self {
        if val == 0 {
            Mpz::new()
        } else {
            Mpz {
                sign: 1,
                mag: alloc::vec![val],
            }
        }
    }

    /// `mpz_set_sll (dest, val)` (numeric.c:198, COB_EXPERIMENTAL): set to a signed 64-bit host integer
    /// — magnitude `|val|` in one limb, `_mp_size` carrying the sign of `val`.
    pub fn mpz_set_sll(val: i64) -> Self {
        let mag = (val as i128).unsigned_abs() as u64;
        if mag == 0 {
            Mpz::new()
        } else {
            Mpz {
                sign: if val < 0 { -1 } else { 1 },
                mag: alloc::vec![mag],
            }
        }
    }

    /// `mpz_get_ull (src)` (numeric.c:216, COB_EXPERIMENTAL): the low 64-bit limb of the magnitude
    /// (`_mp_d[0]`), or 0 when the value is zero — wrapping past 64 bits exactly as the C does.
    pub fn mpz_get_ull(&self) -> u64 {
        self.mag.first().copied().unwrap_or(0)
    }

    /// `mpz_get_sll (src)` (numeric.c:236, COB_EXPERIMENTAL): reconstruct a signed 64-bit host integer
    /// from the low limb and the sign. Mirrors the C bit-for-bit: positive yields `vtmp & COB_MAX_LL`,
    /// negative yields `~((vtmp - 1) & COB_MAX_LL)`, with `COB_MAX_LL == i64::MAX`.
    pub fn mpz_get_sll(&self) -> i64 {
        if self.sign == 0 {
            return 0;
        }
        let vtmp = self.mag.first().copied().unwrap_or(0);
        if self.sign > 0 {
            (vtmp as i64) & i64::MAX
        } else {
            !(((vtmp as i64).wrapping_sub(1)) & i64::MAX)
        }
    }

    fn trim(mag: &mut Vec<u64>) {
        while mag.last() == Some(&0) {
            mag.pop();
        }
    }
    fn norm(sign: i8, mut mag: Vec<u64>) -> Self {
        Self::trim(&mut mag);
        if mag.is_empty() {
            Mpz { sign: 0, mag }
        } else {
            Mpz { sign, mag }
        }
    }

    // ---- set / get scalars (mpz_set_ui/si/ull/sll, mpz_get_ui/si/ull/sll) ----

    /// `mpz_set_ui`.
    pub fn set_ui(&mut self, v: u64) {
        *self = Self::from_u64(v);
    }
    /// `mpz_set_si`.
    pub fn set_si(&mut self, v: i64) {
        *self = Self::from_i64(v);
    }
    pub fn from_u64(v: u64) -> Self {
        if v == 0 {
            Mpz::new()
        } else {
            Mpz {
                sign: 1,
                mag: alloc::vec![v],
            }
        }
    }
    pub fn from_i64(v: i64) -> Self {
        if v == 0 {
            Mpz::new()
        } else if v > 0 {
            Mpz {
                sign: 1,
                mag: alloc::vec![v as u64],
            }
        } else {
            Mpz {
                sign: -1,
                mag: alloc::vec![(v as i128).unsigned_abs() as u64],
            }
        }
    }
    pub fn from_u128(v: u128) -> Self {
        Self::norm(
            if v == 0 { 0 } else { 1 },
            alloc::vec![v as u64, (v >> 64) as u64],
        )
    }
    pub fn from_i128(v: i128) -> Self {
        let neg = v < 0;
        let u = v.unsigned_abs();
        let m = Self::from_u128(u);
        if neg {
            Mpz { sign: -1, ..m }
        } else {
            m
        }
    }
    /// `mpz_get_ui`: the low 64 bits of the **absolute value** (GMP ignores the sign).
    pub fn get_ui(&self) -> u64 {
        self.mag.first().copied().unwrap_or(0)
    }
    /// `mpz_get_si`: low bits with sign, saturating like GMP's documented behavior for in-range values.
    pub fn get_si(&self) -> i64 {
        let lo = self.get_ui();
        if self.sign < 0 {
            (lo as i64).wrapping_neg()
        } else {
            lo as i64
        }
    }
    /// `mpz_fits_ulong_p` (treating ulong as u64): non-negative and a single limb.
    pub fn fits_ulong(&self) -> bool {
        self.sign >= 0 && self.mag.len() <= 1
    }

    // ---- sign / compare (mpz_sgn, mpz_cmp, mpz_cmpabs, mpz_abs, mpz_neg) ----

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
    fn cmp_mag(a: &[u64], b: &[u64]) -> Ordering {
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
    /// `mpz_cmpabs`: compare absolute values.
    pub fn cmpabs(&self, other: &Mpz) -> Ordering {
        Self::cmp_mag(&self.mag, &other.mag)
    }
    /// `mpz_cmp`: signed compare.
    #[allow(clippy::should_implement_trait)]
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

    // ---- magnitude add/sub helpers ----

    fn mag_add(a: &[u64], b: &[u64]) -> Vec<u64> {
        let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        let mut out = Vec::with_capacity(long.len() + 1);
        let mut carry = 0u128;
        for i in 0..long.len() {
            let mut cur = long[i] as u128 + carry;
            if i < short.len() {
                cur += short[i] as u128;
            }
            out.push(cur as u64);
            carry = cur >> 64;
        }
        if carry != 0 {
            out.push(carry as u64);
        }
        out
    }
    /// `a - b` for `a >= b` (magnitudes).
    fn mag_sub(a: &[u64], b: &[u64]) -> Vec<u64> {
        let mut out = Vec::with_capacity(a.len());
        let mut borrow = 0i128;
        for i in 0..a.len() {
            let bi = if i < b.len() { b[i] as i128 } else { 0 };
            let mut cur = a[i] as i128 - bi - borrow;
            if cur < 0 {
                cur += 1i128 << 64;
                borrow = 1;
            } else {
                borrow = 0;
            }
            out.push(cur as u64);
        }
        Self::trim(&mut out);
        out
    }
    fn add_signed(asign: i8, amag: &[u64], bsign: i8, bmag: &[u64]) -> Mpz {
        if asign == 0 {
            return Self::norm(bsign, bmag.to_vec());
        }
        if bsign == 0 {
            return Self::norm(asign, amag.to_vec());
        }
        if asign == bsign {
            Self::norm(asign, Self::mag_add(amag, bmag))
        } else {
            match Self::cmp_mag(amag, bmag) {
                Ordering::Equal => Mpz::new(),
                Ordering::Greater => Self::norm(asign, Self::mag_sub(amag, bmag)),
                Ordering::Less => Self::norm(bsign, Self::mag_sub(bmag, amag)),
            }
        }
    }

    /// `mpz_add`.
    pub fn add(&self, other: &Mpz) -> Mpz {
        Self::add_signed(self.sign, &self.mag, other.sign, &other.mag)
    }
    /// `mpz_sub`.
    pub fn sub(&self, other: &Mpz) -> Mpz {
        Self::add_signed(self.sign, &self.mag, -other.sign, &other.mag)
    }
    /// `mpz_add_ui`.
    pub fn add_ui(&self, v: u64) -> Mpz {
        self.add(&Mpz::from_u64(v))
    }
    /// `mpz_sub_ui`.
    pub fn sub_ui(&self, v: u64) -> Mpz {
        self.sub(&Mpz::from_u64(v))
    }

    // ---- multiply (mpz_mul, mpz_mul_ui, mpz_mul_2exp, mpz_ui_pow_ui) ----

    fn mag_mul(a: &[u64], b: &[u64]) -> Vec<u64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let mut out = alloc::vec![0u64; a.len() + b.len()];
        for (i, &ai) in a.iter().enumerate() {
            let mut carry = 0u128;
            for (j, &bj) in b.iter().enumerate() {
                let cur = out[i + j] as u128 + ai as u128 * bj as u128 + carry;
                out[i + j] = cur as u64;
                carry = cur >> 64;
            }
            out[i + b.len()] = out[i + b.len()].wrapping_add(carry as u64);
        }
        Self::trim(&mut out);
        out
    }
    /// `mpz_mul`.
    pub fn mul(&self, other: &Mpz) -> Mpz {
        Self::norm(self.sign * other.sign, Self::mag_mul(&self.mag, &other.mag))
    }
    /// `mpz_mul_ui`.
    pub fn mul_ui(&self, v: u64) -> Mpz {
        self.mul(&Mpz::from_u64(v))
    }
    /// `mpz_mul_2exp`: `self << bits`.
    pub fn mul_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let limb_shift = (bits / 64) as usize;
        let bit_shift = bits % 64;
        let mut m = alloc::vec![0u64; limb_shift];
        if bit_shift == 0 {
            m.extend_from_slice(&self.mag);
        } else {
            let mut carry = 0u64;
            for &l in &self.mag {
                m.push((l << bit_shift) | carry);
                carry = l >> (64 - bit_shift);
            }
            if carry != 0 {
                m.push(carry);
            }
        }
        Self::norm(self.sign, m)
    }
    /// `mpz_ui_pow_ui`: `base^exp`.
    pub fn ui_pow_ui(base: u64, exp: u32) -> Mpz {
        let mut r = Mpz::from_u64(1);
        let b = Mpz::from_u64(base);
        for _ in 0..exp {
            r = r.mul(&b);
        }
        r
    }

    /// `mpz_pow_ui (r, base, n)`: `self ^ n` by repeated squaring.
    pub fn pow_ui(&self, n: u32) -> Mpz {
        let mut r = Mpz::from_u64(1);
        let mut base = self.clone();
        let mut e = n;
        while e != 0 {
            if e & 1 == 1 {
                r = r.mul(&base);
            }
            e >>= 1;
            if e != 0 {
                base = base.mul(&base);
            }
        }
        r
    }

    // ---- divide (truncating: mpz_tdiv_q/r/q_ui/ui, fdiv_r/q_2exp) ----

    /// Divide magnitude by a single u64, returning `(quotient_mag, remainder)`.
    fn mag_divmod_u64(a: &[u64], d: u64) -> (Vec<u64>, u64) {
        let mut q = alloc::vec![0u64; a.len()];
        let mut rem: u128 = 0;
        for i in (0..a.len()).rev() {
            let cur = (rem << 64) | a[i] as u128;
            q[i] = (cur / d as u128) as u64;
            rem = cur % d as u128;
        }
        Self::trim(&mut q);
        (q, rem as u64)
    }
    /// `mpz_tdiv_q_ui`: truncated quotient by a u64. Returns the quotient.
    pub fn tdiv_q_ui(&self, d: u64) -> Mpz {
        let (q, _) = Self::mag_divmod_u64(&self.mag, d);
        Self::norm(self.sign, q)
    }
    /// `mpz_tdiv_ui`: the absolute remainder mod `d` (GMP returns |r|).
    pub fn tdiv_ui(&self, d: u64) -> u64 {
        Self::mag_divmod_u64(&self.mag, d).1
    }
    /// `mpz_divisible_ui_p`.
    pub fn divisible_ui(&self, d: u64) -> bool {
        self.sign == 0 || Self::mag_divmod_u64(&self.mag, d).1 == 0
    }
    /// Full truncated division `self / d`, returning `(quotient, remainder)` with the remainder
    /// taking the sign of the dividend (`mpz_tdiv_q` + `mpz_tdiv_r`).
    pub fn tdiv_qr(&self, d: &Mpz) -> (Mpz, Mpz) {
        debug_assert!(d.sign != 0);
        if Self::cmp_mag(&self.mag, &d.mag) == Ordering::Less {
            return (Mpz::new(), self.clone());
        }
        let (qmag, rmag) = Self::mag_divmod(&self.mag, &d.mag);
        let q = Self::norm(self.sign * d.sign, qmag);
        let r = Self::norm(self.sign, rmag);
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
    /// `mpz_fdiv_r_2exp`: the low `bits` bits (floor remainder by 2^bits). For non-negative values
    /// this is a bit mask; libcob only calls it on non-negative magnitudes here.
    pub fn fdiv_r_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let limbs = (bits / 64) as usize;
        let rem_bits = bits % 64;
        let mut m: Vec<u64> = self.mag.iter().take(limbs + 1).copied().collect();
        if rem_bits != 0 && m.len() > limbs {
            if let Some(top) = m.get_mut(limbs) {
                *top &= (1u64 << rem_bits) - 1;
            }
        } else if rem_bits == 0 {
            m.truncate(limbs);
        }
        Self::norm(self.sign, m)
    }
    /// `mpz_fdiv_q_2exp`: `self >> bits` (floor, but non-negative here).
    pub fn fdiv_q_2exp(&self, bits: u32) -> Mpz {
        if self.sign == 0 {
            return Mpz::new();
        }
        let limb_shift = (bits / 64) as usize;
        let bit_shift = bits % 64;
        if limb_shift >= self.mag.len() {
            return Mpz::new();
        }
        let mut m: Vec<u64> = self.mag[limb_shift..].to_vec();
        if bit_shift != 0 {
            let mut carry = 0u64;
            for l in m.iter_mut().rev() {
                let new = (*l >> bit_shift) | carry;
                carry = *l << (64 - bit_shift);
                *l = new;
            }
        }
        Self::norm(self.sign, m)
    }

    /// Schoolbook long division of magnitudes (binary), returning `(quotient, remainder)`.
    fn mag_divmod(a: &[u64], d: &[u64]) -> (Vec<u64>, Vec<u64>) {
        // bit-by-bit; d != 0. a >= d guaranteed by caller for the q!=0 case, but handle generally.
        let nbits = a.len() * 64;
        let mut q = alloc::vec![0u64; a.len()];
        let mut rem = Mpz::new();
        let dm = Mpz::norm(1, d.to_vec());
        for i in (0..nbits).rev() {
            // rem = (rem << 1) | bit_i(a)
            rem = rem.mul_2exp(1);
            let bit = (a[i / 64] >> (i % 64)) & 1;
            if bit != 0 {
                rem = rem.add(&Mpz::from_u64(1));
            }
            if Self::cmp_mag(&rem.mag, d) != Ordering::Less {
                rem = rem.sub(&dm);
                q[i / 64] |= 1u64 << (i % 64);
            }
        }
        Self::trim(&mut q);
        (q, rem.mag)
    }

    // ---- queries / misc ----

    /// `mpz_size`: number of limbs.
    pub fn size(&self) -> usize {
        self.mag.len()
    }
    /// `mpz_sizeinbase(_, 2)`: number of significant bits (1 for zero, like GMP).
    pub fn sizeinbase2(&self) -> usize {
        match self.mag.last() {
            None => 1,
            Some(&top) => (self.mag.len() - 1) * 64 + (64 - top.leading_zeros() as usize),
        }
    }
    /// `mpz_sqrt`: the floor of the integer square root, `floor(sqrt(self))`, for `self >= 0` (0 for a
    /// non-positive value). Newton iteration on integers — converges to the exact floor.
    pub fn isqrt(&self) -> Mpz {
        if self.sgn() <= 0 {
            return Mpz::new();
        }
        // initial over-estimate: 2^ceil(bits/2) >= sqrt(self)
        let mut x = Mpz::from_u64(1).mul_2exp(self.sizeinbase2().div_ceil(2) as u32);
        loop {
            // y = floor((x + floor(self/x)) / 2)
            let y = x.add(&self.tdiv_q(&x)).fdiv_q_2exp(1);
            if y.cmp(&x) != Ordering::Less {
                return x;
            }
            x = y;
        }
    }

    /// `mpz_remove(_, _, 10)`: divide out all factors of ten, returning the count removed.
    pub fn remove_pow10(&mut self) -> u32 {
        if self.sign == 0 {
            return 0;
        }
        let mut count = 0;
        loop {
            let (q, r) = Self::mag_divmod_u64(&self.mag, 10);
            if r != 0 {
                break;
            }
            self.mag = q;
            count += 1;
        }
        Self::trim(&mut self.mag);
        if self.mag.is_empty() {
            self.sign = 0;
        }
        count
    }
    /// `mpz_com`: one's complement (`-self - 1`).
    pub fn com(&self) -> Mpz {
        self.add(&Mpz::from_u64(1)).into_neg()
    }
    fn into_neg(mut self) -> Mpz {
        self.sign = -self.sign;
        self
    }

    /// The value as `i128` if it fits (≤ 2 limbs and in range), else `None`.
    pub fn to_i128(&self) -> Option<i128> {
        if self.mag.len() > 2 {
            return None;
        }
        let u = self.mag.first().copied().unwrap_or(0) as u128
            | ((self.mag.get(1).copied().unwrap_or(0) as u128) << 64);
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

    /// `mpz_get_str(_, 10, _)`: decimal string.
    pub fn to_decimal_string(&self) -> String {
        if self.sign == 0 {
            return "0".to_string();
        }
        let mut m = self.mag.clone();
        let mut chunks: Vec<u64> = Vec::new();
        while !m.is_empty() {
            let (q, r) = Self::mag_divmod_u64(&m, 1_000_000_000_000_000_000);
            chunks.push(r);
            m = q;
        }
        let mut s = String::new();
        if self.sign < 0 {
            s.push('-');
        }
        for (i, c) in chunks.iter().rev().enumerate() {
            if i == 0 {
                s.push_str(&c.to_string());
            } else {
                s.push_str(&alloc::format!("{c:018}"));
            }
        }
        s
    }
    /// `mpz_set_str(_, _, 10)`: parse a decimal string (optional leading sign).
    pub fn from_decimal_string(s: &str) -> Mpz {
        let s = s.trim();
        let (neg, digits) = match s.strip_prefix('-') {
            Some(d) => (true, d),
            None => (false, s.strip_prefix('+').unwrap_or(s)),
        };
        let mut r = Mpz::new();
        let bytes = digits.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let end = (i + 18).min(bytes.len());
            let mut chunk: u64 = 0;
            let mut p = 10u64.pow((end - i) as u32 - 1);
            for &b in &bytes[i..end] {
                if b.is_ascii_digit() {
                    chunk += (b - b'0') as u64 * p;
                    p /= 10;
                }
            }
            let scale = Mpz::from_u64(10u64.pow((end - i) as u32));
            r = r.mul(&scale).add(&Mpz::from_u64(chunk));
            i = end;
        }
        if neg && r.sign != 0 {
            r.sign = -1;
        }
        r
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(m: &Mpz) -> String {
        m.to_decimal_string()
    }

    #[test]
    fn basic_arith() {
        let a = Mpz::from_decimal_string("123456789012345678901234567890");
        let b = Mpz::from_decimal_string("987654321098765432109876543210");
        assert_eq!(s(&a.add(&b)), "1111111110111111111011111111100");
        assert_eq!(s(&b.sub(&a)), "864197532086419753208641975320");
        assert_eq!(
            s(&a.mul(&b)),
            "121932631137021795226185032733622923332237463801111263526900"
        );
        assert_eq!(a.cmp(&b), Ordering::Less);
        assert_eq!(a.cmpabs(&b), Ordering::Less);
    }

    #[test]
    fn signs_and_zero() {
        let a = Mpz::from_i64(-5);
        let b = Mpz::from_i64(5);
        assert_eq!(a.add(&b), Mpz::new());
        assert_eq!(a.sgn(), -1);
        assert_eq!(s(&a.mul(&b)), "-25");
        assert_eq!(a.get_ui(), 5); // |value| low bits
        assert_eq!(a.get_si(), -5);
    }

    #[test]
    fn division_truncating() {
        let a = Mpz::from_decimal_string("-1000000000000000000000007");
        let d = Mpz::from_u64(1000);
        let (q, r) = a.tdiv_qr(&d);
        assert_eq!(s(&q), "-1000000000000000000000"); // toward zero
        assert_eq!(s(&r), "-7"); // remainder takes dividend sign
        assert_eq!(a.tdiv_ui(1000), 7); // |remainder|
    }

    #[test]
    fn shifts_and_remove() {
        let mut x = Mpz::from_decimal_string("123000");
        assert_eq!(x.remove_pow10(), 3);
        assert_eq!(s(&x), "123");
        let y = Mpz::from_u64(1).mul_2exp(100);
        assert_eq!(s(&y), "1267650600228229401496703205376");
        assert_eq!(y.fdiv_q_2exp(100), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(0b1011).fdiv_r_2exp(2), Mpz::from_u64(0b11));
        assert_eq!(
            Mpz::ui_pow_ui(10, 20),
            Mpz::from_decimal_string("100000000000000000000")
        );
    }

    #[test]
    fn sizeinbase_and_str_roundtrip() {
        assert_eq!(Mpz::from_u64(255).sizeinbase2(), 8);
        assert_eq!(Mpz::new().sizeinbase2(), 1);
        for v in ["0", "-1", "42", "-1000000000000000000000000000000001"] {
            assert_eq!(s(&Mpz::from_decimal_string(v)), v);
        }
    }

    #[test]
    fn sll_ull_roundtrip() {
        // mpz_set_sll / mpz_get_sll: signed 64-bit round-trip incl. sign + zero.
        for &v in &[0i64, 1, -1, 42, -42, i64::MAX, i64::MIN + 1] {
            assert_eq!(Mpz::mpz_set_sll(v).mpz_get_sll(), v, "sll round-trip {v}");
        }
        // mpz_get_ull: low 64-bit magnitude limb.
        assert_eq!(Mpz::mpz_set_sll(0).mpz_get_ull(), 0);
        assert_eq!(Mpz::mpz_set_sll(123456789).mpz_get_ull(), 123456789);
        // get_ull ignores the sign (magnitude limb only).
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
        // zero +/- zero
        assert_eq!(z.add(&z), z);
        assert_eq!(z.sub(&z), z);
        assert_eq!(z.mul(&z), z);
        assert!(z.divisible_ui(7));
        assert_eq!(z.tdiv_ui(7), 0);
        assert_eq!(z.tdiv_q_ui(7), z);
        // zero shifts
        assert_eq!(z.mul_2exp(100), z);
        assert_eq!(z.fdiv_q_2exp(100), z);
        assert_eq!(z.fdiv_r_2exp(100), z);
        // zero sqrt
        assert_eq!(z.isqrt(), z);
        // zero one's complement
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
        // Neg trait
        assert_eq!(s(&(-&a)), "123");
        assert_eq!(s(&(-a)), "123");
        // positive abs is no-op
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
        // u64 -> Mpz -> back
        for v in &[0u64, 1, 255, u64::MAX / 2, u64::MAX] {
            let m = Mpz::from(*v);
            assert_eq!(Mpz::from_u64(*v), m);
            assert_eq!(m.get_ui(), *v);
            if *v <= i64::MAX as u64 {
                assert_eq!(m.get_si(), *v as i64);
            }
            assert!(m.fits_ulong());
        }
        // i64 -> Mpz -> i128 round-trip
        for v in &[0i64, 1, -1, i64::MAX, i64::MIN, 1234567890123456789] {
            let m = Mpz::from(*v);
            assert_eq!(m.to_i128(), Some(*v as i128));
        }
        // u128 -> Mpz -> string
        let big = 0xdeadbeef_cafebabe_12345678_9abcdef0u128;
        let m = Mpz::from(big);
        assert_eq!(m.size(), 2);
        assert_eq!(m.to_i128(), None); // doesn't fit in i128
                                       // i128
        let v = i128::MAX;
        let m = Mpz::from(v);
        assert_eq!(m.to_i128(), Some(v));
        let v = i128::MIN;
        let m = Mpz::from(v);
        assert_eq!(m.to_i128(), Some(v));
    }

    #[test]
    fn operator_overloads() {
        let a = Mpz::from_decimal_string("100000000000000000000");
        let b = Mpz::from_decimal_string("1");
        // Add by ref/mixed/owned combos
        assert_eq!(&a + &b, Mpz::from_decimal_string("100000000000000000001"));
        assert_eq!(
            &a + b.clone(),
            Mpz::from_decimal_string("100000000000000000001")
        );
        assert_eq!(
            a.clone() + &b,
            Mpz::from_decimal_string("100000000000000000001")
        );
        assert_eq!(
            a.clone() + b.clone(),
            Mpz::from_decimal_string("100000000000000000001")
        );
        // Sub
        assert_eq!(&a - &b, Mpz::from_decimal_string("99999999999999999999"));
        assert_eq!(
            &a - b.clone(),
            Mpz::from_decimal_string("99999999999999999999")
        );
        assert_eq!(
            a.clone() - &b,
            Mpz::from_decimal_string("99999999999999999999")
        );
        assert_eq!(
            a.clone() - b,
            Mpz::from_decimal_string("99999999999999999999")
        );
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
            let m = Mpz::from_decimal_string(src);
            assert_eq!(alloc::format!("{m}"), *expected);
        }
    }

    #[test]
    fn large_multiplication() {
        // Multiply two large numbers
        let a = Mpz::from_decimal_string("123456789012345678901234567890");
        let b = Mpz::from_decimal_string("987654321098765432109876543210");
        assert_eq!(
            s(&a.mul(&b)),
            "121932631137021795226185032733622923332237463801111263526900"
        );
        // squaring: (10^30)^2 = 10^60 = "1" + 60 zeros -> 61 digits
        let ten_30 = Mpz::ui_pow_ui(10, 30);
        let sq = ten_30.mul(&ten_30);
        let expected = "1".to_string() + &"0".repeat(60);
        assert_eq!(s(&sq), expected);
    }

    #[test]
    fn power_tests() {
        // 2^100
        let p = Mpz::ui_pow_ui(2, 100);
        assert_eq!(s(&p), "1267650600228229401496703205376");
        // 3^10 = 59049
        let p = Mpz::ui_pow_ui(3, 10);
        assert_eq!(p, Mpz::from_u64(59049));
        // pow_ui via repeated squaring
        let base = Mpz::from_u64(5);
        let p = base.pow_ui(15);
        assert_eq!(s(&p), "30517578125");
        // 0^0 = 1 (convention)
        assert_eq!(Mpz::ui_pow_ui(0, 0), Mpz::from_u64(1));
    }

    #[test]
    fn division_edge_cases() {
        // self / 1 = self
        let a = Mpz::from_decimal_string("-12345678901234567890");
        let one = Mpz::from_u64(1);
        let (q, r) = a.tdiv_qr(&one);
        assert_eq!(q, a);
        assert_eq!(s(&r), "0");
        // self / self = 1
        let (q, r) = a.tdiv_qr(&a);
        assert_eq!(s(&q), "1");
        assert_eq!(s(&r), "0");
        // negative dividend, positive divisor
        let (q, r) = a.tdiv_qr(&Mpz::from_u64(2));
        assert_eq!(s(&q), "-6172839450617283945");
        assert_eq!(s(&r), "0");
        // divide by small
        assert_eq!(s(&Mpz::from_u64(100).tdiv_q_ui(3)), "33");
        assert_eq!(Mpz::from_u64(100).tdiv_ui(3), 1);
        assert!(!Mpz::from_u64(100).divisible_ui(3));
        assert!(Mpz::from_u64(100).divisible_ui(5));
        assert!(Mpz::new().divisible_ui(42));
    }

    #[test]
    fn isqrt_tests() {
        // perfect squares
        assert_eq!(Mpz::from_u64(0).isqrt(), Mpz::from_u64(0));
        assert_eq!(Mpz::from_u64(1).isqrt(), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(4).isqrt(), Mpz::from_u64(2));
        assert_eq!(Mpz::from_u64(9).isqrt(), Mpz::from_u64(3));
        assert_eq!(Mpz::from_u64(144).isqrt(), Mpz::from_u64(12));
        // non-perfect squares (floor)
        assert_eq!(Mpz::from_u64(2).isqrt(), Mpz::from_u64(1));
        assert_eq!(Mpz::from_u64(8).isqrt(), Mpz::from_u64(2));
        assert_eq!(Mpz::from_u64(15).isqrt(), Mpz::from_u64(3));
        assert_eq!(Mpz::from_u64(99).isqrt(), Mpz::from_u64(9));
        // negative -> 0
        assert_eq!(Mpz::from_i64(-1).isqrt(), Mpz::new());
        // large
        let big = Mpz::from_decimal_string("15241578750190521027815549000000000000000");
        let r = big.isqrt();
        // 12345678901234567890^2 = 15241578750190521027815549000000000000000... let's check
        // self * self should be <= big, (self+1)^2 > big
        let r_plus_1 = r.add(&Mpz::from_u64(1));
        assert!(r.mul(&r).cmp(&big) != Ordering::Greater);
        assert!(r_plus_1.mul(&r_plus_1).cmp(&big) == Ordering::Greater);
    }

    #[test]
    fn com_tests() {
        assert_eq!(Mpz::new().com(), Mpz::from_i64(-1));
        assert_eq!(Mpz::from_u64(5).com(), Mpz::from_i64(-6));
        assert_eq!(Mpz::from_i64(-5).com(), Mpz::from_i64(4));
    }

    #[test]
    fn from_decimal_string_edge_cases() {
        // positive with explicit plus sign
        assert_eq!(s(&Mpz::from_decimal_string("+42")), "42");
        // leading/trailing whitespace
        assert_eq!(s(&Mpz::from_decimal_string("  42  ")), "42");
        // very small
        assert_eq!(s(&Mpz::from_decimal_string("1")), "1");
        // 18 digits (one chunk)
        let s18 = "123456789012345678";
        assert_eq!(s(&Mpz::from_decimal_string(s18)), s18);
        // 19 digits (two chunks)
        let s19 = "1234567890123456789";
        assert_eq!(s(&Mpz::from_decimal_string(s19)), s19);
        // 36 digits (two 18-digit chunks)
        let s36 = "123456789012345678901234567890123456";
        assert_eq!(s(&Mpz::from_decimal_string(s36)), s36);
        // negative large
        let neg = "-999999999999999999999999999999999999999999999999999999";
        assert_eq!(s(&Mpz::from_decimal_string(neg)), neg);
    }

    #[test]
    fn arithmetic_with_carry_and_borrow() {
        // Addition that triggers a carry across limbs
        let a = Mpz::from_u128(u64::MAX as u128);
        let b = Mpz::from_u128(1);
        let c = a.add(&b);
        assert_eq!(c, Mpz::from_u128((u64::MAX as u128) + 1));
        assert_eq!(c.size(), 2);
        // Subtraction with borrow
        let big = Mpz::from_u128(1u128 << 64);
        let one = Mpz::from_u64(1);
        let d = big.sub(&one);
        assert_eq!(d, Mpz::from_u128(u64::MAX as u128));
        assert_eq!(d.size(), 1);
        // 2-limb + 2-limb with carry to 3rd limb: u128::MAX + u128::MAX = 2^129 - 2
        let e = Mpz::from_u128(u128::MAX);
        let f = Mpz::from_u128(u128::MAX);
        let g = e.add(&f);
        assert_eq!(g.size(), 3);
        assert_eq!(s(&g), "680564733841876926926749214863536422910");
    }

    #[test]
    fn shift_operations() {
        // shift by zero is identity
        assert_eq!(Mpz::from_u64(42).mul_2exp(0), Mpz::from_u64(42));
        // shift by exact limb boundary
        assert_eq!(Mpz::from_u64(1).mul_2exp(64), Mpz::from_u128(1u128 << 64));
        // shift by non-limb boundary
        assert_eq!(Mpz::from_u64(1).mul_2exp(63), Mpz::from_u64(1u64 << 63));
        // finite field divide q_2exp
        let x = Mpz::from_u128((1u128 << 100) - 1);
        assert_eq!(x.fdiv_q_2exp(99), Mpz::from_u64(1));
        assert_eq!(x.fdiv_r_2exp(99), Mpz::from_u128((1u128 << 99) - 1));
        // divmod by power of 2 beyond size => zero
        assert_eq!(Mpz::from_u64(1).fdiv_q_2exp(128), Mpz::new());
    }

    #[test]
    fn from_i64_negative_boundaries() {
        assert_eq!(Mpz::from_i64(i64::MIN).to_i128(), Some(i64::MIN as i128));
        assert_eq!(Mpz::from_i64(i64::MIN).get_si(), i64::MIN);
        assert_eq!(Mpz::from_i64(i64::MIN).sgn(), -1);
        assert!(Mpz::from_i64(i64::MIN).fits_ulong() == false);
    }

    #[test]
    fn cmp_and_cmpabs() {
        let a = Mpz::from_i64(-50);
        let b = Mpz::from_i64(30);
        let c = Mpz::from_i64(-30);
        // cmpabs: |-50| > |30|
        assert_eq!(a.cmpabs(&b), Ordering::Greater);
        assert_eq!(b.cmpabs(&c), Ordering::Equal);
        // signed: -50 < 30
        assert_eq!(a.cmp(&b), Ordering::Less);
        // -30 < 30
        assert_eq!(c.cmp(&b), Ordering::Less);
        // -50 < -30
        assert_eq!(a.cmp(&c), Ordering::Less);
    }

    #[test]
    fn set_ull_and_mpz_set_sll_consistency() {
        for v in &[0u64, 1, 42, u64::MAX] {
            let m1 = Mpz::set_ull(*v);
            let m2 = Mpz::from_u64(*v);
            assert_eq!(m1, m2);
        }
        for v in &[0i64, 1, -1, 42, -42, i64::MAX, i64::MIN + 1] {
            let m1 = Mpz::mpz_set_sll(*v);
            let m2 = Mpz::from_i64(*v);
            assert_eq!(m1, m2);
        }
    }

    #[test]
    fn sizeinbase2() {
        // 2^63 (63 zero bits + 1 one bit = 64 bits)
        assert_eq!(Mpz::from_u64(1u64 << 63).sizeinbase2(), 64);
        // 2^64 - 1 (64 bits)
        assert_eq!(Mpz::from_u64(u64::MAX).sizeinbase2(), 64);
        // 2^64 => needs 2 limbs
        let two_64 = Mpz::from_u128(1u128 << 64);
        assert_eq!(two_64.sizeinbase2(), 65);
        // zero is 1 bit per GMP convention
        assert_eq!(Mpz::new().sizeinbase2(), 1);
    }

    #[test]
    fn remove_pow10_edge_cases() {
        // no trailing zeros
        let mut x = Mpz::from_decimal_string("12345");
        assert_eq!(x.remove_pow10(), 0);
        assert_eq!(s(&x), "12345");
        // remove all (value becomes 0)
        let mut y = Mpz::from_u64(0);
        assert_eq!(y.remove_pow10(), 0);
        // negative with trailing zeros
        let mut z = Mpz::from_i64(-3000);
        assert_eq!(z.remove_pow10(), 3);
        assert_eq!(s(&z), "-3");
    }

    #[test]
    fn to_i128_limits() {
        // fits: 2 limb positive
        let v = Mpz::from_i128(i128::MAX);
        assert_eq!(v.to_i128(), Some(i128::MAX));
        // fits: 2 limb negative
        let v = Mpz::from_i128(i128::MIN);
        assert_eq!(v.to_i128(), Some(i128::MIN));
        // too big (3 limbs)
        let big3 = Mpz::from_decimal_string("340282366920938463463374607431768211456"); // 2^128
        assert_eq!(big3.to_i128(), None);
        // negative too big
        let neg_big3 = Mpz::from_decimal_string("-340282366920938463463374607431768211456");
        assert_eq!(neg_big3.to_i128(), None);
    }
}
