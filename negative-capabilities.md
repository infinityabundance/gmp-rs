# gmp-rs: Forensic Negative-Capabilities Audit vs GNU GMP 6.3.0 `mpz_*`

> **Date:** 2026-07-20
> **gmp-rs version:** 0.1.0
> **GMP reference:** 6.3.0 (released 2023)
> **Constraint envelope:** `no_unsafe` ✓, `no_std` ✓, `no_alloc` ❌ (see §N.1)
>
> This document exhaustively catalogues every deviation, omission, semantic quirk,
> correctness edge-case, performance pathology, and architectural limitation of
> gmp-rs relative to the real GNU GMP `mpz_*` integer API. It is organised by
> GMP category, with a concluding synthesis and risk register.
>
> **Convention:** Entries are tagged with severity:
> - **[MISSING]** — function or variant absent from gmp-rs
> - **[SEMANTIC]** — different behaviour for the same-named operation
> - **[CORRECTNESS]** — provably wrong in a documented edge-case
> - **[PERF]** — algorithmically or asymptotically inferior
> - **[SAFETY]** — missing precondition check that could lead to UB / panic
> - **[API]** — ergonomic or structural difference that changes how callers must interact
> - **[INFO]** — informational, not a defect

---

## Table of Contents

1.  [Representation & Type System](#1-representation--type-system)
2.  [Initialization Functions (§3.3 GMP manual)](#2-initialization-functions)
3.  [Assignment Functions (§3.4)](#3-assignment-functions)
4.  [Combined Init + Assignment (§3.5)](#4-combined-init--assignment)
5.  [Conversion Functions (§3.6)](#5-conversion-functions)
6.  [Arithmetic Functions (§3.7)](#6-arithmetic-functions)
7.  [Division Functions (§3.8)](#7-division-functions)
8.  [Exponentiation Functions (§3.9)](#8-exponentiation-functions)
9.  [Root Extraction Functions (§3.10)](#9-root-extraction-functions)
10. [Number Theoretic Functions (§3.11)](#10-number-theoretic-functions)
11. [Comparison Functions (§3.12)](#11-comparison-functions)
12. [Logical & Bit Manipulation Functions (§3.13)](#12-logical--bit-manipulation-functions)
13. [I/O Functions (§3.14)](#13-io-functions)
14. [Random Number Functions (§3.15)](#14-random-number-functions)
15. [Integer Import/Export (§3.16)](#15-integer-importexport)
16. [Miscellaneous Functions (§3.17)](#16-miscellaneous-functions)
17. [Low-Level / Limb-Access Functions (§3.18)](#17-low-level--limb-access-functions)
18. [Synthesis & Risk Register](#18-synthesis--risk-register)
19. [The `no_alloc` Constraint](#19-the-no_alloc-constraint)

---

## 1. Representation & Type System

### 1.1 Limb width is hardcoded to 64 bits

| Aspect | GMP | gmp-rs |
|--------|-----|--------|
| Limb type | `mp_limb_t` — `unsigned long` (32 or 64 bit, platform dependent) | `u64` — unconditionally 64-bit |
| Limbs per `mpz_t` | `_mp_size` can be any `mp_size_t` (signed) | `mag: Vec<u64>` — `usize` length |

- **[INFO]** gmp-rs will not compile for 32-bit targets. GMP supports 32-bit limbs natively.
- **[API]** Callers on 32-bit hosts cannot use gmp-rs without a 64-bit CPU or emulation layer.

### 1.2 Sign encoding

| Aspect | GMP | gmp-rs |
|--------|-----|--------|
| Zero | `_mp_size == 0` | `sign == 0 && mag.is_empty()` |
| Positive | `_mp_size > 0` (limb count) | `sign == 1 && !mag.is_empty()` |
| Negative | `_mp_size < 0` (negated limb count) | `sign == -1 && !mag.is_empty()` |

- **[API]** GMP uses the size field both as sign *and* limb count, enabling the invariant `abs(_mp_size) == number_of_limbs`. gmp-rs stores them separately, doubling the state that must be kept coherent.
- **[CORRECTNESS]** gmp-rs's invariant is `sign == 0 iff mag.is_empty()`. A bug that sets sign to non-zero with an empty mag, or vice versa, would produce a corrupted value. The `norm()` constructor enforces this, but there is no runtime invariant checker.

### 1.3 No `mp_bitcnt_t` abstraction

- **[MISSING]** GMP uses `mp_bitcnt_t` (unsigned long) for all bit counts. gmp-rs uses `u32` for `mul_2exp`/`fdiv_*_2exp` and `usize` for `sizeinbase2`. This limits max shift distance to `u32::MAX` (~4 billion). GMP's `mp_bitcnt_t` can be wider on 64-bit platforms (typically 64 bits).

---

## 2. Initialization Functions

GMP provides **6 functions** in this category. gmp-rs provides **1**.

| GMP function | Equivalent in gmp-rs | Status |
|---|---|---|
| `mpz_init(x)` | `Mpz::new()` | ✅ |
| `mpz_inits(x, ...)` | — | **[MISSING]** Vararg init of multiple values |
| `mpz_init2(x, n)` | — | **[MISSING]** Pre-allocate `n` bits of storage |
| `mpz_clear(x)` | — (no-op, Rust RAII) | **[API]** GMP requires explicit deallocation; Rust's Drop handles it. No equivalent needed. |
| `mpz_clears(x, ...)` | — | **[MISSING]** Vararg clear (redundant in Rust) |
| `mpz_realloc2(x, n)` | — | **[MISSING]** Resize existing allocation to `n` bits |

### 2.1 `mpz_init2` — preallocation hint

- **[PERF]** GMP's `mpz_init2` lets callers pre-allocate storage for `n` bits, avoiding repeated reallocations as the value grows. gmp-rs's `new()` starts with a zero-capacity `Vec`, so every growth triggers a reallocation + copy. For values that grow incrementally (e.g. `from_decimal_string` parsing 18-digit chunks), this causes `O(log n)` reallocations. GMP can eliminate all reallocations with one `init2` call.
- **[API]** Callers who know an upper bound on precision cannot express that hint.

### 2.2 `mpz_realloc2`

- **[MISSING]** GMP can shrink or grow an existing allocation. gmp-rs relies entirely on `Vec`'s own `reallocate` strategy (typically double-when-full). No explicit resize/compact API exists.

---

## 3. Assignment Functions

| GMP function | Equivalent in gmp-rs | Status |
|---|---|---|
| `mpz_set(rop, op)` | `Clone::clone` | **[API]** not a method, but `clone()` works |
| `mpz_set_ui(rop, ulong)` | `set_ui(&mut self, u64)` | ✅ |
| `mpz_set_si(rop, slong)` | `set_si(&mut self, i64)` | ✅ |
| `mpz_set_d(rop, double)` | — | **[MISSING]** |
| `mpz_set_q(rop, mpq_t)` | — | **[MISSING]** |
| `mpz_set_f(rop, mpf_t)` | — | **[MISSING]** |
| `mpz_set_str(rop, str, base)` | `from_decimal_string(s)` | **[SEMANTIC]** base-10 only, no error return |
| `mpz_swap(rop1, rop2)` | — | **[MISSING]** |

### 3.1 `mpz_set_str` — base flexibility

- **[MISSING]** GMP accepts bases 2–62 (and -36 to -2 for lowercase-only). gmp-rs only accepts decimal strings (implied base 10).
- **[SEMANTIC]** GMP returns `-1` on parse error (empty string, non-digit in base, etc.). gmp-rs silently returns zero for malformed input and silently skips non-digit characters. See §5.1.

### 3.2 `mpz_set_d` — double assignment

- **[MISSING]** GMP can construct an `mpz` from a `f64` (truncating toward zero). gmp-rs has no float conversion.

### 3.3 `mpz_swap`

- **[MISSING]** GMP swaps two `mpz_t` values in constant time by swapping internal pointers. In gmp-rs, `core::mem::swap(&mut a, &mut b)` works thanks to `Mpz` being owned, but this is not exposed as a method.

---

## 4. Combined Init + Assignment

| GMP function | Equivalent in gmp-rs | Status |
|---|---|---|
| `mpz_init_set(rop, op)` | `op.clone()` | **[API]** |
| `mpz_init_set_ui(rop, ulong)` | `Mpz::from_u64(v)` | ✅ (as a constructor) |
| `mpz_init_set_si(rop, slong)` | `Mpz::from_i64(v)` | ✅ (as a constructor) |
| `mpz_init_set_d(rop, double)` | — | **[MISSING]** |
| `mpz_init_set_str(rop, str, base)` | `Mpz::from_decimal_string(s)` | **[SEMANTIC]** base-10 only, no error return |

---

## 5. Conversion Functions

### 5.1 String parsing: `from_decimal_string`

```
gmp-rs:   fn from_decimal_string(s: &str) -> Mpz
GMP:      int mpz_set_str(mpz_t rop, const char *str, int base)
```

- **[CORRECTNESS]** `from_decimal_string("")` silently returns zero. GMP returns `-1` (error).
- **[CORRECTNESS]** `from_decimal_string("-")` silently returns zero. GMP returns `-1` (error).
- **[CORRECTNESS]** `from_decimal_string("+")` silently returns zero. GMP returns `-1` (error).
- **[CORRECTNESS]** `from_decimal_string("12a34")` silently returns 1234 (non-digit characters are skipped). GMP returns `-1` for invalid characters. **This is the most dangerous parsing bug** — a caller typing `"123O"` (letter O instead of zero) would silently get 123 instead of an error.
- **[CORRECTNESS]** `from_decimal_string("  -42")` with leading whitespace: `trim()` handles this. But `from_decimal_string("-  42")` with whitespace between sign and digits: `strip_prefix('-')` returns `Some("  42")`, so it would try to parse `"  42"`. `trim()` is not called again. The leading space in the digit portion causes `bytes[i]` not to be `is_ascii_digit()` and is skipped during chunk parsing. So `"-  42"` parses as `"-  42"` → trimmed to `"-  42"` → neg=true, digits = `"  42"` → spaces skipped, `42` parsed. Result: `-42`. **This actually works but silently ignores the space.**
- **[CORRECTNESS]** `from_decimal_string("--42")`: `strip_prefix('-')` returns `Some("-42")`. Then parsing `"-42"` → `strip_prefix('-')` returns `Some("42")` → neg=true, digits=`"42"` → result=-42. So `"--42"` = `-42`. This is questionable but not dangerous.
- **[CORRECTNESS]** `from_decimal_string("+-42")`: `strip_prefix('-')` returns `None` (starts with '+'). Then `strip_prefix('+')` returns `Some("-42")`. neg=false, digits = `"-42"`. Parsing "-42": the '-' is not ascii_digit so it's skipped during chunk parsing. Chunks: start=0, end=3, chunk="-42" parsed digit-by-digit: '-' skipped (not digit), '4' → chunk=4, '2' → chunk=42. Result: 42. But should have been an error. **BUG**: silently returns wrong sign.

### 5.2 `to_decimal_string`

| Aspect | GMP `mpz_get_str` | gmp-rs `to_decimal_string` |
|--------|--------------------|----------------------------|
| Base | 2–62 (or -36 to -2) | 10 only |
| Buffer | Caller-provided or malloc'd | Always allocates new `String` |
| Error | Returns NULL on allocation failure | Infallible (panics on OOM) |

- **[MISSING]** No support for bases other than 10.
- **[MISSING]** No `mpz_get_d` (to `f64`) or `mpz_get_d_2exp` (frexp-style).

### 5.3 Missing `get_d` / `get_d_2exp`

- **[MISSING]** GMP can convert an `mpz` to a `double` (with optional exponent extraction). This is essential for interoperability with floating-point code. gmp-rs has no float conversion in either direction.

### 5.4 Native integer boundary issues

| Function | gmp-rs returns | GMP equivalent | Notes |
|---|---|---|---|
| `get_ui()` | Low 64-bit limb of |value| | `mpz_get_ui`: low `unsigned long` bits | ✅ matches |
| `get_si()` | `(low_limb as i64).wrapping_neg()` if negative | `mpz_get_si`: truncated to `signed long` | ✅ matches GMP doc |
| `to_i128()` | `Option<i128>`, None if > 2 limbs | No direct equivalent | **[API]** extension |

- **[INFO]** `to_i128` is a gmp-rs extension with no GMP equivalent. It's useful but its 2-limb limit is arbitrary and undocumented from a GMP perspective.

---

## 6. Arithmetic Functions

### 6.1 Implemented vs Missing

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_add(rop, op1, op2)` | `add(&self, other) -> Mpz` | ✅ (returns new value) |
| `mpz_add_ui(rop, op, ulong)` | `add_ui(&self, u64) -> Mpz` | ✅ (returns new value) |
| `mpz_sub(rop, op1, op2)` | `sub(&self, other) -> Mpz` | ✅ |
| `mpz_sub_ui(rop, op, ulong)` | `sub_ui(&self, u64) -> Mpz` | ✅ |
| `mpz_ui_sub(rop, ulong, op)` | — | **[MISSING]** |
| `mpz_mul(rop, op1, op2)` | `mul(&self, other) -> Mpz` | ✅ |
| `mpz_mul_si(rop, op, slong)` | — | **[MISSING]** |
| `mpz_mul_ui(rop, op, ulong)` | `mul_ui(&self, u64) -> Mpz` | ✅ |
| `mpz_addmul(rop, op1, op2)` | — | **[MISSING]** rop += op1 × op2 |
| `mpz_addmul_ui(rop, op, ulong)` | — | **[MISSING]** |
| `mpz_submul(rop, op1, op2)` | — | **[MISSING]** rop -= op1 × op2 |
| `mpz_submul_ui(rop, op, ulong)` | — | **[MISSING]** |
| `mpz_mul_2exp(rop, op, bits)` | `mul_2exp(&self, u32) -> Mpz` | ✅ |
| `mpz_neg(rop, op)` | `neg(&mut self)` (in-place only) | **[SEMANTIC]** no separate dest/src variant |
| `mpz_abs(rop, op)` | `abs(&mut self)` (in-place only) | **[SEMANTIC]** no separate dest/src variant |

### 6.2 `mpz_ui_sub` — unsigned-int-minus-mpz

- **[MISSING]** `mpz_ui_sub(rop, ulong, mpz)` computes `ulong - mpz` without constructing a temporary `Mpz` from the ulong. gmp-rs callers must write `Mpz::from_u64(v).sub(&x)`, which allocates.

### 6.3 `mpz_mul_si` — signed-long multiply

- **[MISSING]** GMP has both `mul_ui` (unsigned) and `mul_si` (signed). gmp-rs only has `mul_ui`. Multiplying by a negative i64 forces callers to construct a temporary `Mpz::from_i64(v)`.

### 6.4 `mpz_addmul` / `mpz_submul` — fused multiply-add/sub

- **[MISSING]** These are GMP's fused multiply-add/sub operations. `mpz_addmul(rop, a, b)` sets `rop += a * b`. Without them, callers must write `let tmp = a.mul(&b); rop = rop.add(&tmp);` which allocates an intermediate `Mpz`. GMP can fuse the operation and allocate only the final result.

### 6.5 `neg` / `abs` — no non-mutating variant

- **[API]** gmp-rs's `neg(&mut self)` and `abs(&mut self)` mutate in place. GMP's `mpz_neg(rop, op)` and `mpz_abs(rop, op)` take separate destination and source, allowing `rop = -op` without cloning `op`. Callers who want `let y = -x` without mutating `x` must write `let y = x.clone(); y.neg();` or use the `Neg` trait (`-&x`), which clones.

---

## 7. Division Functions

This is the **largest gap category** in gmp-rs. GMP provides **30 division-related functions** across 6 subgroups. gmp-rs provides **9**, and several of those have critical issues.

### 7.1 Coverage matrix

| GMP subgroup | Total GMP funcs | gmp-rs provides | Gap |
|---|---|---|---|
| Ceiling division (`mpz_cdiv_*`) | 9 | 0 | **[MISSING]** (all 9) |
| Floor division (`mpz_fdiv_*`) | 9 | 2 (`fdiv_q_2exp`, `fdiv_r_2exp`) | **[MISSING]** (all general floor variants) |
| Truncate division (`mpz_tdiv_*`) | 9 | 7 | Partial (missing `tdiv_q_2exp`, `tdiv_r_2exp`, `tdiv_r_ui`, `tdiv_qr_ui`) |
| Modulo (`mpz_mod`, `mpz_mod_ui`) | 2 | 0 | **[MISSING]** |
| Exact division (`mpz_divexact`, `mpz_divexact_ui`) | 2 | 0 | **[MISSING]** |
| Divisibility/congruence tests | 6 | 1 (`divisible_ui`) | **[MISSING]** (5 missing) |

### 7.2 No ceiling division

- **[MISSING]** GMP's `mpz_cdiv_q`, `mpz_cdiv_r`, `mpz_cdiv_qr` implement ceiling division (quotient rounds toward +∞). Required for any application that needs a rounding-up semantic for negative dividends.

### 7.3 No general floor division

- **[MISSING]** GMP's `mpz_fdiv_q`, `mpz_fdiv_r`, `mpz_fdiv_qr` implement floor division (quotient rounds toward −∞). Only the power-of-two variants (`fdiv_q_2exp`, `fdiv_r_2exp`) are implemented.

### 7.4 No modulo operation

- **[MISSING]** GMP's `mpz_mod(r, n, d)` returns the *non-negative* remainder (always `0 ≤ r < |d|`). This is distinct from `mpz_tdiv_r` (which can be negative). Without it, callers must manually adjust the result of `tdiv_r`. GMP also has `mpz_mod_ui` returning `unsigned long`.

### 7.5 No exact division

- **[MISSING]** GMP's `mpz_divexact(q, n, d)` and `mpz_divexact_ui(q, n, ulong)` assume `d` divides `n` exactly and use faster algorithms (linear time) without computing a remainder. Useful for `remove_pow10` and similar exact-factor operations.

### 7.6 Missing `tdiv_q_2exp` / `tdiv_r_2exp`

- **[MISSING]** GMP provides truncating versions of right-shift and low-bits extraction. gmp-rs only has the floor variants (`fdiv_q_2exp`, `fdiv_r_2exp`). For non-negative values these are identical, but for negative values the floor and truncate semantics differ.

### 7.7 Missing `tdiv_r_ui`

- **[MISSING]** GMP's `mpz_tdiv_r_ui(rop, n, d)` computes the remainder as an `mpz_t` (not a scalar). gmp-rs only has `tdiv_ui` which returns `u64`. The `Mpz`-valued remainder variant is absent.

### 7.8 Division algorithm: bit-by-bit (`mag_divmod`)

```rust
fn mag_divmod(a: &[u64], d: &[u64]) -> (Vec<u64>, Vec<u64>) {
    // Schoolbook long division of magnitudes (binary)
    let nbits = a.len() * 64;  // ← iterates over EVERY bit
    ...
    for i in (0..nbits).rev() {
        rem = rem.mul_2exp(1);                 // clones + allocates
        let bit = ...;                          // bit extraction
        if bit != 0 {
            rem = rem.add(&Mpz::from_u64(1));   // clones + allocates
        }
        if Self::cmp_mag(&rem.mag, d) != Ordering::Less {
            rem = rem.sub(&dm);                  // clones + allocates
            q[i / 64] |= 1u64 << (i % 64);
        }
    }
    ...
}
```

- **[PERF — CRITICAL]** This is **bit-by-bit long division**: `O(n²)` in the *number of bits* for each `tdiv_qr` call (where `n` is the bit width of the dividend). For a 256-bit dividend, this is 256 iterations. For a 1000-bit dividend, 1000 iterations. Each iteration allocates multiple `Mpz` values:
  - `rem.mul_2exp(1)` allocates a new `Vec`
  - `Mpz::from_u64(1)` allocates a new `Vec`
  - `rem.add(...)` allocates a new `Vec`
  - `rem.mag.clone()` inside `cmp_mag` (actually, `cmp_mag` takes slices so no clone — correct)
  - `rem.sub(&dm)` allocates a new `Vec`  
  - `dm = Mpz::norm(1, d.to_vec())` allocates once (but is hoistable — currently allocated once per `mag_divmod` call, which is fine)

  **Total allocations per `tdiv_qr` on a 256-bit dividend**: ~2,048 `Vec` allocations for the inner loop operations alone. For a 1000-bit dividend: ~8,000 allocations.

- **[PERF — CRITICAL]** GMP uses **Knuth's Algorithm D (schoolbook division)** for small sizes, transitioning to **Barrett division** or **Newton division** for larger operands. These are `O(n²)` in *limbs* (not bits) for schoolbook, and `O(n log n)` for Newton-based methods. gmp-rs's bit-at-a-time approach is slower by a factor equal to the limb bit-width (64× worse asymptotically in the constant factor, same asymptotic complexity).

- **[PERF — CRITICAL]** Every call to `mul_2exp(1)`, `add`, `sub`, `from_u64(1)` inside the loop reconstructs a full `Mpz` value. The bit-by-bit algorithm should be rewritten as a bit-by-bit scan with native u64 words, only extending `rem` when there's a borrow/carry. The current implementation's allocation rate is catastrophic for large values.

### 7.9 Division-by-zero safety

- **[SAFETY]** `tdiv_qr` has `debug_assert!(d.sign != 0)` — only checked in debug builds. In release mode, `mag_divmod` is called with an empty divisor (`d.mag` is empty for zero). The outer loop runs (`nbits = a.len() * 64`), `dm = Mpz::norm(1, d.to_vec())` produces `Mpz { sign: 0, mag: vec![] }` (zero). Then `rem.sub(&dm)` does nothing (subtracting zero). The loop runs to completion producing garbage quotient bits. **Division by zero in release builds produces a silently wrong result.**
- **[SAFETY]** `mag_divmod_u64(a, 0)` causes a division-by-zero panic in `(cur / d as u128)`. At least this panics rather than silently producing garbage.
- **[SAFETY]** `tdiv_q_ui`, `tdiv_ui`, `divisible_ui` all call `mag_divmod_u64` without checking `d != 0`.

---

## 8. Exponentiation Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_pow_ui(rop, base, exp)` | `pow_ui(&self, n)` | ✅ (returns new value) |
| `mpz_ui_pow_ui(rop, base, exp)` | `ui_pow_ui(base: u64, exp: u32)` | ✅ |
| `mpz_powm(rop, base, exp, mod)` | — | **[MISSING]** |
| `mpz_powm_ui(rop, base, ulong, mod)` | — | **[MISSING]** |
| `mpz_powm_sec(rop, base, exp, mod)` | — | **[MISSING]** |

### 8.1 No modular exponentiation

- **[MISSING]** `mpz_powm` (RSA-style modular exponentiation) is entirely absent. This is the most requested bignum operation for cryptography. GMP uses the k-ary sliding window method for `powm`, with a side-channel-resistant variant (`powm_sec`).

### 8.2 `ui_pow_ui` uses naive O(n) loop

```rust
pub fn ui_pow_ui(base: u64, exp: u32) -> Mpz {
    let mut r = Mpz::from_u64(1);
    let b = Mpz::from_u64(base);
    for _ in 0..exp {
        r = r.mul(&b);
    }
    r
}
```

- **[PERF — CRITICAL]** This is `O(exp)` multiplications. `pow_ui` (on `&self`) correctly uses exponentiation by squaring (`O(log exp)`), but `ui_pow_ui` is naive. `Mpz::ui_pow_ui(2, 10_000)` would do 10,000 full-precision multiplications, each on growing integers. GMP's `mpz_ui_pow_ui` uses the same squaring algorithm as `mpz_pow_ui`.

### 8.3 `pow_ui` correctness

- **[INFO]** gmp-rs's `pow_ui` uses exponentiation by squaring. `0.pow_ui(0)` returns 1 (matching GMP's documented convention).
- **[CORRECTNESS]** gmp-rs's `pow_ui` clones `self` at the start. For very large `self` with large `exp`, this clone is a one-time cost. Correct.

---

## 9. Root Extraction Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_root(rop, op, n)` | — | **[MISSING]** General n-th root |
| `mpz_rootrem(root, rem, u, n)` | — | **[MISSING]** n-th root with remainder |
| `mpz_sqrt(rop, op)` | `isqrt(&self) -> Mpz` | ✅ (returns new value) |
| `mpz_sqrtrem(rop1, rop2, op)` | — | **[MISSING]** sqrt with remainder |
| `mpz_perfect_power_p(op)` | — | **[MISSING]** |
| `mpz_perfect_square_p(op)` | — | **[MISSING]** |

### 9.1 `isqrt` performance

- **[PERF]** `isqrt` uses Newton iteration on integer `Mpz` values. Each iteration calls `self.tdiv_q(&x)`, which invokes the bit-by-bit division algorithm (`mag_divmod`). For a Newton iteration count of ~`log2(bits)`, and each division being `O(bits²)`, the total cost is `O(bits² log bits)`. GMP's `mpz_sqrt` uses Newton with much faster division (multi-limb schoolbook or Barrett), and early-terminates iterations as precision converges.

### 9.2 No `perfect_square_p` / `perfect_power_p`

- **[MISSING]** GMP can test whether an integer is a perfect square or perfect power without computing the full root (uses fast checks). Without these, callers must call `isqrt` and square check manually.

---

## 10. Number Theoretic Functions

**This entire category is missing from gmp-rs.**

| GMP function | Count | Status |
|---|---|---|
| `mpz_probab_prime_p` | 1 | **[MISSING]** |
| `mpz_nextprime` / `mpz_prevprime` | 2 | **[MISSING]** |
| `mpz_gcd` / `mpz_gcd_ui` / `mpz_gcdext` | 3 | **[MISSING]** |
| `mpz_lcm` / `mpz_lcm_ui` | 2 | **[MISSING]** |
| `mpz_invert` | 1 | **[MISSING]** |
| `mpz_jacobi` / `mpz_legendre` / `mpz_kronecker` variants | 6 | **[MISSING]** |
| `mpz_remove` | 1 | **[SEMANTIC]** `remove_pow10` is a hardcoded special case for radix 10 only |
| `mpz_fac_ui` / `mpz_2fac_ui` / `mpz_mfac_uiui` | 3 | **[MISSING]** |
| `mpz_primorial_ui` | 1 | **[MISSING]** |
| `mpz_bin_ui` / `mpz_bin_uiui` | 2 | **[MISSING]** |
| `mpz_fib_ui` / `mpz_fib2_ui` | 2 | **[MISSING]** |
| `mpz_lucnum_ui` / `mpz_lucnum2_ui` | 2 | **[MISSING]** |

### 10.1 `remove_pow10` — special case of `mpz_remove`

- **[SEMANTIC]** GMP's `mpz_remove(rop, op, f)` removes *all* occurrences of an arbitrary factor `f`, not just 10. gmp-rs's `remove_pow10` is hardcoded to divide by 10 repeatedly. For internal use in porting arithmetic this may be sufficient, but it's not a general `mpz_remove` equivalent.

---

## 11. Comparison Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_cmp(op1, op2)` | `cmp(&self, other)` via `impl Ord` | ✅ |
| `mpz_cmp_d(op, double)` | — | **[MISSING]** |
| `mpz_cmp_si(op, slong)` | — | **[MISSING]** |
| `mpz_cmp_ui(op, ulong)` | — | **[MISSING]** |
| `mpz_cmpabs(op1, op2)` | `cmpabs(&self, other)` | ✅ |
| `mpz_cmpabs_d(op, double)` | — | **[MISSING]** |
| `mpz_cmpabs_ui(op, ulong)` | — | **[MISSING]** |
| `mpz_sgn(op)` | `sgn(&self) -> i32` | ✅ |

### 11.1 Missing scalar comparisons

- **[MISSING]** GMP provides `mpz_cmp_si` and `mpz_cmp_ui` for comparing directly against native integers without constructing a temporary `Mpz`. gmp-rs callers must write `self.cmp(&Mpz::from_i64(v))`, which allocates.

---

## 12. Logical & Bit Manipulation Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_and(rop, op1, op2)` | — | **[MISSING]** |
| `mpz_ior(rop, op1, op2)` | — | **[MISSING]** |
| `mpz_xor(rop, op1, op2)` | — | **[MISSING]** |
| `mpz_com(rop, op)` | `com(&self) -> Mpz` | ✅ |
| `mpz_popcount(op)` | — | **[MISSING]** |
| `mpz_hamdist(op1, op2)` | — | **[MISSING]** |
| `mpz_scan0(op, start)` | — | **[MISSING]** |
| `mpz_scan1(op, start)` | — | **[MISSING]** |
| `mpz_setbit(rop, bit)` | — | **[MISSING]** |
| `mpz_clrbit(rop, bit)` | — | **[MISSING]** |
| `mpz_combit(rop, bit)` | — | **[MISSING]** |
| `mpz_tstbit(op, bit)` | — | **[MISSING]** |

### 12.1 Missing bitwise operations

- **[MISSING — CRITICAL FOR SAFETY]** GMP's bitwise operations (`mpz_and`, `mpz_ior`, `mpz_xor`) operate on *two's complement* representations, not sign-magnitude. The sign-magnitude representation used by gmp-rs means that a naive bitwise and on `vec[u64]` limbs would produce an incorrect result for negative values. To correctly implement these operations, gmp-rs would need to convert to two's complement internally, perform the operation, and convert back to sign-magnitude — or implement them entirely at the semiotic level using arithmetic identities (`x & y = (x | y) - (x ^ y)` etc.). This is a significant correctness hazard for anyone attempting to add bitwise ops without understanding the two's-complement requirement.

### 12.2 `com()` correctness

- **[CORRECTNESS]** gmp-rs's `com()` computes `add(&Mpz::from_u64(1)).into_neg()` which gives `- (x + 1)`. For sign-magnitude, `-x - 1` is NOT the same as two's-complement bitwise NOT for all values, but it IS identical to the mathematical result of `mpz_com`. GMP's documentation states: `mpz_com` computes the one's complement, i.e., `~op`. For `mpz_t` values which use a sign-magnitude representation internally, GMP defines `~op = -op - 1`. So gmp-rs's implementation is correct by definition.

### 12.3 No `tstbit` / `setbit` / `clrbit` / `combit`

- **[MISSING]** Without `tstbit`, callers cannot query individual bits without extracting limbs manually. Without `setbit`/`clrbit`/`combit`, callers cannot modify bits without decomposing and recomposing the value.

### 12.4 No `popcount` / `hamdist`

- **[MISSING]** GMP can count 1-bits (popcount) and compute Hamming distance (popcount of XOR). Useful in combinatorics and cryptography.

---

## 13. I/O Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_out_str(FILE*, base, op)` | `Display::fmt` (implicitly base 10) | **[SEMANTIC]** |
| `mpz_inp_str(rop, FILE*, base)` | — | **[MISSING]** |
| `mpz_out_raw(FILE*, op)` | — | **[MISSING]** |
| `mpz_inp_raw(rop, FILE*)` | — | **[MISSING]** |

- **[MISSING]** No raw binary serialization. GMP's `mpz_out_raw`/`mpz_inp_raw` produce a portable binary representation (4-byte limb count in big-endian, followed by limbs in little-endian). Without this, there is no standardized way to serialize an `Mpz` to bytes and deserialize unambiguously on any platform.

---

## 14. Random Number Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_urandomb(rop, state, n)` | — | **[MISSING]** |
| `mpz_urandomm(rop, state, n)` | — | **[MISSING]** |
| `mpz_rrandomb(rop, state, n)` | — | **[MISSING]** |
| `mpz_random(rop, max_size)` | — | **[MISSING]** (obsolete in GMP) |
| `mpz_random2(rop, max_size)` | — | **[MISSING]** (obsolete in GMP) |

- **[MISSING — COMPLETE]** The entire random number category is absent. For a `no_std` crate, random number generation would need an external RNG provider. GMP bundles its own (`gmp_randstate_t` backed by Mersenne Twister or a hardware source).

---

## 15. Integer Import/Export

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_import(rop, count, order, size, endian, nails, op)` | — | **[MISSING]** |
| `mpz_export(rop, countp, order, size, endian, nails, op)` | — | **[MISSING]** |

- **[MISSING — COMPLETE]** GMP's `mpz_import`/`mpz_export` are the canonical way to convert between opaque byte buffers and `mpz_t` values. They handle endianness, word size, and "nail" bits. Without them, callers must manually construct `Mpz` from bytes by iterating limbs — error-prone and non-portable.

---

## 16. Miscellaneous Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_fits_ulong_p(op)` | `fits_ulong(&self) -> bool` | ✅ |
| `mpz_fits_slong_p(op)` | — | **[MISSING]** |
| `mpz_fits_uint_p(op)` | — | **[MISSING]** |
| `mpz_fits_sint_p(op)` | — | **[MISSING]** |
| `mpz_fits_ushort_p(op)` | — | **[MISSING]** |
| `mpz_fits_sshort_p(op)` | — | **[MISSING]** |
| `mpz_odd_p(op)` | — | **[MISSING]** (trivial: `mag[0] & 1`) |
| `mpz_even_p(op)` | — | **[MISSING]** (trivial: `mag[0] & 1 == 0`) |
| `mpz_sizeinbase(op, base)` | `sizeinbase2(&self) -> usize` | **[SEMANTIC]** base-2 only, returns bit count |

### 16.1 `sizeinbase2` vs `mpz_sizeinbase`

- **[SEMANTIC]** GMP's `mpz_sizeinbase(op, base)` returns an *upper bound* on the number of characters needed for the given base. For base 2, `mpz_sizeinbase(op, 2)` returns the number of significant bits (same as gmp-rs's `sizeinbase2`). But GMP also supports any base 2–62, and importantly, for bases that are powers of two, it computes a tighter bound. gmp-rs only supports base 2.
- **[CORRECTNESS]** GMP's `mpz_sizeinbase(0, 2)` returns 1 (one bit/character for "0"). gmp-rs matches this.

### 16.2 Missing `odd_p` / `even_p`

- **[MISSING]** Trivial to implement (`self.mag.first().map(|l| l & 1).unwrap_or(0) != 0`), but absent.

---

## 17. Low-Level / Limb-Access Functions

| GMP function | gmp-rs | Status |
|---|---|---|
| `mpz_size(op)` | `size(&self) -> usize` | ✅ (returns limb count) |
| `mpz_getlimbn(op, n)` | — | **[MISSING]** |
| `mpz_limbs_read(x)` | — | **[MISSING]** |
| `mpz_limbs_write(x, n)` | — | **[MISSING]** |
| `mpz_limbs_modify(x, n)` | — | **[MISSING]** |
| `mpz_limbs_finish(x, s)` | — | **[MISSING]** |
| `mpz_roinit_n(x, xp, xs)` | — | **[MISSING]** |
| `MPZ_ROINIT_N(xp, xs)` (macro) | — | **[MISSING]** |

### 17.1 No limb-level access

- **[MISSING]** GMP's limb access functions allow callers to read/write individual limbs directly. gmp-rs keeps `mag` private with no public accessors. Callers who need the raw limbs (for interop with C or hand-optimized algorithms) cannot access them.

---

## 18. Synthesis & Risk Register

### 18.1 Coverage summary

| Category | GMP count | gmp-rs count | Coverage |
|---|---|---|---|
| Initialization | 6 | 1 | 17% |
| Assignment | 8 | 4 | 50% |
| Combined Init+Assign | 5 | 4 | 80% |
| Conversion | 5 | 4 | 80% |
| Arithmetic | 14 | 10 | 71% |
| Division | 30 | 9 | 30% |
| Exponentiation | 5 | 2 | 40% |
| Root Extraction | 6 | 1 | 17% |
| Number Theoretic | 26 | 0 | 0% |
| Comparison | 8 | 3 | 38% |
| Logical / Bit | 12 | 1 | 8% |
| I/O | 4 | 0 | 0% |
| Random Numbers | 5 | 0 | 0% |
| Import/Export | 2 | 0 | 0% |
| Miscellaneous | 9 | 2 | 22% |
| Low-Level | 9 | 1 | 11% |
| **Total** | **~154** | **42** | **~27%** |

### 18.2 Risk register (severity: CRITICAL / HIGH / MEDIUM / LOW)

| ID | Severity | Category | Issue |
|---|---|---|---|
| R01 | **CRITICAL** | Division §7.8 | `mag_divmod` is bit-by-bit `O(n²)` in bits, allocating ~8 Vecs per bit. A 256-bit division does ~2,000 allocations. Production use for >128-bit values will be catastrophically slow. |
| R02 | **CRITICAL** | Division §7.9 | Division by zero is not checked in release builds (`debug_assert!` only). Silent wrong result instead of panic/error. |
| R03 | **CRITICAL** | Parsing §5.1 | `from_decimal_string` silently skips non-digit characters (`"12O34"` → `1234`) and silently ignores malformed inputs (`""`, `"-"`, `"+"` → `0`). This is a data corruption risk. |
| R04 | **HIGH** | Exponentiation §8.2 | `ui_pow_ui` uses `O(exp)` multiplications instead of `O(log exp)`. `ui_pow_ui(2, 10_000)` will take seconds or minutes. |
| R05 | **HIGH** | Division §7.3, §7.4 | No floor division or modulo operation for general (non-power-of-2) operands. Many numeric algorithms require these. |
| R06 | **HIGH** | Alloc §N.1 | The crate is NOT `no_alloc`. Every arithmetic operation allocates. Systems that cannot tolerate allocation (embedded, real-time, safety-critical) cannot use gmp-rs. |
| R07 | **HIGH** | Number Theory §10 | Entire category missing. No GCD, no primality testing, no factorial, no binomial, no Fibonacci. |
| R08 | **MEDIUM** | Bitwise §12 | All bitwise ops except `com` are missing. Sign-magnitude representation makes correct implementation non-trivial. |
| R09 | **MEDIUM** | I/O §13 | No raw binary serialization. No way to portably store/transmit an `Mpz`. |
| R10 | **MEDIUM** | Comparison §11 | No scalar comparisons (`cmp_si`, `cmp_ui`). Every comparison against a native integer allocates a temporary `Mpz`. |
| R11 | **MEDIUM** | Arithmetic §6.6 | `add_ui`/`sub_ui`/`mul_ui` allocate a temporary `Mpz` from the scalar. GMP processes the single-limb scalar inline. |
| R12 | **LOW** | Init §2.1 | No preallocation hint (`mpz_init2`). Growing values reallocate repeatedly. |
| R13 | **LOW** | Root §9.1 | `isqrt` uses Newton + slow division. No `sqrtrem` or `root` variants. |
| R14 | **LOW** | Misc §16 | Missing `odd_p`/`even_p` checks. |
| R15 | **LOW** | Repr §1.3 | Bit shifts limited to `u32`. `mp_bitcnt_t` supports platform-native width. |
| R16 | **LOW** | Repr §1.2 | No invariant checker. Corrupted state (sign != 0 with empty mag) would cause silent UB. |

---

## 19. The `no_alloc` Constraint

### 19.1 Current status: NOT `no_alloc`

The current implementation uses `Vec<u64>` as its core storage. Every operation that produces a value allocates:

```rust
pub struct Mpz {
    sign: i8,
    mag: Vec<u64>,  // <-- heap allocation
}
```

**Allocation points per operation:**

| Operation | Allocations |
|---|---|
| `from_u64(v)` | 1 (`vec![v]`) |
| `add(&a, &b)` | 1 (`mag_add` output Vec) |
| `sub(&a, &b)` | 1 (`mag_sub` output Vec) |
| `mul(&a, &b)` | 1 (`mag_mul` output Vec) |
| `tdiv_qr(&a, &d)` | ~`2 * bitlen(a)` + 2 (bit-by-bit loop) |
| `to_decimal_string()` | 1 + chunk Vec |
| `from_decimal_string(s)` | 1 per 18-digit chunk |
| `clone()` | 1 (full `Vec` clone) |

Even `add_ui(5)` allocates: it constructs `Mpz::from_u64(5)` (1 alloc) then calls `add` (1 alloc on the result), totaling 2 allocations for a single-limb addition.

### 19.2 What true `no_alloc` would require

True `no_alloc` means no heap allocation, no `Vec`, no `Box`, no `String`. This implies:

1. **Fixed-capacity limb storage**, e.g.:
   ```rust
   use heapless::Vec;
   const MAX_LIMBS: usize = 64; // 4096-bit precision
   pub struct Mpz {
       sign: i8,
       mag: Vec<u64, MAX_LIMBS>,
   }
   ```

2. **Fallible arithmetic** — every operation that might exceed `MAX_LIMBS` must return `Result`:
   ```rust
   pub fn add(&self, other: &Mpz) -> Result<Mpz, CapacityError>;
   pub fn mul(&self, other: &Mpz) -> Result<Mpz, CapacityError>;
   ```

3. **No `String` in the public API** — `to_decimal_string` would need a caller-provided buffer:
   ```rust
   pub fn to_decimal_string(&self, buffer: &mut [u8]) -> Result<&str, BufferTooSmall>;
   ```

4. **No allocation in internal algorithms** — the bit-by-bit division loop, which currently allocates on every iteration, would need to work entirely on fixed scratch arrays.

5. **Dependency on `heapless` or `arrayvec`** — `heapless::Vec` is the standard `no_alloc` dynamic-array type, but it adds the crate as a dependency.

### 19.3 Design tradeoffs

| Approach | Pros | Cons |
|---|---|---|
| **Current: `Vec<u64>` (external alloc)** | Unlimited precision, simple code, rich API | Heaps allocations, not `no_alloc` |
| **Fixed-capacity: `heapless::Vec<u64, N>`** | True `no_alloc`, deterministic, embedded-friendly | Fixed max precision, fallible API, must choose `N` up front |
| **Stack-optimized small-vec** (small `[u64; 2]` inline, spill to alloc) | Zero alloc for small values, unlimited precision | Complex, still allocates for large values |

### 19.4 Recommendation

For a crate that claims `no_alloc`, the `Vec<u64>` representation is a **hard violation**. Either:

- **(a)** Remove `no_alloc` from the description and document it as `alloc`-dependent (current reality).
- **(b)** Refactor to `heapless::Vec<u64, N>` with `N = 64` (4096 bits, covering all practical cases for the intended decimal port use case) and make all operations fallible.
- **(c)** Adopt a small-vector optimisation strategy where single-limb values live on the stack and multi-limb values spill to `alloc::Box<[u64; N]>` (intermediate — still allocates, but avoids per-op allocation for single-limb values).

The safest path for "serious critical safety infrastructure" is **(b)**, but it requires a major rearchitecture of the entire crate.

---

## Appendix A: gmp-rs-Specific Extensions (not in GMP)

These functions exist in gmp-rs with no GMP `mpz_*` equivalent:

| Function | Purpose | Dependence on gmp-rs internals |
|---|---|---|
| `from_u128(v)` / `from_i128(v)` | Construct from Rust native 128-bit types | No |
| `to_i128(&self) -> Option<i128>` | Convert to i128 when it fits | No |
| `from_decimal_string(s)` → `Mpz` | Parse decimal string (no error return) | No |
| `to_decimal_string(&self) -> String` | Format as decimal string | Yes (`String`) |
| `remove_pow10(&mut self) -> u32` | Factor out powers of 10 | No |
| `set_ull(val)` / `mpz_get_ull` | 64-bit extension of `mpz_set_ui`/`mpz_get_ui` | No |
| `mpz_set_sll(val)` / `mpz_get_sll` | 64-bit signed extension | No |
| `fits_ulong(&self) -> bool` | "Fits in unsigned long" predicate | No |
| `sizeinbase2(&self) -> usize` | Bit count (GMP calls `mpz_sizeinbase(op, 2)`) | No |
| `Neg` trait, operator overloads, `From` impls | Rust ergonomics | No |

## Appendix B: GMP Functions That Are Not Recommended for Port (Obsolete/internal)

GMP marks the following as obsolete or internal. They are excluded from the gap analysis:

- `mpz_random` / `mpz_random2` — obsolete (replaced by `mpz_urandomb`/`mpz_urandomm`)
- `mpz_array_init` — obsolete, do not use
- `mpz_oddfac_1` — internal, undocumented
- `_mpz_realloc` — internal, use `mpz_realloc2`

## Appendix C: Rust Trait Implementations vs GMP Idiom

GMP uses a **destination-first** calling convention:
```c
mpz_add(rop, op1, op2);   // rop = op1 + op2
```

gmp-rs uses a **return-value** calling convention (matching Rust idioms):
```rust
let rop = op1.add(op2);        // returns new value
let rop = &op1 + &op2;         // operator overload
```

This means every arithmetic operation in gmp-rs **allocates new storage for the result**, whereas GMP can reuse the destination's existing allocation. A GMP caller who writes:
```c
mpz_add(x, x, y);   // x += y, reuses x's storage
```
would write the equivalent in gmp-rs as:
```rust
x = x.add(&y);      // old x dropped, new x allocated
```
or, using the `AddAssign` trait (which gmp-rs does NOT implement):
```rust
x += &y;            // would reuse x's allocation if AddAssign were implemented
```

- **[MISSING]** No `AddAssign`, `SubAssign`, `MulAssign`, etc. traits. These would allow in-place mutation without always allocating a new result.
- **[PERF]** Every arithmetic expression like `a + b + c` allocates at least 2 `Vec`s. GMP with destination reuse needs at most 1 allocation for the final result.
