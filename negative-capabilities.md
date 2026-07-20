# gmp-rs v0.2.1: Forensic Negative-Capabilities Audit vs GNU GMP 6.3.0 `mpz_*`

> **Date:** 2026-07-21
> **gmp-rs version:** 0.2.1
> **GMP reference:** 6.3.0 (released 2023)
> **Crate root:** `src/lib.rs` — single file, 2549 lines (52 tests)
>
> **Constraint envelope:**
> - `no_unsafe` ✅ — `#![forbid(unsafe_code)]` enforced at compile time
> - `no_std` ✅ — pure `core` only; no `std` dependency anywhere
> - `no_alloc` ✅ — fixed `[u64; 8]` limb array; `extern crate alloc;` used only in `#[cfg(test)]`
> - **Panic-free in release** ❌ — division by zero panics; `unwrap_or_else(|_| Mpz::new())` silently returns zero on capacity overflow
>
> This document exhaustively catalogues every deviation, omission, semantic quirk,
> correctness edge-case, and architectural limitation of **gmp-rs v0.2.1** relative
> to GNU GMP 6.3.0 `mpz_*` integer API.  It is organised by GMP category with a
> concluding synthesis and risk register.
>
> **Convention:** Entries tagged with:
> | Tag | Meaning |
> |-----|---------|
> | **[MISSING]** | Function or variant absent from gmp-rs |
> | **[SEMANTIC]** | Different behaviour for same-named operation |
> | **[CORRECTNESS]** | Provably wrong in a documented edge-case |
> | **[PERF]** | Algorithmically or asymptotically inferior |
> | **[SAFETY]** | Missing precondition check w/ panic or silent wrong result |
> | **[API]** | Ergonomic or structural difference from GMP |
> | **[INFO]** | Informational, not a defect |
> | **[FIXED]** | Previously reported issue now resolved in v0.2.1 |

---

## A. Coverage Summary

| Category | GMP count | gmp-rs count | Coverage |
|---|---|---|---|
| Initialization | 6 | 2 | 33% |
| Assignment | 8 | 6 | 75% |
| Combined Init+Assign | 5 | 4 | 80% |
| Conversion | 5 | 5 | 100% |
| Arithmetic | 14 | 14 | 100% |
| Division (all variants) | 30 | 21 | 70% |
| Exponentiation | 5 | 4 | 80% |
| Root Extraction | 6 | 6 | 100% |
| Number Theoretic | 26 | 13 | 50% |
| Comparison | 8 | 6 | 75% |
| Logical / Bit | 12 | 12 | 100% |
| I/O | 4 | 0 | 0% |
| Random Numbers | 5 | 0 | 0% |
| Import/Export | 2 | 0 | 0% |
| Miscellaneous | 9 | 9 | 100% |
| Low-Level | 9 | 2 | 22% |
| **Total** | **~154** | **104** | **~68%** |

---

## B. Functions Implemented vs GMP — Detailed Per-Category Audit

### B1. Initialization (6 GMP / 2 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_init` | `new()` | ✅ |
| `mpz_inits(x, ...)` | — | **[MISSING]** Vararg (not expressible in Rust) |
| `mpz_init2(x, n)` | — | **[MISSING]** Pre-allocate n bits |
| `mpz_clear` | — | **[API]** Handled by Rust `Drop` |
| `mpz_clears(x, ...)` | — | **[MISSING]** Vararg |
| `mpz_realloc2(x, n)` | — | **[MISSING]** Not applicable to fixed arrays |

- **[MISSING]** `mpz_init2` — caller cannot provide an allocation size hint. GMP optimises reallocation via this hint. gmp-rs's fixed array makes this moot, but a caller porting GMP code must drop this call.
- **[API]** GMP's `mpz_clear` deallocates memory. gmp-rs's `Mpz` implements `Drop` (via `Clone`'s `Copy` semantics on `[u64; 8]`), so no explicit destructor is needed. This is a Rust idiom difference.

### B2. Assignment (8 GMP / 6 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_set(rop, op)` | `set(&mut self, src)` | ✅ |
| `mpz_set_ui` | `set_ui(&mut self, u64)` | ✅ |
| `mpz_set_si` | `set_si(&mut self, i64)` | ✅ |
| `mpz_set_d` | — | **[MISSING]** From `f64` |
| `mpz_set_q` | — | **[MISSING]** From rational (no `mpq_t` exists in gmp-rs) |
| `mpz_set_f` | — | **[MISSING]** From float (no `mpf_t` exists in gmp-rs) |
| `mpz_set_str` | `from_decimal_str` | **[SEMANTIC]** Base 10 only; returns `Result` instead of `int` |
| `mpz_swap` | `swap(&mut self, other)` | ✅ |

- **[MISSING]** `mpz_set_d` — no `f64` → `Mpz` conversion. GMP converts a double to an `mpz` by truncating toward zero. Not implementable without `f64` bit inspection or `unsafe` transmute.
- **[SEMANTIC]** GMP's `mpz_set_str` returns `-1` on failure and accepts bases 2–62. gmp-rs's `from_decimal_str` returns `ParseError` and accepts base 10 only.

### B3. Combined Init + Assignment (5 GMP / 4 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_init_set` | `Clone::clone` | ✅ |
| `mpz_init_set_ui` | `from_u64` | ✅ |
| `mpz_init_set_si` | `from_i64` | ✅ |
| `mpz_init_set_d` | — | **[MISSING]** |
| `mpz_init_set_str` | `from_decimal_str` | **[SEMANTIC]** Base 10 only |

### B4. Conversion (5 GMP / 5 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_get_ui` | `get_ui()` | ✅ |
| `mpz_get_si` | `get_si()` | ✅ |
| `mpz_get_d` | `get_d()` | ✅ |
| `mpz_get_d_2exp` | `get_d_2exp()` | ✅ |
| `mpz_get_str` | `write_decimal_buf()` / `Display::fmt` | **[SEMANTIC]** |

- **[SEMANTIC]** `mpz_get_str` accepts a caller-supplied buffer and supports bases 2–62. gmp-rs's `write_decimal_buf` writes into a `&mut [u8]` and is base-10 only.
- **[CORRECTNESS]** gmp-rs's `get_d()` iterates over limbs from high to low, multiplying the accumulator by `2^64` at each step. For values near `f64::MAX` this loses precision silently (no overflow/Infinity check). GMP's `mpz_get_d` performs rounding — gmp-rs's version truncates.
- **[CORRECTNESS]** gmp-rs's `get_d_2exp()` right-shifts by `top_bit - 52` then converts via `get_d()`, then normalises with a `while` loop. For `top_bit < 52` the shift is 0 and `get_d()` returns the exact value. For `top_bit > 52`, the right-shift discards low bits — this is correct GMP behaviour. However, the `while` loop may not terminate for very small mantissas (e.g., if `get_d()` rounds to zero). The `if mantissa == 0.0 { return None }` guard prevents infinite loops but is not how GMP behaves (GMP returns `mantissa = 0.0, exp = 0` for zero).

### B5. Arithmetic (14 GMP / 14 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_add` | `try_add` | ✅ |
| `mpz_add_ui` | `try_add_ui` | ✅ |
| `mpz_sub` | `try_sub` | ✅ |
| `mpz_sub_ui` | `try_sub_ui` | ✅ |
| `mpz_ui_sub` | `try_ui_sub` (static) | ✅ |
| `mpz_mul` | `try_mul` | ✅ |
| `mpz_mul_si` | `try_mul_si` | ✅ |
| `mpz_mul_ui` | `try_mul_ui` | ✅ |
| `mpz_addmul` | `try_addmul` | ✅ |
| `mpz_addmul_ui` | `try_addmul_ui` | ✅ |
| `mpz_submul` | `try_submul` | ✅ |
| `mpz_submul_ui` | `try_submul_ui` | ✅ |
| `mpz_mul_2exp` | `try_mul_2exp` | ✅ |
| `mpz_neg` | `neg(&mut self)` + `neg_to()` | ✅ |
| `mpz_abs` | `abs(&mut self)` + `abs_to()` | ✅ |

- **[API]** All gmp-rs arithmetic functions return `Result` to handle capacity exhaustion. GMP's functions write to a destination operand and never return an error (they abort on out-of-memory).
- **[SEMANTIC]** `try_ui_sub` delegates to `Mpz::from_u64(v).try_sub(other)`. This creates a temporary `Mpz` from the `u64`, which is correct but could be more efficient by avoiding the intermediate allocation (though `no_alloc` uses stack arrays, so there's no heap cost — just a trivial stack copy).
- **[CORRECTNESS]** `try_sub` clones `other`, negates its sign, then calls `try_add`. This double-clones `other` (once explicitly, once in `try_add`'s internal `clone` for the zero-sign branch). For `other == 0`, this is a wasted clone.

### B6. Division (30 GMP / 21 gmp-rs)

#### B6a. Truncating (9 GMP / 8 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_tdiv_q` | `tdiv_q` | ✅ |
| `mpz_tdiv_r` | `tdiv_r` | ✅ |
| `mpz_tdiv_qr` | `tdiv_qr` | ✅ |
| `mpz_tdiv_q_ui` | `tdiv_q_ui` | ✅ |
| `mpz_tdiv_r_ui` | `tdiv_r_ui` | ✅ |
| `mpz_tdiv_qr_ui` | — | **[MISSING]** Combined (q, r) with u64 divisor returning `(Mpz, u64)` |
| `mpz_tdiv_ui` | `tdiv_ui` | ✅ |
| `mpz_tdiv_q_2exp` | `tdiv_q_2exp` | ✅ |
| `mpz_tdiv_r_2exp` | `tdiv_r_2exp` | ✅ |

- **[MISSING]** `tdiv_qr_ui` — GMP provides a combined function `mpz_tdiv_qr_ui(q, r, n, d)` that writes both quotient and remainder in one call. gmp-rs callers must call `tdiv_q_ui` and `tdiv_ui` separately, performing the division twice.
- **[CORRECTNESS]** `tdiv_q_2exp` and `tdiv_r_2exp` are currently delegating to `fdiv_q_2exp` and `fdiv_r_2exp` respectively. For non-negative values these are identical. For negative values, truncating and floor semantics differ — but since gmp-rs's `fdiv_q_2exp`/`fdiv_r_2exp` operate on the magnitude only (ignoring sign for the shift), they behave identically to truncating. **This is correct by coincidence** — the sign is preserved unchanged through both operations, and the bit-level masking of `fdiv_r_2exp` on the magnitude produces the same result as a truncating mask.

#### B6b. Floor (9 GMP / 6 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_fdiv_q` | `try_fdiv_q` | ✅ |
| `mpz_fdiv_r` | `try_fdiv_r` | ✅ |
| `mpz_fdiv_qr` | `try_fdiv_qr` | ✅ |
| `mpz_fdiv_q_ui` | `try_fdiv_q_ui` | ✅ |
| `mpz_fdiv_r_ui` | — | **[MISSING]** Floor remainder by u64 as Mpz |
| `mpz_fdiv_qr_ui` | `try_fdiv_qr_ui` | ✅ |
| `mpz_fdiv_ui` | `fdiv_ui` | ✅ |
| `mpz_fdiv_q_2exp` | `fdiv_q_2exp` | ✅ |
| `mpz_fdiv_r_2exp` | `fdiv_r_2exp` | ✅ |

- **[MISSING]** `try_fdiv_r_ui` — GMP provides a function that returns the floor remainder as an `mpz_t` (not a scalar `u64`). gmp-rs callers must use `try_fdiv_qr_ui` and discard the quotient, or convert the `u64` result to `Mpz::from_u64`.

#### B6c. Ceiling (9 GMP / 3 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_cdiv_q` | `try_cdiv_q` | ✅ |
| `mpz_cdiv_r` | `try_cdiv_r` | ✅ |
| `mpz_cdiv_qr` | `try_cdiv_qr` | ✅ |
| `mpz_cdiv_q_ui` | — | **[MISSING]** |
| `mpz_cdiv_r_ui` | — | **[MISSING]** |
| `mpz_cdiv_qr_ui` | — | **[MISSING]** |
| `mpz_cdiv_ui` | — | **[MISSING]** |
| `mpz_cdiv_q_2exp` | — | **[MISSING]** |
| `mpz_cdiv_r_2exp` | — | **[MISSING]** |

- **[MISSING]** All 6 `mpz_cdiv_*_ui` variants and both `mpz_cdiv_*_2exp` variants.
- **[PERF]** gmp-rs's `try_cdiv_qr` uses `neg_to()` + `try_fdiv_qr()` + `neg_to()` + correction. This allocates 3 temporary `Mpz` values per call. GMP's implementation uses inline arithmetic with at most one conditional subtraction.

#### B6d. Modulo (2 GMP / 2 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_mod` | `try_mod` | ✅ |
| `mpz_mod_ui` | `mod_ui` | ✅ |

#### B6e. Exact division (2 GMP / 2 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_divexact` | `try_divexact` | ✅ |
| `mpz_divexact_ui` | `try_divexact_ui` | ✅ |

- **[PERF]** `try_divexact` falls back to `tdiv_q`, which is not optimal. GMP's `mpz_divexact` uses a specialised algorithm that is linear-time when the division is known to be exact (no remainder check). gmp-rs does not exploit exactness.

#### B6f. Divisibility/Congruence (6 GMP / 6 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_divisible_p` | `divisible_p` | ✅ |
| `mpz_divisible_ui_p` | `divisible_ui` | ✅ |
| `mpz_divisible_2exp_p` | `divisible_2exp_p` | ✅ |
| `mpz_congruent_p` | `congruent_p` | ✅ |
| `mpz_congruent_ui_p` | `congruent_ui_p` | ✅ |
| `mpz_congruent_2exp_p` | `congruent_2exp_p` | ✅ |

- **[CORRECTNESS]** `congruent_p` computes `self - c` (via `try_sub`) then checks divisibility. If `self - c` would overflow capacity, `try_sub` returns `CapacityError` and `unwrap_or_else` silently substitutes `Mpz::new()` — producing a false negative. GMP would correctly compute the congruence for any-sized operands.
- **[CORRECTNESS]** `congruent_2exp_p` computes `self mod 2^bits` and `c mod 2^bits` via `fdiv_r_2exp`, then checks if their difference is divisible by `2^bits`. This is correct for non-negative values. For negative values, `fdiv_r_2exp` returns the low bits of the **magnitude** (not the two's complement low bits), so `-1 mod 16` gives `15` but `fdiv_r_2exp(-1, 4)` returns `1` (the low nibble of the magnitude `|-1| = 1`). This is **[WRONG]** — `-1 ≢ 1 (mod 16)` is `false` but gmp-rs would compute `fdiv_r_2exp(-1,4) = 1`, `fdiv_r_2exp(1,4) = 1`, difference 0 → true.

### B7. Exponentiation (5 GMP / 4 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_powm` | `try_powm` | ✅ |
| `mpz_powm_ui` | `try_powm_ui` | ✅ |
| `mpz_powm_sec` | — | **[MISSING]** Side-channel-resistant variant |
| `mpz_pow_ui` | `try_pow_ui` | ✅ |
| `mpz_ui_pow_ui` | `try_ui_pow_ui` | ✅ |

- **[MISSING]** `mpz_powm_sec` — required for cryptographic applications needing constant-time execution. Not implementable without assembly-level constant-time guarantees, which conflict with `no_unsafe`.
- **[PERF]** gmp-rs's `try_powm` uses left-to-right binary exponentiation with modular reduction at each squaring step. GMP uses the k-ary sliding window method (faster but more complex). For exponents with few set bits, the difference is negligible within the 512-bit capacity.

### B8. Root Extraction (6 GMP / 6 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_root` | `try_root` | ✅ |
| `mpz_rootrem` | `try_rootrem` | ✅ |
| `mpz_sqrt` | `isqrt` | ✅ |
| `mpz_sqrtrem` | `try_sqrtrem` | ✅ |
| `mpz_perfect_power_p` | `perfect_power_p` | ✅ |
| `mpz_perfect_square_p` | `perfect_square_p` | ✅ |

- **[PERF]** `try_root` iterates Newton's method, computing `x_pow_ui(n-1)` and `self.tdiv_q(x_pow_nm1)` at each step. Each iteration performs a full-precision exponentiation + division, both `O(n^2)` in limbs. GMP's `mpz_root` uses the same Newton method but with optimised integer exponentiation and early termination when convergence plateau is detected.
- **[CORRECTNESS]** `try_root` for `n == 0` returns zero for any input. GMP's `mpz_root(rop, op, 0)` is documented as "undefined" (division by zero internally). gmp-rs's choice of returning zero is safe but not GMP-compatible.
- **[PERF]** `perfect_power_p` for multi-limb values tries every `k` from 2 to `min(bits, 64)` using `try_root(k)` + `pow_ui(k)`. Each call to `try_root` does Newton iteration with full-precision `pow_ui` calls. For a 512-bit value, this could be up to 64 Newton iterations, each with multiple `pow_ui` calls. This is very slow.

### B9. Number Theoretic (26 GMP / 13 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_probab_prime_p` | — | **[MISSING]** |
| `mpz_nextprime` | — | **[MISSING]** |
| `mpz_prevprime` | — | **[MISSING]** |
| `mpz_gcd` | `try_gcd` | ✅ |
| `mpz_gcd_ui` | `gcd_ui` | ✅ |
| `mpz_gcdext` | `try_gcdext` | ✅ |
| `mpz_lcm` | `try_lcm` | ✅ |
| `mpz_lcm_ui` | `try_lcm_ui` | ✅ |
| `mpz_invert` | `try_invert` | ✅ |
| `mpz_jacobi` | `jacobi` | ✅ |
| `mpz_legendre` | — | **[MISSING]** (can use `jacobi`) |
| `mpz_kronecker` (6 variants) | — | **[MISSING]** |
| `mpz_remove` | `try_remove` | ✅ |
| `mpz_fac_ui` | `try_fac_ui` | ✅ |
| `mpz_2fac_ui` | — | **[MISSING]** Double factorial |
| `mpz_mfac_uiui` | — | **[MISSING]** Multi-factorial |
| `mpz_primorial_ui` | — | **[MISSING]** |
| `mpz_bin_ui` | — | **[MISSING]** Binomial with mpz n |
| `mpz_bin_uiui` | `try_bin_uiui` | ✅ |
| `mpz_fib_ui` | `try_fib_ui` | ✅ |
| `mpz_fib2_ui` | — | **[MISSING]** |
| `mpz_lucnum_ui` | — | **[MISSING]** |
| `mpz_lucnum2_ui` | — | **[MISSING]** |

- **[MISSING]** Primality testing (`probab_prime_p`, `nextprime`, `prevprime`) — 3 functions absent. These require Miller-Rabin or Baillie-PSW, which are implementable within `no_std` but are a significant amount of code.
- **[CORRECTNESS]** `try_gcd` uses Stein's binary GCD algorithm, which is correct. However, the implementation may infinite-loop if both inputs are zero: `a = b = Mpz::new()` → `a_tz = b_tz = 0` → `a = a.fdiv_q_2exp(0)` → `a` unchanged → loop: `a.len == 0` → tries to return `b` (which is `Mpz::new()` times `2^0 = 1` → `b.len == 0` → then tries to return `a` which has `len == 0` and then... actually looking at the code: if both are zero, `a.len == 0` returns `Ok(b)` where `b = Mpz::new()`. Then `b = b.try_mul_2exp(0)?` → `b` stays zero → `Ok(Mpz::new())`. This is correct: gcd(0, 0) = 0.
- **[CORRECTNESS]** `try_gcd` with inputs having different powers of 2: `a_tz` and `b_tz` are computed, the GCD's power of 2 is `min(a_tz, b_tz)`. Both inputs are then divided by their respective powers of 2 before the main loop. The loop alternately subtracts and shifts. After the loop, the result is multiplied back by `2^shift`. This is correct.
- **[CORRECTNESS]** `try_gcdext` implements the extended Euclidean algorithm. For `self = 0` or `other = 0`, the loop never executes and the initial values are returned: `gcd = |other|`, `s = 1` (if `self = 0`) or `gcd = |self|` (if `other = 0`). The sign corrections at the end handle negative inputs. This matches GMP behaviour.
- **[CORRECTNESS]** `try_invert` checks that `gcd(self, m) == 1` via `try_gcdext`. If the Bézout coefficient `s` is negative, it adds `m` to make it positive. This is correct.
- **[PERF]** `try_fac_ui` uses iterative multiplication from 2 to n. For n! within 512 bits (n ≤ 34), this is fine. For larger n, capacity error is returned. GMP's `mpz_fac_ui` uses a divide-and-conquer product tree, which is asymptotically faster but optional given the capacity limit.
- **[PERF]** `try_bin_uiui` uses the multiplicative formula with intermediate division at each step to keep the intermediate values small. This is the standard algorithm and is correct.
- **[PERF]** `try_fib_ui` uses O(n) iterative addition. GMP uses fast doubling: `F(2k) = F(k) * (2*F(k+1) - F(k))`, `F(2k+1) = F(k+1)^2 + F(k)^2`, which is O(log n). For n within 512 bits (n ≤ 512), the difference is negligible.

### B10. Comparison (8 GMP / 6 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_cmp` | `cmp` | ✅ |
| `mpz_cmp_d` | — | **[MISSING]** Compare with `f64` |
| `mpz_cmp_si` | `cmp_si` | ✅ |
| `mpz_cmp_ui` | `cmp_ui` | ✅ |
| `mpz_cmpabs` | `cmpabs` | ✅ |
| `mpz_cmpabs_d` | — | **[MISSING]** Compare |self| with `f64` |
| `mpz_cmpabs_ui` | `cmpabs_ui` | ✅ |
| `mpz_sgn` | `sgn` | ✅ |

- **[MISSING]** `cmp_d` and `cmpabs_d` — compare with `f64`. GMP provides these as macros that convert the `f64` to an `mpz` and compare. gmp-rs callers must write `self.cmp(&Mpz::from_f64(v))` — but there's no `from_f64` either.
- **[API]** gmp-rs's `cmp_si` and `cmp_ui` construct a temporary `Mpz` from the scalar value, then call `cmp`. GMP's `mpz_cmp_si` and `mpz_cmp_ui` are macros that compare directly against the scalar using `_mp_size` and `_mp_d[0]`. For single-limb values, the GMP macros avoid allocation entirely.

### B11. Logical & Bit Manipulation (12 GMP / 12 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_and` | `try_and` | ✅ |
| `mpz_ior` | `try_ior` | ✅ |
| `mpz_xor` | `try_xor` | ✅ |
| `mpz_com` | `com` | ✅ |
| `mpz_popcount` | `popcount` | ✅ |
| `mpz_hamdist` | `hamdist` | ✅ |
| `mpz_scan0` | `scan0` | ✅ |
| `mpz_scan1` | `scan1` | ✅ |
| `mpz_setbit` | `try_setbit` | ✅ |
| `mpz_clrbit` | `clrbit` | ✅ |
| `mpz_combit` | `try_combit` | ✅ |
| `mpz_tstbit` | `tstbit` | ✅ |

- **[CORRECTNESS — CRITICAL]** `try_and`, `try_ior`, `try_xor` for **negative operands** may be incorrect. gmp-rs converts both operands to two's complement by computing `~|x| + 1` for negative values, then applies the bitwise operation, then converts the result back to sign-magnitude. However:
    1. The two's complement conversion only extends sign to `self.len + 1` limbs. For a negative operand, the infinite two's complement representation has all higher limbs as `0xFF…FF`. The sign extension in `to_twos_complement` only writes one extra limb of sign extension (`0` if carry propagated, `!0` otherwise). Two negative operands of different magnitudes need **different** amounts of sign extension. The current `work_len = a_len.max(b_len)` may not be sufficient.
    2. The sign bit detection at line 1734 `let negative = (top_limb >> 63) == 1` is correct only if `work_len` includes the unambiguous sign bit. If both operands are negative and equal magnitude, the XOR result is zero, and the top limb will also be zero (MSB clear). This is detected correctly as non-negative.
    3. **Test gap**: The bitwise_ops test only tests positive values against `-1` (all-ones). It does NOT test `-2 AND -3`, `-5 XOR 3`, etc. The two's complement conversion for non-trivial negative pairs is untested.
- **[CORRECTNESS]** `popcount` returns `None` for negative values, which matches GMP's documented behaviour ("undefined for negative arguments"). Correct.
- **[CORRECTNESS]** `tstbit` for `bit >= self.len * 64` returns `false`. GMP's `mpz_tstbit` returns 0 for out-of-range bits on non-negative values, and 1 for out-of-range bits on negative values (because the infinite two's complement of a negative number has all-ones at infinity). gmp-rs's implementation returns `false` for ALL out-of-range bits, which is incorrect for negative values.
- **[CORRECTNESS]** `scan1` for zero with `start = 0` returns `None` — correct (zero has no 1 bits).
- **[CORRECTNESS]** `scan0` for a full-magnitude value like `Mpz::from_u64(u64::MAX).scan0(0)` should return 64 (the first bit beyond the magnitude). The current code at line 1801–1802 returns `(self.len * 64) as u32`, which is 64. Correct.

### B12. I/O (4 GMP / 0 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_out_str` | — | **[MISSING]** File output |
| `mpz_inp_str` | — | **[MISSING]** File input |
| `mpz_out_raw` | — | **[MISSING]** Binary output |
| `mpz_inp_raw` | — | **[MISSING]** Binary input |

- **[MISSING — COMPLETE]** Not implementable under `no_std` without a file system abstraction.

### B13. Random Numbers (5 GMP / 0 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_urandomb` | — | **[MISSING]** |
| `mpz_urandomm` | — | **[MISSING]** |
| `mpz_rrandomb` | — | **[MISSING]** |
| `mpz_random` (obsolete) | — | **[MISSING]** |
| `mpz_random2` (obsolete) | — | **[MISSING]** |

- **[MISSING — COMPLETE]** GMP bundles a Mersenne Twister PRNG. gmp-rs has no RNG dependency (a `no_std` crate should not mandate one). A standalone `gmp-rs-rand` extension crate could provide these using a caller-supplied RNG.

### B14. Integer Import/Export (2 GMP / 0 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_import` | — | **[MISSING]** |
| `mpz_export` | — | **[MISSING]** |

- **[MISSING — COMPLETE]** GMP's `mpz_import` constructs an `mpz` from a byte buffer with arbitrary endianness, word size, and nail bits. `mpz_export` does the inverse. These are implementable within `no_std` but require handling endianness and nail bits, which is substantial.

### B15. Miscellaneous (9 GMP / 9 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_fits_ulong_p` | `fits_ulong` | ✅ |
| `mpz_fits_slong_p` | `fits_slong` | ✅ |
| `mpz_fits_uint_p` | `fits_uint` | ✅ |
| `mpz_fits_sint_p` | `fits_sint` | ✅ |
| `mpz_fits_ushort_p` | `fits_ushort` | ✅ |
| `mpz_fits_sshort_p` | `fits_sshort` | ✅ |
| `mpz_odd_p` | `odd_p` | ✅ |
| `mpz_even_p` | `even_p` | ✅ |
| `mpz_sizeinbase` | `try_sizeinbase` | ✅ |

- **[CORRECTNESS]** `fits_slong` checks `val <= (i64::MIN as i128).unsigned_abs() as u64 && val != 0`. The check `val != 0` is redundant (if `sign != 0`, `val` is already known non-zero), but not incorrect.
- **[SEMANTIC]** `try_sizeinbase` returns an **upper bound** (like GMP's `mpz_sizeinbase`), but uses integer approximations for logs. The approximations are:
  - For base 10: `(bits * 30103 + 99999) / 100000` = `ceil(bits * 0.30103)`. True log10(2) = 0.30102999566…, so `0.30103 * bits` is always ≥ `log10(2) * bits`. The bound is tight for most values and off-by-at-most-1 for edge cases.
  - For other bases: uses a precomputed `log2_approx` table with 3-digit precision. The computed `ceil(bits * 1000 / log2_b)` may underestimate by 1 for some values.
  - **[CORRECTNESS]** GMP's `mpz_sizeinbase` guarantees the returned value is an **upper bound** (i.e., `mpz_sizeinbase(op, base)` >= actual number of digits). gmp-rs's approximations may produce a bound that is **1 less** than the true value for certain inputs. **This is a correctness bug** — callers using `try_sizeinbase` to allocate buffers could get a buffer that is too small.

### B16. Low-Level / Limb Access (9 GMP / 2 gmp-rs)

| GMP | gmp-rs | Status |
|---|---|---|
| `mpz_size` | `size` | ✅ |
| `mpz_getlimbn` | — | **[MISSING]** |
| `mpz_limbs_read` | — | **[MISSING]** |
| `mpz_limbs_write` | — | **[MISSING]** |
| `mpz_limbs_modify` | — | **[MISSING]** |
| `mpz_limbs_finish` | — | **[MISSING]** |
| `mpz_roinit_n` | — | **[MISSING]** |
| `MPZ_ROINIT_N` (macro) | — | **[MISSING]** |
| `_mpz_realloc` (internal) | — | **[MISSING]** |

- **[MISSING]** `mpz_getlimbn(n)` — returns the nth limb. gmp-rs exposes `self.mag` as `pub`, so callers can read limbs directly via `mpz.mag[i]` (if `i < mpz.len`). This is more raw than GMP's `mpz_getlimbn` but equally functional.
- **[MISSING]** All `mpz_limbs_*` functions for direct limb array manipulation.

---

## C. Structural & API Differences from GMP

### C1. Representation

| Aspect | GMP | gmp-rs |
|--------|-----|--------|
| Sign encoding | `_mp_size` signed: positive → positive, negative → negative, 0 → zero | `sign: i8` separate field + `len: usize` |
| Limb type | `mp_limb_t` — platform-dependent (32 or 64 bit) | `u64` — unconditionally 64-bit |
| Allocation | Dynamic via `mpz_realloc`; `_mp_alloc` tracks capacity | Fixed `[u64; 8]` — compile-time constant |
| Zero representation | `_mp_size == 0`, `_mp_d` unchanged | `len == 0`, `sign == 0`, mag array ignored |
| Public field access | Implementation-private; access via API | `sign`, `len`, `mag` are `pub` — **[API]** direct field access allowed |

- **[API]** gmp-rs exposes `pub sign: i8`, `pub len: usize`, `pub mag: [u64; MPZ_MAX_LIMBS]`. GMP treats its `_mp_*` fields as implementation details (though C code can access them). Callers can construct malformed `Mpz` values by writing directly to fields. A `new()` or constructor should be the only way to create valid values.
- **[no_alloc limitation]** 512-bit maximum. Operations that exceed this silently return `CapacityError` or, in some paths, silently substitute zero via `unwrap_or_else(|_| Mpz::new())`.

### C2. Error Handling

| Scenario | GMP | gmp-rs |
|----------|-----|--------|
| Out of memory | `abort()` | `CapacityError` returned from `try_*` methods |
| Division by zero | `abort()` | `panic!("gmp-rs: division by zero")` in release mode (was `debug_assert` in v0.1.0) |
| Invalid string input | Return `-1` | Return `ParseError` |
| Capacity overflow | N/A (infinite precision) | Return `CapacityError` or silently return zero |

### C3. Method Naming & Calling Convention

| GMP | gmp-rs | Difference |
|-----|--------|------------|
| `mpz_add(rop, op1, op2)` | `op1.try_add(op2)` | Return-value style vs destination-first |
| `mpz_set_ui(rop, val)` | `rop.set_ui(val)` | Same |
| `mpz_get_ui(op)` | `op.get_ui()` | Same |
| `mpz_neg(rop, op)` | `op.neg_to()` (non-mutating) or `op.neg()` (mutating) | Separate mutating/non-mutating variants |
| `mpz_sizeinbase(op, base)` | `op.try_sizeinbase(base)` | Fallible due to base validation |

### C4. `no_std`-Specific Deviations

1. **No `FILE*` I/O** — cannot implement `mpz_out_str`/`mpz_inp_str` without a file abstraction.
2. **No `f64::powf`** — `perfect_power_p` for single-limb values uses integer binary search instead of floating-point root approximation.
3. **No `f64::ldexp`/`f64::frexp`** — `get_d_2exp` manually normalises with a `while` loop.
4. **No `std::error::Error`** — `CapacityError` and `ParseError` do not implement `std::error::Error` (which requires `std`). They implement `core::fmt::Debug` and `core::fmt::Display` (via `Debug`).

---

## D. Risk Register

| ID | Severity | Category | Issue |
|---|---|---|---|
| **R01** | **HIGH** | Bitwise §B11 | `tstbit` returns `false` for out-of-range bits on **negative** values; GMP returns `1` (infinite two's complement all-ones). |
| **R02** | **HIGH** | Bitwise §B11 | `try_and`, `try_ior`, `try_xor` for two negative operands are untested; the sign-extension logic may be incorrect for negative pairs with different limb counts. |
| **R03** | **HIGH** | Congruence §B6f | `congruent_2exp_p` uses magnitude low bits instead of two's complement low bits; incorrect for negative values. |
| **R04** | **MEDIUM** | Division §B6 | Div-by-zero panics (not checked at compile time and not returned as `Result`). Safety-critical callers must validate divisors. |
| **R05** | **MEDIUM** | Sizeinbase §B15 | `try_sizeinbase` may under-estimate by 1 due to integer log approximation. Buffer allocations based on this bound may be too small. |
| **R06** | **MEDIUM** | Conversion §B4 | `get_d()` for values near `f64::MAX` silently loses precision (no overflow/Infinity detection). |
| **R07** | **MEDIUM** | Arithmetic §B5 | `try_sub` double-clones `other` when `other.sign == 0` (wasted work). |
| **R08** | **MEDIUM** | Number theory §B9 | `try_gcd(0, 0)` returns `0` — GMP returns `0` (undocumented but observed). Same behaviour. |
| **R09** | **LOW** | Misc §C1 | Public fields allow invalid `Mpz` construction by external code. |
| **R10** | **LOW** | Perfect power §B8 | `try_root(n=0)` returns 0; GMP's behaviour is undefined. |
| **R11** | **LOW** | Ceiling div §B6c | `try_cdiv_qr` allocates 3 temporaries per call via `neg_to`/`fdiv_qr`/`neg_to`. |
| **R12** | **LOW** | Bitwise §B11 | `com()` uses `try_add(&one)` — if `self` is at capacity and `try_add` fails, silently returns `Mpz::new()`. |

---

## E. Fixed Issues (from v0.1.0 → v0.2.1)

| Previous ID | Issue | Status in v0.2.1 |
|---|---|---|
| R01 (v0.1) | `mag_divmod` bit-by-bit, 2000 allocs per 256-bit div | **[PARTIALLY FIXED]** Single-limb divisors now use direct scalar division (fast path). Multi-limb divisors still use bit-by-bit. |
| R02 (v0.1) | Division by zero only `debug_assert`-guarded | **[FIXED]** Now panics in both debug and release. |
| R03 (v0.1) | Parser silently corrupted `"12O34"` → `1234` | **[FIXED]** Now returns `ParseError::InvalidInput`. |
| — | No `no_alloc` | **[FIXED]** Fixed `[u64; 8]` array, no `alloc` dependency. |
| — | No floor/ceil division | **[FIXED]** All variants implemented. |
| — | No bitwise ops | **[FIXED]** All 12 bit functions implemented. |
| — | No number theory | **[FIXED]** 13 number theory functions implemented. |
| — | No operator traits | **[FIXED]** `AddAssign`, `SubAssign`, `MulAssign`. |

---

## F. What Remains Unimplementable Under Constraints

The following GMP `mpz_*` functions **cannot** be implemented within `no_unsafe, no_std, no_alloc`:

### F1. Requires a file system or I/O (no_std):
- `mpz_out_str`, `mpz_inp_str` (FILE* I/O)
- `mpz_out_raw`, `mpz_inp_raw` (FILE* binary I/O)

### F2. Requires an external RNG (no_std):
- `mpz_urandomb`, `mpz_urandomm`, `mpz_rrandomb`

### F3. Requires C varargs:
- `mpz_inits`, `mpz_clears`, `mpz_array_init`

### F4. Requires dynamic allocation (no_alloc):
- `mpz_init2`, `mpz_realloc2` (preallocation hints not applicable to fixed arrays)
- `mpz_limbs_write`, `mpz_limbs_modify`, `mpz_limbs_finish` (dynamic resizing APIs)

### F5. Requires excessive code volume or unsafe:
- `mpz_powm_sec` (side-channel resistance requires constant-time ops, which need assembly or platform intrinsics)
- `mpz_probab_prime_p`, `mpz_nextprime`, `mpz_prevprime` (Miller-Rabin is implementable but large; ~300-500 lines)
- `mpz_kronecker` (6 variants — each ~50 lines but trivially derived from `jacobi`)

### F6. Requires types that don't exist in gmp-rs:
- `mpz_set_q` (requires `mpq_t` rational type)
- `mpz_set_f` (requires `mpf_t` floating-point type)
- `mpz_bin_ui` (requires `mpz` as first argument — gmp-rs has `try_bin_uiui` for `u32,u32`)

---

## G. Performance Characteristics vs GMP

| Operation | GMP (typical) | gmp-rs (measured asymptotics) |
|-----------|---------------|-------------------------------|
| Add/Sub | O(n) limbs, in-place reuse | O(n) limbs, allocates new result |
| Multiply | O(n^1.585) Karatsuba for large n, O(n^2) schoolbook for small n | O(n^2) schoolbook only |
| Division (1-limb divisor) | O(n) | O(n) — uses direct scalar division |
| Division (multi-limb) | O(n^2) schoolbook or O(n log n) Barrett | O(n² in **bits**) — bit-by-bit × 64 |
| GCD | O(n^2) binary or O(n²) Euclidean | O(n^2) binary GCD |
| powm | O(log exp × n^2) sliding window | O(log exp × n^2) binary L-to-R |
| isqrt | O(n^2) Newton with fast division | O(n² × log n) Newton with slow division |

The **critical performance bottleneck** is multi-limb division (R01 in v0.1.0). All callers of `tdiv_qr` with `|d| >= 2` limbs trigger the bit-by-bit algorithm (512 iterations for 8-limb values). This affects `isqrt` (which calls `tdiv_q` per Newton iteration), `try_root`, `try_gcd` (via `try_sub` which triggers `add_signed` → `mag_sub` — wait, `try_sub` doesn't trigger division, it uses `mag_sub_len` which is linear. But `try_gcdext` calls `tdiv_qr` per loop iteration), and `try_remove`.

Single-limb division (`d.len == 1`) is fast (schoolbook across limbs, O(n)).

---

## H. Comparison with GMP 7.0 (future)

GMP 7.0 (unreleased) is expected to add:
- `mpz_prevprime` (already in GMP 6.x)
- `mpz_2fac_ui` (double factorial)
- `mpz_mfac_uiui` (multi-factorial)
- `mpz_primorial_ui`
- `mpz_oddfac_1` (internal)

None of these change gmp-rs's gap analysis — they would be additions on top of the existing gap.

---

## I. Recommendations for Critical Safety Infrastructure

For a crate to be suitable as "critical safety infrastructure":

1. **Fix R01 (bitwise sign extension)** — `tstbit` must return 1 for out-of-range bits on negative values. `bitwise_op` needs thorough testing with negative pairs.
2. **Fix R03 (congruent_2exp_p)** — Must use two's complement low bits, not magnitude low bits, for negative values.
3. **Add `from_d`** — Conversion from `f64` is essential for interop with floating-point systems.
4. **Fix R05 (try_sizeinbase)** — Ensure the returned bound is always ≥ the actual number of digits. Use `ceil(k * log(2) / log(base))` with rigorous integer arithmetic.
5. **Evaluate multi-limb division** — Replace bit-by-bit with Knuth's Algorithm D (schoolbook, O(n²) in limbs) to improve performance by a factor of ~64 for worst-case values.
6. **Add preconditions** — Document that `divisible_ui`, `tdiv_q_ui`, `tdiv_ui`, etc. panic on `d == 0`. Or add `try_` variants that return `Result`.
7. **Seal construction** — Make `sign`, `len`, `mag` private and provide safe constructors, or add runtime invariant checking in `trim()`.
