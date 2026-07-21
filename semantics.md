# gmp-rs Formal Semantics

> **Version:** 0.3.0
> **Date:** 2026-07-21
> **Constraint envelope:** `no_unsafe` ✓, `no_std` ✓, `no_alloc` ✓

## 1. Mathematical Model

### 1.1 The `Mpz` Type

An `Mpz` value represents the integer:

```
val(Mpz{sign, len, mag}) = sign × Σ(mag[i] × 2^(64×i) for i in 0..len)
```

with the **invariant**: `(len = 0) ⇔ (sign = 0)`.

That is:
- **Zero** is represented uniquely as `(sign=0, len=0)` (the `mag` array is ignored).
- **Non-zero** values have `sign ∈ {-1, 1}`, `len ∈ {1, …, MPZ_MAX_LIMBS}`,
  `mag[len-1] ≠ 0` (no trailing zero limb).
- **Overflow limit:** `|val| < 2^(MPZ_MAX_LIMBS × 64)` — the value must fit
  in `MPZ_MAX_LIMBS` 64-bit limbs.

### 1.2 Capacity

The representable integer set is:

```
ℤ_cap = { n ∈ ℤ : |n| < 2^512 }   (where 512 = MPZ_MAX_LIMBS × 64)
```

All arithmetic operations are **partial functions** on ℤ × ℤ → ℤ_cap.
If the exact mathematical result falls outside ℤ_cap, the operation
returns `Err(CapacityError)`.

### 1.3 Zero Representation

The zero value is unique — there is no negative zero:

```
∀ m : Mpz,  val(m) = 0  ⇔  m.sign = 0 ∧ m.len = 0
```

---

## 2. Operation Semantics

### 2.1 Construction

| Function | Semantics | Domain |
|----------|-----------|--------|
| `new() → Mpz` | Returns `z` where `val(z) = 0` | Total |
| `from_u64(v) → Mpz` | Returns `z` where `val(z) = v` | Total (v ∈ [0, 2^64)) |
| `from_i64(v) → Mpz` | Returns `z` where `val(z) = v` | Total (v ∈ [-2^63, 2^63)) |
| `from_u128(v) → Mpz` | Returns `z` where `val(z) = v` | Total (v ∈ [0, 2^128)) |
| `from_i128(v) → Mpz` | Returns `z` where `val(z) = v` | Total (v ∈ [-2^127, 2^127)) |
| `from_d(v) → Result(Mpz, CapacityError)` | Returns `z` where `val(z) = trunc(v)` (truncation toward zero) | Partial: `Err` if `\|v\| ≥ 2^512` or v is NaN |
| `from_decimal_str(s) → Result(Mpz, ParseError)` | Returns `z` where `val(z) = decimal_parse(s)` | Partial: `Err(InvalidInput)` if s contains non-digit chars; `Err(CapacityOverflow)` if result exceeds ℤ_cap |

### 2.2 Conversion

| Function | Semantics |
|----------|-----------|
| `get_ui(m) → u64` | `val(m) mod 2^64` (low 64 bits of absolute value) |
| `get_si(m) → i64` | Low 64 bits of `val(m)` as signed integer |
| `get_d(m) → f64` | `round_to_f64(val(m), toward_zero)` (may lose precision) |
| `get_d_2exp(m) → Option(f64, i64)` | Returns `(mantissa, exponent)` such that `val(m) = mantissa × 2^exponent` with 0.5 ≤ \|mantissa\| < 1. Returns `None` for zero. |
| `to_i128(m) → Option(i128)` | `val(m)` if `\|val(m)\| ≤ 2^127`, else `None` |
| `write_decimal_buf(m, buf) → usize` | Writes decimal representation of `val(m)` into `buf`, returns bytes written |

### 2.3 Sign and Comparison

| Function | Semantics |
|----------|-----------|
| `sgn(m) → i32` | `sign(val(m))` — returns -1, 0, or 1 |
| `cmp(a, b) → Ordering` | `compare(val(a), val(b))` — standard total order on ℤ |
| `cmpabs(a, b) → Ordering` | `compare(\|val(a)\|, \|val(b)\|)` |
| `cmp_si(a, v) → Ordering` | `compare(val(a), v)` |
| `cmp_ui(a, v) → Ordering` | `compare(val(a), v)` |
| `cmp_d(a, v) → Ordering` | `compare(val(a), trunc(v))` |
| `odd_p(m) → bool` | `val(m) mod 2 ≠ 0` |
| `even_p(m) → bool` | `val(m) mod 2 = 0` |

### 2.4 Arithmetic

| Function | Semantics |
|----------|-----------|
| `try_add(a, b) → Ok(z)` iff `val(z) = val(a) + val(b)` and `val(z) ∈ ℤ_cap` |
| `try_sub(a, b) → Ok(z)` iff `val(z) = val(a) - val(b)` and `val(z) ∈ ℤ_cap` |
| `try_mul(a, b) → Ok(z)` iff `val(z) = val(a) × val(b)` and `val(z) ∈ ℤ_cap` |
| `try_mul_2exp(a, k) → Ok(z)` iff `val(z) = val(a) × 2^k` and `val(z) ∈ ℤ_cap` |
| `try_pow_ui(a, n) → Ok(z)` iff `val(z) = val(a)^n` and `val(z) ∈ ℤ_cap` |
| `try_powm(a, e, m) → Ok(z)` iff `val(z) = val(a)^val(e) mod val(m)` and `val(z) ∈ ℤ_cap` |
| `try_addmul(a, b, c)` sets `val(a) := val(a) + val(b) × val(c)` |
| `try_submul(a, b, c)` sets `val(a) := val(a) - val(b) × val(c)` |

**CapacityError convention:** For `try_add`, `try_sub`, `try_mul`, the error is
returned when the exact result exceeds ℤ_cap. For `try_pow_ui` and `try_powm`,
intermediate values may also exceed ℤ_cap during computation.

### 2.5 Division Semantics

GMP defines three rounding modes for integer division. gmp-rs implements all three:

| Rounding | Condition | Quotient q | Remainder r |
|----------|-----------|------------|-------------|
| **Truncate** (toward zero) | `a = q × d + r` | `q = trunc(a/d)` | `0 ≤ \|r\| < \|d\|`, `sign(r) = sign(a)` |
| **Floor** (toward −∞) | `a = q × d + r` | `q = floor(a/d)` | `0 ≤ r < \|d\|` |
| **Ceiling** (toward +∞) | `a = q × d + r` | `q = ceil(a/d)` | `-\|d\| < r ≤ 0` with `sign(r) ≠ sign(d)` |

**Identity:** For all rounding modes: `a = q × d + r`.

**Special cases:**
- Division by zero: **panics** (this is a precondition violation; callers must validate).
- `0 / d = 0` with remainder `0` for all `d ≠ 0`.
- `a / 1 = a` with remainder `0`.
- `a / a = 1` with remainder `0` (for `a ≠ 0`).

**Floor and modulo relationship:** `try_mod(a, d)` returns the floor remainder,
so `try_mod(a, d) = try_fdiv_r(a, d)`.

| Division function | Rounding | Return type |
|---|---|---|
| `tdiv_qr` | Truncate | `(Mpz, Mpz)` — infallible for valid inputs |
| `tdiv_q` | Truncate | `Mpz` |
| `tdiv_r` | Truncate | `Mpz` |
| `try_fdiv_qr` | Floor | `Result<(Mpz, Mpz), CapacityError>` |
| `try_cdiv_qr` | Ceiling | `Result<(Mpz, Mpz), CapacityError>` |
| `try_mod` | Floor (non-negative remainder) | `Result<Mpz, CapacityError>` |
| `tdiv_q_ui` | Truncate by u64 | `Mpz` |
| `tdiv_ui` | Truncate remainder as u64 | `u64` |
| `fdiv_ui` | Floor remainder as u64 | `u64` |
| `cdiv_ui` | Ceiling remainder as u64 | `u64` |

### 2.6 Division by Power of 2

| Function | Semantics |
|----------|-----------|
| `fdiv_q_2exp(a, k)` | `floor(val(a) / 2^k)` |
| `fdiv_r_2exp(a, k)` | `val(a) mod 2^k` (floor, non-negative) |
| `tdiv_q_2exp(a, k)` | `trunc(val(a) / 2^k)` (same as floor for non-negative `a`) |
| `tdiv_r_2exp(a, k)` | `trunc(val(a) mod 2^k)` (same as floor for non-negative `a`) |
| `cdiv_q_2exp(a, k)` | `ceil(val(a) / 2^k)` |
| `cdiv_r_2exp(a, k)` | `val(a) - ceil(val(a) / 2^k) × 2^k` |

### 2.7 Divisibility and Congruence

| Function | Semantics |
|----------|-----------|
| `divisible_p(a, d)` | `d ≠ 0 ∧ d \| a` |
| `divisible_ui(a, d)` | `d ≠ 0 ∧ d \| a` |
| `divisible_2exp_p(a, k)` | `2^k \| a` |
| `congruent_p(a, c, d)` | `a ≡ c (mod d)`, i.e. `d \| (a - c)` |
| `congruent_ui_p(a, c, d)` | `a ≡ c (mod d)` with `c, d ∈ ℕ` |
| `congruent_2exp_p(a, c, k)` | `a ≡ c (mod 2^k)`, using two's complement representation |

**Note:** `congruent_2exp_p` for negative values interprets `a` and `c` in two's
complement with infinite sign extension. This matches GMP's behaviour where
`-1 ≡ 2^k - 1 (mod 2^k)` for all `k ≥ 0`.

### 2.8 Bitwise Operations

Binary bitwise operations (`and`, `ior`, `xor`) are defined on the **infinite two's
complement** representation:

1. Convert both operands to two's complement with infinite sign extension
   (positive → infinite zero extension, negative → infinite one extension).
2. Apply the bitwise operation pointwise across all bits.
3. Convert the result back to sign-magnitude.

| Function | Semantics |
|----------|-----------|
| `try_and(a, b)` | Bitwise AND of two's complement representations |
| `try_ior(a, b)` | Bitwise inclusive OR |
| `try_xor(a, b)` | Bitwise exclusive OR |
| `com(a)` | `-a - 1` (one's complement, equivalently bitwise NOT of two's complement) |
| `popcount(a)` | Number of 1 bits in two's complement representation. Returns `None` for negative `a`. |
| `hamdist(a, b)` | Number of bits differing between `a` and `b` in two's complement. Returns `None` for mixed signs. |
| `tstbit(a, k)` | `k`-th bit of two's complement representation (0-indexed). For negative `a` and `k` ≥ bitlen, returns `1`. |
| `scan0(a, start)` | Position of first 0 bit at or after `start` |
| `scan1(a, start)` | Position of first 1 bit at or after `start`, or `None` |
| `try_setbit(a, k)` | Set bit `k` to 1 |
| `clrbit(a, k)` | Set bit `k` to 0 |
| `try_combit(a, k)` | Toggle bit `k` |

### 2.9 Root Extraction

| Function | Semantics |
|----------|-----------|
| `isqrt(a)` | `floor(sqrt(a))` for `a ≥ 0`; returns `0` for `a ≤ 0` |
| `try_sqrtrem(a)` | `(floor(sqrt(a)), a - floor(sqrt(a))²)` for `a ≥ 0` |
| `try_root(a, n)` | `floor(a^(1/n))` for `a ≥ 0, n ≥ 1`; returns `0` for `a ≤ 0, n = 0` |
| `try_rootrem(a, n)` | `(floor(a^(1/n)), a - floor(a^(1/n))^n)` |
| `perfect_square_p(a)` | `∃z : z² = a` |
| `perfect_power_p(a)` | `∃z ∈ ℤ, ∃k ≥ 2 : z^k = a` |

### 2.10 Number Theoretic Functions

| Function | Semantics |
|----------|-----------|
| `try_gcd(a, b)` | greatest common divisor of `a` and `b`, always non-negative |
| `gcd_ui(a, v)` | `gcd(a, v)` as `u64` |
| `try_gcdext(a, b)` | `(g, s, t)` where `g = gcd(a, b) = a·s + b·t` |
| `try_lcm(a, b)` | least common multiple, always non-negative |
| `try_lcm_ui(a, v)` | `lcm(a, v)` |
| `try_invert(a, m)` | `x` where `a·x ≡ 1 (mod m)`, or `Err(())` if not invertible |
| `jacobi(a, n)` | Jacobi symbol `(a/n)` for odd `n > 0` |
| `try_fac_ui(n)` | `n!` |
| `try_bin_uiui(n, k)` | binomial coefficient `C(n, k)` |
| `try_fib_ui(n)` | `F_n` (Fibonacci: F_0 = 0, F_1 = 1) |
| `try_fib2_ui(n)` | `(F_n, F_{n-1})` |
| `try_lucnum_ui(n)` | `L_n` (Lucas: L_0 = 2, L_1 = 1) |
| `try_lucnum2_ui(n)` | `(L_n, L_{n-1})` |
| `try_probab_prime_p(n, reps)` | Miller-Rabin: returns 2 (prime), 1 (probable prime), 0 (composite) |
| `try_remove(a, f)` | `(a / f^cnt, cnt)` where `cnt` = multiplicity of factor `f` in `a` |

### 2.11 Constant-Time Operations

When compiled with `const_time` feature:

| Function | Guarantee |
|----------|-----------|
| `ct_add(a, b)` | Same as `try_add` but iteration over all `MPZ_MAX_LIMBS` limbs. Execution time is independent of the values (but may leak whether overflow occurred via `Err`/`Ok`). |
| `ct_sub(a, b)` | Same guarantee as `ct_add`. |
| `ct_cmp(a, b)` | Returns `Ordering` without branching on value bits. Execution time is independent of values. |
| `ct_select(a, b, bit)` | Returns `a` if `bit = 0`, `b` if `bit = 1` using bitwise selection. Constant-time in the values of `a` and `b`. |

**Limitation:** `ct_add` and `ct_sub` may still leak the CapacityError distinction
(overflow vs no overflow) via timing of the `Err` vs `Ok` path.

---

## 3. Error Model

### 3.1 Error Types

```
CapacityError  — The exact mathematical result exceeds the fixed 512-bit capacity.
                 The operation was not performed; the state is unchanged.

ParseError::InvalidInput — The input string could not be parsed as a valid
                           decimal integer with optional sign.

ParseError::CapacityOverflow — The input string parsed correctly but the value
                               exceeds the 512-bit capacity.
```

### 3.2 Panic Conditions

The following operations panic on precondition violation:

| Operation | Panic condition |
|-----------|----------------|
| `tdiv_qr`, `tdiv_q`, `tdiv_r` | Division by zero (`d = 0`) |
| `try_divexact` | Division by zero (inherits from `tdiv_q`) |
| All `mag_divmod_u64` paths | Division by zero (`d = 0`) |

**Safety-critical callers must validate that divisors are non-zero before calling**
these functions.

### 3.3 Silent Degradation

The following code paths silently substitute zero on capacity overflow:

- `AddAssign`, `SubAssign`, `MulAssign` — `unwrap_or_else(|_| Mpz::new())`
- `congruent_p` — `unwrap_or_else(|_| Mpz::new())` on the difference

This is a **known limitation** — callers using the `try_*` variants get explicit
error handling.

---

## 4. Representation Invariants

### 4.1 Internal State

An `Mpz` value always satisfies:

```
0 ≤ len ≤ MPZ_MAX_LIMBS
sign ∈ {-1, 0, 1}
(len = 0) ⇔ (sign = 0)
∀i ∈ [len, MPZ_MAX_LIMBS): mag[i] = 0   (outside the public limbs() view)
len > 0 ⇒ mag[len-1] ≠ 0                 (no trailing zero limb)
```

### 4.2 Serialisation

`write_decimal_buf` writes the decimal representation without leading zeros.
`from_decimal_str` accepts an optional leading `+` or `-` followed by ASCII
decimal digits. The round-trip property holds:

```
∀m : Mpz,  m ≠ −0   (no negative zero)
∀m : Mpz,  from_decimal_str(write_decimal_buf(m)) = Ok(m)
∀s valid decimal string:  write_decimal_buf(from_decimal_str(s)) = s (canonical)
```

---

## 5. GMP Compatibility

gmp-rs implements ~105 of ~154 GMP `mpz_*` functions (~68% coverage). The
remaining functions are intentionally absent due to:

- **Requires `std` (file I/O):** `mpz_out_str`, `mpz_inp_str`, `mpz_out_raw`, `mpz_inp_raw`
- **Requires unsafe (pointer manipulation):** `mpz_limbs_read`, `mpz_limbs_write`, `mpz_roinit_n`
- **Requires dynamic allocation (no_alloc):** `mpz_init2`, `mpz_realloc2`
- **Requires types not in crate:** `mpz_set_q` (needs `mpq_t`), `mpz_set_f` (needs `mpf_t`)
- **Not implemented (code volume):** `mpz_nextprime`, `mpz_prevprime`, `mpz_mfac_uiui`
- **C-specific (not applicable):** `mpz_inits`, `mpz_clears` (vararg), `mpz_clear` (Rust Drop)

See `capability_map()` for the complete per-function status table.
