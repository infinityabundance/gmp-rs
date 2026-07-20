# gmp-rs

[![crates.io](https://img.shields.io/crates/v/gmp-rs.svg)](https://crates.io/crates/gmp-rs)
[![Documentation](https://docs.rs/gmp-rs/badge.svg)](https://docs.rs/gmp-rs)

A **no-unsafe**, **`no_std`** pure-Rust arbitrary-precision signed integer library (`Mpz`), faithful to the
GMP `mpz_*` operations. Built for use as the integer foundation in higher-level decimal/numeric ports.

## Features

- **Zero `unsafe` code** — `#![forbid(unsafe_code)]` enforced at compile time.
- **`no_std`** — only depends on the `alloc` crate; no standard library required.
- **Sign–magnitude representation** — little-endian base-2⁶⁴ limbs, no trailing zero limb.
- **Arithmetic:** addition, subtraction, multiplication, truncated division, shifts, powers, integer square root.
- **Conversions:** to/from `u64`, `i64`, `u128`, `i128`, and decimal strings.
- **Faithful to GMP semantics** — matches the `mpz_*` surface used by GnuCOBOL's numeric port.

## Quick start

```rust
use gmp_rs::Mpz;

let a = Mpz::from_decimal_string("123456789012345678901234567890");
let b = Mpz::from_decimal_string("987654321098765432109876543210");

let sum = &a + &b;
let product = a.mul(&b);

println!("{}", sum);       // 1111111110111111111011111111100
println!("{}", product);   // 1219326311370217952261850327336...
```

## License

Licensed under the **GNU Lesser General Public License v3.0 or later** ([LGPL-3.0-or-later](COPYING.LESSER)).
