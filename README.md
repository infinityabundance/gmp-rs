# gmp-rs

[![crates.io](https://img.shields.io/crates/v/gmp-rs.svg)](https://crates.io/crates/gmp-rs)
[![Documentation](https://docs.rs/gmp-rs/badge.svg)](https://docs.rs/gmp-rs)

A **no-unsafe**, **`no_std`**, **`no_alloc`** pure-Rust arbitrary-precision signed integer library
(`Mpz`), faithful to the GMP `mpz_*` operations.

## Guarantees

- **Zero `unsafe` code** — `#![forbid(unsafe_code)]` enforced at compile time.
- **`no_std`** — no standard library dependency.
- **`no_alloc`** — zero heap allocations. Fixed-capacity limb storage (`[u64; 8]`, 512 bits/~154 decimal digits).
  Operations that would exceed capacity return `Err(CapacityError)`.

## Quick start

```rust
use gmp_rs::Mpz;

let a = Mpz::from_decimal_str("123456789012345678901234567890").unwrap();
let b = Mpz::from_decimal_str("987654321098765432109876543210").unwrap();

let sum = a.try_add(&b).unwrap();
let product = a.try_mul(&b).unwrap();

println!("{}", sum);       // 1111111110111111111011111111100
println!("{}", product);   // 1219326311370217952261850327336...
```

## License

Licensed under the **GNU Lesser General Public License v3.0 or later** ([LGPL-3.0-or-later](COPYING.LESSER)).
