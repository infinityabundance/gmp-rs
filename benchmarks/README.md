# Benchmarks

## Rust benchmarks (criterion)

```sh
# gmp-rs only (no external dependencies)
cargo bench --bench arithmetic

# gmp-rs vs raw GMP C API via FFI (requires libgmp-dev)
cargo bench --features gmp_cross_check --bench gmp_comparison
```

## C benchmarks (raw GMP)

For a true apples-to-apples comparison against GMP's C API:

```sh
# Compile the C benchmark (requires libgmp-dev)
gcc -O2 -o benchmarks/gmp_c_bench benchmarks/gmp_c_bench.c -lgmp

# Run it
./benchmarks/gmp_c_bench
```

### Expected results

The C benchmark measures raw GMP throughput with no overhead from string conversion or Rust FFI. The Rust `gmp_comparison` benchmark measures gmp-rs throughput via its native API. The difference is:

- **gmp-rs overhead**: gmp-rs constructs `Mpz` values from scratch for each operation (allocates + frees the stack array), while the C benchmark reuses pre-initialised `mpz_t` values.
- **FFI overhead**: The Rust FFI benchmark converts `Mpz` to strings, passes them through C FFI, and GMP parses them back — this is intentionally a worst-case comparison. The `arithmetic` benchmark (gmp-rs only) is the appropriate measure for real gmp-rs usage.

For embedded/safety-critical use, the absolute performance of gmp-rs matters more than the comparison to GMP, since GMP cannot be used in those environments.
