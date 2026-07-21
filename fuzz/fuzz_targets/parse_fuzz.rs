#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Feed arbitrary bytes into from_decimal_str.
    // Since from_decimal_str rejects non-digit chars, this primarily
    // tests the error path and the digit-only path.
    if let Ok(s) = core::str::from_utf8(data) {
        let _ = gmp_rs::Mpz::from_decimal_str(s);
    }

    // Also test from_d with arbitrary bit patterns (NaN, inf, etc.)
    if data.len() >= 8 {
        let bits = u64::from_ne_bytes([
            data[0], data[1], data[2], data[3],
            data[4], data[5], data[6], data[7],
        ]);
        let f = f64::from_bits(bits);
        let _ = gmp_rs::Mpz::from_d(f);
    }

    // Test try_import with arbitrary byte sequences
    if data.len() >= 4 {
        let _ = gmp_rs::Mpz::try_import(
            1, gmp_rs::Endian::Little, 1, gmp_rs::Endian::Little, data,
        );
        let _ = gmp_rs::Mpz::try_import(
            1, gmp_rs::Endian::Big, 1, gmp_rs::Endian::Big, data,
        );
    }
});
