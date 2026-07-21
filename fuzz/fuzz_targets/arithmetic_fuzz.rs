#![no_main]

use libfuzzer_sys::fuzz_target;

// Fuzz arithmetic operations.
fuzz_target!(|data: &[u8]| {
    if data.len() < 16 {
        return;
    }

    // Construct two random 64-bit values from the fuzz input.
    let a_val = u64::from_ne_bytes([
        data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
    ]);
    let b_val = u64::from_ne_bytes([
        data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15],
    ]);

    let a = gmp_rs::Mpz::from_u64(a_val);
    let b = gmp_rs::Mpz::from_u64(b_val);

    // ── add ──
    let _sum = a.try_add(&b);

    // ── sub ──
    let _diff = a.try_sub(&b);

    // ── mul ──
    let _prod = a.try_mul(&b);

    // ── tdiv_qr (only if b != 0) ──
    if b_val != 0 {
        let (_q, _r) = a.tdiv_qr(&b);
    }

    // ── floor division ──
    if b_val != 0 {
        let _ = a.try_fdiv_qr(&b);
    }

    // ── ceiling division ──
    if b_val != 0 {
        let _ = a.try_cdiv_qr(&b);
    }

    // ── mod ──
    if b_val != 0 {
        let _ = a.try_mod(&b);
    }

    // ── bitwise ──
    let _ = a.try_and(&b);
    let _ = a.try_ior(&b);
    let _ = a.try_xor(&b);

    // ── shift ──
    if data.len() > 16 {
        let shift = data[16] as u32;
        let _ = a.try_mul_2exp(shift);
    }

    // ── comparison ──
    let _ = a.cmp(&b);
    let _ = a.cmpabs(&b);

    // ── gcd ──
    let _ = a.try_gcd(&b);

    // ── lcm ──
    let _ = a.try_lcm(&b);

    // ── pow (small exponents only) ──
    if data.len() > 17 {
        let exp = data[17] as u32 % 10;
        let _ = a.try_pow_ui(exp);
    }

    // ── isqrt (non-negative only) ──
    let _ = a.isqrt();
    let a_pos = gmp_rs::Mpz::from_u64(a_val | 1); // ensure non-zero for division
    let _ = a_pos.try_root(data[17] as u32 % 10 + 1);

    // ── from_d ──
    let f = f64::from_bits(u64::from_ne_bytes([
        data[0] ^ data[8],
        data[1] ^ data[9],
        data[2] ^ data[10],
        data[3] ^ data[11],
        data[4] ^ data[12],
        data[5] ^ data[13],
        data[6] ^ data[14],
        data[7] ^ data[15],
    ]));
    let _ = gmp_rs::Mpz::from_d(f);
});
