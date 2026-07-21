#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    if data.len() < 16 {
        return;
    }

    let a_val = u64::from_ne_bytes([
        data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7],
    ]);
    let b_val = u64::from_ne_bytes([
        data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15],
    ]);

    let a = gmp_rs::Mpz::from_u64(a_val);
    let b = gmp_rs::Mpz::from_u64(b_val);

    // Bitwise ops
    let _ = a.try_and(&b);
    let _ = a.try_ior(&b);
    let _ = a.try_xor(&b);
    let _ = a.com();
    let _ = a.popcount();
    let _ = a.hamdist(&b);

    // Scan
    if data.len() > 16 {
        let start = data[16] as u32;
        let _ = a.scan0(start);
        let _ = a.scan1(start);
    }

    // tstbit
    for bit in 0..64 {
        let _ = a.tstbit(bit);
    }

    // setbit, clrbit, combit (on a clone to avoid mutating original)
    let mut c = a.clone();
    for bit in [0u32, 1, 7, 63, 64, 127] {
        let _ = c.try_setbit(bit);
        c.clrbit(bit);
        let _ = c.try_combit(bit);
    }

    // tdiv_2exp variants
    if data.len() > 17 {
        let k = data[17] as u32;
        let _ = a.tdiv_q_2exp(k);
        let _ = a.tdiv_r_2exp(k);
        let _ = a.fdiv_q_2exp(k);
        let _ = a.fdiv_r_2exp(k);
        let _ = a.cdiv_q_2exp(k);
        let _ = a.cdiv_r_2exp(k);
    }
});
