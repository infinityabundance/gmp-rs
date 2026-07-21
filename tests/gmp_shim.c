#include <gmp.h>
#include <string.h>

void gmp_add(const char *a_str, const char *b_str, char *out_buf, size_t out_len) {
    mpz_t a, b, r;
    mpz_init(a); mpz_init(b); mpz_init(r);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);
    mpz_add(r, a, b);
    mpz_get_str(out_buf, 10, r);
    mpz_clear(a); mpz_clear(b); mpz_clear(r);
}

void gmp_sub(const char *a_str, const char *b_str, char *out_buf, size_t out_len) {
    mpz_t a, b, r;
    mpz_init(a); mpz_init(b); mpz_init(r);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);
    mpz_sub(r, a, b);
    mpz_get_str(out_buf, 10, r);
    mpz_clear(a); mpz_clear(b); mpz_clear(r);
}

void gmp_mul(const char *a_str, const char *b_str, char *out_buf, size_t out_len) {
    mpz_t a, b, r;
    mpz_init(a); mpz_init(b); mpz_init(r);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);
    mpz_mul(r, a, b);
    mpz_get_str(out_buf, 10, r);
    mpz_clear(a); mpz_clear(b); mpz_clear(r);
}

void gmp_cmp(const char *a_str, const char *b_str) {
    mpz_t a, b;
    mpz_init(a); mpz_init(b);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);
    int cmp = mpz_cmp(a, b);
    mpz_clear(a); mpz_clear(b);
    // Return via integer; handled by FFI
}

void gmp_bits(const char *a_str, unsigned long *bits_out) {
    mpz_t a;
    mpz_init(a);
    mpz_set_str(a, a_str, 10);
    *bits_out = mpz_sizeinbase(a, 2);
    mpz_clear(a);
}

void gmp_tdiv_qr(const char *a_str, const char *b_str,
                 char *q_buf, size_t q_len,
                 char *r_buf, size_t r_len) {
    mpz_t a, b, q, r;
    mpz_init(a); mpz_init(b); mpz_init(q); mpz_init(r);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);
    mpz_tdiv_qr(q, r, a, b);
    mpz_get_str(q_buf, 10, q);
    mpz_get_str(r_buf, 10, r);
    mpz_clear(a); mpz_clear(b); mpz_clear(q); mpz_clear(r);
}
