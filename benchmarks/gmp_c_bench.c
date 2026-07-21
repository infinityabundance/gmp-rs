/**
 * benchmarks/gmp_c_bench.c
 *
 * Standalone C benchmark for raw GMP performance.
 * Compile:  gcc -O2 -o gmp_c_bench gmp_c_bench.c -lgmp
 * Run:      ./gmp_c_bench
 * Compare:  cargo bench --features gmp_cross_check --bench gmp_comparison
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <time.h>

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static const char *OP_SINGLE_A = "1000000000000";
static const char *OP_SINGLE_B = "2000000000000";
static const char *OP_FOUR_STR =
    "6277101735386680763835789423207666416102355444464034512895";
static const char *OP_EIGHT_STR =
    "134078079299425970995740249982058461274793658205923933777235614437"
    "217640300735469768018742981669034276900318581864860508537538828119"
    "46569946433649006084095";

typedef void (*bench_fn)(mpz_t, const mpz_t, const mpz_t);

static double run_bench(bench_fn fn, const char *label,
                        const char *a_str, const char *b_str,
                        int iterations) {
    mpz_t a, b, r;
    mpz_init(a); mpz_init(b); mpz_init(r);
    mpz_set_str(a, a_str, 10);
    mpz_set_str(b, b_str, 10);

    double start = now_sec();
    for (int i = 0; i < iterations; i++) { fn(r, a, b); }
    double elapsed = now_sec() - start;

    double per_op_ns = (elapsed / iterations) * 1e9;
    printf("%-40s  %8d  %10.1f ns/op\n", label, iterations, per_op_ns);
    mpz_clear(a); mpz_clear(b); mpz_clear(r);
    return per_op_ns;
}

static void op_add(mpz_t r, const mpz_t a, const mpz_t b) { mpz_add(r, a, b); }
static void op_sub(mpz_t r, const mpz_t a, const mpz_t b) { mpz_sub(r, a, b); }
static void op_mul(mpz_t r, const mpz_t a, const mpz_t b) { mpz_mul(r, a, b); }

int main(void) {
    printf("\n=== GMP C API benchmarks ===\n\n");
    printf("%-40s  %9s  %12s\n", "Operation", "Iters", "ns/op");
    printf("------------------------------------------------\n");
    int m = 1000000, l = 100000;
    run_bench(op_add, "add (single limb)", OP_SINGLE_A, OP_SINGLE_B, m);
    run_bench(op_add, "add (eight limbs)", OP_EIGHT_STR, OP_EIGHT_STR, m);
    run_bench(op_sub, "sub (single limb)", OP_SINGLE_A, OP_SINGLE_B, m);
    run_bench(op_sub, "sub (eight limbs)", OP_EIGHT_STR, OP_EIGHT_STR, m);
    run_bench(op_mul, "mul (single limb)", OP_SINGLE_A, OP_SINGLE_B, m);
    run_bench(op_mul, "mul (four × eight)", OP_FOUR_STR, OP_EIGHT_STR, l);
    printf("\n");
    return 0;
}
