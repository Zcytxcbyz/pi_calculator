#include "pi.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#define ENABLE_CACHE 1 // Set to 1 to enable cache for factorials and powers
#define USE_DEBUG 0 // Set to 1 to enable debug output
#define OMP_SCHEDULE schedule(guided) // OpenMP schedule type

static mpz_t CONST_X_BASE; // CONST_X_BASE = -262537412640768000
static mpz_t CONST_L_K; // CONST_L_K = 545140134
static mpz_t CONST_L_ADD; // CONST_L_ADD = 13591409

// Initialize constants (executed before entering the parallel region for the first time)
static void init_constants() {
    #pragma omp single
    {
        mpz_init_set_str(CONST_X_BASE, "-262537412640768000", 10);
        mpz_init_set_ui(CONST_L_K, 545140134);
        mpz_init_set_ui(CONST_L_ADD, 13591409);
    }
}

// Chudnovsky algorithm calculates PI
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads) {
    // Set sufficient precision
    mpf_set_default_prec((digits + 2) * log2(10));

    mpf_t C, S, temp;

    // Initialize variable
    mpf_inits(C, temp, NULL);
    mpf_init_set_ui(S, 0);

    // Constant C = 426880 * sqrt(10005)
    mpf_set_ui(C, 426880);
    mpf_sqrt_ui(temp, 10005);
    mpf_mul(C, C, temp);

    // Calculate the required number of iterations (empirical formula)
    unsigned long iterations = (digits / 14) + 1;

    // Set the number of OpenMP threads
    omp_set_num_threads(num_threads);

    // Initialize global constants
    init_constants();

    #if USE_DEBUG
    int cache_hit_count = 0;
    #endif

    #pragma omp parallel private(temp) shared(S, CONST_X_BASE, CONST_L_K, CONST_L_ADD)
    {
        mpf_t S_private, term, temp_f;
        mpz_t temp, M, L, X, K, k_fact, three_k_fact, six_k_fact;

        mpf_init_set_ui(S_private, 0);
        mpf_inits(term, temp_f, NULL);
        mpz_inits(temp, M, L, X, K, k_fact, three_k_fact, six_k_fact, NULL);

        #if ENABLE_CACHE
        // Thread private cache
        struct {
            unsigned long k_M, K_X;
            mpz_t k_fact, three_k_fact, six_k_fact, X;
        } cache;

        cache.k_M = -2;
        cache.K_X = -2;
        mpz_inits(cache.k_fact, cache.three_k_fact, cache.six_k_fact, cache.X);
        #endif

        #pragma omp for OMP_SCHEDULE
        for (unsigned long k = 0; k < iterations; k++) {
            mpz_set_ui(K, k);

            // Calculate M = (6k)! / ((3k)! * (k!)^3)
            #if ENABLE_CACHE
            if (k == 0 || k - 1 != cache.k_M) {
            #endif
                // Calculate factorials
                mpz_fac_ui(six_k_fact, 6 * k);
                mpz_fac_ui(three_k_fact, 3 * k);
                mpz_fac_ui(k_fact, k);
            #if ENABLE_CACHE
            } else {
                // Recursive calculation of factorial
                // k! = (k-1)! * k
                // (3k)! = (3k-1)!*(3k)
                // (6k)! = (6k-1)!*(6k)
                mpz_mul_ui(six_k_fact, cache.six_k_fact, 6 * k);
                mpz_mul_ui(three_k_fact, cache.three_k_fact, 3 * k);
                mpz_mul_ui(k_fact, cache.k_fact, k);
                #if USE_DEBUG
                #pragma omp atomic
                ++cache_hit_count;
                #endif
            }
            #endif

            mpz_pow_ui(temp, k_fact, 3);
            mpz_mul(temp, temp, three_k_fact);
            mpz_divexact(M, six_k_fact, temp);

            // Calculate L = 545140134k + 13591409
            mpz_mul_ui(temp, CONST_L_K, k);
            mpz_add(L, temp, CONST_L_ADD);

            // Calculate X = (-262537412640768000)^k
            if (k == 0) {
                // K = 0 -> X = 1
                mpz_set_ui(X, 1);
            #if ENABLE_CACHE
            } else if (k - 1 == cache.K_X) {
                // Recursive calculation power
                // (-262537412640768000)^k = (-262537412640768000)^(k-1)*(-262537412640768000)
                mpz_mul(X, cache.X, CONST_X_BASE);
                #if USE_DEBUG
                #pragma omp atomic
                ++cache_hit_count;
                #endif
            #endif
            } else {
                // Calculate power
                mpz_pow_ui(X, CONST_X_BASE, k);
            }

            // Calculate the current item: term = M * L / X
            mpf_set_z(term, M);
            mpf_set_z(temp_f, L);
            mpf_mul(term, term, temp_f);
            mpf_set_z(temp_f, X);
            mpf_div(term, term, temp_f);

            // Accumulate to thread private variables
            mpf_add(S_private, S_private, term);

            #if ENABLE_CACHE
            // Cache variables
            if (k > cache.k_M || k == 0) {
                cache.k_M = k;
                mpz_set(cache.six_k_fact, six_k_fact);
                mpz_set(cache.three_k_fact, three_k_fact);
                mpz_set(cache.k_fact, k_fact);
            }

            if (k > cache.K_X || k == 0) {
                cache.K_X = k;
                mpz_set(cache.X, X);
            }
            #endif

        }

        // Merge thread private variables into global variables
        #pragma omp critical
        {
            mpf_add(S, S, S_private);
        }

        // Clean up thread variables
        mpf_clears(S_private, term, temp_f, NULL);
        mpz_clears(temp, M, L, X, K, k_fact, three_k_fact, six_k_fact, NULL);

        #if ENABLE_CACHE
        // Clean up cache
        mpz_clears(cache.k_fact, cache.three_k_fact, cache.six_k_fact, cache.X, NULL);
        #endif

    }

    // Calculate PI = C / S
    mpf_div(pi, C, S);

    // Clean up global constants
    mpz_clears(CONST_X_BASE, CONST_L_K, CONST_L_ADD, NULL);

    // Clean up variables
    mpf_clears(C, S, temp, NULL);

    #if USE_DEBUG
    printf("Cache hit count: %d\n", cache_hit_count);
    #endif
}

void write_pi_to_file(const mpf_t pi, unsigned long digits, const char* filename, double computation_time) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Failed to open file");
        return;
    }

    fprintf(file, "Pi calculated to %lu digits\n", digits);
    fprintf(file, "Computation time: %.2f seconds\n", computation_time);
    fprintf(file, "3.\n");

    mp_exp_t exp;
    // Obtain the string representation of PI
    char* pi_str = mpf_get_str(NULL, &exp, 10, digits + 2, pi);
    if (!pi_str) {
        fprintf(stderr, "Failed to convert pi to string\n");
        fclose(file);
        return;
    }

    // Write file, 100 bits per line
    char buffer[1024];
    size_t buffer_index = 0;

    for (unsigned long i = 1; i < digits + 1; i++) {
        buffer[buffer_index++] = pi_str[i];
        if (i % 100 == 0 || i == digits) {
            buffer[buffer_index++] = '\n';
        } else if (i % 10 == 0) {
            buffer[buffer_index++] = ' ';
        }

        if (buffer_index >= sizeof(buffer) - 2 || i == digits) {
            buffer[buffer_index] = '\0';
            fputs(buffer, file);
            buffer_index = 0;
        }
    }

    free(pi_str);
    fclose(file);
}
