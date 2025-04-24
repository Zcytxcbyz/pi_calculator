#include "pi.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#define USE_DEBUG 0 // Set to 1 to enable debug output
#define OUTPUT_FORMAT 1 // Set to 1 to enable output formatting
#define ENABLE_CACHE 1 // Set to 1 to enable cache for factorials and powers
#define OMP_SCHEDULE guided // OpenMP schedule type

#if USE_DEBUG
int cache_hit_count = 0;
#endif

static mpz_t CONST_X_BASE; // CONST_X_BASE = -262537412640768000
static mpz_t CONST_L_K; // CONST_L_K = 545140134
static mpz_t CONST_L_ADD; // CONST_L_ADD = 13591409

// Type definition for thread private variables
typedef struct {
    mpf_t S, term, temp_f;
    mpz_t temp, M, L, X, K, k_fact, three_k_fact, six_k_fact;
} ThreadVariables;

#if ENABLE_CACHE
// Cache for factorials and powers
typedef struct {
    unsigned long k_M, K_X;
    mpz_t k_fact, three_k_fact, six_k_fact, X;
} ThreadCache;
#endif

// Initialize constants (executed before entering the parallel region for the first time)
static void init_constants() {
    #pragma omp single
    {
        mpz_init_set_str(CONST_X_BASE, "-262537412640768000", 10);
        mpz_init_set_ui(CONST_L_K, 545140134);
        mpz_init_set_ui(CONST_L_ADD, 13591409);
    }
}

// Clean up constants (executed after exiting the parallel region)
static void clean_constants() {
    mpz_clears(CONST_X_BASE, CONST_L_K, CONST_L_ADD, NULL);
}

// Initialize thread variables
void init_thread_variables(ThreadVariables* var) {
    mpf_init_set_ui(var->S, 0);
    mpf_inits(var->term, var->temp_f, NULL);
    mpz_inits(var->temp, var->M, var->L, var->X, var->K, var->k_fact, var->three_k_fact, var->six_k_fact, NULL);
}

// Clean up thread variables
void clean_thread_variables(ThreadVariables* var) {
    mpf_clears(var->S, var->term, var->temp_f, NULL);
    mpz_clears(var->temp, var->M, var->L, var->X, var->K, var->k_fact, var->three_k_fact, var->six_k_fact, NULL);
}

#if ENABLE_CACHE
// Initialize thread cache
void init_thread_cache(ThreadCache* cache) {
    cache->k_M = -2;
    cache->K_X = -2;
    mpz_inits(cache->k_fact, cache->three_k_fact, cache->six_k_fact, cache->X, NULL);
}

// Clean up thread cache
void clean_thread_cache(ThreadCache* cache) {
    mpz_clears(cache->k_fact, cache->three_k_fact, cache->six_k_fact, cache->X, NULL);
}

// Set cache variables
void set_cache(unsigned long k, ThreadCache* cache, ThreadVariables* var) {
    if (k > cache->k_M || k == 0) {
        cache->k_M = k;
        mpz_set(cache->six_k_fact, var->six_k_fact);
        mpz_set(cache->three_k_fact, var->three_k_fact);
        mpz_set(cache->k_fact, var->k_fact);
    }

    if (k > cache->K_X || k == 0) {
        cache->K_X = k;
        mpz_set(cache->X, var->X);
    }
}
#endif

// Calculate M = (6k)! / ((3k)! * (k!)^3)
#if ENABLE_CACHE
void calculate_M(unsigned long k, ThreadVariables* var, ThreadCache* cache) {
#else
void calculate_M(unsigned long k, ThreadVariables* var) {
#endif
// (3(k-1))! = (3k-1)!*(3k)
    #if ENABLE_CACHE
    if (k != 0 && k - 1 == cache->k_M) {
        // Recursive calculation of factorial
        // k! = (k-1)!*k
        // (3k)! = (3k-1)!*(3k)
        // (6k)! = (6k-1)!*(6k)
        mpz_mul_ui(var->six_k_fact, cache->six_k_fact, 6 * k);

        // Calculate 3k! = 3(k-1)! * [3(k-1)+1, 3(k-1)+2, ..., 3k]
        unsigned long prev_3k = 3 * (k - 1);
        mpz_set(var->three_k_fact, cache->three_k_fact);
        for (unsigned long i = prev_3k + 1; i <= 3 * k; i++) {
            mpz_mul_ui(var->three_k_fact, var->three_k_fact, i);
        }

        // Calculate 6k! = 6(k-1)! * [6(k-1)+1, ..., 6k]
        unsigned long prev_6k = 6 * (k - 1);
        mpz_set(var->six_k_fact, cache->six_k_fact);
        for (unsigned long i = prev_6k + 1; i <= 6 * k; i++) {
            mpz_mul_ui(var->six_k_fact, var->six_k_fact, i);
        }

        #if USE_DEBUG
        #pragma omp atomic
        ++cache_hit_count;
        #endif
    } else {
    #endif
        // Calculate factorials
        mpz_fac_ui(var->six_k_fact, 6 * k);
        mpz_fac_ui(var->three_k_fact, 3 * k);
        mpz_fac_ui(var->k_fact, k);
    #if ENABLE_CACHE
    }
    #endif

    mpz_pow_ui(var->temp, var->k_fact, 3);
    mpz_mul(var->temp, var->temp, var->three_k_fact);
    mpz_divexact(var->M, var->six_k_fact, var->temp);
}

// Calculate L = 545140134k + 13591409
void calculate_L(unsigned long k, ThreadVariables* var) {
    mpz_mul_ui(var->temp, CONST_L_K, k);
    mpz_add(var->L, var->temp, CONST_L_ADD);
}

// Calculate X = (-262537412640768000)^k
#if ENABLE_CACHE
void calculate_X(unsigned long k, ThreadVariables* var, ThreadCache* cache) {
#else
void calculate_X(unsigned long k, ThreadVariables* var) {
#endif
    // Calculate X = (-262537412640768000)^k
    if (k == 0) {
        // K = 0 -> X = 1
        mpz_set_ui(var->X, 1);
    #if ENABLE_CACHE
    } else if (k - 1 == cache->K_X) {
        // Recursive calculation power
        // (-262537412640768000)^k = (-262537412640768000)^(k-1)*(-262537412640768000)
        mpz_mul(var->X, cache->X, CONST_X_BASE);
        #if USE_DEBUG
        #pragma omp atomic
        ++cache_hit_count;
        #endif
    #endif
    } else {
        // Calculate power
        mpz_pow_ui(var->X, CONST_X_BASE, k);
    }
}

// Calculate the current item: term = M * L / X
void calculate_term(unsigned long k, ThreadVariables* var) {
    mpf_set_z(var->term, var->M);
    mpf_set_z(var->temp_f, var->L);
    mpf_mul(var->term, var->term, var->temp_f);
    mpf_set_z(var->temp_f, var->X);
    mpf_div(var->term, var->term, var->temp_f);
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

    #pragma omp parallel shared(S, CONST_X_BASE, CONST_L_K, CONST_L_ADD)
    {

        ThreadVariables var; // Thread private variables
        init_thread_variables(&var); // Initialize thread variables

        #if ENABLE_CACHE
        ThreadCache cache; // Thread var cache
        init_thread_cache(&cache); // Initialize thread cache
        #endif

        #pragma omp for schedule(OMP_SCHEDULE)
        for (unsigned long k = 0; k < iterations; k++) {
            mpz_set_ui(var.K, k);

            // Calculate M = (6k)! / ((3k)! * (k!)^3)
            #if ENABLE_CACHE
            calculate_M(k, &var, &cache);
            #else
            calculate_M(k, &var);
            #endif

            // Calculate L = 545140134k + 13591409
            calculate_L(k, &var);

            // Calculate X = (-262537412640768000)^k
            #if ENABLE_CACHE
            calculate_X(k, &var, &cache);
            #else
            calculate_X(k, &var);
            #endif

            // Calculate the current item: term = M * L / X
            calculate_term(k, &var);

            // Accumulate to thread private variables
            mpf_add(var.S, var.S, var.term);

            #if ENABLE_CACHE
            // Set cache variables
            set_cache(k, &cache, &var);
            #endif

        }

        // Merge thread private variables into global variables
        #pragma omp critical
        {
            mpf_add(S, S, var.S);
        }

        clean_thread_variables(&var); // Clean up thread variables

        #if ENABLE_CACHE
        clean_thread_cache(&cache); // Clean up thread cache
        #endif

    }

    // Calculate PI = C / S
    mpf_div(pi, C, S);

    // Clean up global constants
    clean_constants();

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

    fprintf(file, "Pi calculated to %lu digits. ", digits);
    fprintf(file, "Computation time: %.2f seconds.\n\n", computation_time);
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
        buffer[buffer_index++] = pi_str[i]; // Copy the digit

        #if OUTPUT_FORMAT
        // Add formatting for output
        if (i % 100 == 0 || i == digits) {
            buffer[buffer_index++] = '\n';
        } else if (i % 10 == 0) {
            buffer[buffer_index++] = ' ';
        }
        #endif

        // Flush the buffer if it reaches a certain size or if it's the last digit
        if (buffer_index >= sizeof(buffer) - 2 || i == digits) {
            buffer[buffer_index] = '\0';
            fputs(buffer, file);
            buffer_index = 0;
        }
    }

    free(pi_str);
    fclose(file);
}
