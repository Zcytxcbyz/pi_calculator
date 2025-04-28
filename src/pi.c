#include "pi.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

#ifdef DEBUG
unsigned long cache_hit_count = 0; // Count cache hits
#endif

static mpz_t CONST_X_BASE;  // CONST_X_BASE = -262537412640768000
static mpz_t CONST_L_K;     // CONST_L_K = 545140134
static mpz_t CONST_L_ADD;   // CONST_L_ADD = 13591409

// Type definition for thread private variables
typedef struct {
    mpf_t S, term, temp_f;
    mpz_t temp, M, L, X, K, k_fact, three_k_fact, six_k_fact;
} ThreadVariables;

#ifdef ENABLE_CACHE
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

#ifdef ENABLE_CACHE
// Initialize thread cache
void init_thread_cache(ThreadCache* cache) {
    cache->k_M = 0;
    cache->K_X = 0;
    mpz_inits(cache->k_fact, cache->three_k_fact, cache->six_k_fact, cache->X, NULL);
    mpz_set_ui(cache->k_fact, 1);         // k = 0 -> k! = 1
    mpz_set_ui(cache->three_k_fact, 1);   // k = 0 -> (3k)! = 1
    mpz_set_ui(cache->six_k_fact, 1);     // k = 0 -> (6k)! = 1
    mpz_set_ui(cache->X, 1);              // k = 0 -> (-262537412640768000)^k = 1
}

// Clean up thread cache
void clean_thread_cache(ThreadCache* cache) {
    mpz_clears(cache->k_fact, cache->three_k_fact, cache->six_k_fact, cache->X, NULL);
}

// Set cache variables
void set_cache(unsigned long k, ThreadCache* cache, ThreadVariables* var) {
    if (k > cache->k_M) {
        cache->k_M = k;
        mpz_set(cache->six_k_fact, var->six_k_fact);
        mpz_set(cache->three_k_fact, var->three_k_fact);
        mpz_set(cache->k_fact, var->k_fact);
    }

    if (k > cache->K_X) {
        cache->K_X = k;
        mpz_set(cache->X, var->X);
    }
}
#endif

// Calculate M = (6k)! / ((3k)! * (k!)^3)
#ifdef ENABLE_CACHE
void calculate_M(unsigned long k, ThreadVariables* var, ThreadCache* cache) {
#else
void calculate_M(unsigned long k, ThreadVariables* var) {
#endif
    if (k == 0) {
        // k = 0 -> k! = 1
        // k = 0 -> (3k)! = 1
        // k = 0 -> (6k)! = 1
        mpz_set_ui(var->k_fact, 1);
        mpz_set_ui(var->three_k_fact, 1);
        mpz_set_ui(var->six_k_fact, 1);
    #ifdef ENABLE_CACHE
    } else if (k - 1 == cache->k_M && k > 1) {
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

        #ifdef DEBUG
        #pragma omp atomic
        ++cache_hit_count;
        #endif
    #endif
    } else {
        // Calculate factorials
        mpz_fac_ui(var->six_k_fact, 6 * k);
        mpz_fac_ui(var->three_k_fact, 3 * k);
        mpz_fac_ui(var->k_fact, k);
    }

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
#ifdef ENABLE_CACHE
void calculate_X(unsigned long k, ThreadVariables* var, ThreadCache* cache) {
#else
void calculate_X(unsigned long k, ThreadVariables* var) {
#endif
    // Calculate X = (-262537412640768000)^k
    if (k == 0) {
        // K = 0 -> X = 1
        mpz_set_ui(var->X, 1);
    #ifdef ENABLE_CACHE
    } else if (k - 1 == cache->K_X && k > 1) {
        // Recursive calculation power
        // (-262537412640768000)^k = (-262537412640768000)^(k-1)*(-262537412640768000)
        mpz_mul(var->X, cache->X, CONST_X_BASE);
        #ifdef DEBUG
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
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads, const char* omp_schedule, int chunk_size) {
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

    // Set OpenMP schedule type and chunk size
    omp_sched_t schedule_type;
    if (strcmp(omp_schedule, "static") == 0) {
        schedule_type = omp_sched_static;
    } else if (strcmp(omp_schedule, "dynamic") == 0) {
        schedule_type = omp_sched_dynamic;
    } else { // Default to "guided"
        schedule_type = omp_sched_guided;
    }
    omp_set_schedule(schedule_type, chunk_size);

    #ifdef DEBUG
    const char* debug_schedule_name;
    omp_sched_t debug_schedule_id;
    int debug_chunk_size;
    omp_get_schedule(&debug_schedule_id, &debug_chunk_size);
    switch (debug_schedule_id) {
        case omp_sched_static:
        debug_schedule_name = "static";
            break;
        case omp_sched_dynamic:
        debug_schedule_name = "dynamic";
            break;
        case omp_sched_guided:
        debug_schedule_name = "guided";
            break;
        default:
        debug_schedule_name = 0;
    }
    printf("OpenMP schedule type: %s, chunk size: %d\n", debug_schedule_name, debug_chunk_size);
    #endif

    // Initialize global constants
    init_constants();

    #pragma omp parallel shared(S, CONST_X_BASE, CONST_L_K, CONST_L_ADD)
    {

        ThreadVariables var; // Thread private variables
        init_thread_variables(&var); // Initialize thread variables

        #ifdef ENABLE_CACHE
        ThreadCache cache; // Thread var cache
        init_thread_cache(&cache); // Initialize thread cache
        #endif

        #pragma omp for schedule(runtime)
        for (unsigned long k = 0; k < iterations; k++) {
            mpz_set_ui(var.K, k);

            // Calculate M = (6k)! / ((3k)! * (k!)^3)
            #ifdef ENABLE_CACHE
            calculate_M(k, &var, &cache);
            #else
            calculate_M(k, &var);
            #endif

            // Calculate L = 545140134k + 13591409
            calculate_L(k, &var);

            // Calculate X = (-262537412640768000)^k
            #ifdef ENABLE_CACHE
            calculate_X(k, &var, &cache);
            #else
            calculate_X(k, &var);
            #endif

            // Calculate the current item: term = M * L / X
            calculate_term(k, &var);

            // Accumulate to thread private variables
            mpf_add(var.S, var.S, var.term);

            #ifdef ENABLE_CACHE
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

        #ifdef ENABLE_CACHE
        clean_thread_cache(&cache); // Clean up thread cache
        #endif

    }

    // Calculate PI = C / S
    mpf_div(pi, C, S);

    // Clean up global constants
    clean_constants();

    // Clean up variables
    mpf_clears(C, S, temp, NULL);

    #ifdef DEBUG
    printf("Cache hit count: %lu\n", cache_hit_count);
    printf("Cache hit ratio: %.2f%%\n", (double)cache_hit_count / (iterations * 2) * 100);
    #endif
}

void write_pi_to_file(const mpf_t pi, unsigned long digits, const char* filename, double computation_time, bool format_output, size_t buffer_size) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        perror("Failed to open file");
        return;
    }

    // Write header
    fprintf(file, "Pi calculated to %lu digits. ", digits);
    fprintf(file, "Computation time: %.2f seconds.\n\n", computation_time);

    // Write the first two digits of PI
    fprintf(file, "3.\n");

    mp_exp_t exp;
    // Obtain the string representation of PI
    char* pi_str = mpf_get_str(NULL, &exp, 10, digits + 2, pi);
    if (!pi_str) {
        fprintf(stderr, "Failed to convert pi to string\n");
        fclose(file);
        return;
    }

    // Adjust the exponent to get the correct number of digits
    if(exp != 1) {
        fprintf(stderr, "Unexpected exponent value: %ld\n", exp);
        free(pi_str);
        fclose(file);
        return;
    }

    #ifdef DEBUG
    int flush_count = 0; // Count the number of times the buffer is flushed
    #endif

    // Allocate buffer dynamically based on the specified size
    char* buffer = (char*)malloc(buffer_size);
    size_t buffer_index = 0;

    if(format_output) {
        // Block write optimization: every 100 characters on a line, add spaces every 10 characters
        for (unsigned long line_start = 1; line_start <= digits; line_start += 100) {
            unsigned long line_end = line_start + 100;
            if (line_end > digits + 1) line_end = digits + 1;

            // Process blocks of every 10 characters
            for (unsigned long block_start = line_start; block_start < line_end; block_start += 10) {
                    unsigned long block_end = block_start + 10;
                if (block_end > line_end) block_end = line_end;

                // Copy character block
                size_t block_length = block_end - block_start;
                if (buffer_index + block_length >= buffer_size) {
                    // Buffer full, write to file
                    fwrite(buffer, sizeof(char), buffer_index, file);
                    buffer_index = 0;

                    #ifdef DEBUG
                    ++flush_count;
                    #endif
                }
                memcpy(buffer + buffer_index, pi_str + block_start, block_length);
                buffer_index += block_length;

                // Add spaces (do not add to the last block)
                if (block_end < line_end) {
                    if (buffer_index > buffer_size - 1) {
                        // Buffer full, write to file
                        fwrite(buffer, sizeof(char), buffer_index, file);
                        buffer_index = 0;

                        #ifdef DEBUG
                        ++flush_count;
                        #endif
                    }
                    buffer[buffer_index++] = ' ';
                }
            }

            // Add line breaks (do not add to the last line)
            if (line_end <= digits) {
                if (buffer_index > buffer_size - 1) {
                    // Buffer full, write to file
                    fwrite(buffer, sizeof(char), buffer_index, file);
                    buffer_index = 0;

                    #ifdef DEBUG
                    ++flush_count;
                    #endif
                }
                buffer[buffer_index++] = '\n';
            }
        }

        // Write remaining buffer to file
        if (buffer_index > 0) {
            fwrite(buffer, sizeof(char), buffer_index, file);

            #ifdef DEBUG
            ++flush_count;
            #endif
        }

    } else {
        // Write directly without formatting.
        size_t offset = 1; // Skip '3.'
        unsigned long remaining = digits;
        while (remaining > 0) {
            // Calculate the chunk size
            size_t chunk = (remaining > buffer_size) ? buffer_size : remaining;
            if (offset + chunk > strlen(pi_str)) chunk = strlen(pi_str) - offset;

            // Copy the chunk to the buffer
            memcpy(buffer, pi_str + offset, chunk);
            fwrite(buffer, sizeof(char), chunk, file);
            offset += chunk;
            remaining -= chunk;

            #ifdef DEBUG
            ++flush_count;
            #endif
        }
    }

    #ifdef DEBUG
    printf("Buffer flush count: %d\n", flush_count); // Number of times the buffer was flushed
    #endif

    free(buffer);
    free(pi_str);
    fclose(file);
}
