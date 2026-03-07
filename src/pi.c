#include "pi.h"
#include "checkpoint.h"
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
    #ifdef ENABLE_BLOCK_FACTORIAL
    // Block factorial variables
    mpz_t block_prod; // Block product for block factorial
    unsigned long block_size; // Block size for block factorial
    #endif
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
    mpz_init_set_str(CONST_X_BASE, "-262537412640768000", 10);
    mpz_init_set_ui(CONST_L_K, 545140134);
    mpz_init_set_ui(CONST_L_ADD, 13591409);
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
    #ifdef ENABLE_BLOCK_FACTORIAL
    mpz_init(var->block_prod); // Initialize block product
    #endif
}

// Clean up thread variables
void clean_thread_variables(ThreadVariables* var) {
    mpf_clears(var->S, var->term, var->temp_f, NULL);
    mpz_clears(var->temp, var->M, var->L, var->X, var->K, var->k_fact, var->three_k_fact, var->six_k_fact, NULL);
    #ifdef ENABLE_BLOCK_FACTORIAL
    mpz_clear(var->block_prod); // Clean up block product
    #endif
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

#ifdef ENABLE_BLOCK_FACTORIAL
// Block factorial calculation
void block_factorial(unsigned long start, unsigned long end, unsigned long block_size, mpz_t block_prod, mpz_t fact) {
    for (unsigned long i = start; i <= end; i += block_size) {
        unsigned long block_end = (i + block_size <= end) ? i + block_size : end + 1;
        mpz_set_ui(block_prod, 1);
        for (unsigned long j = i; j < block_end; j++) {
            mpz_mul_ui(block_prod, block_prod, j);
        }
        mpz_mul(fact, fact, block_prod);
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
    } else if (k - 1 == cache->k_M) {
        // Recursive calculation of factorial
        // k! = (k-1)!*k
        // (3k)! = (3k-1)!*(3k)
        // (6k)! = (6k-1)!*(6k)

        // Calculate k! = (k-1)! * k
        mpz_mul_ui(var->k_fact, cache->k_fact, k);

        // Calculate 3k! = 3(k-1)! * [3(k-1)+1, 3(k-1)+2, ..., 3k]
        unsigned long prev_3k = 3 * (k - 1);
        mpz_set(var->three_k_fact, cache->three_k_fact);
        #ifdef ENABLE_BLOCK_FACTORIAL
        block_factorial(prev_3k + 1, 3 * k, var->block_size, var->block_prod, var->three_k_fact);
        #else
        mpz_set(var->three_k_fact, cache->three_k_fact);
        for (unsigned long i = prev_3k + 1; i <= 3 * k; i++) {
            mpz_mul_ui(var->three_k_fact, var->three_k_fact, i);
        }
        #endif

        // Calculate 6k! = 6(k-1)! * [6(k-1)+1, ..., 6k]
        unsigned long prev_6k = 6 * (k - 1);
        mpz_set(var->six_k_fact, cache->six_k_fact);
        #ifdef ENABLE_BLOCK_FACTORIAL
        block_factorial(prev_6k + 1, 6 * k, var->block_size, var->block_prod, var->six_k_fact);
        #else
        for (unsigned long i = prev_6k + 1; i <= 6 * k; i++) {
            mpz_mul_ui(var->six_k_fact, var->six_k_fact, i);
        }
        #endif

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
    } else if (k - 1 == cache->K_X) {
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
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads, const char* omp_schedule,
    int chunk_size VAR_BLOCK_SIZE, bool show_progress, int progress_freq) {
    /* completed_count Used solely for progress display;
     * Does not increment if progress is disabled, avoiding atomic operation overhead */
    unsigned long long completed_count = 0;

    // Thread Count Legitimacy Verification
    if (num_threads <= 0) {
        fprintf(stderr, "Warning: invalid thread count (%d), using 1 thread.\n", num_threads);
        num_threads = 1;
    }

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

    // Array Block Reduction: Assigns each thread an independent segment and slot
    int max_threads = num_threads;  // Number of threads actually used (specified by the user)
    mpf_t* thread_S = (mpf_t*) malloc(max_threads * sizeof(mpf_t));
    if (!thread_S) {
        fprintf(stderr, "Error: Failed to allocate thread_S array\n");
        clean_constants();
        mpf_clears(C, S, temp, NULL);
        exit(1);
    }
    for (int i = 0; i < max_threads; i++) {
        mpf_init_set_ui(thread_S[i], 0);
    }

    #pragma omp parallel shared(S, thread_S, CONST_X_BASE, CONST_L_K, CONST_L_ADD, completed_count, show_progress, progress_freq)
    {
        int tid = omp_get_thread_num();  // Get the current thread ID
        ThreadVariables var; // Thread private variables
        init_thread_variables(&var); // Initialize thread variables
        #ifdef ENABLE_BLOCK_FACTORIAL
        var.block_size = block_size; // Set block size for block factorial
        #endif

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

            // Progress display: Atomic counter update (only when progress is enabled)
            if (show_progress) {
                #pragma omp atomic
                ++completed_count;

                if (completed_count % progress_freq == 0 && tid == 0) {
                    // Use the thread ID to determine execution on the main thread
                    fprintf(stderr, "\rProgress: %.2f%%", (double) completed_count / iterations * 100);
                    fflush(stderr);
                }
            }

        }

        /* This code is deprecated
        // Merge thread private variables into global variables
        #pragma omp critical
        {
            mpf_add(S, S, var.S);
        }
        */


        // Store the segment of this thread into the corresponding array slot
        mpf_set(thread_S[tid], var.S);

        clean_thread_variables(&var); // Clean up thread variables

        #ifdef ENABLE_CACHE
        clean_thread_cache(&cache); // Clean up thread cache
        #endif

    }

    // The main thread merges all parts
    for (int i = 0; i < max_threads; i++) {
        mpf_add(S, S, thread_S[i]);
        mpf_clear(thread_S[i]);
    }
    free(thread_S);

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

// Write the PI value to file
void write_pi_to_file(const mpf_t pi, unsigned long digits, const char* filename, double computation_time,
    bool format_output, size_t buffer_size, bool raw_output) {
    FILE* file = fopen(filename, raw_output ? "wb" : "w");
    if (!file) {
        perror("Failed to open file");
        return;
    }
    write_pi_to_stream(pi, digits, file, computation_time, format_output, buffer_size, raw_output);
    fclose(file);
}

// Write the PI value to stream
void write_pi_to_stream(const mpf_t pi, unsigned long digits, FILE* stream, double computation_time,
    bool format_output, size_t buffer_size, bool raw_output) {
    // Write header only if not in raw mode
    if (!raw_output) {
        fprintf(stream, "Pi calculated to %lu digits. ", digits);
        fprintf(stream, "Computation time: %.2f seconds.\n\n", computation_time);
    }

    mp_exp_t exp;
    // Obtain the string representation of PI
    char* pi_str = mpf_get_str(NULL, &exp, 10, digits + 2, pi);
    if (!pi_str) {
        fprintf(stderr, "Failed to convert pi to string\n");
        return;
    }

    // Adjust the exponent to get the correct number of digits
    if(exp != 1) {
        fprintf(stderr, "Unexpected exponent value: %ld\n", exp);
        free(pi_str);
        return;
    }

    // In raw mode: write "3." + digits (no extra newline after "3.")
    fprintf(stream,"3.");
    if (!raw_output) {
        fprintf(stream, "\n");
    }

    #ifdef DEBUG
    int flush_count = 0; // Count the number of times the buffer is flushed
    #endif

    // Allocate buffer dynamically based on the specified size
    char* buffer = (char*)malloc(buffer_size);
    if (!buffer) {
        perror("malloc failed");
        free(pi_str);
        return;
    }
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
                    fwrite(buffer, sizeof(char), buffer_index, stream);
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
                        fwrite(buffer, sizeof(char), buffer_index, stream);
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
                    fwrite(buffer, sizeof(char), buffer_index, stream);
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
            fwrite(buffer, sizeof(char), buffer_index, stream);

            #ifdef DEBUG
            ++flush_count;
            #endif
        }

    } else {
        // Write directly without formatting.
        size_t offset = 1; // Skip '3.'
        unsigned long remaining = digits;
        size_t pi_strlen = strlen(pi_str);
        while (remaining > 0) {
            // Calculate the chunk size
            size_t chunk = (remaining > buffer_size) ? buffer_size : remaining;
            if (offset + chunk > pi_strlen) chunk = pi_strlen - offset;

            // Copy the chunk to the buffer
            memcpy(buffer, pi_str + offset, chunk);
            fwrite(buffer, sizeof(char), chunk, stream);
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
}


/* ---------------------------------------------------------------------------------------------------------------------
 * calculate_pi_checkpoint: A π calculation function with checkpointing capability.
 *
 * This function is algorithmically identical to calculate_pi but incorporates a checkpointing mechanism:
 * - Divides the total iterations into blocks (size specified by checkpoint_freq). After each block completes,
 * saves the current progress (global state and iteration index) to a file.
 - Attempts to load existing checkpoint files at startup to resume computation from the interrupted point.
 *
* This function is implemented separately from the original calculate_pi to preserve the simplicity and purity of the original function,
* preventing checkpoint-related logic (file I/O, blocked loops, state saving, etc.) from intruding into the core computation flow.
*
 * All core computation units (e.g., calculate_M, calculate_L, calculate_X, calculate_term) reuse the original implementations,
 * ensuring algorithm results are identical to the non-checkpoint version.
 *
 * When users enable checkpointing via the command-line argument --checkpoint-enable, the main function calls this function.
 */
void calculate_pi_checkpoint(mpf_t pi, unsigned long digits, int num_threads, const char* omp_schedule,
    int chunk_size VAR_BLOCK_SIZE, bool show_progress, int progress_freq, bool quiet_flag,
    unsigned long checkpoint_freq, const char* checkpoint_file) {

    if (checkpoint_freq < 100 && !quiet_flag) {
        fprintf(stderr, "Warning: checkpoint_freq is very small (%lu). This may cause frequent I/O and significantly slow down the calculation.\n", checkpoint_freq);
    }

    unsigned long long completed_count = 0;  // Progress counter (used only when show_progress is enabled)

    if (num_threads <= 0) {
        fprintf(stderr, "Warning: invalid thread count (%d), using 1 thread.\n", num_threads);
        num_threads = 1;
    }

    mpf_set_default_prec((digits + 2) * log2(10));

    mpf_t C, temp;
    mpf_inits(C, temp, NULL);
    mpf_set_ui(C, 426880);
    mpf_sqrt_ui(temp, 10005);
    mpf_mul(C, C, temp);

    unsigned long iterations = (digits / 14) + 1;

    // Checkpoint Recovery
    unsigned long start_k = 0;
    mpf_t global_S;
    mpf_init_set_ui(global_S, 0);

    uint32_t saved_threads = 0, saved_flags = 0, current_flags=0;

    // Check compilation option flags
    #ifdef ENABLE_CACHE
    current_flags |= CHECKPOINT_FLAG_CACHE;
    #endif
    #ifdef ENABLE_BLOCK_FACTORIAL
    current_flags |= CHECKPOINT_FLAG_BLOCK_FACTORIAL;
    #endif

    int ret = load_checkpoint(checkpoint_file, &start_k, global_S, digits,
                              &saved_threads, &saved_flags, quiet_flag);
    if (ret == 0) {
        if (!quiet_flag) {
            printf("Resuming from iteration %lu (%.2f%%)\n",
                   start_k, (double)start_k / iterations * 100);
        }
        // Optional warning: Thread count or compilation options do not match
        if (saved_threads != (uint32_t)num_threads && !quiet_flag) {
            fprintf(stderr, "Warning: thread count mismatch (saved=%u, current=%d). Performance may be degraded due to load imbalance and cache inefficiency.\n",
                    saved_threads, num_threads);
        }

        if (saved_flags != current_flags && !quiet_flag) {
            fprintf(stderr, "Warning: compilation flags mismatch (saved=0x%x, current=0x%x). Performance may be affected.\n",
                    saved_flags, current_flags);
        }
    } else if (ret == -1) {
        if (!quiet_flag) printf("No checkpoint found, starting from 0.\n");
    } else {
        if (!quiet_flag) fprintf(stderr, "Warning: Checkpoint file invalid, starting from 0.\n");
        start_k = 0;
        mpf_set_ui(global_S, 0);
    }
    // ------------------------------------------------

    omp_set_num_threads(num_threads);
    omp_sched_t schedule_type;
    if (strcmp(omp_schedule, "static") == 0)        schedule_type = omp_sched_static;
    else if (strcmp(omp_schedule, "dynamic") == 0)  schedule_type = omp_sched_dynamic;
    else                                            schedule_type = omp_sched_guided;
    omp_set_schedule(schedule_type, chunk_size);

    init_constants();

    int max_threads = num_threads;
    mpf_t* thread_S = (mpf_t*) malloc(max_threads * sizeof(mpf_t));
    if (!thread_S) {
        fprintf(stderr, "Error: Failed to allocate thread_S array\n");
        clean_constants();
        mpf_clears(C, global_S, temp, NULL);
        exit(1);
    }
    for (int i = 0; i < max_threads; i++) {
        mpf_init_set_ui(thread_S[i], 0);
    }

    // ------------------ Block loop ------------------
    unsigned long current_k = start_k;
    while (current_k < iterations) {
        unsigned long block_end = current_k + checkpoint_freq;
        if (block_end > iterations) block_end = iterations;

        // Reset thread segments and arrays (zeroed before each block begins)
        for (int i = 0; i < max_threads; i++) {
            mpf_set_ui(thread_S[i], 0);
        }

        #pragma omp parallel shared(thread_S, completed_count, show_progress, progress_freq)
        {
            int tid = omp_get_thread_num();
            ThreadVariables var;
            init_thread_variables(&var);
            #ifdef ENABLE_BLOCK_FACTORIAL
            var.block_size = block_size;  // Note: block_size must be passed as a parameter.
            #endif
            #ifdef ENABLE_CACHE
            ThreadCache cache;
            init_thread_cache(&cache);
            #endif

            #pragma omp for schedule(runtime)
            for (unsigned long k = current_k; k < block_end; k++) {
                mpz_set_ui(var.K, k);
                #ifdef ENABLE_CACHE
                calculate_M(k, &var, &cache);
                #else
                calculate_M(k, &var);
                #endif
                calculate_L(k, &var);
                #ifdef ENABLE_CACHE
                calculate_X(k, &var, &cache);
                #else
                calculate_X(k, &var);
                #endif
                calculate_term(k, &var);
                mpf_add(var.S, var.S, var.term);

                #ifdef ENABLE_CACHE
                set_cache(k, &cache, &var);
                #endif
                if (show_progress) {
                    #pragma omp atomic
                    ++completed_count;
                    if (completed_count % progress_freq == 0 && tid == 0) {
                        fprintf(stderr, "\rProgress: %.2f%%",
                                (double)(completed_count) / iterations * 100);
                        fflush(stderr);
                    }
                }
            }

            // Add the portion and cumulative value of this thread to `thread_S[tid]`.
            mpf_add(thread_S[tid], thread_S[tid], var.S);

            clean_thread_variables(&var);
            #ifdef ENABLE_CACHE
            clean_thread_cache(&cache);
            #endif
        } // End of parallel section

        // Accumulate the contributions of each thread in the current block to global_S.
        for (int i = 0; i < max_threads; i++) {
            mpf_add(global_S, global_S, thread_S[i]);
        }

        current_k = block_end;

        // Save checkpoint

        if (save_checkpoint(checkpoint_file, current_k, global_S, digits,
                            (uint32_t)num_threads, current_flags, quiet_flag) != 0) {
            if (!quiet_flag) fprintf(stderr, "Warning: Failed to save checkpoint\n");
        } else if (show_progress && !quiet_flag) {
            fprintf(stderr, "\nCheckpoint saved at iteration %lu\n", current_k);
        }
    }
    // ------------------------------------------------

    // Calculate the final PI = C / global_S
    mpf_div(pi, C, global_S);

    clean_constants();
    for (int i = 0; i < max_threads; i++) mpf_clear(thread_S[i]);
    free(thread_S);
    mpf_clears(C, global_S, temp, NULL);
}

// This function is a merge of the differences between calculate_pi and calculate_pi_checkpoint to facilitate subsequent refactoring.
// Chudnovsky algorithm calculates PI
void calculate_pi_checkpoint_diff(mpf_t pi, unsigned long digits, int num_threads, const char *omp_schedule,
                                  int chunk_size VAR_BLOCK_SIZE, bool show_progress, int progress_freq, bool quiet_flag,
                                  bool enable_checkpoint, unsigned long checkpoint_freq, const char *checkpoint_file) {
    /* completed_count Used solely for progress display;
     * Does not increment if progress is disabled, avoiding atomic operation overhead */
    unsigned long long completed_count = 0;

    // Thread Count Legitimacy Verification
    if (num_threads <= 0) {
        fprintf(stderr, "Warning: invalid thread count (%d), using 1 thread.\n", num_threads);
        num_threads = 1;
    }

    // Set sufficient precision
    mpf_set_default_prec((digits + 2) * log2(10));

    // Original S changed to global_S
    mpf_t C, global_S, temp;

    // Initialize variable
    mpf_inits(C, temp, NULL);
    mpf_init_set_ui(global_S, 0);

    // Constant C = 426880 * sqrt(10005)
    mpf_set_ui(C, 426880);
    mpf_sqrt_ui(temp, 10005);
    mpf_mul(C, C, temp);

    // Calculate the required number of iterations (empirical formula)
    unsigned long iterations = (digits / 14) + 1;

    // ------------------ Checkpoint code begins ---------------------
    // Checkpoint Recovery
    unsigned long start_k = 0;
    uint32_t saved_threads = 0, saved_flags = 0, current_flags = 0;

    if (enable_checkpoint) {
        if (checkpoint_freq < 100 && !quiet_flag) {
            fprintf(
                stderr,
                "Warning: checkpoint_freq is very small (%lu). This may cause frequent I/O and significantly slow down the calculation.\n",
                checkpoint_freq);
        }

        // Check compilation option flags
        #ifdef ENABLE_CACHE
        current_flags |= CHECKPOINT_FLAG_CACHE;
        #endif
        #ifdef ENABLE_BLOCK_FACTORIAL
        current_flags |= CHECKPOINT_FLAG_BLOCK_FACTORIAL;
        #endif

        int ret = load_checkpoint(checkpoint_file, &start_k, global_S, digits,
                                  &saved_threads, &saved_flags, quiet_flag);
        if (ret == 0) {
            if (!quiet_flag) {
                printf("Resuming from iteration %lu (%.2f%%)\n",
                       start_k, (double) start_k / iterations * 100);
            }
            // Optional warning: Thread count or compilation options do not match
            if (saved_threads != (uint32_t) num_threads && !quiet_flag) {
                fprintf(
                    stderr,
                    "Warning: thread count mismatch (saved=%u, current=%d). Performance may be degraded due to load imbalance and cache inefficiency.\n",
                    saved_threads, num_threads);
            }

            if (saved_flags != current_flags && !quiet_flag) {
                fprintf(
                    stderr,
                    "Warning: compilation flags mismatch (saved=0x%x, current=0x%x). Performance may be affected.\n",
                    saved_flags, current_flags);
            }
        } else if (ret == -1) {
            if (!quiet_flag) printf("No checkpoint found, starting from 0.\n");
        } else {
            if (!quiet_flag) fprintf(stderr, "Warning: Checkpoint file invalid, starting from 0.\n");
            start_k = 0;
            mpf_set_ui(global_S, 0);
        }
    }
    // ------------------ Checkpoint code ends   ---------------------

    // Set the number of OpenMP threads
    omp_set_num_threads(num_threads);

    // Set OpenMP schedule type and chunk size
    omp_sched_t schedule_type;
    if (strcmp(omp_schedule, "static") == 0) {
        schedule_type = omp_sched_static;
    } else if (strcmp(omp_schedule, "dynamic") == 0) {
        schedule_type = omp_sched_dynamic;
    } else {
        // Default to "guided"
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

    // Array Block Reduction: Assigns each thread an independent segment and slot
    int max_threads = num_threads; // Number of threads actually used (specified by the user)
    mpf_t* thread_S = (mpf_t*) malloc(max_threads * sizeof(mpf_t));
    if (!thread_S) {
        fprintf(stderr, "Error: Failed to allocate thread_S array\n");
        clean_constants();
        mpf_clears(C, global_S, temp, NULL);
        exit(1);
    }
    for (int i = 0; i < max_threads; i++) {
        mpf_init_set_ui(thread_S[i], 0);
    }

    // ------------------ Checkpoint code begins ---------------------
    unsigned long current_k = enable_checkpoint ? start_k : 0;
    while (current_k < iterations) {
        unsigned long block_end = iterations;
        if (enable_checkpoint) {
            block_end = current_k + checkpoint_freq;
            if (block_end > iterations) block_end = iterations;

            // Reset thread segments and arrays (zeroed before each block begins)
            for (int i = 0; i < max_threads; i++) {
                mpf_set_ui(thread_S[i], 0);
            }
        }

        // ------------------ Checkpoint code ends   ---------------------

        #pragma omp parallel shared(global_S, thread_S, CONST_X_BASE, CONST_L_K, CONST_L_ADD, completed_count, show_progress, progress_freq)
        {
            int tid = omp_get_thread_num(); // Get the current thread ID
            ThreadVariables var; // Thread private variables
            init_thread_variables(&var); // Initialize thread variables
            #ifdef ENABLE_BLOCK_FACTORIAL
            var.block_size = block_size; // Set block size for block factorial
            #endif

            #ifdef ENABLE_CACHE
            ThreadCache cache; // Thread var cache
            init_thread_cache(&cache); // Initialize thread cache
            #endif

            // calculate_pi: for (unsigned long k = 0; k < iterations; k++) {
            #pragma omp for schedule(runtime)
            for (unsigned long k = enable_checkpoint ? current_k : 0; k < block_end; k++) {
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

                // Progress display: Atomic counter update (only when progress is enabled)
                if (show_progress) {
                    #pragma omp atomic
                    ++completed_count;

                    if (completed_count % progress_freq == 0 && tid == 0) {
                        // Use the thread ID to determine execution on the main thread
                        fprintf(stderr, "\rProgress: %.2f%%", (double) completed_count / iterations * 100);
                        fflush(stderr);
                    }
                }
            }

            /* This code is deprecated
            // Merge thread private variables into global variables
            #pragma omp critical
            {
                mpf_add(S, S, var.S);
            }
            */


            // Store the segment of this thread into the corresponding array slot
            mpf_set(thread_S[tid], var.S);

            clean_thread_variables(&var); // Clean up thread variables

            #ifdef ENABLE_CACHE
            clean_thread_cache(&cache); // Clean up thread cache
            #endif
        } // End of parallel section

        /* calculate_pi:
            // The main thread merges all parts
            for (int i = 0; i < max_threads; i++) {
                mpf_add(global_S, global_S, thread_S[i]);
                mpf_clear(thread_S[i]);
            }
            free(thread_S);
        */

        // The main thread merges all parts
        for (int i = 0; i < max_threads; i++) {
            mpf_add(global_S, global_S, thread_S[i]);
        }

        // ------------------ Checkpoint code begins ---------------------
        if (enable_checkpoint) {
            current_k = block_end;

            // Save checkpoint

            if (save_checkpoint(checkpoint_file, current_k, global_S, digits,
                                (uint32_t) num_threads, current_flags, quiet_flag) != 0) {
                if (!quiet_flag) fprintf(stderr, "Warning: Failed to save checkpoint\n");
            } else if (show_progress && !quiet_flag) {
                fprintf(stderr, "\nCheckpoint saved at iteration %lu\n", current_k);
            }
        } else {
            break;
        }
    } // End of while section
    // ------------------ Checkpoint code ends   ---------------------

    // Calculate PI = C / S
    mpf_div(pi, C, global_S);

    // Clean up global constants
    clean_constants();

    // Clean the thread_S array
    for (int i = 0; i < max_threads; i++) {
        mpf_clear(thread_S[i]);
    }
    free(thread_S);

    // Clean up variables
    mpf_clears(C, global_S, temp, NULL);

    #ifdef DEBUG
    printf("Cache hit count: %lu\n", cache_hit_count);
    printf("Cache hit ratio: %.2f%%\n", (double) cache_hit_count / (iterations * 2) * 100);
    #endif
}
