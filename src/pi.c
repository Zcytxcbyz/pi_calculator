#include "pi.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Chudnovsky algorithm calculates PI
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads) {
    mpf_set_default_prec((digits + 2) * log2(10)); // Set sufficient precision

    mpf_t C, M, L, X, S, term, temp;
    mpz_t K, six_k, three_k, k_fact, three_k_fact, six_k_fact;

    // Initialize variable
    mpf_inits(C, M, L, X, S, term, temp, NULL);
    mpz_inits(K, six_k, three_k, k_fact, three_k_fact, six_k_fact, NULL);

    // Constant C = 426880 * sqrt(10005)
    mpf_set_ui(C, 426880);
    mpf_sqrt_ui(temp, 10005);
    mpf_mul(C, C, temp);

    mpf_set_ui(S, 0);

    // Calculate the required number of iterations (empirical formula)
    unsigned long iterations = (digits / 14) + 1;

    // Set the number of OpenMP threads
    omp_set_num_threads(num_threads);

    mpf_t S_private;
    mpf_init_set_ui(S_private, 0);

    #pragma omp parallel private(K, six_k, three_k, k_fact, three_k_fact, six_k_fact, term, temp, M, L, X)
    {
        mpf_t S_local, term, temp_f;
        mpz_t temp, M, L, X, K, six_k, three_k, k_fact, three_k_fact, six_k_fact;

        mpf_init(S_local);
        mpf_set_ui(S_local, 0);
        mpf_inits(term, temp_f, NULL);
        mpz_inits(temp, M, L, X, K, six_k, three_k, k_fact, three_k_fact, six_k_fact, NULL);

        #pragma omp for schedule(dynamic, 10)
        for (unsigned long k = 0; k < iterations; k++) {
            mpz_set_ui(K, k);

            // Calculate M = (6k)! / ((3k)! * (k!)^3)
            mpz_fac_ui(six_k_fact, 6 * k);
            mpz_fac_ui(three_k_fact, 3 * k);
            mpz_fac_ui(k_fact, k);

            mpz_pow_ui(temp, k_fact, 3);
            mpz_mul(temp, temp, three_k_fact);
            mpz_divexact(M, six_k_fact, temp);

            // Calculate L = 545140134k + 13591409
            mpz_set_ui(temp, 545140134);
            mpz_mul_ui(temp, temp, k);
            mpz_add_ui(L, temp, 13591409);

            // Calculate X = (-262537412640768000)^k
            if (k == 0) {
                mpz_set_ui(X, 1);
            } else {
                mpz_set_str(X, "-262537412640768000", 10);
                mpz_pow_ui(X, X, k);
            }

            // Calculate the current item: M * L / X
            mpf_set_z(term, M);
            mpf_set_z(temp_f, L);
            mpf_mul(term, term, temp_f);
            mpf_set_z(temp_f, X);
            mpf_div(term, term, temp_f);

            // Accumulate to thread private variables
            mpf_add(S_local, S_local, term);
        }

        // Clean up temporary variables
        mpf_clears(term, temp_f, NULL);
        mpz_clears(temp, M, L, X, K, six_k, three_k, k_fact, three_k_fact, six_k_fact, NULL);

        // Merge thread private variables into global variables
        #pragma omp critical
        {
            mpf_add(S_private, S_private, S_local);
        }

        mpf_clear(S_local);
    }

    // Assign the result to S
    mpf_set(S, S_private);
    mpf_clear(S_private);

    // Calculate PI = C / S
    mpf_div(pi, C, S);

    // Clean up variables
    mpf_clears(C, M, L, X, S, term, temp, NULL);
    mpz_clears(K, six_k, three_k, k_fact, three_k_fact, six_k_fact, NULL);
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
