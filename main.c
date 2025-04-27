#include "pi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

void print_usage(const char* program_name) {
    printf("Usage: %s [options]\n", program_name);
    printf("Options:\n");
    printf("  -d <digits>    Number of digits to calculate (default: 1000)\n");
    printf("  -o <filename>  Output file name (default: pi.txt)\n");
    printf("  -t <threads>   Number of threads to use (default: number of CPU cores)\n");
    printf("  -f             Format output (default: unformatted)\n");
    printf("  -c             Disable output file\n");
    printf("  -b <size>      Set buffer size in bytes (default: 65536)\n");
    printf("  -s <schedule>  Set OpenMP schedule type (static, dynamic, guided) and chunk size (default: guided)\n");
    printf("  -h             Show this help message\n");
}

int main(int argc, char* argv[]) {
    unsigned long digits = 1000;
    char* output_file = "pi.txt";
    int num_threads = omp_get_max_threads();
    bool enable_output = true;     // Default to enabled output
    bool format_output = false;    // Default to unformatted output
    size_t buffer_size = 65536;    // Default buffer size
    char* omp_schedule = "guided"; // Default OpenMP schedule type
    int chunk_size = 1;            // Default chunk size

    // Analyze command-line parameters
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 && i + 1 < argc) {
            digits = strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-f") == 0) {
            format_output = true;
        } else if (strcmp(argv[i], "-c") == 0) {
            enable_output = 0;
        } else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            buffer_size = strtoul(argv[++i], NULL, 10);
            if (buffer_size < 1024) {
                fprintf(stderr, "Buffer size must be at least 1024 bytes.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            char* schedule_arg = argv[++i];
            char* comma_pos = strchr(schedule_arg, ',');
            if (comma_pos) {
                *comma_pos = '\0'; // Split the string into two parts
                omp_schedule = schedule_arg;
                chunk_size = atoi(comma_pos + 1);
                if (chunk_size < 0) {
                    fprintf(stderr, "Invalid chunk size: %d\n", chunk_size);
                    return 1;
                }
            } else {
                omp_schedule = schedule_arg;
                chunk_size = 1; // Default chunk size
            }

            if (strcmp(omp_schedule, "static") != 0 &&
                strcmp(omp_schedule, "dynamic") != 0 &&
                strcmp(omp_schedule, "guided") != 0) {
                fprintf(stderr, "Invalid OpenMP schedule type: %s\n", omp_schedule);
                return 1;
            }
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    printf("Calculating pi to %lu digits using %d threads...\n", digits, num_threads);

    mpf_t pi;
    mpf_init2(pi, (digits + 2) * log2(10)); // Pre allocate sufficient precision

    double start_time = omp_get_wtime();
    calculate_pi(pi, digits, num_threads, omp_schedule, chunk_size);
    double end_time = omp_get_wtime();

    double total_time = end_time - start_time;
    printf("Total time: %.2f seconds\n", total_time);

    if(enable_output) {
        write_pi_to_file(pi, digits, output_file, total_time, format_output, buffer_size);
        printf("Result written to %s\n", output_file);
    }

    mpf_clear(pi);

    return 0;
}
