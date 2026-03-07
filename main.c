#include "pi.h"
#include "pi_ref.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>

void print_usage(const char* program_name) {
    printf("%s Version %s\n", program_name, PROJECT_VERSION);
    printf("Usage: %s [options]\n", program_name);
    printf("Options:\n");
    printf("  -d(--digits) <digits>             Number of digits to calculate (default: 1000)\n");
    printf("  -o(--output) <filename>           Output file name (default: pi.txt)\n");
    printf("  -t(--thread) <threads>            Number of threads to use (default: number of CPU cores)\n");
    printf("  -f(--format)                      Format output (default: unformatted)\n");
    printf("  --disable-output                  Disable output file\n");
    printf("  --buffer-size <size>              Set buffer size in bytes (default: 65536)\n");
    printf("  --schedule <schedule>             Set OpenMP schedule type (static, dynamic, guided) and chunk size (default: guided)\n");
    #ifdef ENABLE_BLOCK_FACTORIAL
    printf("  --block-size <size>               Set block size for factorial calculation (default: 8)\n");
    #endif
    printf("  --raw                             Output raw digits only (no header, no \"3.\" line, no formatting)\n");
    printf("  --quiet                           Suppress all informational output (errors still go to stderr)\n");
    printf("  --stdout                          Write result to standard output instead of a file (overrides -o)\n");
    printf("  --progress                        Show progress during long calculations (disabled by --quiet)\n");
    printf("  --progress-freq <num>             Update progress every <num> iterations (default: 1000, only with --progress)\n");
    printf("  --time-file <filename>            Write computation time to a separate file (even with --quiet)\n");
    printf("  --verify                          Verify first 1000 digits of result against known value (exit code 2 if mismatch)\n");
    printf("  --checkpoint-enable               Enable checkpoint/restart functionality\n");
    printf("  --checkpoint-freq <N>             Save checkpoint every N iterations (default: 1000)\n");
    printf("  --checkpoint-file <filename>      Path to checkpoint file (default: pi_checkpoint.dat)\n");
    printf("  -v(--version)                     Show program version and exit\n");
    printf("  -h(--help)                        Show this help message\n");
}

int main(int argc, char* argv[]) {
    unsigned long digits = 1000;
    char* output_file = "pi.txt";
    int num_threads = omp_get_max_threads();
    bool enable_output = true;                      // Default to enabled output
    bool format_output = false;                     // Default to unformatted output
    size_t buffer_size = 65536;                     // Default buffer size
    char* omp_schedule = "guided";                  // Default OpenMP schedule type
    int chunk_size = 1;                             // Default chunk size
    bool raw_output = false;                        // flag for --raw
    bool quiet_flag = false;                        // flag for --quiet
    bool stdout_flag = false;                       // flag for --stdout
    bool progress_flag = false;                     // flag for --progress
    int progress_freq = 1000;                       // default frequency for progress updates
    char* time_file = NULL;                         // flag for --time-file
    bool verify_flag = false;                       // flag for --verify
    #ifdef ENABLE_BLOCK_FACTORIAL
    unsigned long block_size = 8;                   // Default block size for factorial
    #endif
    bool checkpoint_enable = false;                 // flag for --checkpoint-enable
    unsigned long checkpoint_freq = 1000;           // By default, save every 1000 iterations.
    char* checkpoint_file = "pi_checkpoint.dat";    // Checkpoint File

    // Analyze command-line parameters
    for (int i = 1; i < argc; i++) {
        if ((strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--digits") == 0) && i + 1 < argc) {
            digits = strtoul(argv[++i], NULL, 10);
        } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            output_file = argv[++i];
        } else if ((strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--thread") == 0) && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--format") == 0) {
            format_output = true;
        } else if (strcmp(argv[i], "--disable-output") == 0) {
            enable_output = false;
        } else if (strcmp(argv[i], "--buffer-size") == 0 && i + 1 < argc) {
            buffer_size = strtoul(argv[++i], NULL, 10);
            if (buffer_size < 1024) {
                fprintf(stderr, "Buffer size must be at least 1024 bytes.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--schedule") == 0 && i + 1 < argc) {
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
        #ifdef ENABLE_BLOCK_FACTORIAL
        } else if (strcmp(argv[i], "--block-size") == 0 && i + 1 < argc) {
            block_size = strtoul(argv[++i], NULL, 10);
            if (block_size < 1) {
                fprintf(stderr, "Block size must be at least 1.\n");
                return 1;
            }
        #endif
        } else if (strcmp(argv[i], "--raw") == 0) {
            raw_output = true;
        } else if (strcmp(argv[i], "--quiet") == 0) {
            quiet_flag = true;
        } else if (strcmp(argv[i], "--stdout") == 0) {
            stdout_flag = true;
        } else if (strcmp(argv[i], "--progress") == 0) {
            progress_flag = true;
        } else if (strcmp(argv[i], "--progress-freq") == 0 && i + 1 < argc) {
            int freq = atoi(argv[++i]);
            if (freq <= 0) {
                fprintf(stderr, "Error: progress frequency must be positive.\n");
                return 1;
            }
            progress_freq = freq;
        } else if (strcmp(argv[i], "--time-file") == 0 && i + 1 < argc) {
            time_file = argv[++i];
        } else if (strcmp(argv[i], "--verify") == 0) {
            verify_flag = true;
        } else if (strcmp(argv[i], "--checkpoint-enable") == 0) {
            checkpoint_enable = true;
        } else if (strcmp(argv[i], "--checkpoint-freq") == 0 && i + 1 < argc) {
            checkpoint_freq = strtoul(argv[++i], NULL, 10);
            if (checkpoint_freq == 0) {
                fprintf(stderr, "Error: checkpoint frequency must be positive.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--checkpoint-file") == 0 && i + 1 < argc) {
            checkpoint_file = argv[++i];
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            printf("pi_calculator version %s\n", PROJECT_VERSION);
            return 0;
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    // If both --stdout and -o are specified (and the filename is not the default), issue a warning
    if (stdout_flag && strcmp(output_file, "pi.txt") != 0 && !quiet_flag) {
        fprintf(stderr, "Warning: --stdout overrides -o, ignoring output file.\n");
    }
    // --disable-output takes precedence over --stdout
    if (!enable_output && stdout_flag) {
        if (!quiet_flag) fprintf(stderr, "Warning: --disable-output is set, ignoring --stdout.\n");
        stdout_flag = 0;
    }
    // Warning when outputting large numbers to stdout
    if (stdout_flag && digits > 100000 && !quiet_flag) {
        fprintf(stderr, "Warning: printing %lu digits to stdout may cause terminal slowdown. Consider redirecting to a file.\n", digits);
    }

    if (!quiet_flag) {
        printf("Calculating pi to %lu digits using %d threads...\n", digits, num_threads);
    }

    mpf_t pi;
    mpf_init2(pi, (digits + 2) * log2(10)); // Pre allocate sufficient precision

    double start_time = omp_get_wtime();

    bool show_progress = progress_flag && !quiet_flag;  // Display only when not in silent mode and progress is enabled.

    #ifdef ENABLE_BLOCK_FACTORIAL
    calculate_pi(pi, digits, num_threads, omp_schedule, chunk_size, block_size,
        show_progress, progress_freq, quiet_flag,
        checkpoint_enable, checkpoint_freq, checkpoint_file);
    #else
    calculate_pi(pi, digits, num_threads, omp_schedule, chunk_size,
        show_progress, progress_freq, quiet_flag,
        checkpoint_enable, checkpoint_freq, checkpoint_file);
    #endif

    /* This code is deprecated
    #ifdef ENABLE_BLOCK_FACTORIAL
    calculate_pi(pi, digits, num_threads, omp_schedule, chunk_size, block_size, show_progress, progress_freq);
    #else
    calculate_pi(pi, digits, num_threads, omp_schedule, chunk_size, show_progress, progress_freq);
    #endif
    */

    double end_time = omp_get_wtime();

    double total_time = end_time - start_time;
    if (!quiet_flag) {
        printf("\nTotal time: %.2f seconds\n", total_time);
    }

    if (enable_output) {
        if (stdout_flag) {
            // Output to stdout
            write_pi_to_stream(pi, digits, stdout, total_time, format_output, buffer_size, raw_output);
            if (!quiet_flag) {
                fflush(stdout);
                fprintf(stderr, "\nResult written to stdout\n");
            }
        } else {
            // Output to file
            write_pi_to_file(pi, digits, output_file, total_time, format_output, buffer_size, raw_output);
            if (!quiet_flag) {
                printf("Result written to %s\n", output_file);
            }
        }
    }

    if (time_file) {
        FILE* tf = fopen(time_file, "w");
        if (!tf) {
            perror("Failed to open time file");
        } else {
            fprintf(tf, "Total time: %.2f seconds\n", total_time);
            fclose(tf);
            if (!quiet_flag) {
                printf("Time written to %s\n", time_file);
            }
        }
    }

    // Verification (if requested)
    if (verify_flag && digits >= 1000) {
        mp_exp_t exp;
        char* pi_str = mpf_get_str(NULL, &exp, 10, digits + 2, pi);
        if (!pi_str || exp != 1) {
            fprintf(stderr, "Verification failed: cannot obtain string representation.\n");
        } else {
            // Construct the complete string “3.” + the decimal part
            char computed[1003];
            computed[0] = '3';
            computed[1] = '.';
            strncpy(computed + 2, pi_str + 1, 1000);
            computed[1002] = '\0';

            if (strcmp(computed, KNOWN_PI_1000) == 0) {
                if (!quiet_flag) printf("Verification passed: first 1000 digits match known value.\n");
            } else {
                fprintf(stderr, "Verification FAILED: first 1000 digits do NOT match known value.\n");
                fprintf(stderr, "Computed: %.1000s\n", computed);
                fprintf(stderr, "Expected: %.1000s\n", KNOWN_PI_1000);
                mpf_clear(pi);
                free(pi_str);
                return 2;
            }
            free(pi_str);
        }
    } else if (verify_flag) {
        fprintf(stderr, "Verification requires at least 1000 digits (current: %lu). Skipping.\n", digits);
    }


    mpf_clear(pi);

    return 0;
}
// checkpoint
