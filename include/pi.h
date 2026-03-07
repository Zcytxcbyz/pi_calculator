#ifndef PI_H
#define PI_H

#include <gmp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef ENABLE_BLOCK_FACTORIAL
#define VAR_BLOCK_SIZE , unsigned long block_size
#else
#define VAR_BLOCK_SIZE
#endif

// Calculate PI to the specified number of digits
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads, const char *omp_schedule,
    int chunk_size VAR_BLOCK_SIZE, bool show_progress, int progress_freq, bool quiet_flag,
    bool enable_checkpoint, unsigned long checkpoint_freq, const char *checkpoint_file, bool checkpoint_verbose);

// Write the PI value to file
void write_pi_to_file(const mpf_t pi, unsigned long digits, const char* filename, double computation_time,
    bool format_output, size_t buffer_size, bool raw_output);

// Write the PI value to stream
void write_pi_to_stream(const mpf_t pi, unsigned long digits, FILE* stream, double computation_time,
    bool format_output, size_t buffer_size, bool raw_output);


#ifdef USE_DEPRECATED
/* ---------------------------------------------------------------------------------------------------
 * A π calculation function with checkpointing. Used only when --checkpoint-enable is enabled.
 * Parameters are similar to calculate_pi, but include additional checkpoint_freq and checkpoint_file.
 */
void calculate_pi_checkpoint(mpf_t pi, unsigned long digits, int num_threads, const char* omp_schedule,
    int chunk_size VAR_BLOCK_SIZE, bool show_progress, int progress_freq, bool quiet_flag,
    unsigned long checkpoint_freq, const char* checkpoint_file);
#endif

#endif // PI_H
