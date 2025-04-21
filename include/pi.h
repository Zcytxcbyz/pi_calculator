#ifndef PI_H
#define PI_H

#include <gmp.h>
#include <math.h>
#include <time.h>

// Calculate PI to the specified number of digits
void calculate_pi(mpf_t pi, unsigned long digits, int num_threads);

// Write the PI value to a file
void write_pi_to_file(const mpf_t pi, unsigned long digits, const char* filename, double computation_time);

#endif // PI_H
