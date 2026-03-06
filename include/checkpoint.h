#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>

// Flag Bit Definition
#define CHECKPOINT_FLAG_CACHE           (1U << 0)   // Enable Caching
#define CHECKPOINT_FLAG_BLOCK_FACTORIAL (1U << 1)   // Enable block factorial

// Save checkpoint
int save_checkpoint(const char* filename, unsigned long completed_k, const mpf_t global_S,
    unsigned long digits, uint32_t num_threads, uint32_t flags, bool quiet_flag);

// Load checkpoint
int load_checkpoint(const char* filename, unsigned long* completed_k, mpf_t global_S,
    unsigned long req_digits, uint32_t* num_threads, uint32_t* flags, bool quiet_flag);

#endif
