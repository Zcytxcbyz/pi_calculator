#include "checkpoint.h"
#include <math.h>
#include <string.h>

// Head Structure (Internal Use)
typedef struct {
    char     magic[4];          // "PICK"
    uint8_t  version;           // The current value is 1
    uint8_t  reserved[3];       // Alignment
    uint64_t digits;            // Target digit
    uint32_t num_threads;       // Number of threads
    uint32_t flags;             // Compiler Option Flags
    uint8_t  future[32];        // reserved
} checkpoint_header_t;

// Adjust precision to suit current computational requirements
static int mpf_out_raw(FILE *fp, const mpf_t x) {
    // Write exp
    long exp = x->_mp_exp;
    if (fwrite(&exp, sizeof(exp), 1, fp) != 1)
        return 0;

    // Write the tail number disguised as mpz_t
    mpz_t mantissa;
    mantissa->_mp_alloc = x->_mp_size;
    mantissa->_mp_size = x->_mp_size;
    mantissa->_mp_d = x->_mp_d;

    return mpz_out_raw(fp, mantissa) + sizeof(exp);
}

// Custom mpf_t binary input
static int mpf_inp_raw(FILE *fp, mpf_t x) {
    long exp;
    if (fread(&exp, sizeof(exp), 1, fp) != 1)
        return 0;

    mpz_t mantissa;
    mpz_init(mantissa);
    if (mpz_inp_raw(mantissa, fp) == 0) {
        mpz_clear(mantissa);
        return 0;
    }
    mpf_set_z(x, mantissa);
    x->_mp_exp = exp;
    mpz_clear(mantissa);
    return 1;
}

// Save checkpoint
int save_checkpoint(const char* filename, unsigned long completed_k, const mpf_t global_S,
    unsigned long digits, uint32_t num_threads, uint32_t flags, bool quiet_flag) {
    FILE* fp = fopen(filename, "wb");
    if (!fp) {
        if (!quiet_flag) perror("Warning: Failed to open checkpoint file for writing");
        return -1;
    }

    // Write to header
    checkpoint_header_t header;
    memcpy(header.magic, "PICK", 4);
    header.version = 1;
    memset(header.reserved, 0, sizeof(header.reserved));
    header.digits = digits;
    header.num_threads = num_threads;
    header.flags = flags;
    memset(header.future, 0, sizeof(header.future));

    if (fwrite(&header, sizeof(header), 1, fp) != 1) {
        if (!quiet_flag) perror("Warning: Failed to write checkpoint header");
        fclose(fp);
        return -1;
    }

    // Number of iterations completed
    if (fwrite(&completed_k, sizeof(completed_k), 1, fp) != 1) {
        if (!quiet_flag) perror("Warning: Failed to write completed_k to checkpoint");
        fclose(fp);
        return -1;
    }

    // Write to the global section and (GMP raw format)
    if (mpf_out_raw(fp, global_S) == 0) {
        if (!quiet_flag) fprintf(stderr, "Warning: Failed to write global_S to checkpoint\n");
        fclose(fp);
        return -1;
    }

    fclose(fp);
    return 0;
}

// Load checkpoint
int load_checkpoint(const char* filename, unsigned long* completed_k, mpf_t global_S,
    unsigned long req_digits, uint32_t* num_threads, uint32_t* flags, bool quiet_flag) {
    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        // File not found is not an error; it is handled by the caller.
        return -1;
    }

    // Read Header
    checkpoint_header_t header;
    if (fread(&header, sizeof(header), 1, fp) != 1) {
        if (!quiet_flag) fprintf(stderr, "Warning: Failed to read checkpoint header\n");
        fclose(fp);
        return -2;
    }

    // Verify Magic Number and Version
    if (memcmp(header.magic, "PICK", 4) != 0 || header.version != 1) {
        if (!quiet_flag) fprintf(stderr, "Warning: Invalid checkpoint file (magic/version mismatch)\n");
        fclose(fp);
        return -2;
    }

    // Verify that the number of digits matches
    if (header.digits != req_digits) {
        if (!quiet_flag) {
            fprintf(stderr, "Warning: Checkpoint digits mismatch (saved=%lu, requested=%lu)\n",
                    header.digits, req_digits);
        }
        fclose(fp);
        return -2;
    }

    // Number of output threads and flags
    if (num_threads) *num_threads = header.num_threads;
    if (flags) *flags = header.flags;

    // Read the number of completed iterations
    if (fread(completed_k, sizeof(*completed_k), 1, fp) != 1) {
        if (!quiet_flag) fprintf(stderr, "Warning: Failed to read completed_k from checkpoint\n");
        fclose(fp);
        return -2;
    }

    // Read the global section
    if (mpf_inp_raw(fp, global_S) == 0) {
        if (!quiet_flag) fprintf(stderr, "Warning: Failed to read global_S from checkpoint\n");
        fclose(fp);
        return -2;
    }

    fclose(fp);

    // Adjust precision to suit current computational requirements
    mpf_set_prec(global_S, (req_digits + 2) * log2(10));

    return 0;
}
