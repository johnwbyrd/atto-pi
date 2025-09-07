/*
 * pi.c - Bellard's decimal pi spigot algorithm implementation
 *
 * Author: John Byrd <johnwbyrd at gmail dot com>
 *
 * See README.md for overview of the algorithm and mathematical background
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Arbitrary-precision numbers represented as byte arrays in little-endian
 * format. Lower indices contain fractional part, highest index contains integer
 * part. Precision window [precision_lower, precision_upper] tracks active
 * computation range.
 */
typedef uint8_t *bignum;

/*
 * Lightweight printing functions that avoid printf overhead (~3KB+ on 6502).
 * Critical for fitting within memory constraints of retro systems.
 */
/* Print null-terminated string character by character. */
void print_str(const char *s) {
    while (*s)
        putchar(*s++);
}

/* Print exactly 3 digits with trailing space, used for main digit output. */
void print_three_digit(uint32_t val) {
    putchar('0' + (val / 100) % 10);
    putchar('0' + (val / 10) % 10);
    putchar('0' + val % 10);
    putchar(' ');
}

/* Print unsigned integer without leading zeros. */
void print_uint(uint16_t val) {
    if (val == 0) {
        putchar('0');
        return;
    }
    char buf[6];
    int i = 0;
    do {
        buf[i++] = '0' + (val % 10);
        val /= 10;
    } while (val);
    while (i--)
        putchar(buf[i]);
}

/* Initialize platform-specific optimizations, expanding heap on LLVM-MOS
 * targets. */
void init_platform() {
#ifdef __mos__
    __set_heap_limit(__get_heap_max_safe_size());
#endif
}

/* Global variables for Bellard's pi calculation algorithm.
 * These maintain state across all bignum operations and track the
 * precision boundaries as the calculation progresses. Understanding
 * these variables is crucial for comprehending the algorithm's operation.
 */

/* Divisor used in bignum_div_addsub() for each term of Bellard's formula.
 * Ranges from small values like (10n+1) to large products like (10n+9)*256.
 * The seven terms use different polynomial denominators:
 * - Terms with 10n: (10n+1), (10n+3)*4, (10n+5)*64, (10n+7)*64, (10n+9)*256
 * - Terms with 4n: (4n+1)*8, (4n+3)*256
 * Must be 64-bit to handle the largest denominators without overflow.
 * At 50,000 digits: largest value ≈ (166,670+9)*256 = 42,669,824
 */
uint64_t divisor;

/* Lower index bound for active bignum precision window.
 * Represents the leftmost (least significant) byte still containing
 * meaningful fractional data. As iterations progress, the fractional
 * precision requirements grow, and this bound increments to skip
 * computation on bytes that have become insignificant (leading zeros).
 * Formula: precision_lower = base / 128, where base grows by 159 per iteration.
 * This optimization significantly reduces computation time in later iterations.
 */
int16_t precision_lower;

/* Upper index bound for active bignum precision window.
 * Represents the rightmost (most significant) byte containing non-zero
 * numerator data. The numerator is rescaled by 250/256 each iteration,
 * causing it to shrink and making higher-order bytes become zero.
 * This bound decrements accordingly, reducing unnecessary computation.
 * The precision window [precision_lower, precision_upper] defines the
 * active computation range, optimizing performance as the calculation
 * progresses.
 */
int16_t precision_upper;

/* Size of each bignum array in bytes, calculated from requested digits.
 * Formula: (digits * 415241) / 1000000 + guard_digits
 * The constant 0.415241 = log₂(10)/8, representing bytes needed per decimal
 * digit. For 50,000 digits: ≈ 20,762 + 3 = 20,765 bytes per bignum. This
 * provides sufficient binary precision to accurately represent the requested
 * decimal precision without accumulation errors.
 */
uint16_t bignum_size;

/* Sum accumulator bignum storing the running total of all Bellard formula
 * terms. Structure: [fractional_part...][fractional_part][integer_part]
 * - Integer part (highest byte): Contains 3 digits ready for extraction
 * - Fractional part (lower bytes): Accumulates precision for future iterations
 * After each iteration, integer digits are extracted and masked to zero,
 * then the entire sum is multiplied by 1000 to shift the next 3 digits
 * from fractional to integer position.
 */
bignum sum_accumulator;

/* Numerator bignum used in division operations for each formula term.
 * Initialization: Set to 4 (derived from Bellard's series coefficient)
 * Per-iteration cycle:
 *   1. Divided by each of the 7 term denominators via bignum_div_addsub()
 *   2. Rescaled by factor 250/256 to prevent unbounded growth
 * The 250/256 factor = 1000/1024, compensating for the base-1000 vs base-1024
 * approximation in the decimal adaptation of Bellard's binary series.
 * This rescaling maintains the mathematical precision balance between
 * sum_accumulator and numerator across iterations.
 */
bignum numerator;

/*
 * Bignum arithmetic functions
 * These functions implement arbitrary-precision arithmetic operations
 * on byte arrays representing fixed-point numbers in little-endian format.
 */

/* Initialize a bignum to a small integer value, clearing the fractional part */
void bignum_set(bignum bn, uint32_t value);

/* Perform division and add or subtract the quotient into the sum bignum.
 * This implements the core mathematical operation for each term in Bellard's
 * formula */
void bignum_div_addsub(int is_subtract);

/* Rescale a bignum by the factor 250/256, used to prevent numerator overflow
 * while maintaining precision across iterations */
void bignum_rescale(bignum bn);

/* Zero the integer part of a bignum, removing digits that have been output */
void bignum_mask_digits(bignum bn);

/* Multiply a bignum by 250 with carry propagation.
 * Part of the rescaling operation (250/256 ratio) */
void bignum_multiply_250(bignum bn);

/* Multiply a bignum by 1000 with carry propagation.
 * Shifts the next three decimal digits into the integer position */
void bignum_multiply_1000(bignum bn);

/*
 * Memory management functions
 * Handle allocation and sizing of bignum arrays based on precision requirements
 */

/* Calculate the required size for bignums based on desired digit precision */
uint16_t calculate_bignum_size(uint16_t digits, uint8_t guard);

/* Allocate and initialize both sum and numerator bignums */
void allocate_bignums(uint16_t size);

/* Free allocated bignum memory */
void deallocate_bignums(void);

/*
 * Main algorithm
 */

/* Execute Bellard's formula to calculate pi to the specified precision */
void calculate_pi_bellard(uint16_t digits, uint8_t guard_digits);

/*
 * Initialize a bignum to a small integer value by clearing the fractional
 * part and placing the integer at the highest position. The integer value
 * is written as a 32-bit word at the high end, which may span multiple bytes.
 * This is used to set initial coefficients like 4 for the Bellard formula.
 */
void bignum_set(bignum bn, uint32_t value) {
    /* Clear fractional portion from precision_lower to bignum_size-1 */
    for (uint16_t i = precision_lower; i < bignum_size; i++) {
        bn[i] = 0;
    }
    /* Place integer value at highest position (32-bit write) */
    *(uint32_t *)&bn[bignum_size - 1] = value;
}

/*
 * Core division and accumulation function implementing each term of Bellard's
 * formula. Performs: sum_accumulator ±= (numerator / divisor)
 *
 * Algorithm: Binary long division processing 8 bits at a time
 * 1. For each byte position (most significant to least significant):
 *    - Build remainder by shifting left 8 bits and adding next numerator byte
 *    - Scale divisor by 256 to match remainder magnitude
 *    - Extract quotient bits one by one using comparison and subtraction
 * 2. Add or subtract resulting quotient byte into sum_accumulator with carry
 * propagation
 *
 * The is_subtract parameter controls whether this term is added or subtracted,
 * implementing the alternating signs in Bellard's series.
 */
void bignum_div_addsub(int is_subtract) {
    uint64_t remainder = 0;              /* Running remainder for division */
    uint64_t divisor_original = divisor; /* Preserve original divisor */

    /* Process bignum from most significant to least significant byte */
    for (int16_t i = precision_upper; i >= precision_lower; i--) {
        /* Build up remainder: shift left 8 bits and add next byte */
        remainder = remainder * 256 + numerator[i];
        uint16_t quotient_byte = 0; /* Quotient byte being constructed */
        divisor *= 256;             /* Scale divisor to match remainder */

        /* Binary division: extract 8 bits one at a time */
        for (uint8_t j = 0; j <= 7; j++) {
            quotient_byte *= 2; /* Shift quotient bit left */
            divisor /= 2;       /* Halve divisor for this bit */
            if (remainder >= divisor) {
                remainder -= divisor; /* Remainder -= divisor */
                quotient_byte += 1;   /* Set quotient bit */
            }
        }

        /* Add or subtract quotient byte into sum with carry propagation */
        uint16_t sum_idx = i;
        do {
            int32_t tmp = is_subtract
                              ? (sum_accumulator[sum_idx] - quotient_byte)
                              : (sum_accumulator[sum_idx] + quotient_byte);
            sum_accumulator[sum_idx++] = tmp & 255; /* Store low byte */
            /* Carry/borrow continues if result outside byte range */
            quotient_byte = (tmp >= 0 && tmp <= 255) ? 0 : 1;
        } while (quotient_byte);
    }

    /* Restore original divisor for next term */
    divisor = divisor_original;
}

/*
 * Rescale bignum by factor 250/256 = 1000/1024 to prevent numerator growth.
 * This implements the radix correction factor in Bellard's decimal adaptation.
 *
 * Two-step process:
 * 1. Multiply by 250 using carry propagation
 * 2. Divide by 256 by shifting all bytes down one position (equivalent to >> 8)
 *
 * The 250/256 ratio compensates for the approximation 10³ ≈ 2¹⁰ used in
 * adapting Bellard's binary series for decimal digit extraction.
 */
void bignum_rescale(bignum bn) {
    bignum_multiply_250(bn); /* Multiply by 250 */
    /* Divide by 256: shift all bytes down one position */
    for (uint16_t i = precision_lower; i < bignum_size; i++) {
        bn[i] = bn[i + 1];
    }
    bn[bignum_size] = 0; /* Clear vacated high byte */
}

/*
 * Zero the integer part of bignum after extracting output digits.
 * This removes the digits that have been output, leaving only the
 * fractional part for the next iteration. Clears both the integer
 * byte and one extra byte for safety.
 */
void bignum_mask_digits(bignum bn) {
    uint16_t i = bignum_size - 1;
    bn[i] = 0;     /* Clear integer part */
    bn[i + 1] = 0; /* Clear extra byte for safety */
}

/*
 * Multiply bignum by 250 with carry propagation.
 * This is part of the 250/256 rescaling factor. Each byte is multiplied
 * by 250, with carries propagated to higher bytes. The final carry must
 * be zero to ensure no precision loss (guaranteed by proper sizing).
 */
void bignum_multiply_250(bignum bn) {
    uint32_t temp;
    uint32_t carry = 0;

    /* Process from least to most significant byte */
    for (uint16_t i = precision_lower; i <= bignum_size; i++) {
        temp = (uint32_t)bn[i] * 250 + carry;
        bn[i] = temp & 255; /* Store low byte */
        carry = temp >> 8;  /* Propagate high byte as carry */
    }

    /* Verify no final carry (indicates sufficient precision) */
    if (carry > 0) {
        assert(0 && "bignum_multiply_250 should end with carry = 0");
    }
}

/*
 * Multiply bignum by 1000 to shift next three decimal digits into position.
 * After extracting three digits from the integer part, this operation moves
 * the next three digits from the fractional part into the integer part for
 * the following iteration. Uses carry propagation like multiply_250.
 */
void bignum_multiply_1000(bignum bn) {
    uint32_t temp;
    uint32_t carry = 0;

    /* Process from least to most significant byte */
    for (uint16_t i = precision_lower; i <= bignum_size; i++) {
        temp = (uint32_t)bn[i] * 1000 + carry;
        bn[i] = temp & 255; /* Store low byte */
        carry = temp >> 8;  /* Propagate high byte as carry */
    }

    /* Verify no final carry (indicates sufficient precision) */
    if (carry > 0) {
        assert(0 && "bignum_multiply_1000 should end with carry = 0");
    }
}

/*
 * Main implementation of Bellard's formula for calculating pi.
 * Generates digits using a spigot algorithm that produces three decimal
 * digits per iteration through an alternating series of rational terms.
 */
void calculate_pi_bellard(uint16_t digits, uint8_t guard_digits) {
    /* Initialize memory and precision tracking */
    bignum_size = calculate_bignum_size(digits, guard_digits);
    allocate_bignums(bignum_size);

    /* Set up bignum array bounds: track precision window that shrinks
     * as calculation progresses. Base tracks fractional precision growth
     * in fixed-point arithmetic (scaled by 128 for integer math).
     */
    int32_t base = 0;                  /* Precision tracker (*128) */
    precision_lower = 0;               /* Lower bound starts at 0 */
    precision_upper = bignum_size - 1; /* Upper bound at array end */

    /* Initialize numerator to 4 (coefficient from Bellard's formula) */
    bignum_set(numerator, 4);

    /* Initialize loop counters based on iteration number n:
     * three_n = 3n (0, 3, 6, 9, 12, ...)
     * four_n = 4n (0, 4, 8, 12, 16, ...)
     * ten_n = 10n (0, 10, 20, 30, 40, ...)
     */
    uint16_t three_n = 0; /* 3n counter - tracks total digits produced */
    uint32_t four_n = 0;  /* 4n counter - used in denominators 4n+1, 4n+3 */
    uint64_t ten_n = 0;   /* 10n counter - used in denominators 10n+1, 10n+3,
                             10n+5, 10n+7, 10n+9 */
    uint8_t op = 1;       /* Operation toggle: 1=add, 0=subtract */

    /* Initialize sum accumulator to 4 */
    bignum_set(sum_accumulator, 4);

    uint16_t digit_counter = 0; /* Track digits output for reporting */

    /* Main calculation loop: each iteration produces 3 digits */
    do {
        /* Compute the seven terms of Bellard's decimal spigot formula:
         * Each term has form: coefficient / (polynomial_in_n * scale_factor)
         * The alternating signs and specific denominators derive from Bellard's
         * adaptation of his binary BBP formula for decimal digit extraction.
         */

        /* Term 1: -32/(10n+1) - but coefficient absorbed into numerator */
        divisor = ten_n + 1;
        if (three_n >
            0) { /* Skip first iteration to avoid double-counting initial 4 */
            bignum_div_addsub(1 - op);
        }

        /* Term 2: -1/((10n+3)*4) = -1/(4*(10n+3)) */
        divisor = (ten_n + 3) * 4;
        bignum_div_addsub(op);

        /* Term 3: +256/((10n+5)*64) = +4/(10n+5) */
        divisor = (ten_n + 5) * 64;
        bignum_div_addsub(op);

        /* Term 4: -64/((10n+7)*64) = -1/(10n+7) */
        divisor = (ten_n + 7) * 64;
        bignum_div_addsub(op);

        /* Term 5: -4/((10n+9)*256) = -1/(64*(10n+9)) */
        divisor = (ten_n + 9) * 256;
        bignum_div_addsub(1 - op);

        /* Term 6: -4/((4n+1)*8) = -1/(2*(4n+1)) */
        divisor = (four_n + 1) * 8;
        bignum_div_addsub(op);

        /* Term 7: +1/((4n+3)*256) */
        divisor = (four_n + 3) * 256;
        bignum_div_addsub(op);

        /* Extract three digits from sum's integer part */
        uint32_t digit_val = *(uint32_t *)&sum_accumulator[bignum_size - 1];

        /* Output digits with appropriate formatting */
        if (three_n > 0) {
            print_three_digit(digit_val);
            digit_counter += 3; /* Each iteration produces 3 digits */
        } else {
            print_uint(digit_val);
            putchar('.');
            putchar('\n');      /* First digit with decimal */
            digit_counter += 1; /* First iteration produces 1 digit */
        }

        /* Prepare for next iteration */
        bignum_mask_digits(sum_accumulator);   /* Remove extracted digits */
        bignum_multiply_1000(sum_accumulator); /* Shift next 3 digits up */
        bignum_rescale(numerator);             /* Apply 250/256 scaling */

        /* Toggle operation sign for alternating series */
        op = 1 - op;

        /* Adjust precision tracking as numerator shrinks from rescaling */
        if (numerator[precision_upper] == 0) {
            precision_upper = precision_upper - 1;
        }

        /* Update loop variables for next iteration.
         * Base grows by 159/128 ≈ 1.2422, empirically determined to match
         * the rate at which precision requirements increase per 3 digits.
         * Each counter advances by its respective step size.
         */
        base += 159;
        precision_lower = base / 128; /* Convert to byte index */
        four_n += 4;                  /* Advance 4n counter: 0→4→8→12... */
        ten_n += 10;                  /* Advance 10n counter: 0→10→20→30... */
        three_n += 3;                 /* Advance 3n counter: 0→3→6→9... */

    } while (digit_counter < digits);

    /* Clean up and report results */
    deallocate_bignums();
    print_str("\nCalculated ");
    print_uint(digit_counter);
    print_str(" digits of pi.\n");
}

/*
 * Each decimal digit requires log2(10) bits of binary precision,
 * approximately 3.321928. For d decimal digits, the total bits needed equal
 * d times log2(10). The bytes needed equal that total divided by 8, or
 * approximately d times 0.415241.
 *
 * Bellard's series converges such that higher digits require more precision.
 * The bignum must maintain enough fractional bits to avoid accumulation
 * errors that compound across iterations. When M, the precision upper bound,
 * drops below L, the precision lower bound, divisions cease, leading to
 * stagnant or zero output.
 *
 * The constant 0.4152410118609203 is precomputed as log2(10) divided by 8.
 * No runtime log2() call is made, which is important for the 6502 target
 * without a math library. Size includes padding for safe 32-bit reads during
 * digit extraction. Memory scales linearly, requiring approximately 20,765
 * bytes for 50,000 digits compared to the old formula's approximately 20,836.
 */
uint16_t calculate_bignum_size(uint16_t digits, uint8_t guard) {
    return (uint16_t)((digits * 415241LL) / 1000000 + guard);
}

/*
 * Allocate memory for both bignum arrays with error checking.
 * Includes padding bytes to ensure safe 32-bit read operations when
 * extracting output digits. Initializes all bytes to zero and exits
 * with error message if allocation fails.
 */
void allocate_bignums(uint16_t size) {
    const uint16_t pad = 4; /* Padding for safe 32-bit reads */

    /* Allocate sum accumulator with error checking */
    sum_accumulator = (bignum)malloc(size + pad);
    if (!sum_accumulator) {
        print_str("Error: Failed to alloc sum_accumulator: ");
        print_uint(size + pad);
        exit(1);
    }

    /* Allocate numerator with error checking */
    numerator = (bignum)malloc(size + pad);
    if (!numerator) {
        print_str("Failed to alloc numerator: ");
        print_uint(size + pad);
        free(sum_accumulator); /* Clean up partial allocation */
        exit(1);
    }

    /* Initialize both arrays to zero */
    for (uint16_t i = 0; i < size + pad; i++) {
        sum_accumulator[i] = 0;
        numerator[i] = 0;
    }
}

/*
 * Free allocated bignum memory.
 * Called at end of calculation to prevent memory leaks.
 */
void deallocate_bignums(void) {
    free(sum_accumulator);
    free(numerator);
}

int main(void) {
    uint16_t digits = 50000;
    uint8_t guard = 3;

    /* Initialize platform-specific heap expansion */
    init_platform();

    print_str("Calculating pi to ");
    print_uint(digits);
    print_str("+ digits...\n");
    calculate_pi_bellard(digits, guard);

    return 0;
}
