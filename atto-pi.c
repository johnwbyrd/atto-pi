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
 * Precision scaling system: MAX_DIGITS_LOG10 determines maximum supported
 * digit count as 10^MAX_DIGITS_LOG10. Types are automatically selected to
 * handle the specified order of magnitude without overflow.
 */
#ifndef MAX_DIGITS_LOG10
#define MAX_DIGITS_LOG10 5 // Default: up to 10^5 = 100,000 digits
#endif

#if MAX_DIGITS_LOG10 > 9
#error "MAX_DIGITS_LOG10 > 9 (billion digits) requires 128-bit integer support"
#endif

/*
 * Arbitrary-precision numbers represented as byte arrays in little-endian
 * format. Lower indices contain fractional part, highest index contains integer
 * part. Precision window [precision_lower, precision_upper] tracks active
 * computation range.
 */
typedef uint8_t *bignum;

/*
 * Semantic type system based on algorithmic purpose rather than bit width.
 * Types are defined by their mathematical role in Bellard's algorithm,
 * allowing implementers to understand constraints and optimize for their
 * platform.
 */

/* - Digit specification and output */

/* Target and actual digit counts - the fundamental constraint that drives all
 * other requirements. Current: uint16_t allows ~65,535 digits max due to
 * iteration counter limits. Platforms: 8-bit systems practical max ~50,000
 * digits (memory), 16-bit+ can use uint32_t for millions.
 */
#if MAX_DIGITS_LOG10 <= 4
typedef uint16_t digit_count_t;
#elif MAX_DIGITS_LOG10 <= 9
typedef uint32_t digit_count_t;
#else
typedef uint64_t digit_count_t;
#endif

/* 3-digit output progression counter (0,3,6,9,12...) tracking total digits
 * produced. Range: 0 to (digit_count_t * 3). Must handle digit_count * 3
 * without overflow. Current: uint16_t sufficient for current digit_count_t
 * limits.
 */
#if MAX_DIGITS_LOG10 <= 4
typedef uint16_t digit_progression_t;
#elif MAX_DIGITS_LOG10 <= 9
typedef uint32_t digit_progression_t;
#else
typedef uint64_t digit_progression_t;
#endif

/* Individual 3-digit output values extracted from bignum integer part.
 * Range: 0-999 (three decimal digits). Only needs to hold 10-bit values.
 * Current: uint32_t provides safe headroom and matches *(uint32_t*) extraction.
 */
typedef uint32_t digit_value_t;

/* - Iteration mathematics */

/* 4n iteration counter for Bellard's formula denominators (4n+1), (4n+3).
 * Range: 0 to (digit_count * 4/3), steps of 4. At 50,000 digits: max ~66,664.
 * Used in expressions: (4n+1)*8, (4n+3)*256.
 * Current: uint32_t handles all realistic digit counts.
 */
#if MAX_DIGITS_LOG10 <= 6
typedef uint32_t iter_4n_t;
#else
typedef uint64_t iter_4n_t;
#endif

/* 10n iteration counter for denominators (10n+1), (10n+3), (10n+5), (10n+7),
 * (10n+9). Range: 0 to (digit_count * 10/3), steps of 10. At 50,000 digits: max
 * ~166,660. Used in largest expression: (10n+9)*256 = ~42.7M at max digits.
 * Current: uint64_t prevents overflow in all denominator calculations.
 */
typedef uint64_t iter_10n_t;

/* - Precision window management */

/* Precision window bounds [precision_lower, precision_upper] for bignum
 * optimization. Range: 0 to bignum_size-1, but can be negative during
 * comparisons. precision_lower: grows ~1.24x per iteration. precision_upper:
 * shrinks as numerator rescales. Current: int16_t signed to prevent comparison
 * bugs on 16-bit systems.
 */
#if MAX_DIGITS_LOG10 <= 4
typedef int16_t precision_bound_t;
#elif MAX_DIGITS_LOG10 <= 9
typedef int32_t precision_bound_t;
#else
typedef int64_t precision_bound_t;
#endif

/* Fixed-point precision tracker (scaled by 128 for integer math).
 * Range: 0 to ~(digit_count * 159). Growth: +159 per iteration.
 * Formula: precision_lower = precision_base / 128. This optimization
 * signficantly reduces computation time in later iterations.
 * Current: int32_t handles growth rate across all practical digit counts.
 */
#if MAX_DIGITS_LOG10 <= 7
typedef int32_t precision_base_t;
#else
typedef int64_t precision_base_t;
#endif

/* Guard digits for computational headroom preventing accumulation errors.
 * Range: typically 3-10. More guard digits = more memory but safer computation.
 * Current: uint8_t sufficient for any reasonable guard count.
 */
typedef uint8_t guard_count_t;

/* - Bignum array management */

/* Total size of each bignum array in bytes.
 * Formula: (digits * 415241) / 1000000 + guard_digits.
 * Range: ~418 bytes at 1K digits to ~20,765 bytes at 50K digits.
 * Current: uint16_t limits to 65,535 bytes - use uint32_t for larger scales.
 */
#if MAX_DIGITS_LOG10 <= 4
typedef uint16_t array_size_t;
#elif MAX_DIGITS_LOG10 <= 9
typedef uint32_t array_size_t;
#else
typedef uint64_t array_size_t;
#endif

/* Array indices for bignum operations.
 * Range: 0 to array_size_t-1. Signed prevents index comparison bugs.
 * Current: int16_t matches precision_bound_t for consistency.
 */
#if MAX_DIGITS_LOG10 <= 4
typedef int16_t array_index_t;
#elif MAX_DIGITS_LOG10 <= 9
typedef int32_t array_index_t;
#else
typedef int64_t array_index_t;
#endif

/* - Division engine arithmetic */

/* Denominators in Bellard's seven-term formula.
 * Range: (10n+1) minimum to (10n+9)*256 maximum. At 50K digits: max ~42.7M.
 * Must handle all denominator expressions: (10n+k)*scale_factor.
 * Current: uint64_t prevents overflow in largest expressions.
 */
typedef uint64_t divisor_t;

/* Running remainder during binary long division.
 * Range: 0 to divisor-1 during computation. Must match divisor_t capacity.
 * Current: uint64_t matches divisor_t requirements.
 */
typedef uint64_t remainder_t;

/* Individual quotient byte constructed during 8-bit division.
 * Range: 0-255 (single byte result from binary long division).
 * Current: uint16_t provides headroom for intermediate calculations.
 */
typedef uint16_t quotient_byte_t;

/* - Carry propagation arithmetic */

/* Carry values during multiply-by-250 and multiply-by-1000 operations.
 * Range: 0 to max(255*250, 255*1000) = 255,000.
 * Must handle: byte * constant + previous_carry.
 * Current: uint32_t handles worst case with headroom.
 */
typedef uint32_t carry_t;

/* Temporary storage for intermediate arithmetic (byte * constant).
 * Range: 0 to 255*1000 = 255,000. Used in: temp = byte * multiplier + carry.
 * Current: uint32_t matches carry_t for consistent intermediate calculations.
 */
typedef uint32_t temp_arith_t;

/* Signed arithmetic for add/subtract operations with borrow/carry.
 * Range: -255 to 510 (byte - byte, or byte + byte + carry).
 * Must handle negative results during subtraction operations.
 * Current: int32_t provides safe range for all bignum arithmetic.
 */
typedef int32_t signed_arith_t;

/* - Algorithm control */

/* Operation toggle for alternating series (1=add, 0=subtract).
 * Range: 0 or 1 only. Controls alternating signs in Bellard's series.
 * Current: uint8_t minimal storage for boolean-like value.
 */
typedef uint8_t operation_toggle_t;

/* Sign control parameter for bignum_div_addsub function.
 * Range: 0 (add) or non-zero (subtract). Treated as boolean in conditionals.
 * Current: int matches function parameter conventions.
 */
typedef int sign_control_t;

/* - Bignum initialization */

/* Small integer coefficients used to initialize bignums.
 * Range: Small positive integers (typically 4 for Bellard's coefficient).
 * Used in bignum_set() to establish initial values.
 * Current: uint32_t allows 32-bit write optimization at bignum high end.
 */
typedef uint32_t bignum_init_value_t;

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
void print_three_digit(digit_value_t val) {
    putchar('0' + (val / 100) % 10);
    putchar('0' + (val / 10) % 10);
    putchar('0' + val % 10);
    putchar(' ');
}

/* Print unsigned integer without leading zeros. */
void print_uint(digit_count_t val) {
    if (val == 0) {
        putchar('0');
        return;
    }
    char
        buf[12]; // Enough for uint32_t max (4,294,967,295 = 10 digits) + safety
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
 * At 50,000 digits: largest value ~= (166,670+9)*256 = 42,671,744
 */
divisor_t divisor;

/* Lower index bound for active bignum precision window.
 * Represents the leftmost (least significant) byte still containing
 * meaningful fractional data. As iterations progress, the fractional
 * precision requirements grow, and this bound increments to skip
 * computation on bytes that have become insignificant (leading zeros).
 * Formula: precision_lower = base / 128, where base grows by 159 per iteration.
 * This optimization significantly reduces computation time in later iterations.
 */
precision_bound_t precision_lower;

/* Upper index bound for active bignum precision window.
 * Represents the rightmost (most significant) byte containing non-zero
 * numerator data. The numerator is rescaled by 250/256 each iteration,
 * causing it to shrink and making higher-order bytes become zero.
 * This bound decrements accordingly, reducing unnecessary computation.
 * The precision window [precision_lower, precision_upper] defines the
 * active computation range, optimizing performance as the calculation
 * progresses.
 */
precision_bound_t precision_upper;

/* Size of each bignum array in bytes, calculated from requested digits.
 * Formula: (digits * 415241) / 1000000 + guard_digits
 * The constant 0.415241 = log_2(10)/8, representing bytes needed per decimal
 * digit. For 50,000 digits: ~= 20,762 + 3 = 20,765 bytes per bignum. This
 * provides sufficient binary precision to accurately represent the requested
 * decimal precision without accumulation errors.
 */
array_size_t bignum_size;

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
 * The 250/256 factor = 1000/1024. This rescaling maintains the mathematical
 * precision balance between sum_accumulator and numerator across iterations.
 */
bignum numerator;

/*
 * Bignum arithmetic functions
 * These functions implement arbitrary-precision arithmetic operations
 * on byte arrays representing fixed-point numbers in little-endian format.
 */

/* Initialize a bignum to a small integer value, clearing the fractional part */
void bignum_set(bignum bn, bignum_init_value_t value);

/* Perform division and add or subtract the quotient into the sum bignum.
 * This implements the core mathematical operation for each term in Bellard's
 * formula */
void bignum_div_addsub(sign_control_t is_subtract);

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
array_size_t calculate_bignum_size(digit_count_t digits, guard_count_t guard);

/* Allocate and initialize both sum and numerator bignums */
void allocate_bignums(array_size_t size);

/* Free allocated bignum memory */
void deallocate_bignums(void);

/*
 * Main algorithm
 */

/* Execute Bellard's formula to calculate pi to the specified precision */
void calculate_pi(digit_count_t digits, guard_count_t guard_digits);

/*
 * Initialize a bignum to a small integer value by clearing the fractional
 * part and placing the integer at the highest position. The integer value
 * is written as a 32-bit word at the high end, which may span multiple bytes.
 * This is used to set initial coefficients like 4 for the Bellard formula.
 */
void bignum_set(bignum bn, bignum_init_value_t value) {
    /* Clear fractional portion from precision_lower to bignum_size-1 */
    for (array_index_t i = precision_lower; i < bignum_size; i++) {
        bn[i] = 0;
    }
    /* Place integer value at highest position (32-bit write) */
    *(bignum_init_value_t *)&bn[bignum_size - 1] = value;
}

/*
 * Core division and accumulation function implementing each term of Bellard's
 * formula. Performs: sum_accumulator +/- = (numerator / divisor)
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
void bignum_div_addsub(sign_control_t is_subtract) {
    remainder_t remainder = 0;            /* Running remainder for division */
    divisor_t divisor_original = divisor; /* Preserve original divisor */

    /* Process bignum from most significant to least significant byte */
    for (array_index_t i = precision_upper; i >= precision_lower; i--) {
#ifdef BRANCH_PREDICTING_CPU
        /* Skip byte entirely if no contribution possible.
         * This optimization helps on modern CPUs with branch prediction
         * but may hurt performance on simple processors like 6502.
         */
        if (remainder == 0 && numerator[i] == 0) {
            continue; /* quotient_byte would be 0, nothing to add/subtract */
        }
#endif

        /* Build up remainder: shift left 8 bits and add next byte */
        remainder = remainder * 256 + numerator[i];
        quotient_byte_t quotient_byte = 0; /* Quotient byte being constructed */
        divisor *= 256; /* Scale divisor to match remainder */

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
        array_index_t sum_idx = i;
        do {
            signed_arith_t tmp =
                is_subtract ? (sum_accumulator[sum_idx] - quotient_byte)
                            : (sum_accumulator[sum_idx] + quotient_byte);
            sum_accumulator[sum_idx++] = tmp & 0xff; /* Store low byte */
            /* Carry/borrow continues if result outside byte range */
            quotient_byte = (tmp >= 0 && tmp <= 0xff) ? 0 : 1;
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
 * The 250/256 ratio compensates for the approximation 10^3 ~= 2^10 used in
 * adapting Bellard's binary series for decimal digit extraction.
 */
void bignum_rescale(bignum bn) {
    bignum_multiply_250(bn); /* Multiply by 250 */
    /* Divide by 256: shift all bytes down one position */
    for (array_index_t i = precision_lower; i < bignum_size; i++) {
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
    array_index_t i = bignum_size - 1;
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
    temp_arith_t temp;
    carry_t carry = 0;

    /* Process from least to most significant byte */
    for (array_index_t i = precision_lower; i <= bignum_size; i++) {
        temp = (temp_arith_t)bn[i] * 250 + carry;
        bn[i] = temp & 0xff; /* Store low byte */
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
    temp_arith_t temp;
    carry_t carry = 0;

    /* Process from least to most significant byte */
    for (array_index_t i = precision_lower; i <= bignum_size; i++) {
        temp = (temp_arith_t)bn[i] * 1000 + carry;
        bn[i] = temp & 0xff; /* Store low byte */
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
void calculate_pi(digit_count_t digits, guard_count_t guard_digits) {
    /* Initialize memory and precision tracking */
    bignum_size = calculate_bignum_size(digits, guard_digits);
    allocate_bignums(bignum_size);

    /* Set up bignum array bounds: track precision window that shrinks
     * as calculation progresses. Base tracks fractional precision growth
     * in fixed-point arithmetic (scaled by 128 for integer math).
     */
    precision_base_t base = 0;         /* Precision tracker (*128) */
    precision_lower = 0;               /* Lower bound starts at 0 */
    precision_upper = bignum_size - 1; /* Upper bound at array end */

    /* Initialize numerator to 4 (coefficient from Bellard's formula) */
    bignum_set(numerator, 4);

    /* Initialize loop counters based on iteration number n:
     * three_n = 3n (0, 3, 6, 9, 12, ...)
     * four_n = 4n (0, 4, 8, 12, 16, ...)
     * ten_n = 10n (0, 10, 20, 30, 40, ...)
     */
    digit_progression_t three_n =
        0;                /* 3n counter - tracks total digits produced */
    iter_4n_t four_n = 0; /* 4n counter - used in denominators 4n+1, 4n+3 */
    iter_10n_t ten_n = 0; /* 10n counter - used in denominators 10n+1, 10n+3,
                            10n+5, 10n+7, 10n+9 */
    operation_toggle_t op = 1; /* Operation toggle: 1=add, 0=subtract */

    /* Initialize sum accumulator to 4 */
    bignum_set(sum_accumulator, 4);

    digit_count_t digit_counter = 0; /* Track digits output for reporting */

    /* Main calculation loop: each iteration produces 3 digits */
    do {
        /* Compute the seven terms of Bellard's decimal spigot formula:
         * Each term has form: coefficient / (polynomial_in_n * scale_factor)
         * The alternating signs and specific denominators derive from Bellard's
         * adaptation of his binary BBP formula for decimal digit extraction.
         */

        divisor = ten_n + 1;
        if (three_n >
            0) { /* Skip first iteration to avoid double-counting initial 4 */
            bignum_div_addsub(1 - op);
        }

        divisor = (ten_n + 3) * 4;
        bignum_div_addsub(op);

        divisor = (ten_n + 5) * 64;
        bignum_div_addsub(op);

        divisor = (ten_n + 7) * 64;
        bignum_div_addsub(op);

        divisor = (ten_n + 9) * 256;
        bignum_div_addsub(1 - op);

        divisor = (four_n + 1) * 8;
        bignum_div_addsub(op);

        divisor = (four_n + 3) * 256;
        bignum_div_addsub(op);

        /* Extract three digits from sum's integer part */
        digit_value_t digit_val =
            *(digit_value_t *)&sum_accumulator[bignum_size - 1];

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
         */
        base += 159;
        precision_lower = base / 128; /* Convert to byte index */
        four_n += 4;                  /* Advance 4n counter: 0->4->8->12... */
        ten_n += 10;  /* Advance 10n counter: 0->10->20->30... */
        three_n += 3; /* Advance 3n counter: 0->3->6->9... */

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
array_size_t calculate_bignum_size(digit_count_t digits, guard_count_t guard) {
    return (array_size_t)((digits * 415241LL) / 1000000 + guard);
}

/*
 * Allocate memory for both bignum arrays with error checking.
 * Includes padding bytes to ensure safe 32-bit read operations when
 * extracting output digits. Initializes all bytes to zero and exits
 * with error message if allocation fails.
 */
void allocate_bignums(array_size_t size) {
    const array_size_t pad = 4; /* Padding for safe 32-bit reads */

    /* Allocate sum accumulator with error checking */
    sum_accumulator = (bignum)malloc(size + pad);
    if (!sum_accumulator) {
        print_str("Failed to alloc sum_accumulator: ");
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
    for (array_index_t i = 0; i < size + pad; i++) {
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
    digit_count_t digits = 1000;
    guard_count_t guard = 3;

    /* Initialize platform-specific heap expansion */
    init_platform();

    print_str("Calculating pi to ");
    print_uint(digits);
    print_str("+ digits...\n");
    calculate_pi(digits, guard);

    return 0;
}
