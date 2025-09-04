/*
 * pi.c
 *
 * Author: John Byrd <johnwbyrd at gmail dot com>
 *
 * This program implements Fabrice Bellard's 2009 formula for computing pi, 
 * which is a spigot algorithm that generates decimal digits sequentially 
 * without storing the entire result. Unlike the more famous Bailey-Borwein-
 * Plouffe formula, Bellard's approach works in base-1000, producing three 
 * decimal digits per iteration rather than hexadecimal digits.
 *
 * The mathematical foundation is the infinite series:
 * pi = (1/2^6) * sum(n=0 to infinity) [(-1)^n/2^(10n)] * [
 *     -2^5/(4n+1) - 1/(4n+3) + 2^8/(10n+1) - 2^6/(10n+3) 
 *     - 2^2/(10n+5) - 2^2/(10n+7) + 1/(10n+9)
 * ]
 *
 * Each iteration evaluates seven rational terms, alternating between addition 
 * and subtraction. The denominators grow as functions of the iteration counter, 
 * requiring arbitrary-precision division to maintain accuracy across thousands 
 * of digits.
 *
 * This particular implementation was inspired by David Banks (hoglet)'s recent
 * work on implementing this spigot on the BBC Micro in BASIC. It was written
 * for the LLVM-MOS C/C++ compiler, but it should work on other reasonably
 * standards compliant compilers as well. 
 *
 * - Arithmetic implementation
 *
 * The core challenge lies in performing high-precision arithmetic on a system 
 * with only 8-bit native operations and 16-bit addressing. The solution uses 
 * "bignums" -- arbitrary-precision numbers represented as arrays of bytes in 
 * little-endian format. Each bignum stores a fixed-point number where the 
 * integer part occupies the highest-indexed byte and the fractional part 
 * extends downward through lower indices.
 *
 * The size of each bignum is calculated as (digits * 5) / 12 + guard, where 
 * the guard digits provide computational headroom to prevent rounding errors 
 * from accumulating. This formula derives from the convergence properties of 
 * Bellard's series and ensures sufficient precision for the requested output.
 *
 * Division operations use binary long division, processing one bit at a time 
 * through an 8-bit window. This approach naturally fits the 6502's byte-
 * oriented architecture while maintaining the precision needed for thousands 
 * of digits. The divisor is repeatedly halved while the remainder is compared 
 * against it, building up the quotient bit by bit.
 *
 * Multiplication by constants (250 and 1000) employs carry propagation across 
 * the entire bignum. Each byte is multiplied by the constant, added to any 
 * carry from the previous position, and the result is split into a stored 
 * byte and a carry for the next position. This process continues until the 
 * carry becomes zero, ensuring no precision is lost.
 *
 * - Algorithmic flow
 *
 * Each iteration begins by computing the seven terms of Bellard's formula 
 * through a sequence of divisions and additions/subtractions into the 
 * accumulator bignum. The divisors follow specific patterns based on the 
 * iteration variables: t increases by 10 each iteration and appears in 
 * expressions like (t+1), (t+3)*4, and so forth, while f increases by 4 
 * and appears in (f+1)*8 and (f+3)*256.
 *
 * After computing all terms, the algorithm extracts three decimal digits 
 * from the integer part of the accumulator. These digits represent the 
 * current contribution to pi's decimal expansion. The extracted digits are 
 * then masked (zeroed) from the accumulator, and the remaining fractional 
 * part is multiplied by 1000 to shift the next three digits into the 
 * integer position for the following iteration.
 *
 * The numerator bignum undergoes rescaling through multiplication by 250/256, 
 * effectively dividing by 256/250 = 1.024. This scaling factor maintains 
 * the precision balance between the accumulator and numerator as the 
 * algorithm progresses. The division by 256 is implemented as a byte shift, 
 * taking advantage of the binary representation.
 *
 * - Precision and limits
 * 
 * This algorithm runs relatively soonish, if 100 digits are requested. It
 * can do 1000, if you're prepared to wait a while.  And it can generate 10000
 * digits of pi on a Commodore 64, but don't expect that result anytime soon,
 * because it seems to require several minutes to generate even 3 digits from 
 * that spigot.
 *
 * The implementation's maximum capacity stems from the integer types chosen 
 * for loop counters. The main iteration counter k uses uint16_t and increments 
 * by 3 each cycle. This creates a hard limit around 65,532 digits, beyond 
 * which k would overflow and wrap to a small positive number, causing an 
 * infinite loop when compared against the target digit count.
 *
 * The secondary counters f and t use larger types to handle their respective 
 * computational demands. Variable f employs uint32_t because it appears in 
 * expressions like (f+3)*256, which would overflow uint16_t after just a 
 * few dozen iterations. At the maximum digit count, f reaches approximately 
 * 87,376, and its largest computation (87,379)*256 produces 22,369,024, 
 * which fits comfortably within uint32_t's range.
 *
 * Variable t uses uint64_t due to its appearance in multiple denominator 
 * calculations. While t itself grows relatively slowly (reaching about 
 * 218,440 at maximum digits), the safety margin provided by 64-bit arithmetic 
 * prevents any possibility of overflow in the various expressions involving t.
 *
 * Memory requirements scale linearly with the requested digit count. For 1000 
 * digits, each bignum requires approximately 417 bytes, totaling 834 bytes 
 * for both the accumulator and numerator. At the theoretical maximum of 
 * 65,532 digits, this grows to about 54,610 bytes total, which approaches 
 * but remains within the 6502's addressing capacity.
 *
 * - Architectural considerations
 *
 * Several design choices specifically accommodate the 6502's limitations and 
 * characteristics. The use of signed int16_t for array indices L and M 
 * prevents subtle bugs that arise from comparing signed loop counters with 
 * unsigned array bounds on 16-bit systems. The little-endian byte ordering 
 * matches the 6502's native format, simplifying multi-byte operations.
 *
 * Memory allocation includes padding beyond the computed bignum size to 
 * ensure safe 32-bit read operations when extracting output digits. This 
 * prevents potential crashes from reading beyond allocated memory boundaries, 
 * a critical consideration on systems without memory protection.
 *
 * The spigot nature of Bellard's algorithm proves particularly well-suited 
 * to resource-constrained systems like the 6502. Rather than computing and 
 * storing all digits before output, each digit can be displayed as soon as 
 * it's calculated, allowing the computation of arbitrarily long pi expansions 
 * without proportional memory growth for the output itself.
 *
 * - Scaling operations in the algorithm
 *
 * Two specific multiplication operations are central to the algorithm's 
 * operation and deserve explanation. The multiply-by-1000 operation extracts 
 * output digits by shifting the fractional part of the accumulator into the 
 * integer position. Since we extract three decimal digits per iteration, 
 * multiplying by 1000 moves the next group of three digits from the fractional 
 * part into the integer part where they can be read and output.
 *
 * The multiply-by-250 operation implements the rescaling factor 250/256 that 
 * maintains precision balance between iterations. This ratio equals 
 * approximately 0.9765625, which slightly reduces the numerator each cycle 
 * to prevent unlimited growth. The division by 256 is implemented as a 
 * byte shift (moving all bytes down one position), making 250/256 an 
 * efficient way to implement the mathematically required scaling factor.
 * This rescaling prevents the numerator from growing without bound while 
 * maintaining sufficient precision for accurate digit extraction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#define BELLARD_TERM_COUNT 7

/* This value needs be overridden from the default malloc() implementation
 * in most cases; this value seems to work up to 10000 digits
 */
size_t __heap_default_limit = 2 << 14;

uint64_t D;
int16_t L;
int16_t M;
uint16_t Big;
uint8_t *SumP;        
uint8_t *NumeratorP;  

/*
 * - Bignum arithmetic functions
 * These functions implement arbitrary-precision arithmetic operations
 * on byte arrays representing fixed-point numbers in little-endian format.
 */

/* Initialize a bignum to a small integer value, clearing the fractional part */
void bignum_set(uint8_t *bignum, uint32_t value);

/* Perform division and add or subtract the quotient into the sum bignum.
 * This implements the core mathematical operation for each term in Bellard's formula */
void bignum_div_addsub(int is_subtract);

/* Rescale a bignum by the factor 250/256, used to prevent numerator overflow
 * while maintaining precision across iterations */
void bignum_rescale(uint8_t *bignum);

/* Zero the integer part of a bignum, removing digits that have been output */
void bignum_mask_digits(uint8_t *bignum);

/* Multiply a bignum by 250 with carry propagation.
 * Part of the rescaling operation (250/256 ratio) */
void bignum_multiply_250(uint8_t *bignum);

/* Multiply a bignum by 1000 with carry propagation.
 * Shifts the next three decimal digits into the integer position */
void bignum_multiply_1000(uint8_t *bignum);

/*
 * - Memory management functions
 * Handle allocation and sizing of bignum arrays based on precision requirements
 */

/* Calculate the required size for bignums based on desired digit precision */
uint16_t calculate_bignum_size(uint16_t digits, uint8_t guard);

/* Allocate and initialize both sum and numerator bignums */
void allocate_bignums(uint16_t size);

/* Free allocated bignum memory */
void deallocate_bignums(void);

/*
 * - Main algorithm function
 */

/* Execute Bellard's formula to calculate pi to the specified precision */
void calculate_pi_bellard(uint16_t digits, uint8_t guard_digits);

/*
 * Initialize a bignum to a small integer value.
 * Clears the fractional part (low indices) and places the integer value
 * at the high end of the array. This is used to set initial coefficients
 * for the Bellard formula terms.
 */
void bignum_set(uint8_t *bignum, uint32_t value) {
    /* Clear fractional portion from L to Big-1 */
    for (uint16_t i = L; i < Big; i++) {
        *(bignum + i) = 0;
    }
    /* Place integer value at highest position (32-bit write) */
    *(uint32_t *)(bignum + Big - 1) = value;
}

/*
 * Core division and accumulation function for Bellard's formula.
 * Divides the numerator bignum by the global divisor D, then adds or subtracts
 * the quotient into the sum bignum. This implements each rational term in the
 * infinite series. Uses binary long division processing 8 bits at a time.
 */
void bignum_div_addsub(int is_subtract) {
    uint32_t T = 0;              /* Running remainder for division */
    uint64_t D_original = D;     /* Preserve original divisor */
    
    /* Process bignum from most significant to least significant byte */
    for (int16_t i = M; i >= L; i--) {
        /* Build up remainder: shift left 8 bits and add next byte */
        T = T * 256 + *(NumeratorP + i);
        uint16_t B = 0;          /* Quotient byte being constructed */
        D *= 256;                /* Scale divisor to match remainder */
        
        /* Binary division: extract 8 bits one at a time */
        for (uint8_t j = 0; j <= 7; j++) {
            B *= 2;              /* Shift quotient bit left */
            D /= 2;              /* Halve divisor for this bit */
            if (T >= D) {
                T -= D;          /* Remainder -= divisor */
                B += 1;          /* Set quotient bit */
            }
        }
        
        /* Add or subtract quotient byte into sum with carry propagation */
        uint8_t *tmps = SumP + i;
        do {
            int32_t tmp = is_subtract ? (*tmps - B) : (*tmps + B);
            *tmps++ = tmp & 255; /* Store low byte */
            B = (tmp >= 0 && tmp <= 255) ? 0 : 1;  /* Carry/borrow */
        } while (B);
    }
    
    /* Restore original divisor for next term */
    D = D_original;
}

/*
 * Rescale bignum by factor 250/256 to prevent numerator growth.
 * First multiplies by 250, then divides by 256 (via byte shift).
 * This maintains precision balance between sum and numerator across iterations.
 */
void bignum_rescale(uint8_t *bignum) {
    bignum_multiply_250(bignum);         /* Multiply by 250 */
    /* Divide by 256: shift all bytes down one position */
    for (uint16_t i = L; i < Big; i++) {
        *(bignum + i) = *(bignum + i + 1);
    }
    *(bignum + Big) = 0;                 /* Clear vacated high byte */
}

/*
 * Zero the integer part of bignum after extracting output digits.
 * This removes the digits that have been output, leaving only the
 * fractional part for the next iteration. Clears both the integer
 * byte and one extra byte for safety.
 */
void bignum_mask_digits(uint8_t *bignum) {
    uint16_t i = Big - 1;
    *(bignum + i) = 0;                   /* Clear integer part */
    *(bignum + i + 1) = 0;               /* Clear extra byte for safety */
}

/*
 * Multiply bignum by 250 with carry propagation.
 * This is part of the 250/256 rescaling factor. Each byte is multiplied
 * by 250, with carries propagated to higher bytes. The final carry must
 * be zero to ensure no precision loss (guaranteed by proper sizing).
 */
void bignum_multiply_250(uint8_t *bignum) {
    uint32_t temp;
    uint32_t carry = 0;
    
    /* Process from least to most significant byte */
    for (uint16_t i = L; i <= Big; i++) {
        temp = (uint32_t)(*(bignum + i)) * 250 + carry;
        *(bignum + i) = temp & 255;      /* Store low byte */
        carry = temp >> 8;               /* Propagate high byte as carry */
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
void bignum_multiply_1000(uint8_t *bignum) {
    uint32_t temp;
    uint32_t carry = 0;
    
    /* Process from least to most significant byte */
    for (uint16_t i = L; i <= Big; i++) {
        temp = (uint32_t)(*(bignum + i)) * 1000 + carry;
        *(bignum + i) = temp & 255;      /* Store low byte */
        carry = temp >> 8;               /* Propagate high byte as carry */
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
    Big = calculate_bignum_size(digits, guard_digits);
    allocate_bignums(Big);
    
    /* Set up bignum array bounds: L tracks leading zeros, M tracks precision */
    float base = 0.0;                    /* Tracks fractional precision needs */
    L = (int16_t) base;                  /* Lower bound starts at 0 */
    M = Big - 1;                         /* Upper bound at end of array */
    
    /* Initialize numerator to 4 (coefficient from Bellard's formula) */
    bignum_set(NumeratorP, 4);
    
    /* Initialize loop counters: k=iteration*3, f=iteration*4, t=iteration*10 */
    uint16_t k = 0;                      /* Main iteration counter */
    uint32_t f = 0;                      /* Secondary counter for f-terms */
    uint64_t t = 0;                      /* Tertiary counter for t-terms */
    uint8_t op = 1;                      /* Operation toggle: 1=add, 0=subtract */
    
    /* Initialize sum accumulator to 4 */
    bignum_set(SumP, 4);
    
    uint16_t digit_counter = 0;          /* Track digits output for reporting */
    uint8_t group_counter = 0;           /* Track groups of 3 digits for line breaks */
    
    /* Main calculation loop: each iteration produces 3 digits */
    do {
        /* Compute the seven terms of Bellard's formula */
        
        /* Term 1: -2^5/(4n+1) = -32/(t+1) */
        D = t + 1;
        if (k > 0) {                     /* Skip first iteration adjustment */
            bignum_div_addsub(1 - op);
        }
        
        /* Term 2: -1/(4n+3) = -1/((t+3)*4) */
        D = (t + 3) * 4;
        bignum_div_addsub(op);
        
        /* Term 3: +2^8/(10n+1) = +256/((t+5)*64) */
        D = (t + 5) * 64;
        bignum_div_addsub(op);
        
        /* Term 4: -2^6/(10n+3) = -64/((t+7)*64) */
        D = (t + 7) * 64;
        bignum_div_addsub(op);
        
        /* Term 5: -2^2/(10n+5) = -4/((t+9)*256) */
        D = (t + 9) * 256;
        bignum_div_addsub(1 - op);
        
        /* Term 6: -2^2/(10n+7) = -4/((f+1)*8) */
        D = (f + 1) * 8;
        bignum_div_addsub(op);
        
        /* Term 7: +1/(10n+9) = +1/((f+3)*256) */
        D = (f + 3) * 256;
        bignum_div_addsub(op);
        
        /* Extract three digits from sum's integer part */
        uint32_t digit_val = *(uint32_t *)(SumP + Big - 1);
         
        /* Output digits with appropriate formatting */
        if (k > 0) {
            printf("%03lu ", (unsigned long)digit_val);
            digit_counter += 3;          /* Each iteration produces 3 digits */
            group_counter++;             /* Count groups of 3 digits */
            if (group_counter == 10) {
#ifdef GROUP_OUTPUTS
                printf("\n");            /* Line break every 10 groups (30 digits) */
#endif // GROUP_OUTPUTS
                group_counter = 0;       /* Reset counter */
            }
        } else {
            printf("%lu.\n", (unsigned long)digit_val);  /* First digit with decimal */
            digit_counter += 1;          /* First iteration produces 1 digit */
            /* Don't count first iteration toward group counter */
        }
        
        /* Prepare for next iteration */
        bignum_mask_digits(SumP);        /* Remove extracted digits */
        bignum_multiply_1000(SumP);      /* Shift next digits into position */
        bignum_rescale(NumeratorP);      /* Apply 250/256 scaling */
        
        /* Toggle operation sign for alternating series */
        op = 1 - op;
        
        /* Adjust precision tracking as numerator shrinks */
        if (*(NumeratorP + M) == 0) {
            M = M - 1;
        }
        
        /* Update loop variables for next iteration */
        base += 3.0 * 106.0 / 256.0;     /* Track fractional precision needs */
        L = (int16_t) base;               /* Update lower bound */
        f = f + 4;                        /* Increment f-counter */
        t = t + 10;                       /* Increment t-counter */
        k = k + 3;                        /* Increment main counter (3 digits) */
        
    } while (digit_counter < digits);
    
    /* Clean up and report results */
    deallocate_bignums();
    printf("\nCalculated %u digits of pi.\n", digit_counter);
}

/*
 * Calculate required bignum size based on digit precision and guard digits.
 * The formula (digits * 5) / 12 + guard derives from the convergence
 * properties of Bellard's series and provides adequate precision for
 * the requested number of output digits.
 */
uint16_t calculate_bignum_size(uint16_t digits, uint8_t guard) {
    return (uint16_t)((digits * 5) / 12 + guard);
}

/*
 * Allocate memory for both bignum arrays with error checking.
 * Includes padding bytes to ensure safe 32-bit read operations when
 * extracting output digits. Initializes all bytes to zero and exits
 * with error message if allocation fails.
 */
void allocate_bignums(uint16_t size) {
    const uint16_t pad = 4;              /* Padding for safe 32-bit reads */
    
    /* Allocate sum accumulator with error checking */
    SumP = (uint8_t *) malloc(size + pad);
    if (!SumP) {
        printf("Error: Failed to allocate %u bytes for SumP\n", size + pad);
        exit(1);
    }
    
    /* Allocate numerator with error checking */
    NumeratorP = (uint8_t *) malloc(size + pad);
    if (!NumeratorP) {
        printf("Error: Failed to allocate %u bytes for NumeratorP\n", size + pad);
        free(SumP);                      /* Clean up partial allocation */
        exit(1);
    }
    
    /* Initialize both arrays to zero */
    for (uint16_t i = 0; i < size + pad; i++) {
        *(SumP + i) = 0;
        *(NumeratorP + i) = 0;
    }
}

/*
 * Free allocated bignum memory.
 * Called at end of calculation to prevent memory leaks.
 */
void deallocate_bignums(void) {
    free(SumP);
    free(NumeratorP);
}

int main(void) {
    uint16_t digits = 100; 
    uint8_t guard = 3;
    
    printf("Calculating pi to %u+ digits...\n", digits);
    calculate_pi_bellard(digits, guard);
    
    return 0;
}
