/*
 * pi.c
 *
 * Author: John Byrd <johnwbyrd at gmail dot com>
 * 
 * See README.md for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#ifdef __llvm_mos__
extern char __heap_start;
#endif

/*
 * Custom print functions to replace printf and keep code size small.
 * These avoid the large formatting overhead of printf (~3,000+ bytes)
 * by using only putchar for output, reducing the binary size significantly.
 */
void print_str(const char *s) {
    while (*s) putchar(*s++);
}

void print_three_digit(uint32_t val) {
    putchar('0' + (val / 100) % 10);
    putchar('0' + (val / 10) % 10);
    putchar('0' + val % 10);
    putchar(' ');
}

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
    while (i--) putchar(buf[i]);
}

#ifdef __llvm_mos__
/* Maximum heap for C64 with BASIC unmapped:
 * Stack at $D000, heap starts after program (~$1CE2).
 * Available: ~44KB. Using 42KB to leave 2KB for stack.
 */
size_t __heap_default_limit = 42 * 1024;
#endif

uint64_t D;
int16_t L;
int16_t M;
uint16_t Big;
uint8_t *SumP;        
uint8_t *NumeratorP;  

/*
 * Bignum arithmetic functions
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
    uint64_t T = 0;              /* Running remainder for division */
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
    int32_t base = 0;                    /* Tracks fractional precision needs (fixed-point, *128) */
    L = 0;                               /* Lower bound starts at 0 */
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
            print_three_digit(digit_val);
            digit_counter += 3;          /* Each iteration produces 3 digits */
        } else {
            print_uint(digit_val); putchar('.'); putchar('\n');  /* First digit with decimal */
            digit_counter += 1;          /* First iteration produces 1 digit */
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
        base += 159;                      /* Track fractional precision needs (159/128 ≈ 1.2421875) */
        L = base / 128;                   /* Update lower bound */
        f = f + 4;                        /* Increment f-counter */
        t = t + 10;                       /* Increment t-counter */
        k = k + 3;                        /* Increment main counter (3 digits) */
        
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
    const uint16_t pad = 4;              /* Padding for safe 32-bit reads */
    
    /* Allocate sum accumulator with error checking */
    SumP = (uint8_t *) malloc(size + pad);
    if (!SumP) {
        print_str("Error: Failed to alloc SumP: ");
        print_uint(size + pad);
        exit(1);
    }
    
    /* Allocate numerator with error checking */
    NumeratorP = (uint8_t *) malloc(size + pad);
    if (!NumeratorP) {
        print_str("Failed to alloc NumeratorP: ");
        print_uint(size + pad);
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
    uint16_t digits = 25000;
    uint8_t guard = 3;

    print_str("Calculating pi to ");
    print_uint(digits);
    print_str("+ digits...\n");
    calculate_pi_bellard(digits, guard);

    return 0;
}
