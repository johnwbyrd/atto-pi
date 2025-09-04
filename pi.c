#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

void bellard(uint32_t digits, uint32_t guard);
void proc_set(uint8_t *bignumP, uint32_t C);
void proc_divaddsub(int C);
void proc_rescale(uint8_t *bignumP);
void proc_mask(uint8_t *bignumP);
void proc_x250(uint8_t *bignumP);
void proc_x1000(uint8_t *bignumP);

// pi calculation using spigotised Bellard formula
// inspired by BBP calculation in youtu.be/SPTzzSuBFlc

// pi to 80 digits is
// 31415926535897932384
// 62643383279502884197
// 16939937510582097494
// 45923078164062862089

// Global Variables are capitalized
// TODO: Re-write to avoid these!

uint64_t D;
int32_t      L;
int32_t      M;
int32_t      Big;
uint8_t  *SumP;
uint8_t  *NumeratorP;

// Debug counter to limit output
int debug_counter = 0;
int digit_counter = 0;

int main(int argc, char *argv[]) {
   uint32_t digits = 0;
   uint32_t guard = 2;
   if (argc >= 2) {
      digits = atoi(argv[1]);
   }
   if (argc >= 3) {
      guard = atoi(argv[2]);
   }
   digits = 150;
   guard = 3; 

   if (digits > 0 && guard >= 0) {
      bellard(digits, guard);
   } else {
      printf("usage: %s ndigits [ guard ] \n", argv[0]);
   }
   return 0;
}


void dump_sum(const char *label) {
   printf("%s", label);
   for (int i = 0; i < Big + 4; i++) {
      printf(" %02x", *(SumP + i));
   }
   printf("\n");
}


void bellard(uint32_t digits, uint32_t guard) {
   // TODO: this isn't the same as the machine code version
   Big = ((digits * 5) / 12) + guard;  // can add 1 or 2 for guard digits

   // pad bignums to ensure 4-byte comparison at the end to see zeros
   int pad = 4;

   SumP = (uint8_t *) malloc(Big + pad); // our bignum series accumulator
   NumeratorP = (uint8_t *) malloc(Big + pad); // bignum numerator

   for (int i = 0; i < Big + pad; i++)
   {
        *(SumP + i) = 0;
        *(NumeratorP + i) = 0;
   }


   float base = 0.0;
   L = (int32_t) base; // our bignums can get shorter as we go
   M = Big - 1;    // leading zeros index for fast forward division

   // our numerator will be rescaled by 250/256 each digit
   proc_set(NumeratorP, 4);

   // the denominator is a (two or three byte) integer
   // and K is the iteration variable
   int32_t k  = 0;
   int32_t f  = 0;      // increments by 4
   uint64_t t  = 0;      // increments by 10
   int32_t op = 1;      // toggle +/- each iteration

   // start producing digits!
   proc_set(SumP, 4);

   do {
      // compute N/K and add/sub the quotient into S
      // we adjust the denominator as needed
      D = t + 1;
      if (k > 0) {
         proc_divaddsub(1 - op);
      }
      D = (t + 3) * 4;
      proc_divaddsub(op);
      D = (t + 5) * 64;
      proc_divaddsub(op);
      D = (t + 7) * 64;
      proc_divaddsub(op);    
      D = (t + 9) * 256;
      proc_divaddsub(1-op);
      D = (f + 1) * 8;
      proc_divaddsub(op);
      D = (f + 3) * 256;
      proc_divaddsub(op);

      // S is updated, now take three digits
      uint32_t digit_val = *(uint32_t *)(SumP + Big - 1);
      
      if (k > 0) {
         printf("%03ld ", digit_val);
         digit_counter++;
      } else {
         printf("%ld.\n", digit_val);
         digit_counter++;
      }

      proc_mask(SumP);  // remove those digits from S
      proc_x1000(SumP);

      // now re-scale the numerator
      proc_rescale(NumeratorP);
      op = 1 - op;
      if (*(NumeratorP + M) == 0) {
         M = M - 1;
      }
      base += 3.0*106.0/256.0;
      L = (int) base;
      f = f + 4;
      t = t + 10;
      k = k + 3;
   } while (k <= digits);

   free(SumP);
   free(NumeratorP);
}

// --------------------------------------------------

void proc_set(uint8_t *bignumP, uint32_t C) {
   // C is a tiny integer - one byte
   // integer end of our bignum is at the far end
   for (int i = L; i < Big; i++) {
      *(bignumP + i) = 0;
   }
   *(uint32_t *)(bignumP + Big - 1) = C;
}

// --------------------------------------------------

// Globals
//   NumeratorP, SumP are bignums
//   SumP is a bignum to add or subtract into
//   D% is an integer, the number to divide by
//   M is the MSB index into NumeratorP
//   L is the LSB index into NumeratorP
// Parameters
//   C% is TRUE for add, FALSE for subtract
// Locals
//   I% byte loop counter (M to L)
//   J% bit loop counter (0 TO 7)
//   T% is a 4-byte temporary numerator
//   B% is the byte being accumulated

// --------------------------------------------------

void proc_divaddsub(int sub) {
   uint64_t T = 0;
   uint64_t D_original = D;
   
   // Debug first few calls only
   if (debug_counter < 4) {
      // printf("divaddsub: D=%llu sub=%d\n", D, sub);
      debug_counter++;
   }
   
   for (int32_t i = M; i >= L; i--) {
      T = T * 256 + *(NumeratorP + i);
      int32_t B = 0;
      D *= 256;
      for (uint8_t j = 0; j <= 7; j++) {
         B *= 2;
         D /= 2;
         if (T >= D) {
            T -= D;
            B += 1;
         }
      }
      
      // This bit differs from the BASIC
      // if (sub) {
      //    *(uint32_t *)(SumP + i) -= B;
      // } else {
      //    *(uint32_t *)(SumP + i) += B;
      // }
      uint8_t *tmps = SumP + i;
      do {
         int32_t tmp = sub ? (*tmps - B) : (*tmps + B);
         *tmps++ = tmp & 255;
         B = (tmp >= 0 && tmp <= 255) ? 0 : 1;
      } while (B);
   }
   
   // Restore D
   D = D_original;
}

// --------------------------------------------------

void proc_rescale(uint8_t *bignumP) {
   // we multiply by 250/256 which is less than one
   proc_x250(bignumP);
   for (int32_t i = L; i < Big; i++) {
      *(bignumP + i) = *(bignumP + i + 1);
   }
   *(bignumP + Big) = 0;
}

// --------------------------------------------------

void proc_mask(uint8_t *bignumP) {
   // set integer part of the bignum to zero
   // assumes bignum has 1 byte before the binary point
   int32_t i = Big - 1;
   *(bignumP + i) = 0;
   *(bignumP + i + 1) = 0; // mask extra byte beyond MSB
}

// --------------------------------------------------

void proc_x250(uint8_t *bignumP) {
   // bignum multiply by small number
   // two byte result from each digit demands a carry
   uint32_t temp;
   uint32_t carry = 0;
   for (int32_t i = L; i <= Big; i++) { // extra element allows one byte of carry to propagate
      temp = (uint32_t)*(bignumP + i) * 250 + carry;
      *(bignumP + i) = temp & 255;
      carry = temp >> 8;
   }
   if (carry > 0) {
      assert("proc_x250 should end with carry = 0");
   }
}

// --------------------------------------------------

void proc_x1000(uint8_t *bignumP) {
   // bignum multiply by small number
   // two byte result from each digit demands a carry
   uint32_t temp;
   uint32_t carry = 0;
   for (int32_t i = L; i <= Big; i++) { // extra element allows one byte of carry to propagate
      temp = (uint32_t)*(bignumP + i) * 1000 + carry;
      *(bignumP + i) = temp & 255;
      carry = temp >> 8;
   }
   if (carry > 0) {
      assert("proc_x1000 should end with carry = 0");
   }
}
