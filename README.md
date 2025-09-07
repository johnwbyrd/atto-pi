# pi.c

**Author**: John Byrd <johnwbyrd at gmail dot com>

This program implements François Bellard's unpublished 1997 decimal spigot algorithm for computing pi, which is a spigot algorithm that generates decimal digits sequentially without storing the entire result. Unlike the more famous [Bailey-Borwein-Plouffe formula](https://observablehq.com/@rreusser/computing-with-the-bailey-borwein-plouffe-formula), this decimal adaptation of Bellard's binary series works in base-1000, producing three decimal digits per iteration rather than hexadecimal digits.

The mathematical foundation is the infinite series:

$$\pi = \sum_{n=0}^{\infty} \left( \frac{-32}{4n+1} - \frac{1}{4n+3} + \frac{256}{10n+1} - \frac{64}{10n+3} - \frac{4}{10n+5} - \frac{4}{4n+1} + \frac{1}{4n+3} \right) \cdot 2^{-10n-6}$$

Each iteration evaluates seven rational terms, alternating between addition and subtraction. The denominators grow as functions of the iteration counter, requiring arbitrary-precision division to maintain accuracy across thousands of digits.

This particular implementation was inspired by [David Banks (hoglet)'s recent work](https://github.com/BigEd/pi-spigot-for-micros) on implementing spigots on the BBC Micro. It was written for the LLVM-MOS C/C++ compiler, but it should work on other reasonably standards compliant compilers as well.

## Arithmetic implementation

The core challenge lies in performing high-precision arithmetic on a system with only 8-bit native operations and 16-bit addressing. The solution uses "bignums" -- arbitrary-precision numbers represented as arrays of bytes in little-endian format. Each bignum stores a fixed-point number where the integer part occupies the highest-indexed byte and the fractional part extends downward through lower indices.

The size of each bignum is calculated as $\text{digits} \times \frac{\log_2(10)}{8} + \text{guard}$, where the constant $\frac{\log_2(10)}{8} \approx 0.415241$ represents the bytes needed per decimal digit, and the guard digits provide computational headroom to prevent rounding errors from accumulating. This formula derives from the convergence properties of Bellard's series and ensures sufficient precision for the requested output.

Division operations use binary long division, processing one bit at a time through an 8-bit window. This approach naturally fits the 6502's byte-oriented architecture while maintaining the precision needed for thousands of digits. The divisor is repeatedly halved while the remainder is compared against it, building up the quotient bit by bit.

Multiplication by constants (250 and 1000) employs carry propagation across the entire bignum. Each byte is multiplied by the constant, added to any carry from the previous position, and the result is split into a stored byte and a carry for the next position. This process continues until the carry becomes zero, ensuring no precision is lost.

## Algorithmic flow

Each iteration begins by computing the seven terms of Bellard's formula through a sequence of divisions and additions/subtractions into the accumulator bignum. The divisors follow specific patterns based on the iteration variables: $t$ increases by 10 each iteration and appears in expressions like $(t+1)$, $(t+3) \times 4$, and so forth, while $f$ increases by 4 and appears in $(f+1) \times 8$ and $(f+3) \times 256$.

After computing all terms, the algorithm extracts three decimal digits from the integer part of the accumulator. These digits represent the current contribution to pi's decimal expansion. The extracted digits are then masked (zeroed) from the accumulator, and the remaining fractional part is multiplied by 1000 to shift the next three digits into the integer position for the following iteration.

The numerator bignum undergoes rescaling through multiplication by $\frac{250}{256}$, effectively dividing by $\frac{256}{250} = 1.024$. This scaling factor maintains the precision balance between the accumulator and numerator as the algorithm progresses. The division by 256 is implemented as a byte shift, taking advantage of the binary representation.

## Precision and limits

This algorithm runs relatively soonish, if 100 digits are requested. It can do 1000, if you're prepared to wait a while. And it can generate 10000 digits of pi on a Commodore 64, but don't expect that result anytime soon, because it seems to require several minutes to generate even 3 digits from that spigot.

The implementation's maximum capacity stems from the integer types chosen for loop counters. The main iteration counter $k$ uses `uint16_t` and increments by 3 each cycle. This creates a hard limit around 65,532 digits, beyond which $k$ would overflow and wrap to a small positive number, causing an infinite loop when compared against the target digit count.

The secondary counters $f$ and $t$ use larger types to handle their respective computational demands. Variable $f$ employs `uint32_t` because it appears in expressions like $(f+3) \times 256$, which would overflow `uint16_t` after just a few dozen iterations. At the maximum digit count, $f$ reaches approximately 87,376, and its largest computation $(87,379) \times 256$ produces 22,369,024, which fits comfortably within `uint32_t`'s range.

Variable $t$ uses `uint64_t` due to its appearance in multiple denominator calculations. While $t$ itself grows relatively slowly (reaching about 218,440 at maximum digits), the safety margin provided by 64-bit arithmetic prevents any possibility of overflow in the various expressions involving $t$.

Memory requirements scale linearly with the requested digit count. For 1000 digits, each bignum requires approximately 417 bytes, totaling 834 bytes for both the accumulator and numerator. At the theoretical maximum of 65,532 digits, this grows to about 54,610 bytes total, which approaches but remains within the 6502's addressing capacity.

## Architectural considerations

Several design choices specifically accommodate the 6502's limitations and characteristics. The use of signed `int16_t` for array indices $L$ and $M$ prevents subtle bugs that arise from comparing signed loop counters with unsigned array bounds on 16-bit systems. The little-endian byte ordering matches the 6502's native format, simplifying multi-byte operations.

Memory allocation includes padding beyond the computed bignum size to ensure safe 32-bit read operations when extracting output digits. This prevents potential crashes from reading beyond allocated memory boundaries, a critical consideration on systems without memory protection.

The spigot nature of Bellard's algorithm proves particularly well-suited to resource-constrained systems like the 6502. Rather than computing and storing all digits before output, each digit can be displayed as soon as it's calculated, allowing the computation of arbitrarily long pi expansions without proportional memory growth for the output itself.

## Scaling operations in the algorithm

Two specific multiplication operations are central to the algorithm's operation and deserve explanation. The multiply-by-1000 operation extracts output digits by shifting the fractional part of the accumulator into the integer position. Since we extract three decimal digits per iteration, multiplying by 1000 moves the next group of three digits from the fractional part into the integer part where they can be read and output.

The multiply-by-250 operation implements the rescaling factor $\frac{250}{256}$ that maintains precision balance between iterations. This ratio equals approximately 0.9765625, which slightly reduces the numerator each cycle to prevent unlimited growth. The division by 256 is implemented as a byte shift (moving all bytes down one position), making $\frac{250}{256}$ an efficient way to implement the mathematically required scaling factor. This rescaling prevents the numerator from growing without bound while maintaining sufficient precision for accurate digit extraction.
