# atto-pi

![Commodore 64 calculating 50000 digits of pi](https://github.com/user-attachments/assets/06a573ca-0f0e-4c05-b425-a32cbac885e8)

Author: John Byrd <johnwbyrd at gmail dot com>

This is a tiny C program for computing pi in small embedded environments. It was written for the LLVM-MOS C/C++ compiler, but it should work on other reasonably standards compliant compilers as well.

When compiled using [LLVM-MOS](https://wwww.llvm-mos.org) for the unexpanded Commodore 64 computer, it can generate up to 50k successive digits of pi.  Total code size is around 8KB, and the remainder of the C64's memory is used for bignums.

This program implements [François Bellard's 1997 decimal spigot algorithm for computing pi](https://bellard.org/pi/pi_bin/pi_bin.html), which is a spigot algorithm that generates decimal digits sequentially without storing the entire result. Unlike the more famous [Bailey-Borwein-Plouffe formula](https://observablehq.com/@rreusser/computing-with-the-bailey-borwein-plouffe-formula), this decimal adaptation of Bellard's binary series works in base-1000, producing three decimal digits per iteration rather than hexadecimal digits.

The mathematical foundation is the infinite series:

$$\pi = \frac{1}{2^{6}} \sum_{n=0}^{\infty} \frac{(-1)^{n}}{2^{10n}} \left[ \frac{-2^{5}}{4n+1} - \frac{1}{4n+3} + \frac{2^{8}}{10n+1} - \frac{2^{6}}{10n+3} - \frac{2^{2}}{10n+5} - \frac{2^{2}}{10n+7} + \frac{1}{10n+9} \right]$$


## Arithmetic implementation

The core challenge lies in performing high-precision arithmetic on a system with only 8-bit native operations and 16-bit addressing. The solution uses "bignums" -- arbitrary-precision numbers represented as arrays of bytes in little-endian format. Each bignum stores a fixed-point number where the integer part occupies the highest-indexed byte and the fractional part extends downward through lower indices.

The size of each bignum is calculated as $\text{digits} \times \frac{\log_2(10)}{8} + \text{guard}$, where the constant $\frac{\log_2(10)}{8} \approx 0.415241$ represents the bytes needed per decimal digit, and the guard digits provide computational headroom to prevent rounding errors from accumulating. This formula derives from the convergence properties of Bellard's series and ensures sufficient precision for the requested output.

Division operations use binary long division, processing one bit at a time through an 8-bit window. This approach naturally fits the 6502's byte-oriented architecture while maintaining the precision needed for thousands of digits. The divisor is repeatedly halved while the remainder is compared against it, building up the quotient bit by bit.

Multiplication by constants (250 and 1000) employs carry propagation across the entire bignum. Each byte is multiplied by the constant, added to any carry from the previous position, and the result is split into a stored byte and a carry for the next position. This process continues until the carry becomes zero, ensuring no precision is lost.

Rescaling operations maintain the mathematical invariant needed for the spigot algorithm. The 250/256 factor specifically addresses the radix conversion error between the binary computation base and the decimal output base. Since 10^3 = 1000 ≈ 1024 = 2^10, each digit extraction introduces a small error that accumulates. The 250/256 rescaling (which equals 1000/1024) exactly compensates for this approximation.

## Algorithmic flow

Each iteration begins by computing the seven terms of Bellard's formula through a sequence of divisions and additions/subtractions into the accumulator bignum. The divisors follow specific patterns based on the iteration variables: $t$ increases by 10 each iteration and appears in expressions like $(t+1)$, $(t+3) \times 4$, and so forth, while $f$ increases by 4 and appears in $(f+1) \times 8$ and $(f+3) \times 256$.

After computing all terms, the algorithm extracts three decimal digits from the integer part of the accumulator. These digits represent the current contribution to pi's decimal expansion. The extracted digits are then masked (zeroed) from the accumulator, and the remaining fractional part is multiplied by 1000 to shift the next three digits into the integer position for the following iteration.

The numerator bignum undergoes rescaling through multiplication by $\frac{250}{256}$, effectively dividing by $\frac{256}{250} = 1.024$. This scaling factor maintains the precision balance between the accumulator and numerator as the algorithm progresses. The division by 256 is implemented as a byte shift, taking advantage of the binary representation.

## Precision and limits

Practical calculation of pi all about limits. It is necessary to track types exactly in order to make sure that C integer overflow does not occur for the digits of pi requested. To this end, every integer type in this program is abstracted semantically.

The algorithm produces three decimal digits per iteration, with the first iteration yielding just the initial digit 3 followed by the decimal point. Each subsequent iteration contributes three more digits to pi's expansion. The performance characteristics vary dramatically across different digit counts - computing 100 digits runs relatively quickly, 1000 digits requires noticeable time, while 10,000 digits on period hardware like a Commodore 64 can take hours due to the quadratic scaling of the arithmetic operations.

Memory consumption follows a predictable linear relationship with the requested digit count. Each bignum requires approximately 0.415241 bytes per decimal digit plus a small guard factor for computational headroom. For 1,000 digits, each bignum occupies roughly 418 bytes, totaling 836 bytes for both the accumulator and numerator. At 10,000 digits, this grows to about 4,155 bytes per bignum.

## Thanks

This particular implementation was inspired by [David Banks (hoglet67)'s recent work](https://github.com/BigEd/pi-spigot-for-micros) on implementing spigots on the BBC Micro.

Thanks also to [mysterymath](https://github.com/mysterymath), christen, and the other people who verified this implementation on real hardware.