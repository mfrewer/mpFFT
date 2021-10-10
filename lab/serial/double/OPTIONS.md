## Serial, double-precision FFT code options

__All codes in this category (serial, double-precision) are of the type:__ decimation-in-time (DIT), split-radix, depth-first-recursive, complex-to-complex, forward direction.

__All timing results shown below were obtained with `gcc -O3 -lm` on the machine:__ Intel(R)Core(TM)i5-6500CPU@<!-- -->3.20GHz, Ubuntu(20.04) 64-bit VMware Workstation.

### \# Input options

__`#define IN 1`__\
Complex random input for both `inr` (real part) and `ini` (imaginary part), uniformly distributed between 0 and 1. This input is used for time measurements.

__`#define IN 2`__\
Complex fixed input: `inr[k]`=1/(1+k<sup>2</sup>), `ini[k]`=1/(1+k<sup>4</sup>). This input is used to check the codes' correctness.

### \# Code options

__`#define CODE 1`__\
__Features:__ SR-2/4, 4m2a, 3-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(N=1024): T ~ 11.4 μs
* FLOPs: mults = 9336, adds = 25488, total = 34824
* rRMSE: ∝ 1e-16

__Comments:__ -

__`#define CODE 2`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(N=1024): T ~ 11.2 μs
* FLOPs: mults = 11376, adds = 26168, total = 37544
* rRMSE: ∝ 1e-16

__Comments:__ Although the complexity (FLOPs) is considerably higher (8% more FLOPs) than in `CODE 1`, due to using only 1 butterfly, there is no time penalty.
On the contrary, it even executes a bit faster than `CODE 1`.
This proves again that for double-precision computations the complexity itself is not the decisive factor for a fast FFT on modern CPUs.
For multiprecision computations, however, this will no longer be the case anymore &mdash; for those, a reduction in complexity is critical for a fast FFT,
in particular the reduction of multiplications.

__`#define CODE 3`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=16, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(N=1024): T ~ 8.9 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ We see a considerable speed-up (20% faster) compared to `CODE 1` and `CODE 2`, due to invoking hard-coded FFTs up to size n=16,
allowing therefore for less recursion calls. The number of recursive function calls is reduced by factor of 4, from `cfc=511` (n=4) down to `cfc=127` (n=16).
Unfortunately, for multiprecision computations this reduction of `cfc` does not play such a decisive role anymore.
Note that the FLOPs differ from `CODE 2`, despite having the same 1-butterfly complexity for n&gt;16.
This difference is due to using algorithms with a different complexity for the FFTs n=8 and n=16, which here in `CODE 3` are both hard-coded,
while in `CODE 2` they are generated. The hard-coded algorithms used are from [Winograd (1978)](https://doi.org/10.1090/S0025-5718-1978-0468306-4).

__`#define CODE 4`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=16, twiddle-index look-up: pass by shift\
__Results:__

* Average time per FFT(N=1024): T ~ 8.7 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ Alternative to `CODE 3`, using a different index generator for the look-up table.
Instead of directly calculating the index from `n` through multiplication, a new additional variable `log2n` is passed to calculate the index through a bit-shift.
This idea is taken from [1,2].

__`#define CODE 5`__\
__Features:__ SR-2/4, conjugate-pair, 4m2a, 3-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(N=1024): T ~ 12.2 μs
* FLOPs: mults = 9336, adds = 25488, total = 34824
* rRMSE: ∝ 1e-16

__Comments:__ This code is to be compared with `CODE 1`, using the same complexity but a different approach,
the conjugate-pair approach [1-3] of the split-radix algorithm. The advantage of this approach is that only 1 complex look-up table is needed
instead of 2 as in the standard split-radix algorithm. Its disadvantage, however, is that a conditional variable `sn` needs to be
passed to account for the cyclic shift in the indices. Since the latter outweighs the former in this implementation,
it results in a less efficient and therefore slower code than `CODE 1`. However, this performance drastically
changes to the better when switching to a different implementation, using only 1 butterfly and a recursion-break of higher order, as done next in `CODE 6`.

__`#define CODE 6`__\
__Features:__ SR-2/4, conjugate-pair, 4m2a, 1-butterfly, recursion-break: n=16, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(N=1024): T ~ 8.1 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ This code is to be compared with `CODE 3`, using the same complexity but the conjugate-pair approach,
gaining a further speed-up of 7%. Contrary to `CODE 5`, the advantage of having only 1 complex look-up table outweighs now its
disadvantage of having to pass a conditional variable.

__`#define CODE 7`__\
__Features:__ SR-2/4, conjugate-pair, 4m2a, 1-butterfly, recursion-break: n=16, twiddle-index look-up: pass by shift\
__Results:__

* Average time per FFT(N=1024): T ~ 7.9 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ This code is to be compared with `CODE 4`, using the same complexity and the same twiddle-index look-up scheme but the conjugate-pair approach,
gaining thus the same speed-up as in `CODE 6`.

### \# Option for the fastest code

__`#define CODE 'F'`__\
This option will select the fastest code in the present category (serial, double-precision), which currently for input size N=1024 is <code>CODE&nbsp;7</code>.
For larger input sizes, the efficiency slowly shifts towards `CODE 6` as the fastest,
but the difference between these two codes always remains moderate within a few percent,
therefore being in the end a matter of taste which code to choose.

## Acronyms
FFT: Fast Fourier Transform\
DFT: Discrete Fourier Transform\
DIT: Decimation in Time Algorithm\
SR-s/r: Split-Radix Algorithm, with radix-r for the even and radix-s for the uneven indices\
CP: Conjugate Pair Algorithm\
4m2a: complex-multiplication-scheme with 4 multiplications and 2 additions\
3m3a: complex-multiplication-scheme with 3 multiplications and 3 additions\
FLOPs: number of floating point operations\
FLOP/s: floating point operations per second\
rRMSE: relative root-mean-square-error

## References
[1] [Blake, 2012: Computing the fast Fourier transform on SIMD microprocessors](https://researchcommons.waikato.ac.nz/handle/10289/6417)\
[2] [Blake et al., 2013: The fastest Fourier transform in the south](https://doi.org/10.1109/TSP.2013.2273199)\
[3] [Johnson &amp; Frigo, 2007: A modified split-radix FFT with reduced arithmetic complexity](https://doi.org/10.1109/TSP.2006.882087)
