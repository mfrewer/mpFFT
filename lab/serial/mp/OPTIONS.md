## Serial, multiprecision FFT code options

__All codes in this category (serial, multiprecision) are of the type:__ decimation-in-time (DIT), split-radix, depth-first-recursive, complex-to-complex, forward direction

__All timing results shown below were obtained with `gcc -O3 -lmpfr -lgmp -lm` on the machine:__ Intel(R)Core(TM)i5-6500CPU@<!-- -->3.20GHz, Ubuntu(20.04) 64-bit VMware Workstation.

### \# Input options

__`#define IN 1`__\
Complex random input for both `inr` (real part) and `ini` (imaginary part), uniformly distributed between 0 and 1. This input is used for time measurements.

__`#define IN 2`__\
Complex fixed input: `inr[k]`=1/(1+k<sup>2</sup>), `ini[k]`=1/(1+k<sup>4</sup>). This input is used to check the codes' correctness.

### \# Code options

__`#define CODE 1`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(P=1024,N=1024): T ~ 3.47 ms
* FLOPs: mults = 11376, adds = 26168, total = 37544
* rRMSE: ∝ 1e-308

__Comments:__ Note that P=1024 is the binary precision, which corresponds to a decimal precision dP=P·log<sub>10</sub>2 ~ 308.

__`#define CODE 2`__\
__Features:__ SR-2/4, 4m2a, 3-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(P=1024,N=1024): T ~ 3.32 ms
* FLOPs: mults = 9336, adds = 25488, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ A speed-up is gained by reducing the number of FLOPs, which is achieved by implementing 3 butterflies instead of only one.

__`#define CODE 3`__\
__Features:__ SR-2/4, 3m3a, 3-butterfly, recursion-break: n=4, twiddle-index look-up: no-pass by multiplication\
__Results:__

* Average time per FFT(P=1024,N=1024): T ~ 3.02 ms
* FLOPs: mults = 7172, adds = 27652, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ A further significant speed-up is gained by reducing the number of multiplications, which is achieved by using the 3m3a complex multiplication scheme instead of the usual 4m2a scheme. Note that although the total number of FLOPs is the same as in `CODE 2`, a lower multiplication count yields a significant speed-up, simply due to that multiplications in multiprecision get more expensive than additions.

__`#define CODE 4`__\
__Features:__ SR-2/4, 3m3a, 3-butterfly, recursion-break: n=4, twiddle-index look-up: pass by shift\
__Results:__

* Average time per FFT(P=1024,N=1024): T ~ 3.04 ms
* FLOPs: mults = 7172, adds = 27652, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ Alternative to `CODE 3`, using a different index generator for the look-up table. Instead of directly calculating the index from `n` through multiplication,
a new additional variable `log2n` is passed to calculate the index through a bit-shift. This idea is taken from [1,2].

### \# Option for the fastest code

__`#define CODE 'F'`__\
This option will select the fastest code in the present category (serial, multiprecision), which currently for binary precision P=1024 and input size N=1024 is <code>CODE&nbsp;3</code>.

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
[2] [Blake et al., 2013: The fastest Fourier transform in the south](https://doi.org/10.1109/TSP.2013.2273199)
