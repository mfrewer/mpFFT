## Parallel, multiprecision FFT code options

__All codes in this category (parallel, multiprecision) are of the type:__ decimation-in-time (DIT), split-radix, depth-first-recursive, complex-to-complex,
forward direction, and twiddle-index-look-up through no-pass-by-multiplication.

__All timing results shown below were obtained with `gcc -O3 -fopenmp -lmpfr -lgmp -lm` on the machine:__ Intel(R)Core(TM)i5-&#8288;6500CPU@<!-- -->3.20GHz,
Ubuntu(20.04) 64-bit VMware Workstation.

### \# Input options

__`#define IN 1`__\
Complex random input for both `inr` (real part) and `ini` (imaginary part), uniformly distributed between 0 and 1. This input is used for time measurements.

__`#define IN 2`__\
Complex fixed input: `inr[k]`=1/(1+k<sup>2</sup>), `ini[k]`=1/(1+k<sup>4</sup>). This input is used to check the codes' correctness.

### \# Code options

__`#define CODE 1`__\
__Features:__ SR-2/4, 3m3a, 3-butterfly, recursion-break: n=4, openmp (omp for): n&gt;4, cores: 4\
__Results:__

* Average time per FFT(Cores=4,P=1024,N=1024): T ~ 1.73 ms
* FLOPs: mults = 7172, adds = 27652, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ This code is to be compared with its corresponding serial implementation (CODE 3). For 4 CPUs, a speed-up of only 1.7x is
gained over its serial version &mdash; see [README](/src/parallel/mp/README.md) for the explanation of this small performance gain.

As for the serial case, here again P=1024 denotes the binary precision, which corresponds to a decimal precision dP=P·log<sub>10</sub>2 ~ 308.

__`#define CODE 2`__\
__Features:__ SR-2/4, 3m3a, 3-butterfly, recursion-break: n=4, openmp (omp for): n&gt;8, cores: 4\
__Results:__

* Average time per FFT(Cores=4,P=1024,N=1024): T ~ 1.58 ms
* FLOPs: mults = 7172, adds = 27652, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ A small speed-up is gained by reducing the overhead of openmp, in deactivating the parallel routine below a minimal FFT-size (here n=16).

__`#define CODE 3`__\
__Features:__ SR-2/4, 3m3a, 3-butterfly, recursion-break: n=4, openmp (omp task): n&gt;8, cores: 4\
__Results:__

* Average time per FFT(Cores=4,P=1024,N=1024): T ~ 1.60 ms
* FLOPs: mults = 7172, adds = 27652, total = 34824
* rRMSE: ∝ 1e-308

__Comments:__ Alternative to `CODE 2`, using the directive `omp task` instead of `omp for`.

### \# Option for the fastest code

__`#define CODE 'F'`__\
This option will select the fastest code in the present category (parallel, multiprecision),
which currently for 4 cores, binary precision P=1024 and input size N=1024 is `CODE 2`.

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
