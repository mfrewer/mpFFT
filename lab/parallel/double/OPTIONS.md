## Parallel, double-precision FFT code options

__All codes in this category (parallel, double-precision) are of the type:__ decimation-in-time (DIT), split-radix, depth-first-recursive, complex-to-complex,
forward direction, and twiddle-index-look-up through no-pass-by-multiplication.

__All timing results shown below were obtained with `gcc -O3 -fopenmp -lm` on the machine:__ Intel(R)Core(TM)i5&#8209;6500CPU@<!-- -->3.20GHz,
Ubuntu(20.04) 64&#8209;bit VMware Workstation.

### \# Input options

__`#define IN 1`__\
Complex random input for both `inr` (real part) and `ini` (imaginary part), uniformly distributed between 0 and 1. This input is used for time measurements.

__`#define IN 2`__\
Complex fixed input: `inr[k]`=1/(1+k<sup>2</sup>), `ini[k]`=1/(1+k<sup>4</sup>). This input is used to check the codes' correctness.

### \# Code options

__`#define CODE 1`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=16, openmp (omp for): n&gt;128, cores: 4\
__Results:__

* Average time per FFT(Cores=4,N=1024): T ~ 38.7 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ This code is to be compared with its corresponding serial implementation (CODE 6).
It is considerably slower than its serial version, showing that for double-precision the openmp-overhead for 4 CPUs is too high for this case to benefit from it.
For multiprecision, however, this will no longer be the case anymore.

__`#define CODE 2`__\
__Features:__ SR-2/4, 4m2a, 1-butterfly, recursion-break: n=16, openmp (omp for): n&gt;128, cores: 2\
__Results:__

* Average time per FFT(Cores=2,N=1024): T ~ 29.3 μs
* FLOPs: mults = 9668, adds = 25828, total = 35496
* rRMSE: ∝ 1e-16

__Comments:__ Same as `CODE 1`, except that it now runs on 2 cores. The speed-up by using less CPUs underlines once again that parallelization
here in this case for double-precision is not an advantage but an obstruction. For multiprecision, however, this will no longer be the case anymore.

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
