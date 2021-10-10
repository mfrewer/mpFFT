## mpFFT

__Multiprecision Fast Fourier Transform__

*by Michael Frewer © 2021*

### Synopsis

mpFFT is an open-source project to implement a high-performance multiprecision Fast Fourier Transform that can compete with non-free software as Mathematica and MATLAB, in both serial and parallel computations.

The need for such a project is clear: The FFT routine is one of the main workhorses in scientific computing and usually runs as a subroutine in a larger program.
Therefore, if a multiprecision computation or result is needed, it is natural to have an efficient and fast implementation that can be easily embedded into an existing project
without requiring a redesign or migration onto a new (mostly non-free) platform.

Many fast open-source FFT packages are available for single- and double-precision, such as the projects from [Ooura](https://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
and [GSL](https://www.gnu.org/software/gsl/doc/html/fft.html) for scalar implementations, or, for the significantly faster vectorised implementations,
the projects [FFTW](https://www.fftw.org/), [FFTE](http://www.ffte.jp/) and  [FFTS](https://github.com/anthonix/ffts), just to name a few popular ones.
For multiprecision computing, however, the need for fast FFT packages has not yet been met and great potential still exists.
Especially since the demand for multiprecision calculations has grown in recent years [[Higham (2017)]](https://doi.org/10.1109/ARITH.2017.24).
One of the issues, for example, is that many scientific computing problems on closer inspection turn out to be hybrid, in that certain parts of the problem can be solved with lower precision,
while other parts need higher precision.

To note is that this project is still in development and far from complete. However, the first implementations available already show promising results:
In serial computation with varying precision size (up to 1 million digits), a more than [2x speed-up](/bench/serial_mp_N1024_time.png) to Mathematica and a [10x speed-up](/bench/serial_mp_N1024_time.png) to MATLAB is achieved,
where for MATLAB an external toolbox had to be used, the for FFT only available non-free multiprecision toolbox [ADVANPIX](https://www.advanpix.com/).
Also for parallel computations with 4 CPU cores, a [10x speed-up](/bench/parallel_mp_C4_N1024_time.png) to MATLAB is achieved for high precision orders &mdash; for Mathematica, unfortunately,
no comparison is possible, as it does not provide any multicore support for its FFT routine.

### Design

`mpFFT.c` is the main program that selects the currently most efficient multiprecision implementations developed in&nbsp;[lab](/lab) for serial and parallel computation. It is the result of an ongoing development
in [lab](/lab) to construct a high-performance multiprecision FFT algorithm through small independent modules, so-called codes `*_mp_fft*.c`, each of which can be compiled and tested separately. This allows for a better and finer control over which implementations and multiprecision arithmetic schemes are performance-hindering and which are performance-enhancing.

There will be a choice between split-radix algorithms of different orders, unscaled and scaled, different hard-coded base-case FFTs to terminate recursion,
two different schemes of complex multiplication (4m2a, 3m3a), different indexing schemes of the twiddle-factor lookup tables, and finally the choice of choosing either 1, 2 or 3 butterfly operations. All these options will first be studied in double-precision, with the most efficient ones then serving as templates for the multiprecision implementations.

For parallel computations there will be an additional choice between [OpenMP](https://www.openmp.org/) and a self-designed threadpool using the `gcc` library `pthread.h`.

### Implementation

All codes in this project will be:

* __only 1D:__ The reason for this is to first find the most efficient FFT algorithm and implementation that suits multiprecision computation best, and this is naturally to be done in 1D.
Once it is found, it can then form the basis for all higher-dimensional FFTs, since they are all just iterated 1D FFTs.
To achieve efficient implementations, however, the problem then inevitably shifts to the parallelization issue where the main performance bottleneck in higher dimensions is the communication among the CPUs.
Nevertheless, a fast 3D FFT relies on a fast 1D FFT.

* __only scalar:__ Vectorised instruction sets such as SIMD, which are relevant for single- and double-precision, lose their importance when calculating with multiprecision numbers going far beyond 64- or 128-bit,
since the requested jobs are too large to fit into the CPU caches. The arithmetic efficiency of this project fully relies on the implementation of the underlying GNU library [MPFR](https://www.mpfr.org/),
which will be used throughout to perform all multiprecision computations.

* __only 2<sup>n</sup>:__ Input signal lengths of power 2 will be considered only.

* __only DIT:__ The FFT approach throughout this project will be based on the decimation-in-time (DIT) decomposition,
where the indices of the input sequence are separated into even and odd classes. The dual approach, the decimation-in-frequency (DIF) where the output sequence is divided, will not be considered.

* __only split-radix:__ The reason for choosing the split-radix algorithm is the advantage of having low complexity, since it aims to compute the FFT with the least number of multiplications.
For single- and double-precision computations the complexity-issue is not so much of a concern, but becomes highly relevant for multiprecision computations,
as multiplications become increasingly more expensive than additions the higher the order of precision gets.
Different classes and variations of the split-radix algorithm will be implemented and tested, to find for multiprecision the most efficient one.

* __only recursive:__ The reason for choosing a recursive rather than an iterative scheme is the advantage of memory locality [[Frigo&nbsp;&amp;&nbsp;Johnson&nbsp;(1998)]](https://doi.org/10.1109/ICASSP.1998.681704),
which is a critical component for fast multiprecision computations. In addition, a recursive first-depth scheme eliminates the need for the computationally expensive bit-reversal permutation in the indices,
since the recursion already implicitly performs the permutation.

  Overall, the recursive DIT-split-radix implementation here follows the notation and derivation as concisely
  and directly given in [Wikipedia: Split-radix FFT algorithm](https://en.wikipedia.org/wiki/Split-radix_FFT_algorithm).

* __only complex-split format:__ The real and imaginary part of complex data will be stored in separate arrays. The interleaved format,
where the real and imaginary part are stored adjacently through a single array in memory, is not used.

* __only out-of-place:__ Two different arrays `in` and `out` are used for input and output. The input is not overwritten by the output as the program executes.
Due to the complex-split format used, the input and output are thus controlled by four separate arrays: `inr`, `ini` for the real and imaginary input, and `outr`, `outi` for the real and imaginary output.

* __only complex input:__ Complex-to-complex FFTs will be considered only, where both input arrays `inr` and `ini` are non-zero.
To generate efficient real-to-complex FFTs, where `ini=0`, is, in the setting of this project, straightforward to achieve:
Due to the DIT-split-radix and complex-split format used, the input array `ini` needs only to be removed from the code,
while operations for the complex output array X<sub>k</sub>:=`outr[k]+I*outi[k]` can be easily reduced into half,
by recognizing the redundancy in the computation of the elements X<sub>k+2N/4</sub> = (X<sub>(N/4&#8209;k)+N/4</sub>)\* and X<sub>k+3N/4</sub> = (X<sub>N/4-k</sub>)\*,
i.e., X<sub>k+2N/4</sub> and X<sub>k+3N/4</sub> can be determined from X<sub>k+1N/4</sub> and X<sub>k+0N/4</sub>, respectively,
due to the Hermitian symmetry of the discrete Fourier Transform (DFT) for real input: X<sub>k</sub> = (X<sub>N-k</sub>)\*, where 0≤k&lt;N/4, and N the input length and&nbsp;*&nbsp;symbolizing the complex conjugate.
This holds both for the classical and conjugate-pair split-radix algorithm.

  Note that while the decimation-in-time (DIT) approach is the natural choice for the real-to-complex FFT,
  the dual decimation-in-frequency (DIF) approach is the natural choice for the inverse complex-to-real FFT [[Sorensen et al. (1987)]](https://doi.org/10.1109/TASSP.1987.1165220).

* __only forward direction:__ All FFTs are performed here only in forward direction, with the weight factor (root-of-unity) defined as <nobr>ω<sub>N</sub>=*e*<sup>-&#8288;2πi/N</sup>.</nobr>
The inverse FFT can be obtained straightforwardly from the forward FFT without any additional computational cost,
when efficiently implemented as presented in [Duhamel et al. (1988)](https://doi.org/10.1109/29.1519). The inverse FFT will be used herein for error analysis.

* __arbitrary precision:__ No restrictions will be imposed on the order of precision; it can be set arbitrarily.
The upper (binary) precision limit is dictated through `long int MPFR_PREC_MAX`, and on a 64-bit machine given by: (2<sup>64-&#8288;1</sup>-&#8288;1)-&#8288;2<sup>8</sup>,
which is about 9·10<sup>18</sup>. However, since MPFR needs to increase the precision internally, in order to provide accurate results and correct rounding,
it is not recommended to set the precision to any value near `MPFR_PREC_MAX`.

### Compiling

All builds used `gcc` (v9.3.0), with full optimization flag <code>&#8209;O3</code>,
and static linking <code>&#8209;lmpfr&nbsp;&#8209;lgmp&nbsp;&#8209;lm</code> to the GNU libraries [MPFR&nbsp;(v4.1.0)](https://www.mpfr.org/), [GMP&nbsp;(v6.2.1)](https://gmplib.org/)
and the standard C math-library. For parallel computations the extra flags <code>&#8209;pthread&nbsp;&#8209;fopenmp</code> need to be included.

Make sure that MPFR is installed 'thread-safe' in order to run correctly and reliably in parallel. This is done by setting the option <code>&#8209;&#8209;enable&#8209;thread&#8209;safe</code> during installation.

In the current stage, the project's main code can be compiled and run with\
<code>gcc -O3 mpFFT.c -fopenmp -lmpfr -lgmp -lm -o mpFFT && ./mpFFT</code>

To run a code in [lab](/lab), it is divided into several categories

* [serial & double-precision](/lab/serial/double)
* [serial & multi-precision](/lab/serial/mp)
* [parallel & double-precision](/lab/parallel/double)
* [parallel & multi-precision](/lab/parallel/mp)

each containing a `*_main.c` file which can be compiled independently from all other ones along with the project's header file `mpFFT.h`.
Inside each category different implementations and inputs can be chosen and selected in `*_main.c` through the macros `#define CODE` and <code>#define&nbsp;IN</code>.
For more details, please see the accompanying OPTIONS- and README-files.

### Documentation

For further information and more details on particular implementations and their benchmarks,
please see the other README-files in each of the different categories of this project.
A detailed LaTeX-written documentation on all theoretical aspects and mathematical underpinnings of this project will follow in due course.

### License

The content of this project itself is licensed under the [Creative Commons Attribution 4.0 International (CC BY 4.0) License](https://creativecommons.org/licenses/by/4.0/),
and all underlying source codes to compile and to run this project are licensed under the [Apache-2.0 License](http://www.apache.org/licenses/LICENSE-2.0).

### Other multiprecision FFT implementations
* mpfft: MPFR FFT radix-2 functions in C: [https://github.com/urrfinjuss/mpfft](https://github.com/urrfinjuss/mpfft)
