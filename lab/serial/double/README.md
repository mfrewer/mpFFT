## Documentation for serial, double-precision FFT codes

The codes presented here in this category (serial, double-precision) serve as a prototype and template for the corresponding multiprecision codes. They are intentionally designed as scalar implementations and do therefore not compete with any vectorised single- or double-precision codes such as FFTW [1] or FFTS [2].

The reason for adhering to scalar codes is that when translating to higher-order floating-point computations, a memory problem naturally and automatically arises, since the requested jobs are too large to fit into the CPU caches. Vectorised instruction sets such as SIMD, which are relevant for single- and double-precision, lose their importance when calculating with multiprecision numbers going far beyond 64- or 128-bit. Nevertheless, the herein presented double-precision codes have been optimized and can compete with fast existing scalar codes, e.g., as presented and developed in [3-4].

Timing and performance benchmarks are given, comparing the currently fastest code obtained herein to the two scalar codes [3] and [4] (see&nbsp;[T<sub>1</sub>](/bench/serial_double_scalar_time.png), [P<sub>1</sub>](/bench/serial_double_scalar_flops.png)), as well as to the much faster vectorized codes from Mathematica and MATLAB (see [T<sub>2</sub>](/bench/serial_double_vector_time.png), [P<sub>2</sub>](/bench/serial_double_vector_flops.png)), where it's to note that [MATLAB's fft-implementation](https://de.mathworks.com/help/matlab/ref/fft.html) is based on FFTW [1]. Also, all benchmarks shown relate solely to the run-time of the FFT-kernel &mdash; any initialization times are not included.

Interesting to see are the different structures in the scalar performance benchmark [P<sub>1</sub>](/bench/serial_double_scalar_flops.png): While mpFFT shows a steady decline in performance with increasing input size throughout, code [3] first gains a maximum and then decreases, while code [4] decreases and then stabilizes.

To note is that by default MATLAB computes in multi-core mode. To force serial (single-core) computations, MATLAB has to be started with the flag <code>&#8209;singleCompThread</code>.

### Implementation details

Presented here are decimation-in-time (DIT) implementations of the split-radix algorithm, based on a recursive first-depth scheme, including both the classical split-radix as well as the conjugate-pair algorithm by using two different versions of indexing the twiddle-factor lookup tables. The implementation follows the notation of [5].

* The reason for choosing a recursive instead of an iterative scheme is the advantage of memory locality [6-7], which will be a crucial component for fast multiprecision computations. In addition, a recursive first-depth scheme eliminates the need for the computationally expensive bit-reversal permutation in the indices, since the recursion already implicitly performs the permutation.

* The reason for choosing the split-radix algorithm is the advantage of having low complexity, since it aims to compute the FFT with the least number of multiplications [8-12]. For single- and double-precision computations the complexity-issue is not so much of a concern, but becomes highly relevant for multiprecision computations, as multiplications become increasingly more expensive than additions the higher the order of precision gets &mdash; while addition scales linearly with the digit number length n, multiplication scales (in the best case) with n·log(n) [13].

  To note is that the current record for the lowest FFT-complexity is still given by [11]. The new record announced in [14] cannot be confirmed, as was already mentioned first in [15].

### Timing, complexity and error measures

1. The FFT kernel will be timed by iterating repeatedly over it and taking the average. This time measurement is controlled by the iteration parameter `M`. To achieve good timing results, `M` should not be chosen too small or too large, otherwise the result will be spoiled by internal system events. The strategy is to start with a small value `M` and then gradually increase it until the timing result for the first time proves to be stationary or stable within a few percent.

    The generation of the time statistics has intentionally not been automated, but has to be done in each case manually by the user himself. This allows for better and finer control and also avoids system dependencies, as each machine behaves differently.

    Next to the run-time of the FFT kernel, its initialization will be timed too. But since in practice all initialization functions for a FFT-routine needs to be called only once, it will also be timed only once. Hence this timing result is not a performance measure, but merely severs as a waiting-time indicator until the FFT-kernel is ready to work and therefore not included in any benchmarks.

2. The FFT kernel is equipped with four global complexity counters, counting for a single FFT the exact number of multiplications `cmul`, additions `cadd`, recursion calls `cfc` and recursion breaks `crb`.

    The measure `cmul + cadd` is the total number of floating-point operations (FLOPs) performed. In ratio with the mean averaged time it gives the performance measure FLOP/s (floating-point operations per second). However, as is usual in the FFT-literature [1-2,7,11], not the exact number of FLOPs is taken, but rather the asymptotic number of floating-point operations 5N·log<sub>2</sub>(N) of the radix-2 Cooley-Tukey algorithm is used as the reference point. Hence, FFT performance is therefore measured in so-called CTGs (Cooley-Tukey gigaflops per second), a pseudo-performance measure computed as 5N·log<sub>2</sub>(N)/t, where t is the run-time in seconds and N the input size. CTGs is a *"rough measure of a particular implementation's efficiency relative to the radix-2 algorithm and the clock speed of the machine"* [2], and thus should be *"viewed as convenient scaling factor rather than as an absolute indicator of CPU performance"* [1].

    Next to the exact total FLOPs performed, also the optimal and currently lowest complexity count of the split-radix algorithm will be displayed for comparison [8-12]. However, as already said and as will be explicitly shown, the total FLOPs in double-precision is not so much the decisive factor to gain better performance. Relevant is rather the overall memory structure and which particular hard-coded short-FFT-algorithms are used to break the recursion and therefore reduce the total number of recursion calls.

   Important to note here is that the complexity count can be switched off in the FFT kernel when performing time measurements. This can be controlled by setting either 1 or 0 (on or off) in the preprocessor macro `#define CC`. If switched off to gain speed, the compiler should nonetheless be flagged to make use of CSE (common subexpression elimination), which for `gcc`, however, is automatically supported with `-O3`, the recommended optimization flag to compile the codes in this project. Otherwise the flag <code>&#8209;fcse&#8209;skip&#8209;blocks</code> has to be given explicitly, to remove the counting code-lines from the kernel during runtime.

3. The error measurements come in two categories:

    The first, slow one is by comparing the FFT-result directly to the defining discrete Fourier transform (DFT). This check only serves as a first orientation whether the FFT-algorithm as a whole was implemented correctly or not.

   The second category provides a detailed error analysis, by transforming the FFT-result back to its input using the inverse FFT, efficiently implemented as proposed in [16]. The absolute error, relative error and the relative root-mean-square-error (rRMSE) between the original (true) and computed (approximated) input is then determined. While the absolute error is a measure of the precision used, the relative error and rRMSE are measures of the accuracy achieved, where it is to note that the relative error gets ill-defined when the true value is close to zero (within the precision used). In such a case, the relative error is to be ignored and the appropriate global measure rRMSE is to be used instead.

    The first column-pair in `data_*_test_error` shows the absolute error, the second pair the relative error and the third pair the global accuracy measure rRMSE.

   To note is that the error measurements do not influence or effect the time measurement of the FFT kernel, since these routines operate separately from the kernel. Nevertheless, they too can be switched off individually, to reduce the overall runtime of the FFT-code chosen.

##
### Definitions

![absolute error](https://latex.codecogs.com/png.image?\dpi{105}%20\text{absolute%20error:%20}%20\Delta%20x[k]:=x_{\text{true}}[k]-x_\text{approx.}[k],%20k=1,\ldots,N)

![relative error](https://latex.codecogs.com/png.image?\dpi{105}%20\text{relative%20error:%20}%20\frac{\Delta%20x[k]}{x_\text{true}[k]})

![rRMSE](https://latex.codecogs.com/png.image?\dpi{105}%20\text{rRMSE:%20}%20\frac{\sqrt{\frac{1}{N}\sum_{k}\big|\Delta%20x[k]\big|^2}}{\sqrt{\frac{1}{N}\sum_k\big|x_\text{true}[k]\big|^2}})

### References
[1] [Frigo &amp; Johnson, 2005: The design and implementation of FFTW3](https://doi.org/10.1109/JPROC.2004.840301)\
[2] [Blake et al., 2013: The fastest Fourier transform in the south](https://doi.org/10.1109/TSP.2013.2273199)\
[3] [GSL - GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/fft.html)\
[4] [Ooura, 2001: General Purpose FFT Package](https://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)\
[5] [Wikipedia: Split-radix FFT algorithm](https://en.wikipedia.org/wiki/Split-radix_FFT_algorithm)\
[6] [Frigo &amp; Johnson, 1998: FFTW: An adaptive software architecture for the FFT](https://doi.org/10.1109/ICASSP.1998.681704)\
[7] [Chellappa et al., 2007: How to write fast numerical code: A small introduction](https://doi.org/10.1007/978-3-540-88643-3_5)\
[8] [Yavne, 1968: An economical method for calculating the discrete Fourier transform](https://doi.org/10.1145/1476589.1476610)\
[9] [Duhamel &amp; Hollmann, 1984: Split-radix FFT algorithm](https://doi.org/10.1049/el:19840012)\
[10] [Sorensen et al., 1986: On computing the split-radix FFT](https://doi.org/10.1109/TASSP.1986.1164804)\
[11] [Johnson &amp; Frigo, 2007: A modified split-radix FFT with fewer arithmetic operations](https://doi.org/10.1109/TSP.2006.882087)\
[12] [Burrus, 2012: Fast Fourier Transforms](https://open.umn.edu/opentextbooks/textbooks/261), [(2nd link)](https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Signal_Processing_and_Modeling/Book%3A_Fast_Fourier_Transforms_(Burrus))\
[13] [Harvey &amp; van der Hoeven, 2021: Integer multiplication in time O(n log n)](https://doi.org/10.4007/annals.2021.193.2.4)\
[14] [Zheng et al., 2014: Scaled radix-2/8 algorithm for efficient computation of length-N=2<sup>m</sup> DFTs](https://doi.org/10.1109/TSP.2014.2310434)\
[15] [Sergeev, 2017: On the real complexity of a complex DFT](https://doi.org/10.1134/S0032946017030103)\
[16] [Duhamel et al., 1988: On computing the inverse DFT](https://doi.org/10.1109/29.1519)
