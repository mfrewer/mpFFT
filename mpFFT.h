
/*
  SPDX-FileCopyrightText: Copyright © 2021 Michael Frewer <frewer.science@gmail.com>
  SPDX-License-Identifier: Apache-2.0
*/

/* All libraries, forward struct and function prototype declarations to run mpFFT with either double- or multi-precision */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

/* serial, double */
unsigned long long int get_time_int(void);
void code_info();
void create_memory();
void lut();
void ic();
void d_sr24_fft(unsigned long,unsigned long,double*,double*,double*,double*);
void d_sr24_fft_v2(unsigned long,unsigned long,unsigned long,double*,double*,double*,double*);
void d_sr24_cp_fft(unsigned long,unsigned long,long int,double*,double*,double*,double*);
void d_sr24_cp_fft_v2(unsigned long,unsigned long,unsigned long,long int,double*,double*,double*,double*);
void fft_complexity();
void check();
void print_fft_result();
void test();
void clear_memory();

/* serial, mp */
unsigned long long int get_time_int(void);
void create_memory();
void lut();
void ic();
void mp_sr24_cp_fft(unsigned long,unsigned long,long int,mpfr_t*,mpfr_t*,mpfr_t*,mpfr_t*);
void mp_sr24_cp_fft_v2(unsigned long,unsigned long,unsigned long,long int,mpfr_t*,mpfr_t*,mpfr_t*,mpfr_t*);
void fft_complexity();
void check();
void print_fft_result_d();
void print_fft_result_q();
void test();
void clear_memory();

/* parallel, double */
unsigned long long int get_time_int(void);
void code_info();
void create_memory();
void lut();
void ic();
void d_sr24_cp_fft(unsigned long,unsigned long,long int,double*,double*,double*,double*);
void fft_complexity();
void check();
void print_fft_result();
void test();
void clear_memory();

/* parallel, mp */
unsigned long long int get_time_int(void);
void code_info();
void create_memory();
void lut();
void ic();
void mp_sr24_cp_fft(unsigned long,unsigned long,long int,mpfr_t*,mpfr_t*,mpfr_t*,mpfr_t*);
void fft_complexity();
void check();
void print_fft_result_d();
void print_fft_result_q();
void test();
void clear_memory();

/* ********************************************************** GLOSSARY **********************************************************

check: first orientation-check if overall code is okay at least to within 1e-10, by comparing to a directly computed dft
       -- this is a slow function and should be switched off (via macro CH) for large input sizes ; use test-function instead
       -- if turned on, this function has no influence on the time measurement of the fft-kernel

code_info: prints all features of the selected code on terminal

create_memory: dynamically allocating memory on the heap
               -- it's a slow function for large arrays and/or high precision

clear_memory: clears all dynamically allocated memory on the heap

d_sr24_fft: double-precision, recursive split-radix-2/4 fft-kernel

d_sr24_fft_v2: alternative version to d_sr24_fft, with twiddle-index look-up through pass-by-shift

d_sr24_cp_fft: double-precision, recursive split-radix-2/4 conjugate-pair fft-kernel

d_sr24_cp_fft_v2: alternative version to d_sr24_cp_fft, with twiddle-index look-up through pass-by-shift

fft_complexity: counts the exact complexity (FLOPs) of the selected code
                -- since this feature is part of the fft-kernel, it has to be switched off (via macro CC) when speed is needed

get_time_int: time measurement function, with time in nanoseconds as integer

ic: initial-condition function, setting the fft input-array

lut: creating the look-up table for the twiddle factors w

mp_sr24_cp_fft: multiprecision, recursive split-radix-2/4 conjugate-pair fft-kernel

mp_sr24_cp_fft_v2: alternative version to mp_sr24_cp_fft, with twiddle-index look-up through pass-by-shift

print_fft_result: prints the fft-result of a double-precision selected code to file

print_fft_result_d: prints the fft-result of a multiprecision selected code to file in double-precision

print_fft_result_q: prints the fft-result of a multiprecision selected code to file in Q-precision, where Q<=P

test: prints the error analysis (via inverse fft) of the selected code to file *_test_error.dat:
      -- 1st column-pair: absolute error, which is a measure of the precision used
      -- 2nd column-pair: relative error, as a measure of accuracy achieved
                          NOTE: this measure is ill-defined when the true-value is close to zero
      -- 3rd column-pair: relative root-mean-square error (rRMSE), as a better/reliable measure of accuracy achieved

********************************************************************************************************************************* */
