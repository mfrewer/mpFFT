
/*
  SPDX-FileCopyrightText: Copyright © 2021 Michael Frewer <frewer.science@gmail.com>
  SPDX-License-Identifier: Apache-2.0
*/

#include "../../../mpFFT.h"

/* ********************************************************** Set Values ******************************************************** */

// set fft-algorithm (see OPTIONS.md for possible choices)
#define CODE 2
//#define CODE 'F'

// set vector-length as power of 2: N=2^n, n>=0
#define N (1UL<<10)

// set number of FFT repetitions M>=1 for time statistics
#define M 1
//#define M 1000

// set input condition (see OPTIONS.md for possible choices)
#define IN 2

// set complexity count for fft-algorithm: 1 (on), 0 (off)
// turn off, when speed is needed (depending on compiler optimizations)
#define CC 1

// set first-orientation check by comparing to DFT (slow): 1 (on), 0 (off)
#define CH 1

// set file printing options: 1 (on), 0 (off)
#define PRA 1 // print fft-result
#define PRB 1 // print fft-test via ifft

/* ****************************************************************************************************************************** */

// read-in code
#if CODE=='F'
#undef CODE
#define CODE 2
#define FFT d_sr24_cp_fft(N,1,0,inr,ini,outr,outi);
#include "parallel_d_fft2.c"
#elif CODE==1
#define FFT d_sr24_cp_fft(N,1,0,inr,ini,outr,outi);
#include "parallel_d_fft1.c"
#elif CODE==2
#define FFT d_sr24_cp_fft(N,1,0,inr,ini,outr,outi);
#include "parallel_d_fft2.c"
#endif

int main()
{
 // disable the buffer for immediate print-out on screen
 setbuf(stdout,NULL);

 // declare time measurement variables
 unsigned long long int start_time,end_time,
                        start_time_init,end_time_init,
                        start_time_fft,end_time_fft,
                        elapsed;

 // set clock for overall time
 start_time=get_time_int();

 // set seed for random number generator
 seed=time(NULL);

 // print code info
 code_info();

 // initialize fft
 start_time_init=get_time_int();
 create_memory();
 lut();
 ic();
 end_time_init=get_time_int(); elapsed=end_time_init-start_time_init;
 printf("\nInitialization time: \t%.2f ms  ( %8.2f \xc2\xb5s)\n",(1e-6)*elapsed,(1e-3)*elapsed);

 // run fft
 start_time_fft=get_time_int();
 for(unsigned long m=1;m<=M;m++) FFT;
 end_time_fft=get_time_int(); elapsed=end_time_fft-start_time_fft;
 printf("\nAverage time per FFT: \t%.2f ms  ( %8.2f \xc2\xb5s),\
 Performance (CTGs): %.2f GFLOP/s\n",(1e-6)*elapsed/M,(1e-3)*elapsed/M,5*N*log2(N)/(elapsed/M));

 // show the performed fft computation complexity
 fft_complexity();

 // aftermath
 check();
 print_fft_result();
 test();

 // exit
 clear_memory();

 end_time=get_time_int(); elapsed=end_time-start_time;
 printf("\nOverall time: \t%.2f s\n\n",(1e-9)*elapsed);

 return 0;
}
