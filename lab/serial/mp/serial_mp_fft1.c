
/*
  SPDX-FileCopyrightText: Copyright © 2021 Michael Frewer <frewer.science@gmail.com>
  SPDX-License-Identifier: Apache-2.0
*/

#ifndef CODE
#define CODE 1
#endif
#ifndef PAR
#define PAR 0
#endif

/* ****************************************************************************************************************************** */

static mpfr_t *LUTr,*LUTi;

static mpfr_t *inr,*ini,*outr,*outi;

static mpfr_t aux,tmp,
              m1r,m1i,m2r,m2i,
              a1r,a1i,a2r,a2i;

unsigned long cmul=0, cadd=0, // complexity counter for multiplication and addition
              cfc=0, crb=0;   // complexity counter for recursive-function calls and its recursion breaks

time_t seed;

/* ****************************************************************************************************************************** */

unsigned long long int get_time_int(void)
// time in nanoseconds as an integer
{
 struct timespec time;
 clock_gettime(CLOCK_MONOTONIC,&time);
 return time.tv_sec*1000000000+time.tv_nsec;
}

/* ****************************************************************************************************************************** */

void code_info()
{
 printf("\nSerial, multi-precision, CODE=%d, N=%lu, IN=%d\n\
 \bFFT-algorithm: DIT, recursive SR-2/4, conjugate-pair, 4m2a, 1-butterfly, recursion-break: n=4, twiddle-index: mul(*)\n",
 CODE,N,IN);
 printf("\nBinary precision set: P=%lu,\
    Decimal precision (digits) set: dP=%lu,\
    Internal decimal precision (digits) used: idP=%zu\n",
    (unsigned long)P,(unsigned long)floor(P*log10(2)),mpfr_get_str_ndigits(10,P));
}

/* ****************************************************************************************************************************** */

void create_memory()
// defining and initializing the globals
{
 LUTr=malloc(N/4*sizeof(mpfr_t)); LUTi=malloc(N/4*sizeof(mpfr_t));
 for(unsigned long k=0;k<N/4;k++) mpfr_inits2(P,LUTr[k],LUTi[k],NULL);

 mpfr_inits2(P,aux,tmp,m1r,m1i,m2r,m2i,a1r,a1i,a2r,a2i,NULL);

 inr =malloc(N*sizeof(mpfr_t)); ini =malloc(N*sizeof(mpfr_t));
 outr=malloc(N*sizeof(mpfr_t)); outi=malloc(N*sizeof(mpfr_t));
 for(unsigned long k=0;k<N;k++) mpfr_inits2(P,inr[k],ini[k],outr[k],outi[k],NULL);
}

/* ****************************************************************************************************************************** */

void lut()
{
 mpfr_const_pi(aux,R);
 mpfr_div_ui(aux,aux,N,R);

 for(unsigned long k=0;k<N/4;k++)
 {
  mpfr_mul_ui(tmp,aux,2*k,R);
  mpfr_cos(LUTr[k],tmp,R); mpfr_sin(LUTi[k],tmp,R); // NOTE: w=wr-Iwi; w*=wr+Iwi
 }
}

/* ****************************************************************************************************************************** */

void ic()
{
 switch(IN)
 {
  case 1: // complex random input, random generator: uniform distribution between 0 and 1
  {
   gmp_randstate_t rstate;
   gmp_randinit_mt(rstate);
   gmp_randseed_ui(rstate,seed);

   for(unsigned long k=0;k<N;k++)
   {
    mpfr_urandom(inr[k],rstate,R); //mpfr_out_str(stdout,10,0,inr[k],R); printf("\n");
    mpfr_urandom(ini[k],rstate,R); //mpfr_out_str(stdout,10,0,ini[k],R); printf("\n");
   }
  }
  break;

  case 2: // fixed complex input: inr[k]=1/(1+k^2); ini[k]=1/(1+k^4)
  for(unsigned long k=0;k<N;k++)
  {
   mpfr_set_ui(aux,k*k,R);
   mpfr_add_ui(inr[k],aux,1,R);
   mpfr_ui_div(inr[k],1,inr[k],R);

   mpfr_mul(aux,aux,aux,R);
   mpfr_add_ui(ini[k],aux,1,R);
   mpfr_ui_div(ini[k],1,ini[k],R);
  }
  break;
 }
}

/* ****************************************************************************************************************************** */

void mp_sr24_cp_fft(unsigned long n, unsigned long s, long int sn, mpfr_t *inr, mpfr_t *ini, mpfr_t *outr, mpfr_t *outi)
{
 if(CC){cmul+=0; cadd+=0; cfc+=1; crb+=0;}

 if(n==1)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  mpfr_set(outr[0],inr[s0],R);

  mpfr_set(outi[0],ini[s0],R);

  if(CC){cmul+=0; cadd+=0; cfc+=0; crb+=1;}
 }
 else if(n==2)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  mpfr_add(outr[0],inr[s0],inr[s],R);
  mpfr_sub(outr[1],inr[s0],inr[s],R);

  mpfr_add(outi[0],ini[s0],ini[s],R);
  mpfr_sub(outi[1],ini[s0],ini[s],R);

  if(CC){cmul+=0; cadd+=4; cfc+=0; crb+=1;}
 }
 else if(n==4)
 {
  unsigned long s0,s2=s<<1,s3=s2+s;
  if(sn<0) s0=N; else s0=0;

  mpfr_add(outr[0],inr[s0],inr[s2],R);
  mpfr_sub(outr[1],inr[s0],inr[s2],R);

  mpfr_add(outi[0],ini[s0],ini[s2],R);
  mpfr_sub(outi[1],ini[s0],ini[s2],R);

  mpfr_add(a1r,inr[s],inr[s3],R);
  mpfr_sub(a2r,inr[s],inr[s3],R);

  mpfr_add(a1i,ini[s],ini[s3],R);
  mpfr_sub(a2i,ini[s],ini[s3],R);

  mpfr_sub(outr[2],outr[0],a1r,R);
  mpfr_add(outr[0],outr[0],a1r,R);
  mpfr_sub(outr[3],outr[1],a2i,R);
  mpfr_add(outr[1],outr[1],a2i,R);

  mpfr_sub(outi[2],outi[0],a1i,R);
  mpfr_add(outi[0],outi[0],a1i,R);
  mpfr_add(outi[3],outi[1],a2r,R);
  mpfr_sub(outi[1],outi[1],a2r,R);

  if(CC){cmul+=0; cadd+=16; cfc+=0; crb+=1;}
 }
 else
 {
  mp_sr24_cp_fft(n/2,2*s,sn,inr,ini,outr,outi);
  mp_sr24_cp_fft(n/4,4*s,+1*s+sn,inr+s,ini+s,outr+2*n/4,outi+2*n/4);
  mp_sr24_cp_fft(n/4,4*s,-1*s+sn,inr-s,ini-s,outr+3*n/4,outi+3*n/4);

  // only 1 general butterfly for all k
  unsigned long k,kn,ks=N/n;
  for(k=0;k<n/4;k++)
  {
   kn=k*ks;

   mpfr_mul(m1r,LUTr[kn],outr[k+2*n/4],R);
   mpfr_mul(m2r,LUTr[kn],outr[k+3*n/4],R);
   mpfr_mul(m1i,LUTr[kn],outi[k+2*n/4],R);
   mpfr_mul(m2i,LUTr[kn],outi[k+3*n/4],R);

   mpfr_mul(tmp,LUTi[kn],outi[k+2*n/4],R);
   mpfr_add(m1r,m1r,tmp,R);
   mpfr_mul(tmp,LUTi[kn],outi[k+3*n/4],R);
   mpfr_sub(m2r,m2r,tmp,R);
   mpfr_mul(tmp,LUTi[kn],outr[k+2*n/4],R);
   mpfr_sub(m1i,m1i,tmp,R);
   mpfr_mul(tmp,LUTi[kn],outr[k+3*n/4],R);
   mpfr_add(m2i,m2i,tmp,R);

   mpfr_add(a1r,m1r,m2r,R);
   mpfr_sub(a2r,m1r,m2r,R);

   mpfr_add(a1i,m1i,m2i,R);
   mpfr_sub(a2i,m1i,m2i,R);

   mpfr_sub(outr[k+2*n/4],outr[k+0*n/4],a1r,R);
   mpfr_add(outr[k+0*n/4],outr[k+0*n/4],a1r,R);
   mpfr_sub(outr[k+3*n/4],outr[k+1*n/4],a2i,R);
   mpfr_add(outr[k+1*n/4],outr[k+1*n/4],a2i,R);

   mpfr_sub(outi[k+2*n/4],outi[k+0*n/4],a1i,R);
   mpfr_add(outi[k+0*n/4],outi[k+0*n/4],a1i,R);
   mpfr_add(outi[k+3*n/4],outi[k+1*n/4],a2r,R);
   mpfr_sub(outi[k+1*n/4],outi[k+1*n/4],a2r,R);

   if(CC){cmul+=8; cadd+=16; cfc+=0; crb+=0;}
  }
 }
}

/* ****************************************************************************************************************************** */

void fft_complexity()
{
 if(CC)
 {
  unsigned long lgN=round(log2(N)),c1FLOPS,c2FLOPS;

  if(N>1)
  {
   c1FLOPS=round(4.*N*lgN-6.*N+8.);
   if((lgN%2)==0) c2FLOPS=round((34./9)*N*lgN-(124./27)*N-2.*lgN-(2./9)*lgN+(16./27)+8.);
   else c2FLOPS=round((34./9)*N*lgN-(124./27)*N-2.*lgN+(2./9)*lgN-(16./27)+8.);
  }
  else {c1FLOPS=0; c2FLOPS=0;}

  printf("\nFFT Computation Complexity: cmul=%lu, cadd=%lu, cfc=%lu, crb=%lu\
  \n\nExact total FLOPs performed: cmul + cadd =\x1b[1m %lu\x1b[0m\
  \nSplit-radix theoretical value:\
  \b\b\x1b[1m %lu\x1b[0m (scales as 4Nlog2(N) for N>>1: Yavne, 1968)\
  \nSplit-radix theoretical (current) lowest value:\
  \b\b\x1b[1m %lu\x1b[0m (scales as (34/9)Nlog2(N) for N>>1: van Buskirk, 2004; Johnson & Frigo, 2007)\n",
  cmul/M,cadd/M,cfc/M,crb/M,(cmul+cadd)/M,c1FLOPS,c2FLOPS);
 }
}

/* ****************************************************************************************************************************** */

void check()
// first orientation-check if overall code is okay at least to within 1e-10, by comparing to a directly computed dft (slow)
{
 if(CH)
 {
  unsigned long j,k,
                c_bad=0, c_ok=0;
  double sumr,sumi,u=2.*M_PI/N;;
  for(k=0;k<N;k++)
  {
   sumr=0.;sumi=0.;
   for(j=0;j<N;j++)
   {
    sumr+=mpfr_get_d(inr[j],R)*cos(u*j*k)+mpfr_get_d(ini[j],R)*sin(u*j*k);
    sumi+=mpfr_get_d(ini[j],R)*cos(u*j*k)-mpfr_get_d(inr[j],R)*sin(u*j*k);
   }
   if(fabs(mpfr_get_d(outr[k],R)-sumr)>(1e-10)) c_bad+=1; else c_ok+=1;
   if(fabs(mpfr_get_d(outi[k],R)-sumi)>(1e-10)) c_bad+=1; else c_ok+=1;
  }
  if(c_bad>0) printf("\n***___BAD___***: %lu\n",c_bad); else printf("\nOK: %lu\n",c_ok);
 }
}

/* ****************************************************************************************************************************** */

void print_fft_result_d()
// print fft-result in double precision
{
 if(PRA)
 {
  FILE *f=fopen("data_serial_mp_fft1_d.dat","w");
  for(unsigned long k=0;k<N;k++)
  {
   fprintf(f,"% -.16e  % -.16ei  |  % -.16e  % -.16ei\n",mpfr_get_d(inr[k],R),mpfr_get_d(ini[k],R),
                                                         mpfr_get_d(outr[k],R),mpfr_get_d(outi[k],R));
  }
  fclose(f);
 }
}

/* ****************************************************************************************************************************** */

void print_fft_result_q()
// print to arbitrary decimal precision (Q<=P)*log(10,2), in the format: %[width].[precision]Re
{
 if(PRB)
 {
  FILE *f=fopen("data_serial_mp_fft1_q.dat","w");
  for(unsigned long k=0;k<N;k++)
  {
   unsigned long dQ=floor(Q*log10(2)),wd=dQ+6,pr=dQ-1;

   mpfr_fprintf(f,"%*.*Re",wd,pr,inr[k]); fprintf(f, "  "); mpfr_fprintf(f,"%*.*Re",wd,pr,ini[k]);
   fprintf(f, "   |  ");
   mpfr_fprintf(f,"%*.*Re",wd,pr,outr[k]); fprintf(f, "  "); mpfr_fprintf(f,"%*.*Re",wd,pr,outi[k]);
   fprintf(f, "\n");
  }
  fclose(f);
 }
}

/* ****************************************************************************************************************************** */

void test()
{
 if(PRC)
 {
  FILE *f1=fopen("data_serial_mp_fft1_test.dat","w"),
       *f2=fopen("data_serial_mp_fft1_test_error.dat","w");

  unsigned long k,
                dQ=floor(Q*log10(2)),wd=dQ+6,pr=dQ-1;

  mpfr_t *aEr=malloc(N*sizeof(mpfr_t)),*aEi=malloc(N*sizeof(mpfr_t)),
         *rEr=malloc(N*sizeof(mpfr_t)),*rEi=malloc(N*sizeof(mpfr_t));
  for(k=0;k<N;k++) mpfr_inits2(P,aEr[k],aEi[k],rEr[k],rEi[k],NULL);

  mpfr_t rRMSEr,rRMSEi,
         sumr,sumi,sumdr,sumdi;
  mpfr_inits2(P,rRMSEr,rRMSEi,sumr,sumi,sumdr,sumdi,NULL);
  mpfr_set_ui(sumr,0,R); mpfr_set_ui(sumi,0,R);
  mpfr_set_ui(sumdr,0,R); mpfr_set_ui(sumdi,0,R);

  // intializing ifft by swapping inr and ini with the current out-values
  for(unsigned long k=0;k<N;k++)
  {
   mpfr_set(inr[k],outi[k],R); mpfr_set(ini[k],outr[k],R);
  }

  // run ifft (up to trivial re-ordering in the real & imaginary part and normalization in the resulting out-values)
  FFT;
  // restoring the orignal in-values
  ic();

  // print (in Q<=P precision) new out-values, normalized by N, which then, when swapped, should equal original in-values
  for(k=0;k<N;k++)
  {
   mpfr_div_ui(outr[k],outr[k],N,R); mpfr_div_ui(outi[k],outi[k],N,R);

   mpfr_fprintf(f1,"%*.*Re",wd,pr,inr[k]); fprintf(f1, "  "); mpfr_fprintf(f1,"%*.*Re",wd,pr,ini[k]);
   fprintf(f1, "   |  ");
   mpfr_fprintf(f1,"%*.*Re",wd,pr,outi[k]); fprintf(f1, "  "); mpfr_fprintf(f1,"%*.*Re",wd,pr,outr[k]);
   fprintf(f1, "\n");

   // get the data for the error measures
   mpfr_sub(aEr[k],inr[k],outi[k],R); mpfr_sub(aEi[k],ini[k],outr[k],R);
   mpfr_div(rEr[k],aEr[k],inr[k],R); mpfr_div(rEi[k],aEi[k],ini[k],R);

   mpfr_mul(aux,inr[k],inr[k],R); mpfr_add(sumr,sumr,aux,R);
   mpfr_mul(aux,ini[k],ini[k],R); mpfr_add(sumi,sumi,aux,R);

   mpfr_mul(aux,aEr[k],aEr[k],R); mpfr_add(sumdr,sumdr,aux,R);
   mpfr_mul(aux,aEi[k],aEi[k],R); mpfr_add(sumdi,sumdi,aux,R);

   mpfr_div(aux,sumdr,sumr,R); mpfr_sqrt(rRMSEr,aux,R);
   mpfr_div(aux,sumdi,sumi,R); mpfr_sqrt(rRMSEi,aux,R);
  }

  // print the error measures between the true (in) and computed (out) values:
  // 1st column-pair: absolute error, which is a measure of the precision used
  // 2nd column-pair: relative error, as a measure of accuracy achieved
  //                  --- NOTE: this measure is ill-defined when the true-value is close to zero
  // 3rd column-pair: relative root-mean-square error (rRMSE), as a better/reliable measure of accuracy achieved
  wd=dQ+6+floor(log10(P*log10(2)))+1; // digit-length of exponent; counting digits of a number x: floor(log10(abs(x)))+1
  for(k=0;k<N;k++)
  {
   mpfr_fprintf(f2,"%*.*Re",wd,pr,aEr[k]); fprintf(f2, "  "); mpfr_fprintf(f2,"%*.*Re",wd,pr,aEi[k]);
   fprintf(f2, "   |  ");
   mpfr_fprintf(f2,"%*.*Re",wd,pr,rEr[k]); fprintf(f2, "  "); mpfr_fprintf(f2,"%*.*Re",wd,pr,rEi[k]);
   fprintf(f2, "   |  ");
   mpfr_fprintf(f2,"%*.*Re",wd,pr,rRMSEr); fprintf(f2, "  "); mpfr_fprintf(f2,"%*.*Re",wd,pr,rRMSEi);
   fprintf(f2, "\n");
  }

  dQ=floor(16*log10(2)),wd=dQ+6,pr=dQ-1;
  mpfr_printf("\nRelative root-mean-square error (rRMSE): %*.*Re (real part), %*.*Re (imaginary part)\n",
               wd,pr,rRMSEr,wd,pr,rRMSEi);

  free(aEr); free(aEi); free(rEr); free(rEi);
  mpfr_clears(rRMSEr,rRMSEi,sumr,sumi,sumdr,sumdi,NULL);

  fclose(f1);
  fclose(f2);
 }
}

/* ****************************************************************************************************************************** */

void clear_memory()
{
 mpfr_clears(aux,tmp,m1r,m1i,m2r,m2i,a1r,a1i,a2r,a2i,NULL);
 free(inr); free(ini); free(outr); free(outi);
 free(LUTr); free(LUTi);
}

/* ****************************************************************************************************************************** */
