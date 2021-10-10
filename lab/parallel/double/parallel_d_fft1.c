
/*
  SPDX-FileCopyrightText: Copyright © 2021 Michael Frewer <frewer.science@gmail.com>
  SPDX-License-Identifier: Apache-2.0
*/

/* ****************************************************************************************************************************** */

static double *LUTr,*LUTi;

static double *inr,*ini,*outr,*outi;

static double w8,
              w16_1,w16_2,w16_3,w16_4;

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
 printf("\nParallel: Cores=4, double-precision, CODE=%d, N=%lu, IN=%d\n\
 \bFFT-algorithm: DIT, recursive SR-2/4, conjugate-pair, 4m2a, 1-butterfly, recursion-break: n=16,\
 openmp (omp for): n>128, twiddle-index: mul(*)\n",CODE,N,IN);
}

/* ****************************************************************************************************************************** */

void create_memory()
{
 LUTr=malloc(N/4*sizeof(double)); LUTi=malloc(N/4*sizeof(double));

 inr =malloc(N*sizeof(double)); ini =malloc(N*sizeof(double));
 outr=malloc(N*sizeof(double)); outi=malloc(N*sizeof(double));
}

/* ****************************************************************************************************************************** */

void lut()
{
 for(unsigned long k=0;k<N/4;k++)
 {
  LUTr[k]=cos(2.*M_PI*k/N); LUTi[k]=sin(2.*M_PI*k/N); // NOTE: w=wr-Iwi; w*=wr+Iwi
 }

 w8=cos(2.*M_PI*1/8);
 w16_1=cos(2.*M_PI*1/16); w16_2=cos(2.*M_PI*1/16)-cos(2.*M_PI*3/16);
 w16_3=cos(2.*M_PI*3/16); w16_4=cos(2.*M_PI*1/16)+cos(2.*M_PI*3/16);
}

/* ****************************************************************************************************************************** */

void ic()
{
 switch(IN)
 {
  case 1: // complex random input between 0 and 1
  {
   srand(seed);
   for(unsigned long k=0;k<N;k++)
   {
    inr[k]=(double)rand()/RAND_MAX; //printf("inr[%lu]=%.15lf\n",k,inr[k]);
    ini[k]=(double)rand()/RAND_MAX; //printf("ini[%lu]=%.15lf\n",k,ini[k]);
   }
  }
  break;

  case 2: // fixed complex input: inr[k]=1/(1+k^2); ini[k]=1/(1+k^4)
  for(unsigned long k=0;k<N;k++)
  {
   double k2=k*k;
   inr[k]=1./(1+k2); ini[k]=1./(1+k2*k2);
  }
  break;
 }
}

/* ****************************************************************************************************************************** */

void d_sr24_cp_fft(unsigned long n, unsigned long s, long int sn, double *inr, double *ini, double *outr, double *outi)
{
 if(CC){cmul+=0; cadd+=0; cfc+=1; crb+=0;}

 if(n==1)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  outr[0]=inr[s0];

  outi[0]=ini[s0];

  if(CC){cmul+=0; cadd+=0; cfc+=0; crb+=1;}
 }
 else if(n==2)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  outr[0]=inr[s0]+inr[1*s];
  outr[1]=inr[s0]-inr[1*s];

  outi[0]=ini[s0]+ini[1*s];
  outi[1]=ini[s0]-ini[1*s];

  if(CC){cmul+=0; cadd+=4; cfc+=0; crb+=1;}
 }
 else if(n==4)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  double a1r=inr[s0]+inr[2*s],
         a2r=inr[s0]-inr[2*s],
         a3r=inr[1*s]+inr[3*s],
         a4r=inr[1*s]-inr[3*s],

         a1i=ini[s0]+ini[2*s],
         a2i=ini[s0]-ini[2*s],
         a3i=ini[1*s]+ini[3*s],
         a4i=ini[1*s]-ini[3*s];

  outr[0]=a1r+a3r;
  outr[1]=a2r+a4i;
  outr[2]=a1r-a3r;
  outr[3]=a2r-a4i;

  outi[0]=a1i+a3i;
  outi[1]=a2i-a4r;
  outi[2]=a1i-a3i;
  outi[3]=a2i+a4r;

  if(CC){cmul+=0; cadd+=16; cfc+=0; crb+=1;}
 }
 else if(n==8)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  double a1r=inr[s0]+inr[4*s],
         a2r=inr[s0]-inr[4*s],
         a3r=inr[2*s]+inr[6*s],
         a4r=inr[2*s]-inr[6*s],
         a5r=inr[1*s]+inr[7*s],
         a6r=inr[1*s]-inr[7*s],
         a7r=inr[3*s]+inr[5*s],
         a8r=inr[3*s]-inr[5*s],

         a1i=ini[s0]+ini[4*s],
         a2i=ini[s0]-ini[4*s],
         a3i=ini[2*s]+ini[6*s],
         a4i=ini[2*s]-ini[6*s],
         a5i=ini[1*s]+ini[7*s],
         a6i=ini[1*s]-ini[7*s],
         a7i=ini[3*s]+ini[5*s],
         a8i=ini[3*s]-ini[5*s],

         m1r=w8*(a6r+a8r),
         m1i=w8*(a6i+a8i),
         m2r=w8*(a5r-a7r),
         m2i=w8*(a5i-a7i),

         b1r=a5r+a7r,
         b1i=a5i+a7i,
         b2r=a1r+a3r,
         b2i=a1i+a3i,
         b3r=a1r-a3r,
         b3i=a1i-a3i,
         b4r=a6r-a8r,
         b4i=a6i-a8i,
         b5r=a2r+m2r,
         b5i=a2i+m2i,
         b6r=a2r-m2r,
         b6i=a2i-m2i,
         b7r=a4r+m1r,
         b7i=a4i+m1i,
         b8r=a4r-m1r,
         b8i=a4i-m1i;

  outr[0]=b2r+b1r;
  outr[1]=b5r+b7i;
  outr[2]=b3r+b4i;
  outr[3]=b6r-b8i;
  outr[4]=b2r-b1r;
  outr[5]=b6r+b8i;
  outr[6]=b3r-b4i;
  outr[7]=b5r-b7i;

  outi[0]=b2i+b1i;
  outi[1]=b5i-b7r;
  outi[2]=b3i-b4r;
  outi[3]=b6i+b8r;
  outi[4]=b2i-b1i;
  outi[5]=b6i-b8r;
  outi[6]=b3i+b4r;
  outi[7]=b5i+b7r;

  if(CC){cmul+=4; cadd+=52; cfc+=0; crb+=1;}
 }
 else if(n==16)
 {
  unsigned long s0;
  if(sn<0) s0=N; else s0=0;

  double a1r=inr[s0]+inr[8*s],
         a2r=inr[s0]-inr[8*s],
         a3r=inr[1*s]+inr[9*s],
         a4r=inr[1*s]-inr[9*s],
         a5r=inr[2*s]+inr[10*s],
         a6r=inr[2*s]-inr[10*s],
         a7r=inr[3*s]+inr[11*s],
         a8r=inr[3*s]-inr[11*s],
         a9r=inr[4*s]+inr[12*s],
         a10r=inr[4*s]-inr[12*s],
         a11r=inr[5*s]+inr[13*s],
         a12r=inr[5*s]-inr[13*s],
         a13r=inr[6*s]+inr[14*s],
         a14r=inr[6*s]-inr[14*s],
         a15r=inr[7*s]+inr[15*s],
         a16r=inr[7*s]-inr[15*s],

         a1i=ini[s0]+ini[8*s],
         a2i=ini[s0]-ini[8*s],
         a3i=ini[1*s]+ini[9*s],
         a4i=ini[1*s]-ini[9*s],
         a5i=ini[2*s]+ini[10*s],
         a6i=ini[2*s]-ini[10*s],
         a7i=ini[3*s]+ini[11*s],
         a8i=ini[3*s]-ini[11*s],
         a9i=ini[4*s]+ini[12*s],
         a10i=-ini[4*s]+ini[12*s],
         a11i=ini[5*s]+ini[13*s],
         a12i=ini[5*s]-ini[13*s],
         a13i=ini[6*s]+ini[14*s],
         a14i=ini[6*s]-ini[14*s],
         a15i=ini[7*s]+ini[15*s],
         a16i=ini[7*s]-ini[15*s],

         b1r=a1r-a9r,
         b2r=a1r+a9r,
         b3r=a5r-a13r,
         b4r=a5r+a13r,
         b5r=b2r+b4r,
         b6r=b2r-b4r,
         b7r=a3r-a11r,
         b8r=a3r+a11r,
         b9r=a7r-a15r,
         b10r=a7r+a15r,
         b11r=b8r+b10r,
         b12r=b8r-b10r,
         b13r=a4r+a16r,
         b14r=a4r-a16r,
         b15r=a8r+a12r,
         b16r=a12r-a8r,

         b1i=a1i-a9i,
         b2i=a1i+a9i,
         b3i=-a5i+a13i,
         b4i=a5i+a13i,
         b5i=b2i+b4i,
         b6i=b2i-b4i,
         b7i=a3i-a11i,
         b8i=a3i+a11i,
         b9i=a7i-a15i,
         b10i=a7i+a15i,
         b11i=b8i+b10i,
         b12i=-b8i+b10i,
         b13i=a4i+a16i,
         b14i=a4i-a16i,
         b15i=a8i+a12i,
         b16i=a12i-a8i,

         m1r=w8*(a6r-a14r),
         m2r=w8*(a6r+a14r),
         m3r=w8*(b7r-b9r),
         m4r=w8*(b7r+b9r),
         m5r=w16_1*(b13r+b15r),
         m6r=-w16_2*b13r,
         m7r=w16_4*b15r,
         m8r=w16_3*(b14r+b16r),
         m9r=w16_4*b14r,
         m10r=-w16_2*b16r,

         m1i=w8*(a6i-a14i),
         m2i=-w8*(a6i+a14i),
         m3i=w8*(b7i-b9i),
         m4i=-w8*(b7i+b9i),
         m5i=-w16_1*(b13i+b15i),
         m6i=w16_2*b13i,
         m7i=-w16_4*b15i,
         m8i=w16_3*(b14i+b16i),
         m9i=w16_4*b14i,
         m10i=-w16_2*b16i,

         c1r=b1r+m3r,
         c2r=b1r-m3r,
         c3r=b3r+m4r,
         c4r=m4r-b3r,
         c5r=a2r+m1r,
         c6r=a2r-m1r,
         c7r=m9r-m8r,
         c8r=m10r-m8r,
         c9r=c5r+c7r,
         c10r=c5r-c7r,
         c11r=c6r+c8r,
         c12r=c6r-c8r,
         c13r=a10r+m2r,
         c14r=a10r-m2r,
         c15r=m5r+m6r,
         c16r=m5r-m7r,
         c17r=c13r+c15r,
         c18r=c13r-c15r,
         c19r=c14r+c16r,
         c20r=c14r-c16r,

         c1i=b1i+m3i,
         c2i=b1i-m3i,
         c3i=b3i+m4i,
         c4i=m4i-b3i,
         c5i=a2i+m1i,
         c6i=a2i-m1i,
         c7i=m9i-m8i,
         c8i=m10i-m8i,
         c9i=c5i+c7i,
         c10i=c5i-c7i,
         c11i=c6i+c8i,
         c12i=c6i-c8i,
         c13i=a10i+m2i,
         c14i=a10i-m2i,
         c15i=m5i+m6i,
         c16i=m5i-m7i,
         c17i=c13i+c15i,
         c18i=c13i-c15i,
         c19i=c14i+c16i,
         c20i=c14i-c16i;

  outr[0]=b5r+b11r;
  outr[1]=c9r-c17i;
  outr[2]=c1r-c3i;
  outr[3]=c12r+c20i;
  outr[4]=b6r-b12i;
  outr[5]=c11r-c19i;
  outr[6]=c2r-c4i;
  outr[7]=c10r+c18i;
  outr[8]=b5r-b11r;
  outr[9]=c10r-c18i;
  outr[10]=c2r+c4i;
  outr[11]=c11r+c19i;
  outr[12]=b6r+b12i;
  outr[13]=c12r-c20i;
  outr[14]=c1r+c3i;
  outr[15]=c9r+c17i;

  outi[0]=b5i+b11i;
  outi[1]=c9i-c17r;
  outi[2]=c1i-c3r;
  outi[3]=c12i+c20r;
  outi[4]=b6i-b12r;
  outi[5]=c11i-c19r;
  outi[6]=c2i-c4r;
  outi[7]=c10i+c18r;
  outi[8]=b5i-b11i;
  outi[9]=c10i-c18r;
  outi[10]=c2i+c4r;
  outi[11]=c11i+c19r;
  outi[12]=b6i+b12r;
  outi[13]=c12i-c20r;
  outi[14]=c1i+c3r;
  outi[15]=c9i+c17r;

  if(CC){cmul+=20; cadd+=148; cfc+=0; crb+=1;}
 }
 else
 {
  d_sr24_cp_fft(n/2,2*s,sn,inr,ini,outr,outi);
  d_sr24_cp_fft(n/4,4*s,+1*s+sn,inr+s,ini+s,outr+2*n/4,outi+2*n/4);
  d_sr24_cp_fft(n/4,4*s,-1*s+sn,inr-s,ini-s,outr+3*n/4,outi+3*n/4);

  // only 1 general butterfly for all k
  unsigned long k,kn,ks=N/n;
  #pragma omp parallel for firstprivate(ks) if(n>128) num_threads(4)
  for(k=0;k<n/4;k++)
  {
   kn=k*ks;

   double U0r=outr[k+0*n/4],
          U1r=outr[k+1*n/4],
          Zpr=outr[k+2*n/4],
          Zmr=outr[k+3*n/4],

          U0i=outi[k+0*n/4],
          U1i=outi[k+1*n/4],
          Zpi=outi[k+2*n/4],
          Zmi=outi[k+3*n/4],

          wr=LUTr[kn], wi=LUTi[kn],

          m1r=wr*Zpr+wi*Zpi, m1i=wr*Zpi-wi*Zpr,
          m2r=wr*Zmr-wi*Zmi, m2i=wr*Zmi+wi*Zmr,

          a1r=m1r+m2r, a1i=m1i+m2i,
          a2r=m1r-m2r, a2i=m1i-m2i;

   outr[k+0*n/4]=U0r+a1r;
   outr[k+2*n/4]=U0r-a1r;
   outr[k+1*n/4]=U1r+a2i;
   outr[k+3*n/4]=U1r-a2i;

   outi[k+0*n/4]=U0i+a1i;
   outi[k+2*n/4]=U0i-a1i;
   outi[k+1*n/4]=U1i-a2r;
   outi[k+3*n/4]=U1i+a2r;

   if(CC){
   #pragma omp critical
   {cmul+=8; cadd+=16; cfc+=0; crb+=0;}}
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
    sumr+=inr[j]*cos(u*j*k)+ini[j]*sin(u*j*k);
    sumi+=ini[j]*cos(u*j*k)-inr[j]*sin(u*j*k);
   }
   if(fabs(outr[k]-sumr)>(1e-10)) c_bad+=1; else c_ok+=1;
   if(fabs(outi[k]-sumi)>(1e-10)) c_bad+=1; else c_ok+=1;
  }
  if(c_bad>0) printf("\n***___BAD___***: %lu\n",c_bad); else printf("\nOK: %lu\n",c_ok);
 }
}

/* ****************************************************************************************************************************** */

void print_fft_result()
{
 if(PRA)
 {
  FILE *f=fopen("data_parallel_d_fft1.dat","w");

  for(unsigned long k=0;k<N;k++)
  fprintf(f,"% -.16e  % -.16ei  |  % -.16e  % -.16ei\n",inr[k],ini[k],outr[k],outi[k]);

  fclose(f);
 }
}

/* ****************************************************************************************************************************** */

void test()
{
 if(PRB)
 {
  FILE *f1=fopen("data_parallel_d_fft1_test.dat","w"),
       *f2=fopen("data_parallel_d_fft1_test_error.dat","w");

  unsigned long k;
  double sumr=0.0, sumi=0.0,
         sumdr=0.0,sumdi=0.0;

  // intializing ifft by swapping inr and ini with the current out-values
  for(k=0;k<N;k++)
  {
   inr[k]=outi[k]; ini[k]=outr[k];
  }

  // run ifft (up to trivial re-ordering in the real & imaginary part and normalization in the resulting out-values)
  FFT;
  // restoring the orignal in-values
  ic();

  // print the new out-values, normalized by N, which then, when swapped, should equal the original in-values
  for(k=0;k<N;k++)
  {
   fprintf(f1,"% -.16e  % -.16ei  |  % -.16e  % -.16ei\n",inr[k],ini[k],outi[k]/N,outr[k]/N);

   // get the data for the relative root-mean-square-error (rRMSE)
   sumr+=inr[k]*inr[k]; sumi+=ini[k]*ini[k];
   sumdr+=(inr[k]-outi[k]/N)*(inr[k]-outi[k]/N); sumdi+=(ini[k]-outr[k]/N)*(ini[k]-outr[k]/N);
  }

  // print the error measures between the true (in) and computed (out) values:
  // 1st column-pair: absolute error, which is a measure of the precision used
  // 2nd column-pair: relative error, as a measure of accuracy achieved
  //                  --- NOTE: this measure is ill-defined when the true-value is close to zero
  // 3rd column-pair: relative root-mean-square error (rRMSE), as a better/reliable measure of accuracy achieved
  for(k=0;k<N;k++)
  {
   fprintf(f2,"% -.16e  % -.16ei  |  % -.16e  % -.16ei  |  % -.16e  % -.16ei\n",
               inr[k]-outi[k]/N,ini[k]-outr[k]/N,
               (inr[k]-outi[k]/N)/inr[k],(ini[k]-outr[k]/N)/ini[k],
               sqrt(sumdr)/sqrt(sumr),sqrt(sumdi)/sqrt(sumi));
  }
  printf("\nRelative root-mean-square error (rRMSE): % -.5e (real part), % -.5e (imaginary part)\n",
          sqrt(sumdr)/sqrt(sumr),sqrt(sumdi)/sqrt(sumi));

  fclose(f1);
  fclose(f2);
 }
}

/* ****************************************************************************************************************************** */

void clear_memory()
{
 free(inr); free(ini); free(outr); free(outi);
 free(LUTr); free(LUTi);
}

/* ****************************************************************************************************************************** */
