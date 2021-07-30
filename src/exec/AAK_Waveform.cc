#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "Globals.h"
#include "GKTrajFast.h"
#include "KSParMap.h"
#include "KSTools.h"
#include "AAK.h"

#include "NRIEKG.h"
#include "NRUtil.h"



//------------ QAAKMOD START -----------------------

// inverse of A
gsl_matrix* GetInverse(gsl_matrix *A)
{
  int n = A->size1;
  int sign = 0;
  gsl_matrix *inverse;

  gsl_matrix *tmpA = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(tmpA, A);
  gsl_permutation *p = gsl_permutation_alloc(n);

  gsl_linalg_LU_decomp(tmpA, p, &sign);
  inverse = gsl_matrix_alloc(n, n);
  gsl_linalg_LU_invert(tmpA, p, inverse);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
  return inverse;
}


double integral_R_Theta(double(* fun)(double,double,double,double,double,double),double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++)
    s+=1./sqrt( (*fun)(EE,LL,QQ,a,m,start+i*h) );
  return(s*h);
}

double integral_J_r(double(* fun)(double,double,double,double,double,double),double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h,Delta_r;
  int i;
  s=0.;
  h=(b-start)/N;
  double M=1.;
  for(i=3;i<N-2;i++){
    Delta_r = (start+i*h)*(start+i*h)-2*M*(start+i*h)+a*a;
    s+=sqrt( (*fun)(EE,LL,QQ,a,m,start+i*h) )/Delta_r;
  }
  return(s*h);
}

double integral_J_theta(double(* fun)(double,double,double,double,double,double),double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  double M=1.;
  for(i=3;i<N-2;i++)
    s+=sqrt( (*fun)(EE,LL,QQ,a,m,start+i*h) );
  return(s*h);
}




double inInterp_q(const double t, double timearr[], double dataarr[],
		       const int N)
{
  int ihi = 1;

  while (timearr[ihi] < t && ihi < N)
    ihi++;

  while (timearr[ihi - 1] > t && ihi > 1)
    ihi--;

  const double frac = (timearr[ihi] - t)/(timearr[ihi] - timearr[ihi - 1]);
  return(frac*dataarr[ihi - 1] + (1. - frac)*dataarr[ihi]);
}

double caculate_integration_2D(double(*f)(double,double,double,double,double,double,double,double),double EE, double LL, double QQ, double a, double m, double B2, double rarray[], double lambda_r[], double thetaarray[], double lambda_theta[], double ax,double bx,double ay,double by,int precision){

  double lenx=bx-ax,leny=by-ay;
        
  int Nx=precision;
  int Ny=precision;
  double interval_x=lenx/precision;
  double interval_y=leny/precision;
  double result=0.;
  for(int i=5;i<Nx-4;i++){
    for(int j=5;j<Ny-4;j++){
      result=result+f(EE,LL,QQ,a,m,B2,inInterp_q(ax+i*interval_x, lambda_r, rarray, precision),inInterp_q(ay+j*interval_y, lambda_theta, thetaarray, precision))*(interval_x*interval_y);
    }
  }
  return result;
}





double Csc(double theta)
{
  return (1./sin(theta));
}

double R_motion(double EE, double LL, double QQ, double a, double m, double r )
{
  return (-((pow(a,2.) - 2.*r + pow(r,2.))*(pow(-(a*EE) + LL,2.) + QQ + pow(m,2.)*pow(r,2.))) + 
  pow(-(a*LL) + EE*(pow(a,2.) + pow(r,2.)),2.) );
}

double Theta_motion(double EE, double LL, double QQ, double a, double m, double theta )
{
  return ( QQ - pow(cos(theta),2.)*(pow(a,2.)*(-pow(EE,2.) + pow(m,2.)) + pow(LL,2.)*pow(Csc(theta),2.)) );
}

double Phi_motion(double EE, double LL, double QQ, double a, double m, double r, double theta )
{
  return ( -((pow(a,2.)*LL)/(pow(a,2.) - 2.*r + pow(r,2.))) + 
  a*EE*(-1. + (pow(a,2.) + pow(r,2.))/(pow(a,2.) - 2.*r + pow(r,2.))) + LL*pow(Csc(theta),2.)
  );
}

double T_motion(double EE, double LL, double QQ, double a, double m, double B2, double r, double theta )
{
  return ( a*LL*(1. - (pow(a,2.) + pow(r,2.))/(pow(a,2.) - 2.*r + pow(r,2.))) + 
  EE*(pow(pow(a,2.) + pow(r,2.),2.)/(pow(a,2.) - 2.*r + pow(r,2.)) - 
  pow(a,2.)*pow(sin(theta),2.)) );
}

double H1_multiply_T(double EE, double LL, double QQ, double a, double m, double B2, double r, double theta )
{
  double M = 1.;
  double r2 = pow(r,2.);
  double a2 = pow(a,2.);
  double sqrt5M_PI = sqrt(5/M_PI);
  double M2 = pow(M,2.);
  double powMr2 = pow(-M + r,2.);
  double costheta = cos(theta);
  double powcostheta2 = pow(costheta,2.);
  double a22rr2 = (a2 - 2.*r + r2);
  double powM4 = pow(M,4.);
  double powM3 = pow(M,3.);
  return ( T_motion( EE, LL,  QQ,  a,  m,B2, r, theta )*(-((2.*(-(a22rr2*(pow(-(a*EE) + LL,2.) + QQ + M2*r2)) + 
           pow(-(a*LL) + EE*(a2 + r2),2.))*
         (r2 + a2*powcostheta2)*
         (-(B2*powM3*sqrt5M_PI*(-1. + 
                 (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
                  (-2.*M*r + r2 + (a2 + M2)*powcostheta2)))/
            (4.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)) + 
           B2*sqrt5M_PI*(-1. + (sqrt(powMr2 + a2*powcostheta2)*
                 (3.*powM4 - 5.*M2*powMr2 + 2.*pow(-M + r,4.) + 
                   (-3.*powM4 + 5.*M2*powMr2 + 
                      a2*(-5.*M2 + 4.*powMr2))*powcostheta2 + 
                   a2*(2.*a2 + 5.*M2)*pow(costheta,4.)))/
               (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,2.5)))))
        /a22rr2 + 
      2.*(r2 + a2*powcostheta2)*
       (-(B2*powM3*sqrt5M_PI*(-1. + 
               (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
                (-2.*M*r + r2 + (a2 + M2)*powcostheta2)))/
          (4.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)) + 
         B2*sqrt5M_PI*(-1. + (sqrt(powMr2 + a2*powcostheta2)*
               (3.*powM4 - 5.*M2*powMr2 + 2.*pow(-M + r,4.) + 
                 (-3.*powM4 + 5.*M2*powMr2 + 
                    a2*(-5.*M2 + 4.*powMr2))*powcostheta2 + 
                 a2*(2.*a2 + 5.*M2)*pow(costheta,4.)))/
             (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,2.5))))*
       (QQ - powcostheta2*(a2*(-pow(EE,2.) + M2) + 
            pow(LL,2.)*pow(Csc(theta),2.))) + 
      a22rr2*pow(-((a2*LL)/a22rr2) + 
         a*EE*(-1. + (a2 + r2)/a22rr2) + 
         LL*pow(Csc(theta),2.),2.)*pow(sin(theta),2.)*
       (-(B2*powM3*sqrt5M_PI*(-1. + 
               (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
                (-2.*M*r + r2 + (a2 + M2)*powcostheta2)))/
          (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)*
            (1. - (2.*M*r)/(r2 + a2*powcostheta2))) + 
         (8.*a2*M2*r2*
            (-(B2*powM3*sqrt5M_PI*
                  (-1. + (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
                     (-2.*M*r + r2 + (a2 + M2)*powcostheta2)))/
               (4.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)) + 
              B2*sqrt5M_PI*(-1. + (sqrt(powMr2 + a2*powcostheta2)*
                    (3.*powM4 - 5.*M2*powMr2 + 2.*pow(-M + r,4.) + 
                      (-3.*powM4 + 5.*M2*powMr2 + 
                         a2*(-5.*M2 + 4.*powMr2))*powcostheta2 + 
                      a2*(2.*a2 + 5.*M2)*pow(costheta,4.)))/
                  (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,2.5))
                 ))*pow(sin(theta),2.))/
          (a22rr2*(r2 + a2*powcostheta2)*
            (-2.*M*r + r2 + a2*powcostheta2))) + 
      (4.*a*M*r*(-(B2*powM3*sqrt5M_PI*
               (-1. + (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
                  (-2.*M*r + r2 + (a2 + M2)*powcostheta2)))/
            (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)) + 
           B2*sqrt5M_PI*(-1. + (sqrt(powMr2 + a2*powcostheta2)*
                 (3.*powM4 - 5.*M2*powMr2 + 2.*pow(-M + r,4.) + 
                   (-3.*powM4 + 5.*M2*powMr2 + 
                      a2*(-5.*M2 + 4.*powMr2))*powcostheta2 + 
                   a2*(2.*a2 + 5.*M2)*pow(costheta,4.)))/
               (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,2.5))))*
         (-((a2*LL)/a22rr2) + 
           a*EE*(-1. + (a2 + r2)/a22rr2) + 
           LL*pow(Csc(theta),2.))*pow(sin(theta),2.)*
         (a*LL*(1. - (a2 + r2)/a22rr2) + 
           EE*(pow(a2 + r2,2.)/a22rr2 - 
              a2*pow(sin(theta),2.))))/(r2 + a2*powcostheta2) - 
      (B2*powM3*sqrt5M_PI*(1. - (2.*M*r)/(r2 + a2*powcostheta2))*
         (-1. + (3.*powcostheta2*(powMr2 + a2*powcostheta2))/
            (-2.*M*r + r2 + (a2 + M2)*powcostheta2))*
         pow(a*LL*(1. - (a2 + r2)/a22rr2) + 
           EE*(pow(a2 + r2,2.)/a22rr2 - 
              a2*pow(sin(theta),2.)),2.))/
       (2.*pow(-2.*M*r + r2 + (a2 + M2)*powcostheta2,1.5)))/
   (2.*pow(r2 + a2*powcostheta2,2.))) );
}


double average_H1(double EE, double LL, double QQ, double e, double p, double iota, double a, double m, double B2){
  //integral(double(* fun)(double),double a,double b, int N)
  double rp = p/(1.+e);
  double ra = p/(1.-e);
  int prec = 100;

  //interpolation of r,theta
  double *lambda_r_pre,*lambda_theta_pre,*rarray_pre,*thetaarray_pre;

  lambda_r_pre=(double*)malloc(prec*sizeof(double));
  lambda_theta_pre=(double*)malloc(prec*sizeof(double));
  rarray_pre=(double*)malloc(prec*sizeof(double));
  thetaarray_pre=(double*)malloc(prec*sizeof(double));


  lambda_r_pre[0]=0.;  //care
  lambda_theta_pre[0]=0.;

  for(int i=0;i<prec;i++) rarray_pre[i] = rp + i*(ra-rp)/(prec-1);
  for(int i=0;i<prec;i++) thetaarray_pre[i] = M_PI/2.-iota + i*(iota)/(prec-1);

  double r_step = (ra-rp)/prec;
  double theta_step = (iota)/prec;

  for(int i=1;i<3;i++){
  lambda_r_pre[i] = lambda_r_pre[i-1] + r_step/sqrt(R_motion(EE,LL,QQ,a,m,rp+3*r_step));
  lambda_theta_pre[i] = lambda_theta_pre[i-1] + theta_step/sqrt(Theta_motion(EE,LL,QQ,a,m,M_PI/2.-iota+3*theta_step));
  }


  for(int i=3;i<prec-2;i++) lambda_r_pre[i] = lambda_r_pre[i-1] + r_step/sqrt(R_motion(EE,LL,QQ,a,m,rp+i*r_step));

  for(int i=3;i<prec-2;i++) lambda_theta_pre[i] = lambda_theta_pre[i-1] + theta_step/sqrt(Theta_motion(EE,LL,QQ,a,m,M_PI/2.-iota+i*theta_step));

  for(int i=prec-2;i<prec;i++){
  lambda_r_pre[i] = lambda_r_pre[i-1] + r_step/sqrt(R_motion(EE,LL,QQ,a,m,rp+(prec-3)*r_step));
  lambda_theta_pre[i] = lambda_theta_pre[i-1] + theta_step/sqrt(Theta_motion(EE,LL,QQ,a,m,M_PI/2.-iota+(prec-3)*theta_step));
  }

  double Gamma_t = caculate_integration_2D(T_motion, EE,  LL,  QQ,  a, m,B2, rarray_pre,  lambda_r_pre, thetaarray_pre,  lambda_theta_pre,  0., lambda_r_pre[prec-1], 0., lambda_theta_pre[prec-1], prec);

  double H1_average = 1./Gamma_t*caculate_integration_2D(H1_multiply_T, EE,  LL,  QQ,  a, m,B2, rarray_pre, lambda_r_pre,  thetaarray_pre, lambda_theta_pre,  0., lambda_r_pre[prec-1], 0., lambda_theta_pre[prec-1], prec);

  free(rarray_pre);
  free(thetaarray_pre);
  free(lambda_r_pre);
  free(lambda_theta_pre);

  return (H1_average);
}

double J_r(double EE, double LL, double QQ, double e, double p, double iota, double a, double m){
  int prec = 500;
  double rp = p/(1.+e);
  double ra = p/(1.-e);
  return(integral_J_r(R_motion,EE,LL,QQ,a,m,rp,ra,prec)/M_PI);
}

double J_theta(double EE, double LL, double QQ, double e, double p, double iota, double a, double m){
  int prec = 500;
  return(integral_J_theta(Theta_motion,EE,LL,QQ,a,m,M_PI/2.-iota,M_PI/2.,prec)*2./M_PI);
}

double J_phi(double EE, double LL, double QQ, double e, double p, double iota, double a, double m){
  return(LL);
}

double J_t(double EE, double LL, double QQ, double e, double p, double iota, double a, double m){
  return(-EE);
}

double dJr_dHH(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++)
      s+=(start+i*h)*(start+i*h)/sqrt( R_motion(EE,LL,QQ,a,m,start+i*h) );
  return(s*h/2./M_PI*2.);
}

double dJr_dQQ(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++)
      s+=1./sqrt( R_motion(EE,LL,QQ,a,m,start+i*h) );
  return(-s*h/4./M_PI*2.);
}

double dJr_dLL(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h,Delta_r,r,M = 1.;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++){
      r = start+i*h;
      Delta_r = (start+i*h)*(start+i*h)-2*M*(start+i*h)+a*a;
      s+=r*(r*LL-2*M*(LL-a*EE))/Delta_r/sqrt( R_motion(EE,LL,QQ,a,m,start+i*h) );
  }
  return(-s*h/2./M_PI*2.);
}

double dJr_dEE(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h,Delta_r,r,M = 1.;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++){
      r = start+i*h;
      Delta_r = (start+i*h)*(start+i*h)-2*M*(start+i*h)+a*a;
      s+=r*(r*EE*(r*r+a*a)-2*M*a*(LL-a*EE))/Delta_r/sqrt( R_motion(EE,LL,QQ,a,m,start+i*h) );
  }
  return(s*h/2./M_PI*2.);
}

double dJtheta_dHH(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++) 
    s+=a*a*cos(start+i*h)*cos(start+i*h)/sqrt( Theta_motion(EE,LL,QQ,a,m,start+i*h) );
  return(s*h/2./M_PI*4.);
}

double dJtheta_dQQ(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++) 
    s+=1/sqrt( Theta_motion(EE,LL,QQ,a,m,start+i*h) );
  return(s*h/4./M_PI*4.);
}

double dJtheta_dLL(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++) 
    s+=LL*cos(start+i*h)*cos(start+i*h)/sin(start+i*h)/sin(start+i*h)/sqrt( Theta_motion(EE,LL,QQ,a,m,start+i*h) );
  return(-s*h/2./M_PI*4.);
}

double dJtheta_dEE(double EE, double LL, double QQ, double a, double m, double start, double b, int N)      
{
  double s,h;
  int i;
  s=0.;
  h=(b-start)/N;
  for(i=3;i<N-2;i++) 
    s+=a*a*EE*cos(start+i*h)*cos(start+i*h)/sqrt( Theta_motion(EE,LL,QQ,a,m,start+i*h) );
  return(s*h/2./M_PI*4.);
}
//------------ QAAKMOD END -----------------------

//------------ Orbital Constants Part Start -----------------------
// This part is mainly refer to EMRI_Kludge_Suite functions to compute E, Lz and Carter constant. For example, RFunc_q is refer to RFunc. 
Real RFunc_q(const Real r, const Real a, const Real z,
		const Real E, const Real Lz, const Real Q)
{
  const Real tmp1 = E*(r*r + a*a) - a*Lz;
  const Real tmp2 = Kerr::Delta(r, a)*(r*r + SQR(Lz - a*E) + Q);
  const Real tmp = SQR(tmp1) - tmp2;

  return(tmp);
}

void func_and_jac_q(Real *x, int n, Real *fvec, Real **fjac, double spin, double cosiota, double p, double ecc)
{

  const Real e = x[1];
  const Real lz = x[2];
  const Real tansqiota=(1. - cosiota)*(1. + cosiota)/(cosiota*cosiota);
  const Real taniota=sqrt(tansqiota);
  const Real tisqr = taniota*taniota;
  const Real cisqr = cosiota*cosiota;
  const Real q = lz*lz*tisqr;
  const Real peri = p/(1. + ecc);
  const Real ap = p/(1. - ecc);

  fvec[1] = RFunc_q(peri, spin, 0, e, lz, q);
  fvec[2] = RFunc_q(ap, spin, 0, e, lz, q);

  const Real aelz = spin*e - lz;
  fjac[1][1] = 2.*peri*(2.*spin*aelz + peri*(spin*spin*e + peri*peri*e));
  fjac[1][2] = -2.*spin*spin*lz*tisqr + 2.*peri*(2.*(lz*tisqr - aelz)
					   - peri*lz/cisqr);
  fjac[2][1] = 2.*ap*(2.*spin*aelz + ap*(spin*spin*e + ap*ap*e));
  fjac[2][2] = -2.*spin*spin*lz*tisqr + 2.*ap*(2.*(lz*tisqr - aelz)
					 - ap*lz/cisqr);
}


void mnewt_q(int ntrial, Real x[], int n, Real tolx, Real tolf, double spin, double cosiota, double pp, double ecc)
{
  int k, i, *indx;
  Real errx, errf, d, *fvec, **fjac, *p;
  
  indx = ivector(1, n);
  p = Realvector(1, n);
  fvec = Realvector(1, n);
  fjac = Realmatrix(1, n, 1, n);
  for (k = 1; k <= ntrial; k++) {
    func_and_jac_q(x, n, fvec, fjac, spin, cosiota, pp, ecc);
    errf = 0.0;
    for (i = 1; i <= n; i++) errf += fabs(fvec[i]);
    if (errf <= tolf) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
    for (i = 1; i <= n; i++) p[i] = -fvec[i];
    ludcmp(fjac, n, indx, &d);
    lubksb(fjac, n, indx, p);
    errx = 0.0;
    for (i = 1; i <= n; i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
    }
    if (errx <= tolx) {
      free_Realmatrix(fjac, 1, n, 1, n);
      free_Realvector(fvec, 1, n);
      free_Realvector(p, 1, n);
      free_ivector(indx, 1, n);
      return;
    }
  }
  free_Realmatrix(fjac, 1, n, 1, n);
  free_Realvector(fvec, 1, n);
  free_Realvector(p, 1, n);
  free_ivector(indx, 1, n);
  return;
}

void Constants_q_once_main(double p_q,double e_q,double cosiota_q, double constants_q[], double spin) 
{

  const Real peri_q = p_q/(1. + e_q);
  const Real semimaj = peri_q/(1. - e_q);
  Real Epro = Kerr::Eeqpro(semimaj, spin);
  Real Eret = Kerr::Eeqret(semimaj, spin);
  if (isnan(Epro) || Epro > 1.-1.e-6) Epro = 1.-1.e-6;
  if (isnan(Eret) || Eret > 1.) Eret = 1.;
  const Real Eguess = 0.5*(cosiota_q + 1.)*Epro + 0.5*(1. - cosiota_q)*Eret;
  const Real Lzguess = cosiota_q*sqrt((1. - e_q)*(1. + e_q)/(2. - 2.*Eguess));
  const double tansqiota=(1. - cosiota_q)*(1. + cosiota_q)/(cosiota_q*cosiota_q);

  const Real tolx = 1.e-14;
  const Real tolf = 1.e-14;
  Real *x; x = Realvector(1, 2);

  x[1] = Eguess; x[2] = Lzguess;

  mnewt_q(200, x, 2, tolx, tolf,spin,cosiota_q,p_q,e_q);
  
  constants_q[1] = x[1]; constants_q[2] = x[2]; constants_q[3] = x[2]*x[2]*tansqiota;
  free_Realvector(x, 1, 2);
}
//------------ Orbital Constants End -----------------------
//------------ QAAKMOD START -----------------------
//------------ calculating dH1d(m,e,p,iota) ------------
void H1_derives (double EE, double LL, double QQ, double e, double p, double iota, double a, double m, double B2, const double d_constants_p[4], const double d_constants_e[4], const double d_constants_iota[4], double result[4], double ratio_diff){
  double fiducial_H1 = average_H1(EE, LL, QQ, e, p, iota, a, m, B2);

  result[0] = (average_H1(EE/m*(m + ratio_diff*m), LL/m*(m + ratio_diff*m), QQ/m/m*(m + ratio_diff*m)*(m + ratio_diff*m), e, p, iota, a, m + ratio_diff*m, B2) - fiducial_H1)/ratio_diff*m;

  result[1] = (average_H1(EE + d_constants_p[1], LL + d_constants_p[2], QQ + d_constants_p[3], e, p + ratio_diff, iota, a, m, B2) - fiducial_H1)/ratio_diff;

  result[2] = (average_H1(EE + d_constants_e[1], LL + d_constants_e[2], QQ + d_constants_e[3], e + ratio_diff, p, iota, a, m, B2) - fiducial_H1)/ratio_diff;

  result[3] = (average_H1(EE + d_constants_iota[1], LL + d_constants_iota[2], QQ + d_constants_iota[3], e, p, iota + ratio_diff, a, m, B2) - fiducial_H1)/ratio_diff;
}
//------------ calculating dJd(m,e,p,iota) ------------
void J_derives (double EE, double LL, double QQ, double e, double p, double iota, double a, double m, double B2, const double d_constants_p[4], const double d_constants_e[4], const double d_constants_iota[4], double result[4][4], double ratio_diff){
  double fiducial_J_r = J_r(EE, LL, QQ, e, p, iota, a, m);
  result[0][0] = (J_r(EE/m*(m + ratio_diff*m), LL/m*(m + ratio_diff*m), QQ/m/m*(m + ratio_diff*m)*(m + ratio_diff*m), e, p, iota, a, m + ratio_diff*m) - fiducial_J_r)/ratio_diff*m;
  result[0][1] = (J_r(EE + d_constants_p[1], LL + d_constants_p[2], QQ + d_constants_p[3], e, p + ratio_diff, iota, a, m) - fiducial_J_r)/ratio_diff;
  result[0][2] = (J_r(EE + d_constants_e[1], LL + d_constants_e[2], QQ + d_constants_e[3], e + ratio_diff, p, iota, a, m) - fiducial_J_r)/ratio_diff;
  result[0][3] = (J_r(EE + d_constants_iota[1], LL + d_constants_iota[2], QQ + d_constants_iota[3], e, p, iota + ratio_diff, a, m) - fiducial_J_r)/ratio_diff;

  double fiducial_J_theta = J_theta(EE, LL, QQ, e, p, iota, a, m);
  result[1][0] = (J_theta(EE/m*(m + ratio_diff*m), LL/m*(m + ratio_diff*m), QQ/m/m*(m + ratio_diff*m)*(m + ratio_diff*m), e, p, iota, a, m + ratio_diff*m) - fiducial_J_theta)/ratio_diff*m;
  result[1][1] = (J_theta(EE + d_constants_p[1], LL + d_constants_p[2], QQ + d_constants_p[3], e, p + ratio_diff, iota, a, m) - fiducial_J_theta)/ratio_diff;
  result[1][2] = (J_theta(EE + d_constants_e[1], LL + d_constants_e[2], QQ + d_constants_e[3], e + ratio_diff, p, iota, a, m) - fiducial_J_theta)/ratio_diff;
  result[1][3] = (J_theta(EE + d_constants_iota[1], LL + d_constants_iota[2], QQ + d_constants_iota[3], e, p, iota + ratio_diff, a, m) - fiducial_J_theta)/ratio_diff;

  double fiducial_J_phi = J_phi(EE, LL, QQ, e, p, iota, a, m);
  result[2][0] = (J_phi(EE/m*(m + ratio_diff*m), LL/m*(m + ratio_diff*m), QQ/m/m*(m + ratio_diff*m)*(m + ratio_diff*m), e, p, iota, a, m + ratio_diff*m) - fiducial_J_phi)/ratio_diff*m;
  result[2][1] = (J_phi(EE + d_constants_p[1], LL + d_constants_p[2], QQ + d_constants_p[3], e, p + ratio_diff, iota, a, m) - fiducial_J_phi)/ratio_diff;
  result[2][2] = (J_phi(EE + d_constants_e[1], LL + d_constants_e[2], QQ + d_constants_e[3], e + ratio_diff, p, iota, a, m) - fiducial_J_phi)/ratio_diff;
  result[2][3] = (J_phi(EE + d_constants_iota[1], LL + d_constants_iota[2], QQ + d_constants_iota[3], e, p, iota + ratio_diff, a, m) - fiducial_J_phi)/ratio_diff;
  double fiducial_J_t = J_t(EE, LL, QQ, e, p, iota, a, m);
  result[3][0] = (J_t(EE/m*(m + ratio_diff*m), LL/m*(m + ratio_diff*m), QQ/m/m*(m + ratio_diff*m)*(m + ratio_diff*m), e, p, iota, a, m + ratio_diff*m) - fiducial_J_t)/ratio_diff*m;
  result[3][1] = (J_t(EE + d_constants_p[1], LL + d_constants_p[2], QQ + d_constants_p[3], e, p + ratio_diff, iota, a, m) - fiducial_J_t)/ratio_diff;
  result[3][2] = (J_t(EE + d_constants_e[1], LL + d_constants_e[2], QQ + d_constants_e[3], e + ratio_diff, p, iota, a, m) - fiducial_J_t)/ratio_diff;
  result[3][3] = (J_t(EE + d_constants_iota[1], LL + d_constants_iota[2], QQ + d_constants_iota[3], e, p, iota + ratio_diff, a, m) - fiducial_J_t)/ratio_diff;
}

void frequency_shift_q(double p,double e,double iota,double a,double B2,double freq_shift[3]){ 
  double constants_q[4]={0.,0.,0.,0.};
  double d_constants_p[4]={0.,0.,0.,0.};
  double d_constants_e[4]={0.,0.,0.,0.};
  double d_constants_iota[4]={0.,0.,0.,0.};
  double dpeiota = 1.e-3;
  double m_dimensionless = 1;//mu/M;//1e-5;
  int prec = 500;
  //Calculating d(E,Lz,Carter Constant: K)/d(p,e,iota).
  Constants_q_once_main(p,e,cos(iota), constants_q, a);
  Constants_q_once_main(p+dpeiota,e,cos(iota), d_constants_p, a);
  Constants_q_once_main(p,e+dpeiota,cos(iota), d_constants_e, a);
  Constants_q_once_main(p,e,cos(iota+dpeiota), d_constants_iota, a);
  for(int i=1;i<=3;i++) d_constants_p[i] = (d_constants_p[i] - constants_q[i])*m_dimensionless;
  for(int i=1;i<=3;i++) d_constants_e[i] = (d_constants_e[i] - constants_q[i])*m_dimensionless;
  for(int i=1;i<=3;i++) d_constants_iota[i] = (d_constants_iota[i] - constants_q[i])*m_dimensionless;
  d_constants_p[3] = (d_constants_p[3])*m_dimensionless;
  d_constants_e[3] = (d_constants_e[3])*m_dimensionless;
  d_constants_iota[3] = (d_constants_iota[3])*m_dimensionless;

  double EE = constants_q[1]*m_dimensionless;
  double LL = constants_q[2]*m_dimensionless;
  double QQ = constants_q[3]*m_dimensionless*m_dimensionless;
  double rp = p/(1.+ e);
  double ra = p/(1.- e);

  double dJrdHH = dJr_dHH(EE, LL, QQ, a, m_dimensionless, rp, ra, prec);
  double dJrdQQ = dJr_dQQ(EE, LL, QQ, a, m_dimensionless, rp, ra, prec);
  double dJrdLL = dJr_dLL(EE, LL, QQ, a, m_dimensionless, rp, ra, prec);
  double dJrdEE = dJr_dEE(EE, LL, QQ, a, m_dimensionless, rp, ra, prec);
  double dJthetadHH = dJtheta_dHH(EE, LL, QQ, a, m_dimensionless, M_PI/2.-iota, M_PI/2., prec);
  double dJthetadQQ = dJtheta_dQQ(EE, LL, QQ, a, m_dimensionless, M_PI/2.-iota, M_PI/2., prec);
  double dJthetadLL = dJtheta_dLL(EE, LL, QQ, a, m_dimensionless, M_PI/2.-iota, M_PI/2., prec);
  double dJthetadEE = dJtheta_dEE(EE, LL, QQ, a, m_dimensionless, M_PI/2.-iota, M_PI/2., prec);

  double omega_denominator = dJrdHH*dJthetadQQ - dJrdQQ*dJthetadHH;
  double omega_r = dJthetadQQ/omega_denominator;
  double omega_theta = -dJrdQQ/omega_denominator;
  double omega_phi = (dJrdQQ*dJthetadLL - dJrdLL*dJthetadQQ)/omega_denominator;
  double Gamma = (dJrdEE*dJthetadQQ - dJrdQQ*dJthetadEE)/omega_denominator;

  //printf("Kerr freq r,th,ph %14.12e %14.12e %14.12e\n",omega_r/Gamma,omega_theta/Gamma,omega_phi/Gamma);
  //printf("E,L,Q %14.12e %14.12e %14.12e\n",constants_q[1],constants_q[2],constants_q[3]);

  double H1_derives_result[4], J_derives_result[4][4],J_inverse[4][4];
  H1_derives (EE, LL, QQ, e, p, iota, a, m_dimensionless, B2, d_constants_p, d_constants_e, d_constants_iota, H1_derives_result, dpeiota);

  J_derives (EE, LL, QQ, e, p, iota, a, m_dimensionless, B2, d_constants_p, d_constants_e, d_constants_iota, J_derives_result, dpeiota);

  //doing inversion
  gsl_matrix * J_in_gsl = gsl_matrix_alloc( 4, 4 );
  gsl_matrix * J_in_gsl_inv;

  for ( int matrix_i = 0; matrix_i < 4; matrix_i++ )
    for ( int matrix_j = 0; matrix_j < 4; matrix_j++ )
      gsl_matrix_set( J_in_gsl, matrix_i, matrix_j, J_derives_result[matrix_i][matrix_j] );

  J_in_gsl_inv = GetInverse(J_in_gsl);

  for ( int print_i = 0; print_i < 4; print_i++ ){        
    for ( int print_j = 0; print_j < 4; print_j++ )  J_inverse[print_i][print_j] = gsl_matrix_get( J_in_gsl_inv, print_i, print_j ) ;
  }

  gsl_matrix_free( J_in_gsl );
  gsl_matrix_free( J_in_gsl_inv );

  //calculating freq
  double frequency_shift[4]={0.,0.,0.,0.};
  for (int print_j = 0; print_j < 4; print_j++) frequency_shift[0] = frequency_shift[0]+ H1_derives_result[print_j]*J_inverse[print_j][0];
  for (int print_j = 0; print_j < 4; print_j++) frequency_shift[1] = frequency_shift[1]+ H1_derives_result[print_j]*J_inverse[print_j][1];
  for (int print_j = 0; print_j < 4; print_j++) frequency_shift[2] = frequency_shift[2]+ H1_derives_result[print_j]*J_inverse[print_j][2];
  for (int print_j = 0; print_j < 4; print_j++) frequency_shift[3] = frequency_shift[3]+ H1_derives_result[print_j]*J_inverse[print_j][3];

  freq_shift[0] = (frequency_shift[0]/Gamma - omega_r*frequency_shift[3]/Gamma/Gamma)/(2.*M_PI);
  freq_shift[1] = (frequency_shift[1]/Gamma - omega_theta*frequency_shift[3]/Gamma/Gamma)/(2.*M_PI);
  freq_shift[2] = (frequency_shift[2]/Gamma - omega_phi*frequency_shift[3]/Gamma/Gamma)/(2.*M_PI);

}
//------------ QAAKMOD END -----------------------

//------------ QAAK_waveform Start -----------------------

void QAAK_main(char AAK_path[100], bool AAK_backint ,bool AAK_LISA ,bool AAK_traj ,bool AAK_SNR ,bool AAK_timing,int AAK_length ,double AAK_dt ,double AAK_p ,double AAK_T ,double AAK_f ,double AAK_T_fit , 
double AAK_mu ,double AAK_M ,double AAK_s ,double AAK_e ,double AAK_iota ,double AAK_gamma ,double AAK_psi , 
double AAK_theta_S ,double AAK_phi_S ,double AAK_theta_K ,double AAK_phi_K ,double AAK_alpha ,double AAK_D, double AAK_q_q){

  clock_t ticks=clock();

  GKTrajFast gktraj3(cos(AAK_iota),AAK_s,AAK_q_q);
  gktraj3.p=AAK_p;
  gktraj3.ecc=AAK_e;
  int maxsteps=100;
  int steps=0;
  double dt_fit=min(AAK_T_fit,AAK_length*AAK_dt/SOLARMASSINSEC/AAK_M/AAK_M*AAK_mu)/(maxsteps-1);
  TrajData *traj3;
  traj3=(TrajData*)malloc((size_t)((maxsteps+1)*sizeof(TrajData)));
  gktraj3.Eccentric(dt_fit,traj3,maxsteps,steps);
  double Omega_t[3],ang[3],map_t[3],e_traj[steps],v_map[steps],M_map[steps],s_map[steps],dt_map[steps];
  double Phi;

  double freq_shift[3];
  double B2=-(AAK_q_q + AAK_s*AAK_s)/sqrt(5./4./M_PI);

  if(fabs(B2)>1.e-16){   //QAAK MOD: If the deviation is too small, we will treat it as 0.
    for(int i=1;i<=steps;i++){
      IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,AAK_s);
      geodesic_t.Frequencies(Omega_t);

      frequency_shift_q(traj3[i].p,traj3[i].ecc,acos(traj3[i].cosiota),AAK_s,B2,freq_shift);
      Omega_t[0]=Omega_t[0]+freq_shift[0];
      Omega_t[1]=Omega_t[1]+freq_shift[1];
      Omega_t[2]=Omega_t[2]+freq_shift[2];
      //printf("freq_shift[012] %14.12e %14.12e %14.12e\n",freq_shift[0],freq_shift[1],freq_shift[2]);
      if(i==1){
      ParAng(ang,AAK_e,AAK_iota,AAK_gamma,AAK_psi,AAK_theta_S,AAK_phi_S,AAK_theta_K,AAK_phi_K,AAK_alpha,geodesic_t.zedminus);
        Phi=ang[0]; // initial mean anomaly
      }
      ParMap(map_t,Omega_t,traj3[i].p,AAK_M,AAK_s,traj3[i].ecc,AAK_iota,AAK_q_q);
      e_traj[i-1]=traj3[i].ecc;
      v_map[i-1]=map_t[0]; // mapped initial velocity in c 
      M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
      s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
      dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*AAK_M*AAK_M/AAK_mu;
      //printf("v, M, s || %lf %lf %lf\n",v_map[i-1],M_map[i-1],s_map[i-1]);
    }
  }else{
    for(int i=1;i<=steps;i++){
      IEKG geodesic_t(traj3[i].p,traj3[i].ecc,traj3[i].cosiota,AAK_s);
      geodesic_t.Frequencies(Omega_t);


      if(i==1){
        ParAng(ang,AAK_e,AAK_iota,AAK_gamma,AAK_psi,AAK_theta_S,AAK_phi_S,AAK_theta_K,AAK_phi_K,AAK_alpha,geodesic_t.zedminus);
        Phi=ang[0]; // initial mean anomaly
      }
      ParMap(map_t,Omega_t,traj3[i].p,AAK_M,AAK_s,traj3[i].ecc,AAK_iota,AAK_q_q);
      e_traj[i-1]=traj3[i].ecc;
      v_map[i-1]=map_t[0]; // mapped initial velocity in c 
      M_map[i-1]=map_t[1]; // mapped BH mass in solar masses
      s_map[i-1]=map_t[2]; // mapped spin parameter a/M = S/M^2
      dt_map[i-1]=traj3[i].t*SOLARMASSINSEC*AAK_M*AAK_M/AAK_mu;
      //printf("v, M, s || %lf %lf %lf\n",v_map[i-1],M_map[i-1],s_map[i-1]);
    }
  }

  double *t,*hI,*hII;
  t=(double*)malloc(AAK_length*sizeof(double));
  hI=(double*)fftw_malloc(AAK_length*sizeof(double));
  hII=(double*)fftw_malloc(AAK_length*sizeof(double));
  GenWave(t,hI,hII,AAK_dt,AAK_length,e_traj,v_map,AAK_M,M_map,AAK_mu,AAK_s,
s_map,AAK_D,AAK_iota,AAK_gamma,Phi,AAK_theta_S,AAK_phi_S,AAK_alpha,AAK_theta_K,
AAK_phi_K,dt_map,steps,AAK_backint,AAK_LISA,false,AAK_q_q);

  ticks=clock()-ticks;
  double secs=((double)ticks)/CLOCKS_PER_SEC;
  // ----------

  // ----- output to file -----
  FILE *file;
  char filename[100];
  strcpy(filename,AAK_path);
  strcat(filename,"_wave.dat");
  if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
  file=fopen(filename,"w");
  for(int i=0;i<AAK_length;i++) fprintf(file,"%14.15e %14.15e %14.15e\n",t[i],hI[i],hII[i]);
  fclose(file);

  if(AAK_traj==true){
    double *pvec,*evec;
    pvec=(double*)malloc(AAK_length*sizeof(double));
    evec=(double*)malloc(AAK_length*sizeof(double));
    GenWave(t,pvec,evec,AAK_dt,AAK_length,e_traj,v_map,AAK_M,M_map,AAK_mu,AAK_s,
s_map,AAK_D,AAK_iota,AAK_gamma,Phi,AAK_theta_S,AAK_phi_S,AAK_alpha,AAK_theta_K,
AAK_phi_K,dt_map,steps,AAK_backint,AAK_LISA,true,AAK_q_q);
    strcpy(filename,AAK_path);
    strcat(filename,"_traj.dat");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    double dt_RR=0.001; // radiation-reaction timestep for downsampling
    int i_RR=(int)(dt_RR*(SOLARMASSINSEC*AAK_M*AAK_M/AAK_mu)/AAK_dt)+1;
    int i_max=0;
    while(pvec[i_max]>0.) i_max++;
    for(int i=0;i<i_max;i++){
      if(i%i_RR==0 || i+i_RR>=i_max){
        IEKG geodesic_t(pvec[i],evec[i],cos(AAK_iota),AAK_s);
        fprintf(file,"%14.12e %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e\n",t[i],pvec[i],evec[i],AAK_iota,geodesic_t.E,geodesic_t.Lz,geodesic_t.Q);
      }
    }
    fclose(file);
    free(pvec);
    free(evec);
  }

  if(AAK_SNR==true||AAK_timing==true){
    strcpy(filename,AAK_path);
    strcat(filename,"_info.txt");
    if(CheckFile(filename)==1) fprintf(stderr,"Output warning: Overwriting %s\n",filename);
    file=fopen(filename,"w");
    if(AAK_SNR==true){
      fftw_plan planhI,planhII;
      planhI=fftw_plan_r2r_1d(AAK_length,hI,hI,FFTW_R2HC,FFTW_ESTIMATE);
      planhII=fftw_plan_r2r_1d(AAK_length,hII,hII,FFTW_R2HC,FFTW_ESTIMATE);
      fftw_execute(planhI);
      fftw_execute(planhII);
      double hI2=InnProd(hI,hI,AAK_dt,AAK_length);
      double hII2=InnProd(hII,hII,AAK_dt,AAK_length);
      fprintf(file,"SNR_I: %.1f\nSNR_II: %.1f\nSNR: %.1f\n\n",sqrt(hI2),sqrt(hII2),sqrt(hI2+hII2));
      fftw_destroy_plan(planhI);
      fftw_destroy_plan(planhII);
      // ----------
    }
    if(AAK_timing==true){
      fprintf(file,"Clock ticks: %.2e\nSeconds: %.2fs\n\n",(double)ticks,secs);
    }
    fclose(file);
  }

  free(traj3);
  fftw_free(hI);
  fftw_free(hII);

//return( hI2+hII2 );

}

//------------ QAAK_waveform End -----------------------


int main(int argc, char *argv[]){

  // ----- interface -----
  if(argc<2){
    fprintf(stderr,"Input error: Specify path to settings/parameters file\n");
    return 0;
  }
  char *AAK_par=argv[1];
  if(CheckFile(AAK_par)==0){
    fprintf(stderr,"Input error: Cannot read settings/parameters file\n");
    return 0;
  }
  // ----------
  // ----- load settings and parameters -----
  SetPar AAK;
  if(LoadSetPar(&AAK,AAK_par)==0) return 0;

  double Q_here = -AAK.s*AAK.s+AAK.q_q; //QAAK MOD: Note that Q_here is quadrupole.

  QAAK_main(AAK.path, AAK.backint ,AAK.LISA ,AAK.traj ,AAK.SNR ,AAK.timing, AAK.length, AAK.dt, AAK.p, AAK.T, AAK.f, AAK.T_fit, AAK.mu, AAK.M, AAK.s, AAK.e, AAK.iota, AAK.gamma, AAK.psi, AAK.theta_S, AAK.phi_S, AAK.theta_K, AAK.phi_K, AAK.alpha, AAK.D, Q_here);

}

