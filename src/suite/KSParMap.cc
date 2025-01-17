#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>

#include "IEKG.h"
#include "KSParMap.h"

//---------- QAAK MOD: We  choose  a  rough  approach  by  replacing  all the 3PN terms quadratic in spin with −Q in the higher order equations. The q_q is quadrupole.

void cross(const gsl_vector *u,const gsl_vector *v,gsl_vector *w){

  gsl_vector_set(w,0,gsl_vector_get(u,1)*gsl_vector_get(v,2)-gsl_vector_get(u,2)*gsl_vector_get(v,1));
  gsl_vector_set(w,1,gsl_vector_get(u,2)*gsl_vector_get(v,0)-gsl_vector_get(u,0)*gsl_vector_get(v,2));
  gsl_vector_set(w,2,gsl_vector_get(u,0)*gsl_vector_get(v,1)-gsl_vector_get(u,1)*gsl_vector_get(v,0));

}

void ParAng(double ang[],double e,double iota,double gamma,double psi,double theta_S,double phi_S,double theta_K,double phi_K,double alpha,double zm){

  ang[0]=atan2(sqrt(1.-e*e)*sin(psi),e+cos(psi))-e*sqrt(1.-e*e)*sin(psi)/(1.+e*cos(psi));

  gsl_vector *pos,*n,*S,*L,*b,*nxS,*LxS,*Sx_nxS,*Lx_LxS;
  gsl_matrix *A,*B;
  pos=gsl_vector_alloc(3);
  n=gsl_vector_alloc(3);
  S=gsl_vector_alloc(3);
  L=gsl_vector_alloc(3);
  b=gsl_vector_calloc(3);
  nxS=gsl_vector_alloc(3);
  LxS=gsl_vector_alloc(3);
  Sx_nxS=gsl_vector_alloc(3);
  Lx_LxS=gsl_vector_alloc(3);
  A=gsl_matrix_alloc(3,3);
  B=gsl_matrix_alloc(3,3);
  gsl_vector_set(n,0,sin(theta_S)*cos(phi_S));
  gsl_vector_set(n,1,sin(theta_S)*sin(phi_S));
  gsl_vector_set(n,2,cos(theta_S));
  gsl_vector_set(S,0,sin(theta_K)*cos(phi_K));
  gsl_vector_set(S,1,sin(theta_K)*sin(phi_K));
  gsl_vector_set(S,2,cos(theta_K));
  gsl_vector_set(L,0,cos(iota)*sin(theta_K)*cos(phi_K)+sin(iota)*(sin(alpha)*sin(phi_K)-cos(alpha)*cos(theta_K)*cos(phi_K)));
  gsl_vector_set(L,1,cos(iota)*sin(theta_K)*sin(phi_K)-sin(iota)*(sin(alpha)*cos(phi_K)+cos(alpha)*cos(theta_K)*sin(phi_K)));
  gsl_vector_set(L,2,cos(iota)*cos(theta_K)+sin(iota)*cos(alpha)*sin(theta_K));
  gsl_vector_set(b,0,cos(gamma+psi));
  gsl_vector_set(b,1,sin(gamma+psi));
  cross(n,S,nxS);
  gsl_blas_dscal(1./gsl_blas_dnrm2(nxS),nxS);
  cross(L,S,LxS);
  gsl_blas_dscal(1./gsl_blas_dnrm2(LxS),LxS);
  cross(S,nxS,Sx_nxS);
  cross(L,LxS,Lx_LxS);
  gsl_matrix_set_col(A,0,nxS);
  gsl_matrix_set_col(A,1,Sx_nxS);
  gsl_matrix_set_col(A,2,S);
  gsl_matrix_set_col(B,0,LxS);
  gsl_matrix_set_col(B,1,Lx_LxS);
  gsl_matrix_set_col(B,2,L);
  gsl_blas_dgemv(CblasNoTrans,1.,B,b,0.,pos);
  gsl_vector_memcpy(b,pos);
  gsl_blas_dgemv(CblasTrans,1.,A,b,0.,pos);

  ang[1]=acos(gsl_vector_get(pos,2)/sqrt(zm));

  ang[2]=atan2(gsl_vector_get(pos,1),gsl_vector_get(pos,0));

  gsl_vector_free(pos);
  gsl_vector_free(n);
  gsl_vector_free(S);
  gsl_vector_free(L);
  gsl_vector_free(b);
  gsl_vector_free(nxS);
  gsl_vector_free(LxS);
  gsl_vector_free(Sx_nxS);
  gsl_vector_free(Lx_LxS);
  gsl_matrix_free(A);
  gsl_matrix_free(B);

}

// ----- 3PN O(e^6) equations (Sago & Fujita, 2015) -----
double dvdt(double v,double e,double Y,double m,double M,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double v6=v3*v3;
  double e2=e*e;
  double e4=e2*e2;
  double e6=e4*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=(32.*pow(1. - e2,1.5)*m*v6*v3*(1. + (7.*e2)/8. + 
       ((-5944. - 7040.*e2 + 8539.*e4)*v2)/2688. - 
       (v3*((-18432. - 55872.*e2 - 7056.*e4 + 49.*e6)*M_PI + 
            48.*(1064. + 1516.*e2 + 475.*e4)*q*Y))/4608. + 
       (v3*v2*((-4791168. - 28113984.*e2 + 12239226.*e4 + 4005097.*e6)*M_PI - 
            144.*(139296. + 58408.*e2 + 249968.*e4 + 106707.*e6)*q*Y))/774144. + 
       (v2*v2*(953325.*e6 + 8.*(34103. + 189.*q2*(-329. + 815.*Y2)) + 
            12.*e2*(-526955. + 126.*q2*(-929. + 1431.*Y2)) + 
            3.*e4*(-1232809. + 63.*q2*(-1051. + 2997.*Y2))))/145152. + 
       v6*(117.72574285227559 - (1712.*EulerGamma)/105. + (16.*PI2)/3. - (331.*q2)/192. - 
          (289.*M_PI*q*Y)/6. + (145759.*q2*Y2)/1344. - (3424.*log(2.))/105. + 
          e2*(764.5906705662063 - (24503.*EulerGamma)/210. + (229.*PI2)/6. + (2129.*q2)/42. - 
             (4225.*M_PI*q*Y)/24. + (27191.*q2*Y2)/224. + (1391.*log(2.))/30. - (78003.*log(3.))/280.) + 
          e4*(723.1604428881717 - (11663.*EulerGamma)/140. + (109.*PI2)/4. - (56239.*q2)/10752. - 
             (17113.*M_PI*q*Y)/192. + (414439.*q2*Y2)/3584. - (418049.*log(2.))/84. + (3042117.*log(3.))/1120.) + 
          e6*(202.0979603432282 - (2461.*EulerGamma)/560. + (23.*PI2)/16. - (3571.*q2)/3584. - 
             (108577.*M_PI*q*Y)/13824. + (41071.*q2*Y2)/1536. + (94138279.*log(2.))/2160. - 
             (42667641.*log(3.))/3584. - (1044921875.*log(5.))/96768.) - 
          (107.*(256. + 1832.*e2 + 1308.*e4 + 69.*e6)*log(v))/1680.)))/(5.*M*M);
  return eq;
}

double dedt(double v,double e,double Y,double m,double M,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double v6=v3*v3;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=-(e*pow(1. - e2,1.5)*m*v6*v2*(2451456. + 975744.*e2 + 
        144.*(-54792. - 18600.*e2 + 22579.*e4)*v2 + 
        84.*v3*((189120. + 286512.*e2 + 24217.*e4)*M_PI - 
           48.*(7032. + 5592.*e2 + 1313.*e4)*q*Y) + 
        3.*v3*v2*((-16885824. - 42474444.*e2 + 5678971.*e4)*M_PI - 
           48.*(106568. + 187862.*e2 + 442811.*e4)*q*Y) + 
        8.*v2*v2*(6.*e2*(-2070667. + 126.*q2*(-2975. + 4009.*Y2)) + 
           8.*(-286397. + 63.*q2*(-3179. + 5869.*Y2)) + 
           e4*(-3506201. + 63.*q2*(-3191. + 9009.*Y2))) + 
        2451456.*v6*(241.02948298890263 - (82283.*EulerGamma)/1995. + (769.*PI2)/57. + 
           (180255.*q2)/8512. - (11809.*M_PI*q*Y)/152. + (598987.*q2*Y2)/8512. - (11021.*log(2.))/285. - 
           (234009.*log(3.))/5320. + e2*(1048.5728116817568 - (297674.*EulerGamma)/1995. + (2782.*PI2)/57. + 
              (536653.*q2)/8512. - (91375.*M_PI*q*Y)/608. + (356845.*q2*Y2)/8512. - 
              (2982946.*log(2.))/1995. + (1638063.*log(3.))/3040.) + 
           e4*(725.9069951397873 - (1147147.*EulerGamma)/15960. + (10721.*PI2)/456. + 
              (56509.*q2)/9728. - (1739605.*M_PI*q*Y)/29184. + (3248951.*q2*Y2)/68096. + 
              (760314287.*log(2.))/47880. - (1022385321.*log(3.))/340480. - (1044921875.*log(5.))/204228.) - 
           (107.*(6152. + 22256.*e2 + 10721.*e4)*log(v))/15960.)))/(120960.*M*M);
  return eq;
}

double dtdm(double v,double e,double Y,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double v4=v2*v2;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=(2.*e4*(240. + v2*(-120. + v*(42.*(-27. + 4.*q2)*v + (-6567. + 1996.*q2)*v3 + 
              48.*q*(8. + 77.*v2)*Y - 4.*q2*v*(90. + 1577.*v2)*Y2))) + 
     e4*e2*(560. + v2*(-360. + 960.*q*v*Y + 8816.*q*v3*Y + 
           v4*(-15565. + 24.*q2*(200. - 629.*Y2)) + 
           v2*(-2742. + 80.*q2*(5. - 11.*Y2)))) - 
     8.*e2*(-48. + v2*(8. - 64.*q*v*Y - 688.*q*v3*Y + 
           2.*v2*(99. + 16.*q2*(-1. + 2.*Y2)) + v4*(1233. + 8.*q2*(-47. + 150.*Y2))))
       + 16.*(16. + v2*(24. + v2*(27.*(2. + 5.*v2) - 48.*q*v*Y + 
              4.*q2*(2. - v2 + (-2. + 3.*v2)*Y2)))))/(256.*v4);
  return eq;
}

double drdm(double v,double e,double Y,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=(16. + 8.*(-3. + e2)*v2 - 16.*(-3. + e2)*q*v3*Y + 
     8.*(33. + 4.*e2 - 3.*e4)*q*v3*v2*Y + 
     v3*v3*(-351. + 132.*q2 + e2*(-135. + 21.*e2 + 5.*e4 + 2.*(7. + e2)*q2) + 
        2.*(-204. + 13.*e2*(-3. + e2))*q2*Y2) + 
     2.*v2*v2*(-45. + 3.*e4 + 4.*q2*(1. - 4.*Y2) + 2.*e2*q2*(1. + Y2)))/(16.*v);
  return eq;
}

double dthetadm(double v,double e,double Y,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=(16. + 8.*(3. + e2)*v2 - 16.*(3. + e2)*q*v3*Y - 
     8.*(3. + e2)*(5. + 3.*e2)*q*v3*v2*Y + 
     v3*v3*(135. - 54.*q2 + e2*(5.*(27. + 9.*e2 + e4) + 2.*(-38. + e2)*q2) + 
        2.*(57. + 90.*e2 + 13.*e4)*q2*Y2) + 
     2.*v2*v2*(27. + 3.*e4 + 2.*q2*(-1. + 7.*Y2) + 2.*e2*(9. + q2*(1. + Y2))))/(16.*v);
  return eq;
}

double dphidm(double v,double e,double Y,double q, double q_q){
  double v2=v*v;
  double v3=v2*v;
  double e2=e*e;
  double e4=e2*e2;
  double Y2=Y*Y;
  double q2=-q_q;
  double eq=(16. + 8.*(3. + e2)*v2 - 16.*q*v3*(-2. + (3. + e2)*Y) - 
     8.*q*v3*v2*(-6. + 15.*Y + 3.*e4*Y + 2.*e2*(-4. + 7.*Y)) + 
     2.*v2*v2*(27. + 3.*e4 + 2.*q2*(-1. + Y)*(1. + 7.*Y) + 2.*e2*(9. + q2*(1. + Y2))) + 
     v3*v3*(5.*e4*e2 + e4*(45. + q2*(2. + 26.*Y2)) + 
        e2*(135. + 4.*q2*(-19. + 5.*Y*(-7. + 9.*Y))) + 3.*(45. + 2.*q2*(-9. + Y*(-6. + 19.*Y)))))/(16.*v);
  return eq;
}

// ----------

int sol_fun(const gsl_vector *x,void *p,gsl_vector *f){

  double Omega_r=((struct sol_par*)p)->Omega_r;
  double Omega_theta=((struct sol_par*)p)->Omega_theta;
  double Omega_phi=((struct sol_par*)p)->Omega_phi;
  double M=((struct sol_par*)p)->M;
  double e=((struct sol_par*)p)->e;
  double iota=((struct sol_par*)p)->iota;
  double q=((struct sol_par*)p)->q;

  const double v_map=gsl_vector_get(x,0);
  const double M_map=gsl_vector_get(x,1);
  const double s_map=gsl_vector_get(x,2);
  const double M_map_INSEC =M_map*SOLARMASSINSEC;

  double Omega_r_map,Omega_theta_map,Omega_phi_map;  //QAAK MOD added
  if(cos(iota)>0){
    Omega_r_map=drdm(v_map,e,cos(iota),s_map,q)/dtdm(v_map,e,cos(iota),s_map,q)/2./M_PI;
    Omega_theta_map=dthetadm(v_map,e,cos(iota),s_map,q)/dtdm(v_map,e,cos(iota),s_map,q)/2./M_PI;
    Omega_phi_map=dphidm(v_map,e,cos(iota),s_map,q)/dtdm(v_map,e,cos(iota),s_map,q)/2./M_PI;
  }
  else{
    Omega_r_map=drdm(v_map,e,-cos(iota),-s_map,q)/dtdm(v_map,e,-cos(iota),-s_map,q)/2./M_PI;
    Omega_theta_map=dthetadm(v_map,e,-cos(iota),-s_map,q)/dtdm(v_map,e,-cos(iota),-s_map,q)/2./M_PI;
    Omega_phi_map=-dphidm(v_map,e,-cos(iota),-s_map,q)/dtdm(v_map,e,-cos(iota),-s_map,q)/2./M_PI;
  }

  const double f0=Omega_r*M_map-Omega_r_map*M;
  const double f1=Omega_theta*M_map-Omega_theta_map*M;
  const double f2=Omega_phi*M_map-Omega_phi_map*M;

  gsl_vector_set(f,0,f0);
  gsl_vector_set(f,1,f1);
  gsl_vector_set(f,2,f2);

  return GSL_SUCCESS;

}

int sol_inv(const gsl_vector *x,void *p,gsl_vector *f){

  double Omega_r_map=((struct sol_par*)p)->Omega_r;
  double Omega_theta_map=((struct sol_par*)p)->Omega_theta;
  double Omega_phi_map=((struct sol_par*)p)->Omega_phi;
  double M_map=((struct sol_par*)p)->M;
  double e=((struct sol_par*)p)->e;
  double iota=((struct sol_par*)p)->iota;

  const double v=gsl_vector_get(x,0);
  const double M=gsl_vector_get(x,1);
  const double s=gsl_vector_get(x,2);

  double Omega[3];
  IEKG geodesic(1./v/v,e,cos(iota),s);
  geodesic.Frequencies(Omega);

  const double f0=Omega[0]*M_map-Omega_r_map*M;
  const double f1=Omega[1]*M_map-Omega_theta_map*M;
  const double f2=Omega[2]*M_map-Omega_phi_map*M;

  gsl_vector_set(f,0,f0);
  gsl_vector_set(f,1,f1);
  gsl_vector_set(f,2,f2);

  return GSL_SUCCESS;

}

void print_state(size_t i, gsl_multiroot_fsolver *sol){

  printf("iter = %zu x = % .3f % .3f % .3f f(x) = % .3e % .3e % .3e\n",i,
    gsl_vector_get(sol->x,0),gsl_vector_get(sol->x,1),gsl_vector_get(sol->x,2),
    gsl_vector_get(sol->f,0),gsl_vector_get(sol->f,1),gsl_vector_get(sol->f,2));

}

void ParMap(double map[],double Omega[],double p,double M,double s,double e,double iota,double q){

  const gsl_multiroot_fsolver_type *sol_typ;
  gsl_multiroot_fsolver *sol;

  int status;
  size_t i=0;

  const size_t n=3;
  struct sol_par par={Omega[0],Omega[1],Omega[2],M,e,iota,q};
  gsl_multiroot_function f={&sol_fun,n,&par};

  double x0[n]={1./sqrt(p),M,s};
  gsl_vector *x=gsl_vector_alloc(n);

  gsl_vector_set(x,0,x0[0]);
  gsl_vector_set(x,1,x0[1]);
  gsl_vector_set(x,2,x0[2]);

  sol_typ=gsl_multiroot_fsolver_hybrids;
  sol=gsl_multiroot_fsolver_alloc(sol_typ,n);
  gsl_multiroot_fsolver_set(sol,&f,x);

  //print_state(i,sol);

  do{
    i++;
    status=gsl_multiroot_fsolver_iterate(sol);
    //print_state(i,sol);
    if(status) break;
    status=gsl_multiroot_test_residual(sol->f,1.e-6);
  }
  while(status==GSL_CONTINUE&&i<1000);

  //printf("status = %s\n",gsl_strerror(status));

  map[0]=gsl_vector_get(sol->x,0);
  map[1]=gsl_vector_get(sol->x,1);
  map[2]=gsl_vector_get(sol->x,2);

  gsl_multiroot_fsolver_free(sol);
  gsl_vector_free(x);

}

void ParInvMap(double map[],double Omega[],double p,double M,double s,double e,double iota){

  const gsl_multiroot_fsolver_type *sol_typ;
  gsl_multiroot_fsolver *sol;

  int status;
  size_t i=0;

  const size_t n=3;
  struct sol_par par={Omega[0],Omega[1],Omega[2],M,e,iota};
  gsl_multiroot_function f={&sol_inv,n,&par};

  double x0[n]={1./sqrt(p),M,s};
  gsl_vector *x=gsl_vector_alloc(n);

  gsl_vector_set(x,0,x0[0]);
  gsl_vector_set(x,1,x0[1]);
  gsl_vector_set(x,2,x0[2]);

  sol_typ=gsl_multiroot_fsolver_hybrids;
  sol=gsl_multiroot_fsolver_alloc(sol_typ,n);
  gsl_multiroot_fsolver_set(sol,&f,x);

  //print_state(i,sol);

  do{
    i++;
    status=gsl_multiroot_fsolver_iterate(sol);
    //print_state(i,sol);
    if(status) break;
    status=gsl_multiroot_test_residual(sol->f,1.e-6);
  }
  while(status==GSL_CONTINUE&&i<1000);

  //printf("status = %s\n",gsl_strerror(status));

  map[0]=gsl_vector_get(sol->x,0);
  map[1]=gsl_vector_get(sol->x,1);
  map[2]=gsl_vector_get(sol->x,2);

  gsl_multiroot_fsolver_free(sol);
  gsl_vector_free(x);

}

void Interp(double *x_in,double *y_in,int n_in,double *x_out,double *y_out,int n_out){

  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline,n_in);
  gsl_spline_init(spline,x_in,y_in,n_in);

  for(int i=0;i<n_out;i++) y_out[i]=gsl_spline_eval(spline,x_out[i],acc);

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

}

void PolyFit(double *coeff,double *x,double *y,int n){

  double chisq;
  gsl_matrix *X,*cov;
  gsl_vector *Y,*c;

  X = gsl_matrix_alloc(n-1,2);
  Y = gsl_vector_alloc(n-1);
  c = gsl_vector_alloc(2);
  cov = gsl_matrix_alloc(2,2);

  for(int i=0;i<n-1;i++){
    gsl_matrix_set(X,i,0,x[i+1]-x[0]);
    gsl_matrix_set(X,i,1,pow(x[i+1]-x[0],2));
    gsl_vector_set(Y,i,y[i+1]-y[0]);
  }

  gsl_multifit_linear_workspace *work=gsl_multifit_linear_alloc(n-1,2);
  gsl_multifit_linear(X,Y,c,cov,&chisq,work);
  gsl_multifit_linear_free(work);

  for(int i=0;i<2;i++){
    coeff[i]=gsl_vector_get(c,i);
  }

  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(c);
  gsl_matrix_free(cov);

}

void RotCoeff(double rot[],double iota,double theta_S,double phi_S,double theta_K,double phi_K,double alpha){

  gsl_vector *n,*L,*S,*nxL,*nxS,*nxL_xn,*nxS_xn;
  n=gsl_vector_alloc(3);
  L=gsl_vector_alloc(3);
  S=gsl_vector_alloc(3);
  nxL=gsl_vector_alloc(3);
  nxS=gsl_vector_alloc(3);

  gsl_vector_set(n,0,sin(theta_S)*cos(phi_S));
  gsl_vector_set(n,1,sin(theta_S)*sin(phi_S));
  gsl_vector_set(n,2,cos(theta_S));
  gsl_vector_set(S,0,sin(theta_K)*cos(phi_K));
  gsl_vector_set(S,1,sin(theta_K)*sin(phi_K));
  gsl_vector_set(S,2,cos(theta_K));
  gsl_vector_set(L,0,cos(iota)*sin(theta_K)*cos(phi_K)+sin(iota)*(sin(alpha)*sin(phi_K)-cos(alpha)*cos(theta_K)*cos(phi_K)));
  gsl_vector_set(L,1,cos(iota)*sin(theta_K)*sin(phi_K)-sin(iota)*(sin(alpha)*cos(phi_K)+cos(alpha)*cos(theta_K)*sin(phi_K)));
  gsl_vector_set(L,2,cos(iota)*cos(theta_K)+sin(iota)*cos(alpha)*sin(theta_K));
  cross(n,L,nxL);
  cross(n,S,nxS);

  double norm=gsl_blas_dnrm2(nxL)*gsl_blas_dnrm2(nxS);
  double dot,cosrot,sinrot;
  gsl_blas_ddot(nxL,nxS,&dot);
  cosrot=dot/norm;
  gsl_blas_ddot(L,nxS,&dot);
  sinrot=dot;
  gsl_blas_ddot(S,nxL,&dot);
  sinrot-=dot;
  sinrot/=norm;

  rot[0]=2.*cosrot*cosrot-1.;
  rot[1]=cosrot*sinrot;
  rot[2]=-rot[1];
  rot[3]=rot[0];

  gsl_vector_free(n);
  gsl_vector_free(L);
  gsl_vector_free(S);
  gsl_vector_free(nxL);
  gsl_vector_free(nxS);

}

void EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, double** q, double** n){
	double Omega = 2.*M_PI*3.1709791983764586e-8;
   	double alpha = Omega*t + kappa0;
	double xi[3];
	xi[0] = lambda0;
	xi[1] = lambda0 + 2.0*M_PI/3.0;
	xi[2] = lambda0 + 4.0*M_PI/3.0;
	double RAU = AUsec;
	R[0] = RAU*cos(alpha);
	R[1] = RAU*sin(alpha);
	R[2] = 0;
	for(int i=0; i<3; i++){ // loop over s/c
    	q[i][0] = ( sin(alpha)*cos(alpha)*sin(xi[i]) - (1.+ sin(alpha)*sin(alpha))*cos(xi[i])  )/(2.0*sqrt(3.));
        q[i][1] = ( sin(alpha)*cos(alpha)*cos(xi[i]) - (1.+ cos(alpha)*cos(alpha))*sin(xi[i])  )/(2.0*sqrt(3.));
        q[i][2] = -0.5*cos(alpha - xi[i]);
	}
	for(int i=0; i<3; i++){ // links: 1st index - s/c, 2nd index - coordinates
        n[0][i] = q[1][i] - q[2][i];
		n[1][i] = q[2][i] - q[0][i];
		n[2][i] = q[0][i] - q[1][i];
	}
}
