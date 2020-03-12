#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#include "IEKG.h"
#include "KSParMap.h"
#include "gpuAAK.h"

#include <cstring>

using namespace std;


// ----- magnitude of azimuthal angular frequency for prograde/retrograde orbits -----
double OmegaPhi(double v, double e, double cosiota, double s, double M){

  double omegaphi;
  if(cosiota>0) omegaphi=dphidm(v,e,cosiota,s)/dtdm(v,e,cosiota,s)/M;
  else omegaphi=dphidm(v,e,-cosiota,-s)/dtdm(v,e,-cosiota,-s)/M;

  return omegaphi;

}
// ----------


void PNevolution(double *t_in, double *e_in, double *v_in, double *M_in, double *S_in,
                double timestep, int vlength, double *par,
                double e_traj[], double v_map[], double M_phys, double M_map[], double S_phys, double S_map[], double dt_map[],
                int steps, int *i_plunge, int *i_buffer, bool backint, double *t_clip, double *interp_timestep, int *init_length){

  double mu=par[0];
  double M0=par[1];
  double S0=par[2];
  double e0=par[3];
  double v0=v_map[0];
  double coslam=cos(par[4]);
  double gim0=par[5];
  double Phi0=par[6];
  double alp0=par[11];

  double edot,vdot,Phidot,gimdot,alpdot;
  double edotp,vdotp,Phidotp,gimdotp,alpdotp;
  double edotm,vdotm,Phidotm,gimdotm,alpdotm;

  double dt_large=dt_map[1];
  double dt_small=dt_map[steps-1]-dt_map[steps-2];
  int j0=(int)((vlength*timestep-dt_map[steps-1])/dt_large+steps+1);
  int j_start=j0;
  int j_min=j_start;
  int j_end=2*j0-1;
  int j_max=j_end;
  int j_max_temp=j_max;
  double t_end=vlength*timestep;

  double *e_AK,*v_AK,*t_fit,*e_fit,*v_fit,*M_fit,*S_fit;
  e_AK=(double*)malloc(2*j0*sizeof(double));
  v_AK=(double*)malloc(2*j0*sizeof(double));
  t_fit=(double*)malloc(2*j0*sizeof(double));
  e_fit=(double*)malloc(2*j0*sizeof(double));
  v_fit=(double*)malloc(2*j0*sizeof(double));
  M_fit=(double*)malloc(2*j0*sizeof(double));
  S_fit=(double*)malloc(2*j0*sizeof(double));

  // ----- evolve AK from t_0 to T_fit -----
  e_AK[j0]=e0;
  v_AK[j0]=v0;
  for(int j=j0;j<j0+steps-1;j++){
    edotm=edot;
    vdotm=vdot;
    edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    if(j==j0){
      edotm=edot;
      vdotm=vdot;
    }
    e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*(dt_map[j-j0+1]-dt_map[j-j0]);
    if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      if(j_max==j_max_temp) j_max=j;
      e_AK[j+1]=e_AK[j_max];
      v_AK[j+1]=v_AK[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- compute fit coefficients -----
  int points=min(steps,j_max-j0+1);
  while(e_traj[points-1]>e_traj[points-2]) points--;
  if(points<3){fprintf(stderr,"Fitting error: t_0 too close to plunge\n"); exit(EXIT_FAILURE);}

  double *e_diff,*v_diff,*e_coeff,*v_coeff,*M_coeff,*S_coeff;
  e_diff=(double*)malloc(points*sizeof(double));
  v_diff=(double*)malloc(points*sizeof(double));
  e_coeff=(double*)malloc(2*sizeof(double));
  v_coeff=(double*)malloc(2*sizeof(double));
  M_coeff=(double*)malloc(2*sizeof(double));
  S_coeff=(double*)malloc(2*sizeof(double));

  for(int i=0;i<points;i++){
    e_diff[i]=e_traj[i]-e_AK[i+j0];
    v_diff[i]=v_map[i]-v_AK[i+j0];
  }

  PolyFit(e_coeff,dt_map,e_diff,points);
  PolyFit(v_coeff,dt_map,v_diff,points);
  PolyFit(M_coeff,dt_map,M_map,points);
  PolyFit(S_coeff,dt_map,S_map,points);
  // ----------

  // ----- evolve AK from T_fit to t_end -----
  for(int j=j0+steps-1;j<j_end;j++){
    edotm=edot;
    vdotm=vdot;
    edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
    e_AK[j+1]=e_AK[j]+(1.5*edot-.5*edotm)*dt_large;
    v_AK[j+1]=v_AK[j]+(1.5*vdot-.5*vdotm)*dt_large;
    if(e_AK[j+1]<0. || v_AK[j+1]<v_AK[j] || isnan(e_AK[j+1]) || isnan(v_AK[j+1])){
      if(j_max==j_max_temp) j_max=j;
      e_AK[j+1]=e_AK[j_max];
      v_AK[j+1]=v_AK[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- fit AK from t_0 to t_end -----
  for(int j=j0;j<=j_end;j++){
  	double dt;
  	if(j<j0+steps) dt=dt_map[j-j0];
    else dt=dt_map[steps-1]+dt_large*(j-j0-steps+1);
  	double dt2=dt*dt;
    t_fit[j]=dt;
    e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
    v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
    M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
    S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    if(e_fit[j]<0. || v_fit[j]<v_fit[max(j0,j-1)]){
      if(j_max==j_max_temp) j_max=j-1;
      e_fit[j]=e_fit[j_max];
      v_fit[j]=v_fit[j_max];
      M_fit[j]=M_fit[j_max];
      S_fit[j]=S_fit[j_max];
    }
  }
  j_max_temp=j_max;
  // ----------

  // ----- check for plunge -----
  int j_RR=(int)((SOLARMASSINSEC*M_phys*SOLARMASSINSEC*M_phys/mu)/dt_large)+1;
  double z1=1.+pow(1.-S_phys*S_phys,1./3.)*(pow(1.+S_phys,1./3.)+pow(1.-S_phys,1./3.));
  double z2=sqrt(3.*S_phys*S_phys+z1*z1);
  double LSO_min=3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2));
  double LSO_max=3.+z2+sqrt((3.-z1)*(3.+z1+2.*z2));

  for(int j=j0;j<=j_max_temp;j++){
    if(j>j_max) break;
    if((1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_min && j_max==j_max_temp){
      j_min=max(j-j_RR,j0);
      j_max=min(j+1,j_max_temp);
    }
    if(j>j0 && (j-j0)%j_RR==0 && (1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.)<LSO_max && j_max==j_max_temp){
      j_min=max(j-j_RR,j0);
      IEKG iekg((1.-e_fit[j]*e_fit[j])*pow(OmegaPhi(v_fit[j],e_fit[j],coslam,S_fit[j],M_fit[j])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j],coslam,S_phys);
      if(iekg.Stable==-1 || iekg.E>1) j_max=min(j+1,j_max_temp);
    }
  }

  while(j_max-j_min>1){
    int j_mid=(j_max+j_min)/2;
    IEKG iekg((1.-e_fit[j_mid]*e_fit[j_mid])*pow(OmegaPhi(v_fit[j_mid],e_fit[j_mid],coslam,S_fit[j_mid],M_fit[j_mid])*M_phys*SOLARMASSINSEC,-2./3.),e_fit[j_mid],coslam,S_phys);
    if(iekg.Stable==-1 || iekg.E>1) j_max=j_mid;
    else j_min=j_mid;
  }
  // ----------

  // ----- evolve and fit AK from t_0 to t_start -----
  	j_end=j_max;
    t_end=min(t_fit[j_end],vlength*timestep);
    j_start=j0+floor((t_end-vlength*timestep)/dt_large);

    //printf("%d %d %d %e %e %e, %e, %e\n", j_end, j_start, j0, t_end, v_AK[j0], e_AK[j0], e_AK[j0+1], e_AK[j0-1]);
    for(int j=j0;j>j_start;j--){
      edotp=edot;
      vdotp=vdot;
      edot=dedt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
      vdot=dvdt(v_AK[j],e_AK[j],coslam,mu,M0,S0);
      if(j==j0){
        edotp=edot;
        vdotp=vdot;
      }
      e_AK[j-1]=e_AK[j]-(1.5*edot-.5*edotp)*dt_large;
      v_AK[j-1]=v_AK[j]-(1.5*vdot-.5*vdotp)*dt_large;
    }

    for(int j=j_start;j<j0;j++){
      double dt_undecayed=dt_large*(j-j0);
      double dt;
      if(j>j0-points) dt=(dt_undecayed+dt_undecayed*(j-j0+points)/points)/2.; // decaying C1 fit
      else dt=-dt_large*points/2.;
      double dt2=dt*dt;
      t_fit[j]=dt_undecayed;
      e_fit[j]=e_AK[j]+e_coeff[0]*dt+e_coeff[1]*dt2;
      v_fit[j]=v_AK[j]+v_coeff[0]*dt+v_coeff[1]*dt2;
      M_fit[j]=M0+(M_coeff[0]*dt+M_coeff[1]*dt2)*SOLARMASSINSEC;
      S_fit[j]=S0+S_coeff[0]*dt+S_coeff[1]*dt2;
    }
  // ----------

  // ----- interpolate AK from t_start to t_end -----
  /*double *t_in,*e_in,*v_in,*M_in,*S_in;
  t_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  e_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  v_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  M_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  S_in=(double*)malloc((j_end-j_start+1)*sizeof(double));
  */
  for(int j=j0;j<=j_end;j++){
    t_in[j-j0]=t_fit[j];
    e_in[j-j0]=e_fit[j];
    v_in[j-j0]=v_fit[j];
    M_in[j-j0]=M_fit[j];
    S_in[j-j0]=S_fit[j];
  }

  int i0=vlength-(int)(t_end/timestep);

  //for (int i=0; i<j_end-j_start+1; i++) printf("%e %e %e %e %e\n", t[i], v[i], e[i], M[i], S[i]);
  // ----------

  // ----- evolve phases from t_0 to t_end -----
  if(j_max==j_max_temp) *i_plunge=vlength-1;
  else *i_plunge=i0+(int)(t_fit[j_min]/timestep);

  double *t_find_plunge = (double*) malloc(1*sizeof(double));
  double *v_find_plunge = (double*) malloc(1*sizeof(double));
  double *e_find_plunge = (double*) malloc(1*sizeof(double));
  double *M_find_plunge = (double*) malloc(1*sizeof(double));
  double *S_find_plunge = (double*) malloc(1*sizeof(double));

  *t_find_plunge = (*i_plunge-i0)*timestep;

  Interp(t_in,e_in,j_end-j0+1,t_find_plunge,e_find_plunge,1);
  Interp(t_in,v_in,j_end-j0+1,t_find_plunge,v_find_plunge,1);
  Interp(t_in,M_in,j_end-j0+1,t_find_plunge,M_find_plunge,1);
  Interp(t_in,S_in,j_end-j0+1,t_find_plunge,S_find_plunge,1);

  //for (int i=0; i<(j_end-j0+1); i++) printf("%e %e %e %e %e\n", t_in[i], v_in[i], e_in[i], M_in[i], S_in[i]);
  //printf("%d %d %d %e %e %e, %e, %e\n", j_end, j_start, j0, t_end, v_AK[j0], e_AK[j0], e_AK[j0+1], e_AK[j0-1]);

  *i_buffer=(int)(10./(OmegaPhi(*v_find_plunge,
                                *e_find_plunge,
                                coslam,
                                *S_find_plunge,*M_find_plunge)/2./M_PI)/timestep)+1; // 10 orbits after plunge

  *t_clip = (*i_plunge + *i_buffer - i0)*timestep;

  if (*t_clip > t_in[j_end-j0]) *t_clip = t_in[j_end-j0];
  *interp_timestep = dt_large;

  *init_length = j_end - j0 + 1;

  *i_plunge = *i_plunge - i0;

  //printf("%e %d %d %d %e %e %e\n", *t_find_plunge, *i_plunge, i0, *i_buffer, timestep, *interp_timestep, *t_clip);

  free(t_find_plunge);
  free(e_find_plunge);
  free(v_find_plunge);
  free(M_find_plunge);
  free(S_find_plunge);



  // we need t,v,M,e,S vectors. We need t_max = timestep*(*i_plunge - i0 + *i_buffer)

  /*
  printf("i0: %d vlength: %d\n", i0, vlength);
  gim[i0]=gim0;
  Phi[i0]=Phi0;
  alp[i0]=alp0;
  for(int i=i0;i<vlength;i++){
    gimdotm=gimdot;
    Phidotm=Phidot;
    alpdotm=alpdot;
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    if(i==i0){
      gimdotm=gimdot;
      Phidotm=Phidot;
      alpdotm=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    gim[i+1]=gim[i]+(1.5*gimdot-.5*gimdotm)*timestep;
    Phi[i+1]=Phi[i]+(1.5*Phidot-.5*Phidotm)*timestep;
    alp[i+1]=alp[i]+(1.5*alpdot-.5*alpdotm)*timestep;
    if(i>*i_plunge+*i_buffer){
      gim[i+1]=gim[*i_plunge];
      Phi[i+1]=Phi[*i_plunge];
      alp[i+1]=alp[*i_plunge];
    }
  }
  nu[vlength]=drdm(v[vlength],e[vlength],coslam,S[vlength])/dtdm(v[vlength],e[vlength],coslam,S[vlength])/(2.*M_PI*M[vlength]);
  gimdotvec[vlength]=(dthetadm(v[vlength],e[vlength],coslam,S[vlength])-drdm(v[vlength],e[vlength],coslam,S[vlength]))/dtdm(v[vlength],e[vlength],coslam,S[vlength])/M[vlength];
  // ----------

  // ----- evolve phases from t_0 to t_start -----
  for(int i=i0;i>0;i--){
    gimdotp=gimdot;
    Phidotp=Phidot;
    alpdotp=alpdot;
    gimdot=(dthetadm(v[i],e[i],coslam,S[i])-drdm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    Phidot=drdm(v[i],e[i],coslam,S[i])/dtdm(v[i],e[i],coslam,S[i])/M[i];
    alpdot=(dphidm(v[i],e[i],coslam,S[i])-dthetadm(v[i],e[i],coslam,S[i]))/dtdm(v[i],e[i],coslam,S[i])/M[i];
    if(i==i0){
      gimdotp=gimdot;
      Phidotp=Phidot;
      alpdotp=alpdot;
    }
    nu[i]=Phidot/2./M_PI;
    gimdotvec[i]=gimdot;
    gim[i-1]=gim[i]-(1.5*gimdot-.5*gimdotp)*timestep;
    Phi[i-1]=Phi[i]-(1.5*Phidot-.5*Phidotp)*timestep;
    alp[i-1]=alp[i]-(1.5*alpdot-.5*alpdotp)*timestep;
  }
  nu[0]=drdm(v[0],e[0],coslam,S[0])/dtdm(v[0],e[0],coslam,S[0])/(2.*M_PI*M[0]);
  gimdotvec[0]=(dthetadm(v[0],e[0],coslam,S[0])-drdm(v[0],e[0],coslam,S[0]))/dtdm(v[0],e[0],coslam,S[0])/M[0];
  // ----------

  */
  free(e_AK);
  free(v_AK);
  free(t_fit);
  free(e_fit);
  free(v_fit);
  free(M_fit);
  free(S_fit);
  free(e_diff);
  free(v_diff);
  free(e_coeff);
  free(v_coeff);
  free(M_coeff);
  free(S_coeff);
  /*free(t_in);
  free(e_in);
  free(v_in);
  free(M_in);
  free(S_in);*/
  return;

}
