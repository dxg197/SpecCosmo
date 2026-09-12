/*@@ Calculates the Tmunu and the RHS of MHD equations @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

#define TINY 1.0e-35
#define LIGHT 2.99792458e8
#define SIE 1.11265006e-17
#define Rho_Min 1.0e-47 // or 1.0e-45
#define Epp_Min 1.0e-13 // or 1.0e-10

void MHD_RHS(CCTK_ARGUMENTS);
void MHD_Tmunu(CCTK_ARGUMENTS);
void MHD_Boundaries(CCTK_ARGUMENTS);
void Reconstruction(CCTK_ARGUMENTS);
void No_Reconstruction(CCTK_ARGUMENTS);
void SpecDeriv_Scalar_Derivative( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void Clean_Read_Scalar(CCTK_ARGUMENTS, CCTK_REAL *scalar); 
void ludcmp( CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *dd );
void lubksb( CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *bb );
int sgn(CCTK_REAL v);
void MC(CCTK_REAL *aaq,CCTK_REAL *bbq,CCTK_REAL *delta);
void monotonize(CCTK_REAL U, CCTK_REAL *Ur, CCTK_REAL *Ul);
CCTK_REAL determinant(CCTK_REAL a[4][4], CCTK_REAL k);
void cofactor(CCTK_REAL num[4][4], CCTK_REAL inverse[4][4], CCTK_REAL f);
void transpose(CCTK_REAL num[4][4], CCTK_REAL fac[4][4], CCTK_REAL inverse[4][4], CCTK_REAL r);
void qmax(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_max);
void qmin(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_min); 

/* Begin MHD_RHS */

void MHD_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
/*  Declare local variables */

  CCTK_INT i,j,k,m,n,p,four,handle;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez, index, ierr;
  CCTK_REAL ss,dtg,divb_temp,rho_star_dot,tau_dot;
  CCTK_REAL S_x_dot,S_y_dot,S_z_dot,beta_2;
  CCTK_REAL bb_x_dot,bb_y_dot,bb_z_dot,conf,sqrtdtg;
  CCTK_REAL gTcrossgNx, gTcrossgNy, gTcrossgNz, aa_ratio;
  CCTK_REAL RelativisticScale, geomCharge, volume;
  CCTK_REAL RelativisticBCorrX, RelativisticBCorrY, RelativisticBCorrZ;
  CCTK_REAL xmin, xmax, ymin, ymax, zmin, zmax, conductivity;
  
/* Declare Arrays */

  CCTK_REAL ***Flux_rho,***Flux_tau,****Flux_SS,****Flux_BB,****bb, **psib;
  CCTK_REAL gg4[4][4][4], bssn_dkij[4][4][4], T_mhd[4][4];
  CCTK_REAL kk3[4][4], igg3[4][4], detg[4], bssn_dk[4];
  CCTK_REAL **gradRho, **gradPress;
  CCTK_REAL lapse[4], bssn_beta[4][4];
  

/* Set up shorthands */

    sizex  = cctk_lsh[0]; 
    sizey  = cctk_lsh[1];
    sizez  = cctk_lsh[2];
    istart = cctk_nghostzones[0];
	jstart = cctk_nghostzones[1];
	kstart = cctk_nghostzones[2];	
	iend   = cctk_lsh[0] - cctk_nghostzones[0];
	jend   = cctk_lsh[1] - cctk_nghostzones[1];
	kend   = cctk_lsh[2] - cctk_nghostzones[2];
    
    four = 4;
    
    psib = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    Flux_rho = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    Flux_tau = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    Flux_SS = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    Flux_BB = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    bb = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    gradRho = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    gradPress = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
                  
    for(m=0; m < 4; m++) {
       psib[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       Flux_rho[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       Flux_tau[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       Flux_SS[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       Flux_BB[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       bb[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       gradRho[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       gradPress[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
      
         for(n=0; n < 4; n++) {
            Flux_rho[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
            Flux_tau[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
            bb[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
            Flux_SS[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
            Flux_BB[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
      
               for(p=0; p < 4; p++) {
                  Flux_SS[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
                  Flux_BB[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
                  bb[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
               }
         }
    }
    
    handle = CCTK_ReductionHandle("maximum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmax, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymax, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmax, 1, CCTK_VarIndex("grid::z"));
	
	handle = CCTK_ReductionHandle("minimum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmin, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymin, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmin, 1, CCTK_VarIndex("grid::z")); 
	
	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&aa_ratio, 1, CCTK_VarIndex("MHD_Analysis::aa_ratio_avg"));
	
	if (aa_ratio == 0.0) {aa_ratio = 1.0;}
	
    volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
	
	if (reconstruction) {
		Reconstruction(CCTK_PASS_CTOC); 
	}
	
	if (!reconstruction) {
		No_Reconstruction(CCTK_PASS_CTOC);
	}
	
	if (clean_scalars) {
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_x); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_y); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_z); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_x); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_y); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_z);
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_xx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_yx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_xz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_xy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_yy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_zy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_zx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_yz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SS_zz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_xy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_yx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_xz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_zx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_yz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_BB_zy); 
	}
	

	for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
/* Keep it Physical */

		if ((rho_star[index] < 0.0) || (isnan(rho_star[index]))) {
		    rho_star[index] = 0.0;
		}
		
		if(isnan(F_rho_x[index])) {F_rho_x[index] = 0.0;} 
		if(isnan(F_rho_y[index])) {F_rho_y[index] = 0.0;} 
		if(isnan(F_rho_z[index])) {F_rho_z[index] = 0.0;}
		if(isnan(F_tau_x[index])) {F_tau_x[index] = 0.0;} 
		if(isnan(F_tau_y[index])) {F_tau_y[index] = 0.0;} 
		if(isnan(F_tau_z[index])) {F_tau_z[index] = 0.0;} 
		if(isnan(F_SS_xx[index])) {F_SS_xx[index] = 0.0;} 
		if(isnan(F_SS_yy[index])) {F_SS_yy[index] = 0.0;} 
		if(isnan(F_SS_zz[index])) {F_SS_zz[index] = 0.0;}
		if(isnan(F_SS_xy[index])) {F_SS_xy[index] = 0.0;} 
		if(isnan(F_SS_xz[index])) {F_SS_xz[index] = 0.0;} 
		if(isnan(F_SS_yz[index])) {F_SS_yz[index] = 0.0;}
		if(isnan(F_SS_yx[index])) {F_SS_yx[index] = 0.0;} 
		if(isnan(F_SS_zx[index])) {F_SS_zx[index] = 0.0;} 
		if(isnan(F_SS_zy[index])) {F_SS_zy[index] = 0.0;}
		if(isnan(F_BB_xy[index])) {F_BB_xy[index] = 0.0;} 
		if(isnan(F_BB_xz[index])) {F_BB_xz[index] = 0.0;} 
		if(isnan(F_BB_yz[index])) {F_BB_yz[index] = 0.0;}
		if(isnan(F_BB_yx[index])) {F_BB_yx[index] = 0.0;} 
		if(isnan(F_BB_zx[index])) {F_BB_zx[index] = 0.0;} 
		if(isnan(F_BB_zy[index])) {F_BB_zy[index] = 0.0;}
			
/* Define Arrays */

		Flux_rho[0][1][index] = F_rho_x[index];
        Flux_rho[0][2][index] = F_rho_y[index];
        Flux_rho[0][3][index] = F_rho_z[index];
        
        Flux_tau[0][1][index] = F_tau_x[index];
        Flux_tau[0][2][index] = F_tau_y[index];
        Flux_tau[0][3][index] = F_tau_z[index];

        Flux_SS[0][1][1][index] = F_SS_xx[index];
        Flux_SS[0][2][1][index] = F_SS_yx[index];
        Flux_SS[0][1][3][index] = F_SS_xz[index];
        Flux_SS[0][1][2][index] = F_SS_xy[index];
        Flux_SS[0][2][2][index] = F_SS_yy[index];
        Flux_SS[0][3][2][index] = F_SS_zy[index];
        Flux_SS[0][3][1][index] = F_SS_zx[index];
        Flux_SS[0][2][3][index] = F_SS_yz[index];
        Flux_SS[0][3][3][index] = F_SS_zz[index];
        
        Flux_BB[0][1][2][index] = F_BB_xy[index];
        Flux_BB[0][2][1][index] = F_BB_yx[index];
        Flux_BB[0][1][3][index] = F_BB_xz[index];
        Flux_BB[0][3][1][index] = F_BB_zx[index];
        Flux_BB[0][2][3][index] = F_BB_yz[index];
        Flux_BB[0][3][2][index] = F_BB_zy[index];
        Flux_BB[0][1][1][index] = 0.0;
        Flux_BB[0][2][2][index] = 0.0;
        Flux_BB[0][3][3][index] = 0.0;
        
        bb[0][0][1][index] = bb_x[index];
        bb[0][0][2][index] = bb_y[index];
        bb[0][0][3][index] = bb_z[index];
        
        psib[0][index] = psib0[index];
     	  
/* Finish Inital Calculations */
		    
		  gradPress[0][index] = pp[index];
		  
		  gradRho[0][index] = rho[index];
        
        }
	}
}

/* Calc derivatives */

  if (CCTK_Equals(diff,"Spectral")) {
 
  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, Flux_rho); 
  
  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, Flux_tau); 
  
  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, Flux_SS);
  
  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, Flux_BB);
   
  SpecDeriv_Vector_Derivative2(CCTK_PASS_CTOC, bb); 
   
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, psib); 
  
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, gradPress); 
  
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, gradRho); 
  
  }
  
  if (CCTK_Equals(diff,"Spectral_1d")) {
 
  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, Flux_rho); 
  
  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, Flux_tau); 
  
  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, Flux_SS);
  
  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, Flux_BB);
   
  SpecDeriv_Vector_Derivative2_1d(CCTK_PASS_CTOC, bb); 
   
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, psib); 
  
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, gradPress); 
  
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, gradRho); 
  
  }
  
  if (CCTK_Equals(diff,"Finite")) {
 
  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, Flux_rho); 
  
  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, Flux_tau); 
  
  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, Flux_SS);
  
  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, Flux_BB);
   
  SpecDeriv_Vector_Derivative2_FD(CCTK_PASS_CTOC, bb); 
   
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, psib); 
  
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, gradPress); 
  
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, gradRho);
  
  }
  
  if (CCTK_Equals(diff,"Finite4")) {
 
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, Flux_rho); 
  
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, Flux_tau); 
  
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, Flux_SS);
  
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, Flux_BB);
   
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, bb); 
   
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, psib); 
  
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, gradPress); 
  
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, gradRho); 
  
  }
  
  if (CCTK_Equals(diff,"SCR3")) {
 
  SpecDeriv_Vector_Derivative_SCR3(CCTK_PASS_CTOC, Flux_rho); 
  
  SpecDeriv_Vector_Derivative_SCR3(CCTK_PASS_CTOC, Flux_tau); 
  
  SpecDeriv_Tensor_Derivative_SCR3(CCTK_PASS_CTOC, Flux_SS);
  
  SpecDeriv_Tensor_Derivative_SCR3(CCTK_PASS_CTOC, Flux_BB);
   
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, bb); 
   
  SpecDeriv_Scalar_Derivative_SCR3(CCTK_PASS_CTOC, psib);
  
  SpecDeriv_Scalar_Derivative_SCR3(CCTK_PASS_CTOC, gradPress); 
  
  SpecDeriv_Scalar_Derivative_SCR3(CCTK_PASS_CTOC, gradRho); 
  
  }
  
/* Here we used the charge per unit mass for the particles that existed in the epoch
relative to the electron's charge per unit mass.  We also used the Plank charge per unit
mass to convert all these into relative dimensionless units.  0.05122124 is the charge 
per unit mass of all matter around the time of the EW Phase Transition relative to the 
electron charge per unit mass.  0.08847884 is the charge per unit mass of all matter 
around the time of the QCD Phase Transition relative to the electron charge per unit mass.  
7.236e21 is the ratio of the electron's charge per unit mass to the Plank charge per
unit mass. 1 + epsilon corrects for relativistic changes to mass density */

if ((CCTK_Equals(Biermann_mode,"EW")) || (CCTK_Equals(Biermann_mode,"none"))) {
	geomCharge = (0.05122124) * charge * (7.236e21/(1.0 + epsilon[index]));
}

if (CCTK_Equals(Biermann_mode,"QCD")) {
	geomCharge = (0.08847884) * charge * (7.236e21/(1.0 + epsilon[index]));
}

if (CCTK_Equals(Biermann_mode,"Nuc")) {
	geomCharge = (0.79209766) * charge * (7.236e21/(1.0 + epsilon[index]));
}

 for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );

/* Initialize arrays */  
 
        conf = exp(4.0*phi[index]);
        
        lapse[0] = alpha[index];
        lapse[1] = alpha_x[index];
        lapse[2] = alpha_y[index];
        lapse[3] = alpha_z[index];
        
        bssn_beta[0][1] = betax[index];
        bssn_beta[0][2] = betay[index];
        bssn_beta[0][3] = betaz[index];
        bssn_beta[1][1] = betax1[index];
        bssn_beta[1][2] = betay1[index];
        bssn_beta[1][3] = betaz1[index];
        bssn_beta[2][1] = betax2[index];
        bssn_beta[2][2] = betay2[index];
        bssn_beta[2][3] = betaz2[index];
        bssn_beta[3][1] = betax3[index];
        bssn_beta[3][2] = betay3[index];
        bssn_beta[3][3] = betaz3[index];
        
        bssn_dk[1] = bssn_dk_x[index];
        bssn_dk[2] = bssn_dk_y[index];
        bssn_dk[3] = bssn_dk_z[index];    
        
        kk3[1][1] = kxx[index];
        kk3[2][2] = kyy[index];
        kk3[3][3] = kzz[index];
        kk3[1][2] = kxy[index];
        kk3[1][3] = kxz[index];
        kk3[2][3] = kyz[index];
        kk3[2][1] = kk3[1][2];
        kk3[3][1] = kk3[1][3];
        kk3[3][2] = kk3[2][3];
        
        bssn_dkij[1][1][1] = bssn_dkij_xxx[index];
        bssn_dkij[1][2][2] = bssn_dkij_xyy[index];
        bssn_dkij[1][3][3] = bssn_dkij_xzz[index];
        bssn_dkij[1][1][2] = bssn_dkij_xxy[index];
        bssn_dkij[1][1][3] = bssn_dkij_xxz[index];
        bssn_dkij[1][2][3] = bssn_dkij_xyz[index];
        bssn_dkij[1][2][1] = bssn_dkij_xxy[index];
        bssn_dkij[1][3][1] = bssn_dkij_xxz[index];
        bssn_dkij[1][3][2] = bssn_dkij_xyz[index];
        
        bssn_dkij[2][1][1] = bssn_dkij_yxx[index];
        bssn_dkij[2][2][2] = bssn_dkij_yyy[index];
        bssn_dkij[2][3][3] = bssn_dkij_yzz[index];
        bssn_dkij[2][1][2] = bssn_dkij_yxy[index];
        bssn_dkij[2][1][3] = bssn_dkij_yxz[index];
        bssn_dkij[2][2][3] = bssn_dkij_yyz[index];
        bssn_dkij[2][2][1] = bssn_dkij_yxy[index];
        bssn_dkij[2][3][1] = bssn_dkij_yxz[index];
        bssn_dkij[2][3][2] = bssn_dkij_yyz[index];
        
        bssn_dkij[3][1][1] = bssn_dkij_zxx[index];
        bssn_dkij[3][2][2] = bssn_dkij_zyy[index];
        bssn_dkij[3][3][3] = bssn_dkij_zzz[index];
        bssn_dkij[3][1][2] = bssn_dkij_zxy[index];
        bssn_dkij[3][1][3] = bssn_dkij_zxz[index];
        bssn_dkij[3][2][3] = bssn_dkij_zyz[index];
        bssn_dkij[3][2][1] = bssn_dkij_zxy[index];
        bssn_dkij[3][3][1] = bssn_dkij_zxz[index];
        bssn_dkij[3][3][2] = bssn_dkij_zyz[index];
  	            
        gg4[0][0][0] = gtt[index];
        gg4[0][0][1] = gtx[index];
        gg4[0][0][2] = gty[index];
        gg4[0][0][3] = gtz[index];
        gg4[0][1][0] = gg4[0][0][1];
        gg4[0][2][0] = gg4[0][0][2];
        gg4[0][3][0] = gg4[0][0][3];
        gg4[0][1][1] = gxx[index];
        gg4[0][2][2] = gyy[index];
        gg4[0][3][3] = gzz[index];
        gg4[0][1][2] = gxy[index];
        gg4[0][1][3] = gxz[index];
        gg4[0][2][3] = gyz[index];
        gg4[0][2][1] = gxy[index];
        gg4[0][3][1] = gxz[index];
        gg4[0][3][2] = gyz[index];

        for(m=1; m < 4; m++) {
           for(n=1; n < 4; n++) {
              for(p=1; p < 4; p++) {
                  gg4[m][n][p] = 4.0*bssn_dk[m]*gg4[0][n][p] + conf*bssn_dkij[m][n][p];
              }
           }
        }
        
        detg[0] = gxx[index] * gyy[index] * gzz[index] 
                + 2.0*gxy[index] * gxz[index] * gyz[index] 
                - gxx[index] * pow(gyz[index],2.0)           
                - gyy[index] * pow(gxz[index],2.0)           
                - gzz[index] * pow(gxy[index],2.0);
        
        for(m=1; m < 4; m++) {  
           detg[m] = gg4[m][1][1]*gg4[0][2][2]*gg4[0][3][3]  +
                     gg4[0][1][1]*gg4[m][2][2]*gg4[0][3][3]  +
                     gg4[0][1][1]*gg4[0][2][2]*gg4[m][3][3]  +
                2.0*(gg4[m][1][2]*gg4[0][1][3]*gg4[0][2][3]  +
                     gg4[0][1][2]*gg4[m][1][3]*gg4[0][2][3]  +
                     gg4[0][1][2]*gg4[0][1][3]*gg4[m][2][3]) -
                     gg4[m][1][1]*gg4[0][2][3]*gg4[0][2][3]  -
                 2.0*gg4[0][1][1]*gg4[m][2][3]*gg4[0][2][3]  - 
                     gg4[m][2][2]*gg4[0][1][3]*gg4[0][1][3]  -
                 2.0*gg4[0][2][2]*gg4[m][1][3]*gg4[0][1][3]  -
                     gg4[m][3][3]*gg4[0][1][2]*gg4[0][1][2]  -
                 2.0*gg4[0][3][3]*gg4[m][1][2]*gg4[0][1][2];
        }    
        
        igg3[1][1] = igxx[index];
        igg3[2][2] = igyy[index];
        igg3[3][3] = igzz[index];
        igg3[1][2] = igxy[index];
        igg3[1][3] = igxz[index];
        igg3[2][3] = igyz[index];
        igg3[2][1] = igg3[1][2];
        igg3[3][1] = igg3[1][3];
        igg3[3][2] = igg3[2][3];
        
        sqrtdtg = sqrtdetg[index];
        
        for(m=1; m < 4; m++) {
           gg4[m][0][0] = -2.0*lapse[0]*lapse[m];
           gg4[m][0][1] = gg4[m][1][1]*bssn_beta[0][1]+gg4[0][1][1]*bssn_beta[m][1]+
           				  gg4[m][1][2]*bssn_beta[0][2]+gg4[0][1][2]*bssn_beta[m][2]+
           				  gg4[m][1][3]*bssn_beta[0][3]+gg4[0][1][3]*bssn_beta[m][3];
           gg4[m][0][2] = gg4[m][1][2]*bssn_beta[0][1]+gg4[0][1][2]*bssn_beta[m][1]+
           				  gg4[m][2][2]*bssn_beta[0][2]+gg4[0][2][2]*bssn_beta[m][2]+
           				  gg4[m][2][3]*bssn_beta[0][3]+gg4[0][2][3]*bssn_beta[m][3];
           gg4[m][0][3] = gg4[m][1][3]*bssn_beta[0][1]+gg4[0][1][3]*bssn_beta[m][1]+
           				  gg4[m][2][3]*bssn_beta[0][2]+gg4[0][2][3]*bssn_beta[m][2]+
           				  gg4[m][3][3]*bssn_beta[0][3]+gg4[0][3][3]*bssn_beta[m][3];
           gg4[m][1][0] = gg4[m][0][1];
           gg4[m][2][0] = gg4[m][0][2];
           gg4[m][3][0] = gg4[m][0][3];
        }
        
        for(m=1; m < 4; m++) {
           for(n=1; n < 4; n++) {
              for(p=1; p < 4; p++) {
                 gg4[m][0][0] += gg4[m][n][p]*bssn_beta[0][n]*bssn_beta[0][p] +
                                 gg4[0][n][p]*bssn_beta[m][n]*bssn_beta[0][p] +
                                 gg4[0][n][p]*bssn_beta[0][n]*bssn_beta[m][p];
              }
           }
        }
        
        T_mhd[0][0] = Ttt_mhd[index];
        T_mhd[0][1] = Ttx_mhd[index];
        T_mhd[0][2] = Tty_mhd[index];
        T_mhd[0][3] = Ttz_mhd[index];
        T_mhd[1][0] = Ttx_mhd[index];
        T_mhd[2][0] = Tty_mhd[index];
        T_mhd[3][0] = Ttz_mhd[index];
        
        T_mhd[1][1] = Txx_mhd[index];
        T_mhd[2][2] = Tyy_mhd[index];
        T_mhd[3][3] = Tzz_mhd[index];
        T_mhd[1][2] = Txy_mhd[index];
        T_mhd[1][3] = Txz_mhd[index];
        T_mhd[2][3] = Tyz_mhd[index];
        T_mhd[2][1] = Txy_mhd[index];
        T_mhd[3][1] = Txz_mhd[index];
        T_mhd[3][2] = Tyz_mhd[index];
        
/* Initialize temp variables */
  
 	ss = 0.0;
	
	for(m=1; m < 4; m++) {
      for(n=1; n < 4; n++) {
	  ss += sqrtdtg*(T_mhd[0][0]*bssn_beta[0][m]*bssn_beta[0][n] + 
	        2.0*T_mhd[0][m]*bssn_beta[0][n] + T_mhd[m][n])*kk3[m][n]; 
	  }
	}
	
	for(m=1; m < 4; m++) {
	 ss += -sqrtdtg*(T_mhd[0][0]*bssn_beta[0][m] + T_mhd[0][m])*lapse[m];
    }
   	    
/* Calculate RHS */
        
        rho_star_dot = 0.0;
        tau_dot = 0.0;
        S_x_dot = 0.0;
        S_y_dot = 0.0;
        S_z_dot = 0.0;
        bb_x_dot = 0.0;
        bb_y_dot = 0.0;
        bb_z_dot = 0.0;
        divb_temp = 0.0;
        
        RelativisticScale = 0.0; 
        
        if ((rho[index] > 1.0e-30) && (geomCharge != 0.0)) { 
        	RelativisticScale = 1.0 / (geomCharge * rho[index] * rho[index]);
        }
        
        gTcrossgNx = (gradPress[2][index] * gradRho[3][index] - gradPress[3][index] * gradRho[2][index]);
        gTcrossgNy = (gradPress[3][index] * gradRho[1][index] - gradPress[1][index] * gradRho[3][index]);
        gTcrossgNz = (gradPress[1][index] * gradRho[2][index] - gradPress[2][index] * gradRho[1][index]);
        
        RelativisticBCorrX = RelativisticScale * gTcrossgNx;
        RelativisticBCorrY = RelativisticScale * gTcrossgNy;
        RelativisticBCorrZ = RelativisticScale * gTcrossgNz;
                
        for(m=1; m < 4; m++) {
        
         rho_star_dot += -Flux_rho[m][m][index];
         
         tau_dot += -Flux_tau[m][m][index];
         
         divb_temp += bb[0][m][m][index];
        }
        
        if (isnan(rho_star_dot)) {
    		CCTK_WARN(1, "NaN in rho_star_dot at point (i,j,k) — check floors");
		} else {
    		rho_star_dt[index] = rho_star_dot;
		}
        
        if (isnan(tau_dot)) {
    		CCTK_WARN(1, "NaN in tau_dot at point (i,j,k) — check floors");
		} else {
    		tau_dt[index] = tau_dot + ss; //*pow(aa_ratio,-1.5);
		}
        
        divb[index] = divb_temp;
        
        psib_dt[index] = -divb_temp*pow(ch,2) - pow((ch/cp),2)*psib[0][index];
        
        for(m=0; m < 4; m++) {
         for(n=0; n < 4; n++) {        
             S_x_dot += 0.5*sqrtdtg*T_mhd[m][n]*gg4[1][m][n];
             S_y_dot += 0.5*sqrtdtg*T_mhd[m][n]*gg4[2][m][n];
             S_z_dot += 0.5*sqrtdtg*T_mhd[m][n]*gg4[3][m][n];
         }
        }
        
        for(m=1; m < 4; m++) {   
          S_x_dot += -Flux_SS[m][m][1][index];
          S_y_dot += -Flux_SS[m][m][2][index];
          S_z_dot += -Flux_SS[m][m][3][index];
        }
        
        	S_x_dt[index] = S_x_dot;
        	S_y_dt[index] = S_y_dot;
        	S_z_dt[index] = S_z_dot;
        	
        if (TrK[index] >= TrK_max) {
		
		rho_star_dt[index] = 0.0;
        tau_dt[index] = 0.0;
        S_x_dt[index] = 0.0;
        S_y_dt[index] = 0.0;
        S_z_dt[index] = 0.0;
        
		}
        
        for(m=1; m < 4; m++) {
          bb_x_dot += -Flux_BB[m][m][1][index] - igg3[1][m]*psib[m][index];
          bb_y_dot += -Flux_BB[m][m][2][index] - igg3[2][m]*psib[m][index];
          bb_z_dot += -Flux_BB[m][m][3][index] - igg3[3][m]*psib[m][index];
        }
        
        if (isnan(RelativisticBCorrX)) {
    		CCTK_WARN(1, "NaN in RelativisticBCorrX at point (i,j,k) — check floors");
		} else {
    		bb_x_dot += RelativisticBCorrX;
		}
        
        if (isnan(RelativisticBCorrY)) {
    		CCTK_WARN(1, "NaN in RelativisticBCorrY at point (i,j,k) — check floors");
		} else {
    		bb_y_dot += RelativisticBCorrY;
		}
		
		if (isnan(RelativisticBCorrZ)) {
    		CCTK_WARN(1, "NaN in RelativisticBCorrZ at point (i,j,k) — check floors");
		} else {
    		bb_z_dot += RelativisticBCorrZ;
		}
		
        
        if (resistance) {
            if (temp_out[index] != 0.0) {
            conductivity = 1.0E11 * pow(( temp_out[index]),2.0);
              for(m=1; m < 4; m++) {
                  for(n=1; n < 4; n++) {
        	  bb_x_dot += (igg3[n][m]*bb[m][n][1][index])/conductivity;
        	  bb_y_dot += (igg3[n][m]*bb[m][n][2][index])/conductivity;
        	  bb_z_dot += (igg3[n][m]*bb[m][n][3][index])/conductivity;
        	     }
        	  }
           }  
        }
        
        if (isnan(bb_x_dot)) {
    		CCTK_WARN(1, "NaN in bb_x_dot at point (i,j,k) — check floors");
		} else {
    		bb_x_dt[index] = bb_x_dot;
		}

        if (isnan(bb_y_dot)) {
    		CCTK_WARN(1, "NaN in bb_y_dot at point (i,j,k) — check floors");
		} else {
    		bb_y_dt[index] = bb_y_dot;
		}

        if (isnan(bb_z_dot)) {
    		CCTK_WARN(1, "NaN in bb_z_dot at point (i,j,k) — check floors");
		} else {
    		bb_z_dt[index] = bb_z_dot;
		}
      
        }
      }
   }
    
    for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           free(Flux_SS[m][n][p]);
           free(Flux_BB[m][n][p]);
           free(bb[m][n][p]);
        }
      }
    }
    
    for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
        free(Flux_SS[m][n]);
        free(Flux_BB[m][n]);
        free(Flux_rho[m][n]);
        free(Flux_tau[m][n]);
        free(bb[m][n]);
      }
    }
    
    for(m=0; m < 4; m++) {
      free(Flux_SS[m]);
      free(Flux_BB[m]);
      free(Flux_rho[m]);
      free(Flux_tau[m]);
      free(bb[m]);
      free(psib[m]);
      free(gradRho[m]);
	  free(gradPress[m]);
    }
   
   free(Flux_rho);
   free(Flux_tau);
   free(Flux_SS);
   free(Flux_BB);
   free(bb);
   free(psib);
   free(gradRho);
   free(gradPress);
}
    
void MHD_Tmunu(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
/* Declare local variables */
  CCTK_INT i,j,k,m,n,q,four, itr, index, handle, ierr;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez;
  CCTK_REAL SS[4],omega,T_mhd[4][4],shift[4],pi, SSI[4], sum;
  CCTK_REAL u[4], beta[4], b[4], bu, u_l[4], b_l[4], Va2, kappa;
  CCTK_REAL ff[4], du0_u[4], dbeta0_u[4], dpp_epsilon, epp, shift2;
  CCTK_REAL err, uu2, hp, dh_epsilon, dh_u[4], dbeta2_u[4], dbeta_u[4][4];
  CCTK_REAL vv3[4], dd, gg[4][4], igg3[4][4], err_p, conf, dtg, sqrtdtg, bssn_dtg;
  CCTK_REAL I_n, deltal, Cs, QQ[4][4], lapse, dQ_epsilon, dQ_u[4], dpp_u[4];
  
/* Declare Arrays */
  CCTK_REAL ***vv, **df, *deltx;
  CCTK_INT *indx;

/* Set up shorthands */

    istart = cctk_nghostzones[0];
    jstart = cctk_nghostzones[1];
    kstart = cctk_nghostzones[2];
  
    iend = cctk_lsh[0] - cctk_nghostzones[0];
    jend = cctk_lsh[1] - cctk_nghostzones[1];
    kend = cctk_lsh[2] - cctk_nghostzones[2];

    sizex  = cctk_lsh[0]; 
    sizey  = cctk_lsh[1];
    sizez  = cctk_lsh[2];
    
    four = 4;
    
    pi = 4.0*atan(1.0);
    
    if (gam > 1.5) {
    	kappa = 1.954180e6;
    } else {
    	kappa = 2.093841e2;
    }
    
 /* Define Arrays */
 
    vv = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    
    for(m=0; m < 4; m++) {
       vv[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
         for(n=0; n < 4; n++) {
            vv[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
         }
    }
  	
  	deltx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*four);
  	indx = (CCTK_INT *)malloc(sizeof(CCTK_INT)*four);
  	df = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    for(m=0; m < 4; m++) {
       df[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*four);
    }
    
    for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
  
             index = CCTK_GFINDEX3D( cctkGH, i, j, k );
             
        if ((rho_star[index] < 0.0) || (isnan(rho_star[index]))) {
		    rho_star[index] = 0.0;
		}
		    
        lapse = alpha[index];
    	
        shift[1] = betax[index];
        shift[2] = betay[index];
        shift[3] = betaz[index];
			
		conf = exp(4.0*phi[index]);
  	        
        gg[0][0] = gtt[index];
        gg[0][1] = gtx[index];
        gg[0][2] = gty[index];
        gg[0][3] = gtz[index];
        gg[1][0] = gg[0][1];
        gg[2][0] = gg[0][2];
        gg[3][0] = gg[0][3];
        gg[1][1] = gxx[index];
        gg[2][2] = gyy[index];
        gg[3][3] = gzz[index];
        gg[1][2] = gxy[index];
        gg[1][3] = gxz[index];
        gg[2][3] = gyz[index];
        gg[2][1] = gxy[index];
        gg[3][1] = gxz[index];
        gg[3][2] = gyz[index];
			
		cofactor(gg,igg3,4);
        
        sqrtdtg = sqrtdetg[index];
        
        beta2[index] = 0.0;
        omega = 0.0;
    
/* Begin Tmunu Calculation */
 	
/* Calculate Temperal part of Tmunu for MHD */

        SS[1] = S_x[index];
  		SS[2] = S_y[index];
  		SS[3] = S_z[index];
  		
/* Calculate B-field */

		b[1] = alpha[index]*bb_x[index]/sqrtdtg;
		b[2] = alpha[index]*bb_y[index]/sqrtdtg;
		b[3] = alpha[index]*bb_z[index]/sqrtdtg;
		
		for (m=1; m<4; m++){
			if(isnan(b[m])) {b[m] = 0.0;}
		}
		
/* Use initial guess to refine solution */

		u[0] = uu_t[index];
        u[1] = uu_x[index];
        u[2] = uu_y[index];
        u[3] = uu_z[index];
        
        epp = epsilon[index];
        
        err = 0.0;
        err_p = 0.0;

/* Begin Newton Raphson */

	if ((rho[index] > 0.0) && (! isnan(rho[index]))) { 

	for(itr=0; itr < 1000; itr++) {
         
/* Re-Initiallize variables */
       
	uu2 = 0.0;
	bu = 0.0;
	beta2[index] = 0.0;
	
	for(m=0; m < 4; m++) {
	  u_l[m] = 0.0;
      b_l[m] = 0.0;
	  du0_u[m] = 0.0;
	  dbeta0_u[m] = 0.0;
	  dbeta2_u[m] = 0.0;
	  indx[m] = 0.0;
	  deltx[m] = 0.0;
	  beta[m] = 0.0;
	  SSI[m] = 0.0;
	}
    
    for(m=0; m < 4; m++) {
	  ff[m] = 0.0;
     for(n=0; n < 4; n++) {
       df[m][n] = 0.0;
	 }
	}
	
/* Calculate initial values for Primative Variables */ 	 

 	for(m=1; m < 4; m++) {
       for(n=1; n < 4; n++) {
 	      u_l[m] += gg[m][n]*u[n];	
 	      b_l[m] += gg[m][n]*b[n];	
 	   }
 	}
 
    for(m=1; m < 4; m++) {
      uu2 += u_l[m]*u[m];  
    }
 	
  	u[0] = sqrt(1.0 + fabs(uu2))/alpha[index];
  	
  	rho[index] = rho_star[index]/(sqrtdtg*u[0]);
  	
  	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) { 
  		pp[index] = (gam-1.0)*rho[index]*epp;
  	} else {
  		if (rho[index] > 0.0) {
  		pp[index] = kappa*pow(rho[index],gam);
  		epp = pp[index]/((gam - 1.0)*rho[index]);
  		} else {
  		pp[index] = 0.0;
  		epp = 0.0;
  		}
  	}
  	
  	if (rho[index] > 0.0) {
    hp = 1.0 + epp + pp[index]/rho[index] + trQ[index]/rho[index];   
    } else {
    hp = 1.0;
    }
  	
/* Other variables  */	
    
    for(m=1; m < 4; m++) {
	 bu += b_l[m]*u[m];
	}

  	beta[0] = bu/(sqrt(4.0*pi)*alpha[index]);

	for(m=1; m < 4; m++) {
  	 beta[m] = (b[m] + bu*u[m])/(sqrt(4.0*pi)*alpha[index]*u[0]);
	}
 
	for(m=0; m < 4; m++) {
       for(n=0; n < 4; n++) {
	     beta2[index] += gg[m][n]*beta[m]*beta[n];
	   }
	}

/* Calculate functions */
        
    if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {        
    ff[0] = alpha[index]*u[0]*rho_star[index] - rho_star[index] + gam*rho_star[index]*epp*alpha[index]*u[0] + 
            alpha[index]*sqrtdtg*beta2[index]*pow(u[0],2.0) + alpha[index]*sqrtdtg*(pp[index]+0.5*beta2[index])*igtt[index] - 
            alpha[index]*sqrtdtg*pow(beta[0],2.0) - tau[index];
    } else {
    if (rho[index] > 0.0) {
    ff[0] = pp[index] - kappa*pow(rho[index],gam);
    } else {
    ff[0] = 0.0;
    }
    }

  	for(m=1; m < 4; m++) {
       for(n=1; n < 4; n++) {
  	     SSI[m] += igg3[m][n]*SS[n]; 
  	   }
  	 }
  	 
  	 if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
  	 for(m=1; m < 4; m++) {
  	   ff[m] += rho_star[index]*hp*u[m] + sqrtdtg*beta2[index]*u[0]*u[m] - sqrtdtg*beta[0]*beta[m] - SSI[m]; 
     }
     } else {
       ff[m] += rho_star[index]*hp*u[m] - SSI[m]; 
     }
                
/* Calculate error */
	
	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {    
		if (tau[index] != 0.0) { 
	  		err = fabs(ff[0]/tau[index]);
		} else {
	  		err = fabs(ff[0]);
		}
	} else {
		if (pp[index] != 0.0) {
			err = fabs(ff[0]/pp[index]);
		} else {
	  		err = fabs(ff[0]);
		}
	}
	
	for(m=1; m < 4; m++) {
	  if (SSI[m] != 0.0) {
	    	err += fabs(ff[m]/SSI[m]);
	  } else {
	  		err += fabs(ff[m]);
	  }
	}
 
  	if ( err < sensitivity )  {
  	 break;
  	}
  	
/* Calculate derivatives */

	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
    	dpp_epsilon = (gam - 1.0)*rho[index];
    } else {
    	dpp_epsilon = 0.0;
    }
    
    for(m=1; m < 4; m++) {
		dpp_u[m] = 0.0;
	}

	for(m=1; m < 4; m++) {
      for(n=1; n < 4; n++) {
 	    du0_u[m] += gg[m][n]*u[n]/(2.0*alpha[index]*sqrt(1.0+fabs(uu2))); 
 	  }
 	}
 	
 	for(m=1; m < 4; m++) {
 	 dbeta0_u[m] = b_l[m]/(sqrt(4.0*pi)*alpha[index]);
 	}

	for(m=1; m < 4; m++) {
      for(n=1; n < 4; n++) {
 	  dbeta_u[m][n] = -(b[m]/(sqrt(4.0*pi)*alpha[index]) + beta[0]*u[m])*du0_u[n]/pow(u[0],2.0) + 
 	                  dbeta0_u[n]*u[m]/u[0];
 	  }
 	}
 	
 	for(m=1; m < 4; m++) {
 	 dbeta_u[m][m] += beta[0]/u[0]; 
 	}
 	
	for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           for(q=1; q < 4; q++) {
 	          dbeta2_u[m] += gg[q][n]*beta[n]*dbeta_u[q][m] + gg[q][n]*beta[q]*dbeta_u[n][m];
 	       }
 	    }
 	 }
 	 
 	 for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
 	          dbeta2_u[m] += gg[0][n]*beta[n]*dbeta0_u[m] + gg[0][n]*beta[0]*dbeta_u[n][m];
 	    }
 	 }
 	 
 	for(m=1; m < 4; m++) {
 	   dbeta2_u[m] += 2.0*gg[0][0]*beta[0]*dbeta0_u[m]; 
 	}
 	
 	if ( add_AV_bulk && (div_vv[index] < 0.0) && (! isnan(Speed_Sound[index])) && (! isnan(div_vv[index])) ) {
 		
 	deltal = fmin(CCTK_DELTA_SPACE(2),fmin(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(0)));
 	
 	dQ_epsilon = (rho_star[index] + sqrtdtg*uu_t[index]*dpp_epsilon)*alpha[index]/(sqrtdtg*pow(uu_t[index],nn)) *
 	             deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*div_vv[index];
 	
 	for(m=1; m < 4; m++) {
 		dQ_u[m] = (sqrtdtg*du0_u[m]*(pp[index] + beta2[index])*alpha[index]/(sqrtdtg*pow(uu_t[index],nn)) -
 	          	  nn*(rho_star[index] + rho_star[index]*epsilon[index] + sqrtdtg*uu_t[index]*(pp[index] + 
    	      	  beta2[index]))*alpha[index]/(sqrtdtg*pow(uu_t[index],nn+1))*du0_u[m])*deltal*(kq*deltal*
    	      	  div_vv[index]-kl*Speed_Sound[index])*div_vv[index];
	}
	
	dpp_epsilon += dQ_epsilon;
	
	for(m=1; m < 4; m++) {
		dpp_u[m] = dQ_u[m];
	}
    
    if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
 		dh_epsilon = gam + dQ_epsilon/rho[index]; 
 	} else {
 		if (rho[index] > 0.0) {
 		dh_epsilon = 1.0 + dQ_epsilon/rho[index]; 
 		} else {
 		dh_epsilon = 1.0;
 		}
 	}
 	
 	for(m=1; m < 4; m++) {
 		dh_u[m] = dQ_u[m]/rho[index];;
 	}
 	
 	} else {
 	
 	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
 		dh_epsilon = gam;
 	} else {
 		dh_epsilon = 1.0;
 	}
 	
 	for(m=1; m < 4; m++) {
 		dh_u[m] = 0.0;
 	}
 	
 	}
 	
/* Calculate the derivatives of functions */

     if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
  	   df[0][0] = gam*alpha[index]*u[0]*rho_star[index] + alpha[index]*sqrtdtg*dpp_epsilon*igtt[index];
  	 } else {
       df[0][0] = 0.0;
     }
  	
  	for(m=1; m < 4; m++) {
  	            
  	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
  	 df[0][m] = alpha[index]*rho_star[index]*du0_u[m] + gam*alpha[index]*epp*rho_star[index]*du0_u[m] +
  	            alpha[index]*sqrtdtg*pow(u[0],2.0)*dbeta2_u[m] + 2.0*alpha[index]*sqrtdtg*beta2[index]*u[0]*du0_u[m] + 
  	            alpha[index]*sqrtdtg*(dpp_u[m] + 0.5*dbeta2_u[m])*igtt[index] - 2.0*alpha[index]*sqrtdtg*beta[0]*dbeta0_u[m];
    } else {
  	 df[0][m] = 0.0;
  	}
  	
   	}
  	
  	for(m=1; m < 4; m++) {
  	 df[m][0] = rho_star[index]*u[m]*dh_epsilon;
    }
    
    if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {  
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
  	       df[m][n] = sqrtdtg*(-beta[0]*dbeta_u[m][n] + beta2[index]*u[m]*du0_u[n] + u[0]*u[m]*dbeta2_u[n] - 
  	                  beta[m]*dbeta0_u[n]) + rho_star[index]*u[m]*dh_u[n]; 
  	    }
  	 }
 
     for(m=1; m < 4; m++) {
         df[m][m] += rho_star[index]*hp + sqrtdtg*u[0]*beta2[index];
     }
     } else {
     for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
  	       df[m][n] = rho_star[index]*u[m]*dh_u[n]; 
  	    }
  	 }
  	 for(m=1; m < 4; m++) {
         df[m][m] += rho_star[index]*hp;
     }
  	 
     }
     
/* Find deltx terms */
	
	for(m=0; m < 4; m++) {
	 deltx[m] = -ff[m];
	}
    
    sum = 0.0;
    for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) {
    		sum += fabs(df[m][n]);    
        }
    }
    
    if (sum == 0.0) {
    	break;
    	}
    
    	ludcmp(df,indx,&dd);    
  	
  		lubksb(df,indx,deltx);
  	
/* Calculate new values for energy density and four-velocity */

	for(m=1; m < 4; m++) {
	   if ((! isnan(deltx[m])) && (rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {
	      u[m] += deltx[m];
	   }
  	}
  	
    if ((! isnan(deltx[0])) && (rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) {
	   epp += deltx[0];
	} 
	
  	err_p = err;

	}
	
	}
	
/* finish Calculation */
	
/* Re-Initiallize variables */

	uu2 = 0.0;
	bu = 0.0;
	beta2[index] = 0.0;	
	
/* Calculate final values for Primative Variables */

	for(m=1; m < 4; m++) {
	   if (fabs(u[m]) >= 1.0) {
	   		u[m] = sgn(u[m])*Speed_Sound[index]; // or 0.0?
	   }
	}
	
	for(m=1; m < 4; m++) {
		if (isnan(u[m])) {u[m] = 0.0;}
	}
 
 	for(m=1; m < 4; m++) {
     for(n=1; n < 4; n++) {
 	  uu2 += gg[m][n]*u[m]*u[n];
 	 }
 	}
 	
  	u[0] = sqrt(1.0 + fabs(uu2))/alpha[index];
  	if (isnan(u[0])) {u[0] = 1.0/alpha[index];}
  	
  	rho[index] = rho_star[index]/(sqrtdtg*u[0]);
  	
  	if ((rho[index] > Rho_Min) && (epp > Epp_Min) && (fabs(uu2) < 1.0)) { 
  		pp[index] = (gam-1.0)*rho[index]*epp;
  	} else {
  		pp[index] = kappa*pow(rho[index],gam);
  		epp = pp[index]/((gam - 1.0)*rho[index]);
  	}
  	
  	if ((rho[index] <= 0.0) || (isnan(rho[index]))) {
		rho[index] = 0.0;
		pp[index] = 0.0;
		epp = 0.0;
		u[0] = 1.0/alpha[index];
		u[1] = 0.0;
		u[2] = 0.0;
		u[3] = 0.0;
		uu2 = 0.0;
	}
	
/* Calculate beta2, omega and velocities */
       	
	for(m=1; m < 4; m++) {
     for(n=1; n < 4; n++) {
	  bu += gg[m][n]*b[m]*u[n];
	 }
	}

  	beta[0] = bu/(sqrt(4.0*pi)*alpha[index]);
  	if (isnan(beta[0])) {beta[0] = 0.0;}

	for(m=1; m < 4; m++) {
	 beta[m] = (b[m] + bu*u[m])/(sqrt(4.0*pi)*alpha[index]*u[0]);
  	 vv3[m] = u[m]/(sqrt(1.0 + fabs(uu2))) - shift[m];
  	 if (isnan(beta[m])) {beta[m] = 0.0;}
  	 if (isnan(vv3[m])) {vv3[m] = 0.0;}
	}
	
	
 
	for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
  	  beta2[index] += gg[m][n]*beta[m]*beta[n];
	 }
	} 
	if (isnan(beta2[index])) {beta2[index] = 0.0;}  
    
    epsilon[index] = epp;
	
	vv_x[index] = vv3[1];
	vv_y[index] = vv3[2];
	vv_z[index] = vv3[3];

	uu_t[index] = u[0];
	uu_x[index] = u[1];
	uu_y[index] = u[2];
	uu_z[index] = u[3];
	
	b_x[index] = b[1];
	b_y[index] = b[2];
	b_z[index] = b[3];
	
	beta_t[index] = beta[0];
	beta_x[index] = beta[1];
	beta_y[index] = beta[2];
	beta_z[index] = beta[3]; 
			
    vv[0][1][index] = vv3[1];
    vv[0][2][index] = vv3[2];
    vv[0][3][index] = vv3[3];
  
            }
 	    }
    } 
    
    if (CCTK_Equals(diff,"Spectral")) {
	
	SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, vv);
	
	}
	
	if (CCTK_Equals(diff,"Spectral_1d")) {
	
	SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, vv);
	
	}
	
	if (CCTK_Equals(diff,"Finite")) {
	
	SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, vv); 
	
	}
	
	if (CCTK_Equals(diff,"Finite4")) {
	
	SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, vv); 
	
	}
	
	if (CCTK_Equals(diff,"SCR3")) {
	
	SpecDeriv_Vector_Derivative_SCR3(CCTK_PASS_CTOC, vv);
	}

 /* Do the */
    
for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
 		
 		div_vv[index] = vv[1][1][index]+vv[2][2][index]+vv[3][3][index];
        
        lapse = alpha[index];
    	
        shift[1] = betax[index];
        shift[2] = betay[index];
        shift[3] = betaz[index];
			
		conf = exp(4.0*phi[index]);
  	        
        gg[0][0] = gtt[index];
        gg[0][1] = gtx[index];
        gg[0][2] = gty[index];
        gg[0][3] = gtz[index];
        gg[1][0] = gg[0][1];
        gg[2][0] = gg[0][2];
        gg[3][0] = gg[0][3];
        gg[1][1] = gxx[index];
        gg[2][2] = gyy[index];
        gg[3][3] = gzz[index];
        gg[1][2] = gxy[index];
        gg[1][3] = gxz[index];
        gg[2][3] = gyz[index];
        gg[2][1] = gxy[index];
        gg[3][1] = gxz[index];
        gg[3][2] = gyz[index];
                
        cofactor(gg,igg3,4);
        
        sqrtdtg = sqrtdetg[index];

/* Calculate Artificial Viscosity */

    deltal = fmin(CCTK_DELTA_SPACE(2),fmin(CCTK_DELTA_SPACE(1),CCTK_DELTA_SPACE(0)));
	
	b[1] = b_x[index];
	b[2] = b_y[index];
	b[3] = b_z[index];
	
	u[0] = uu_t[index];
	u[1] = uu_x[index];
	u[2] = uu_y[index];
	u[3] = uu_z[index];
	
	beta[0] = beta_t[index];
	beta[1] = beta_x[index];
	beta[2] = beta_y[index];
	beta[3] = beta_z[index];
	
	Qxx[index] = 0.0;
    Qyy[index] = 0.0;
    Qzz[index] = 0.0;
    Qxy[index] = 0.0;
    Qxz[index] = 0.0;
    Qyz[index] = 0.0;
    trQ[index] = 0.0;
	
	for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
  	  QQ[m][n] = 0.0;
	 }
	}  
 		
    I_n = (rho_star[index] + rho_star[index]*epsilon[index] + uu_t[index]*sqrtdtg*(pp[index] + 
    	  beta2[index]))*pow(1.0/uu_t[index],nn)/sqrt(dtg);
    
    Va2 = beta2[index]/(rho[index] + rho[index]*epsilon[index] + pp[index] + beta2[index]);
    
    if ((Va2 >= 0.0) && !isnan(Va2)) {
    	Speed_Alfven[index] = pow(Va2,0.5);
    } else {
    	Speed_Alfven[index] = 0.0;
    }
    
    if (rho[index] > 0.0) {
    Cs = gam*pp[index]/(rho[index] + rho[index]*epsilon[index] + pp[index]);
    } else {
    Cs = 0.0;
    }
    
    if ((Cs > 0.0) && !isnan(Cs)) {
    	Speed_Sound[index] = pow(Cs,0.5);
    } else {
    	Speed_Sound[index] = 0.0;
    }
    
    if ( add_AV_bulk && (! isnan(Speed_Sound[index])) && (! isnan(div_vv[index])) && (rho_out_frac[index] >= rho_frac_max)  && (div_vv[index] < 0.0)) { 
   
    //printf("Applying Bulk Artificial Viscosity\n");
     
    trQ[index] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*div_vv[index];
    
    }
    
    if ( add_AV_shear && (! isnan(Speed_Sound[index])) && (! isnan(div_vv[index])) && (rho_out_frac[index] >= rho_frac_max)) { 
          
    QQ[1][1] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*(vv[1][1][index] - div_vv[index]/3.0);
                
    QQ[2][2] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*(vv[2][2][index] - div_vv[index]/3.0);
                 
    QQ[3][3] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*(vv[3][3][index] - div_vv[index]/3.0);
                 
    QQ[1][2] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*0.5*(vv[1][2][index] + vv[2][1][index]);
                 
    QQ[1][3] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*0.5*(vv[1][3][index] + vv[3][1][index]);
                
    QQ[2][3] = I_n*deltal*(kq*deltal*div_vv[index]-kl*Speed_Sound[index])*0.5*(vv[2][3][index] + vv[3][2][index]);
    
    QQ[2][1] = QQ[1][2];
    
    QQ[3][1] = QQ[1][3];
    
    QQ[3][2] = QQ[2][3];
    
    if ( neg_div && (div_vv[index] > 0.0) ) {
    
    for(m=0; m < 4; m++) {
     	for(n=0; n < 4; n++) {
  	  		QQ[m][n] = 0.0;
	 	}
	}  
    
    }
    
    //if (QQ[1][1] != 0.0) {printf("Applying Shear Artificial Viscosity\n");}
    
    for(m=1; m < 4; m++) {
       	Qxx[index] += QQ[1][m]*igg3[m][1];
       	Qyy[index] += QQ[2][m]*igg3[m][2];
       	Qzz[index] += QQ[3][m]*igg3[m][3];
       	Qxy[index] += QQ[1][m]*igg3[m][2];
       	Qxz[index] += QQ[1][m]*igg3[m][3];
       	Qyz[index] += QQ[2][m]*igg3[m][3];
    }
    
    }
	
	omega = rho[index] + pp[index] + rho[index]*epsilon[index] + trQ[index];
		
/* Calculate the Tmunu for MHD */
  		
  		Ttt_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_t[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[0][0] - 
  						 beta[0]*beta[0];
  		
  		Ttx_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_x[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[0][1] - 
  						 beta[0]*beta[1];
  		Tty_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_y[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[0][2] - 
  						 beta[0]*beta[2];
  		Ttz_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_z[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[0][3] - 
  						 beta[0]*beta[3];
  		
        Txx_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_x[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[1][1] - 
                         beta[1]*beta[1] + Qxx[index];
  		Tyy_mhd[index] = (omega+beta2[index])*uu_y[index]*uu_y[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[2][2] - 
  		                 beta[2]*beta[2] + Qyy[index];
  		Tzz_mhd[index] = (omega+beta2[index])*uu_z[index]*uu_z[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[3][3] - 
  		                 beta[3]*beta[3] + Qzz[index];
  		Txy_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_y[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[1][2] -  
  		                 beta[1]*beta[2] + Qxy[index];
  		Txz_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_z[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[1][3] - 
  		                 beta[1]*beta[3] + Qxz[index];
  		Tyz_mhd[index] = (omega+beta2[index])*uu_y[index]*uu_z[index] + (pp[index]+trQ[index]+0.5*beta2[index])*igg3[2][3] - 
  		                 beta[2]*beta[3] + Qyz[index];
  		
  	if (rho[index] <= 0.0) {
  	
    	Ttt_mhd[index] = beta2[index]/pow(alpha[index],2.0) + 0.5*beta2[index]*igg3[0][0];
    	Ttx_mhd[index] = 0.5*beta2[index]*igg3[0][1];
  		Tty_mhd[index] = 0.5*beta2[index]*igg3[0][2];
  		Ttz_mhd[index] = 0.5*beta2[index]*igg3[0][3];
    	Txx_mhd[index] = 0.5*beta2[index]*igg3[1][1] - beta[1]*beta[1];
  		Tyy_mhd[index] = 0.5*beta2[index]*igg3[2][2] - beta[2]*beta[2];
  		Tzz_mhd[index] = 0.5*beta2[index]*igg3[3][3] - beta[3]*beta[3];
  		Txy_mhd[index] = 0.5*beta2[index]*igg3[1][2] - beta[1]*beta[2];
  		Txz_mhd[index] = 0.5*beta2[index]*igg3[1][3] - beta[1]*beta[3];
  		Tyz_mhd[index] = 0.5*beta2[index]*igg3[2][3] - beta[2]*beta[3];
  		
   }
   
   if (isnan(Ttt_mhd[index])) {Ttt_mhd[index] = 0.0;}
   if (isnan(Ttx_mhd[index])) {Ttx_mhd[index] = 0.0;}
   if (isnan(Tty_mhd[index])) {Tty_mhd[index] = 0.0;}
   if (isnan(Ttz_mhd[index])) {Ttz_mhd[index] = 0.0;}
   if (isnan(Txx_mhd[index])) {Txx_mhd[index] = 0.0;}
   if (isnan(Tyy_mhd[index])) {Tyy_mhd[index] = 0.0;}
   if (isnan(Tzz_mhd[index])) {Tzz_mhd[index] = 0.0;}
   if (isnan(Txy_mhd[index])) {Txy_mhd[index] = 0.0;}
   if (isnan(Txz_mhd[index])) {Txz_mhd[index] = 0.0;}
   if (isnan(Tyz_mhd[index])) {Tyz_mhd[index] = 0.0;}
  
            }
 	    }
    }  //printf("Completed Tuv Calculation\n");

    for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
         free(vv[m][n]); 
      }
    }
    
    for(m=0; m < 4; m++) {
      free(vv[m]);
      free(df[m]);
    }
    
    free(vv);
    free(df); 
    free(indx);
    free(deltx);
}

/* Calculate Boundary Conditions */

void MHD_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr=-1;
  int sw[3];

#ifdef DEBUG_MOL
  printf("We've reached the Boundary enforcement; step must be done!\n");
#endif

  /* Set the stencil width */
  sw[0]=1;
  sw[1]=1;
  sw[2]=1;

  /* Uses all default arguments, so invalid table handle -1 can be passed */
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "SpecGRMHD::evolvevars", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }
  
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "SpecGRMHD::calcvars", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }

  return;
}

int sgn (CCTK_REAL v) {
	if (v < 0.0) return -1;
	if (v == 0.0) return 0;
	if (v > 0.0) return 1;
}

void MC(CCTK_REAL *aaq,CCTK_REAL *bbq,CCTK_REAL *delta) {
	
	CCTK_INT m;
	
	for (m = 1; m < 4; m++) {
		if (aaq[m]*bbq[m] <= 0.0) {
		delta[m] = 0.0;
		}
		if (aaq[m]*bbq[m] > 0.0) {
		delta[m] = sgn(0.5*(aaq[m]+bbq[m]))*fmin(0.5*fabs(aaq[m]+bbq[m]),fmin(2.0*fabs(aaq[m]),2.0*fabs(bbq[m])));
		}
	}
}

void qmax(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_max) {

	DECLARE_CCTK_ARGUMENTS
	
	CCTK_INT new_index, m, n, p;
	
	qq_max = qq[index];
  
  for(m=-1; m < 2; m++)
  {
	  for(n=-1; n < 2; n++)
	  {
		  for(p=-1; p < 2; p++)
          {
          
	new_index = index + p + cctkGH->cctk_ash[0] * (n + cctkGH->cctk_ash[1] * m);
	
	if (qq[new_index] > qq_max) {
		qq_max = qq[new_index];
	}
	
          }
      }
  }

}

void qmin(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_min) {

	DECLARE_CCTK_ARGUMENTS
	
	CCTK_INT new_index, m, n, p;
	
	qq_min = qq[index];
  
  for(m=-1; m < 2; m++)
  {
	  for(n=-1; n < 2; n++)
	  {
		  for(p=-1; p < 2; p++)
          {
          
	new_index = index + p + cctkGH->cctk_ash[0] * (n + cctkGH->cctk_ash[1] * m);
	
	if (qq[new_index] < qq_min) {
		qq_min = qq[new_index];
	}
	
          }
      }
  }

}

void monotonize(CCTK_REAL U, CCTK_REAL *Ur, CCTK_REAL *Ul) {

	CCTK_INT m;
	CCTK_REAL dU[4];
	CCTK_REAL mU[4];
	
	for (m = 1; m < 4; m++) {
		dU[m] = Ur[m] - Ul[m];
		mU[m] = 0.5*(Ur[m] + Ul[m]);
	
	if ( (Ur[m] - U)*(U - Ul[m]) <= 0.0 ) {
		Ur[m] = U;
		Ul[m] = U;
	}
	if ( dU[m]*(U-mU[m]) > (1.0/8.0)*pow(dU[m],2.0)) {
		Ul[m] = 3.0*U - 2.0*Ur[m];
	}
	if ( dU[m]*(U-mU[m]) < -(1.0/8.0)*pow(dU[m],2.0)) {
		Ur[m] = 3.0*U - 2.0*Ul[m];
	}
	
	}

}

void Reconstruction(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  // Declare and initialize variables
  CCTK_INT i, j, k, m, n, p, ierr;
  CCTK_INT index,istart,jstart,kstart,iend,jend,kend;
  CCTK_INT uu_index1, vv_index1, uu_index2, vv_index2, uu_index3, vv_index3;
  CCTK_INT uu2_index1, uu2_index2, uu2_index3, vv2_index1, vv2_index2, vv2_index3;
  
  CCTK_REAL pi = 4.0*atan(1.0), Cwave;
  CCTK_REAL Cwave_L[4], Cwave_R[4], Cwave_max[4], Cwave_min[4];
  CCTK_REAL Urho, Urho_uu[4], Urho_vv[4], Urho_uu2[4], Urho_vv2[4], U_rho_L[4], U_rho_R[4];
  CCTK_REAL Utau, Utau_uu[4], Utau_vv[4], Utau_uu2[4], Utau_vv2[4], U_tau_L[4], U_tau_R[4];
  CCTK_REAL US_x, US_x_uu[4], US_x_vv[4], US_x_uu2[4], US_x_vv2[4], U_S_x_L[4], U_S_x_R[4];
  CCTK_REAL US_y, US_y_uu[4], US_y_vv[4], US_y_uu2[4], US_y_vv2[4], U_S_y_L[4], U_S_y_R[4];
  CCTK_REAL US_z, US_z_uu[4], US_z_vv[4], US_z_uu2[4], US_z_vv2[4], U_S_z_L[4], U_S_z_R[4];
  CCTK_REAL Ubb_x, Ubb_x_uu[4], Ubb_x_vv[4], Ubb_x_uu2[4], Ubb_x_vv2[4], U_bb_x_L[4], U_bb_x_R[4];
  CCTK_REAL Ubb_y, Ubb_y_uu[4], Ubb_y_vv[4], Ubb_y_uu2[4], Ubb_y_vv2[4], U_bb_y_L[4], U_bb_y_R[4];
  CCTK_REAL Ubb_z, Ubb_z_uu[4], Ubb_z_vv[4], Ubb_z_uu2[4], Ubb_z_vv2[4], U_bb_z_L[4], U_bb_z_R[4];
  CCTK_REAL F_rho_x_L, F_rho_x_R, F_rho_y_L, F_rho_y_R, F_rho_z_L, F_rho_z_R;
  CCTK_REAL F_tau_x_L, F_tau_x_R, F_tau_y_L, F_tau_y_R, F_tau_z_L, F_tau_z_R;
  CCTK_REAL F_SS_xx_L, F_SS_xx_R, F_SS_xy_L, F_SS_xy_R, F_SS_yx_L, F_SS_yx_R;
  CCTK_REAL F_SS_xz_L, F_SS_xz_R, F_SS_zx_L, F_SS_zx_R, F_SS_yy_L, F_SS_yy_R;
  CCTK_REAL F_SS_yz_L, F_SS_yz_R, F_SS_zy_L, F_SS_zy_R, F_SS_zz_L, F_SS_zz_R;
  CCTK_REAL F_BB_xy_L, F_BB_xy_R, F_BB_xz_L, F_BB_xz_R, F_BB_yz_L, F_BB_yz_R;
  CCTK_REAL F_BB_yx_L, F_BB_yx_R, F_BB_zx_L, F_BB_zx_R, F_BB_zy_L, F_BB_zy_R;
  CCTK_REAL aaq[4], bbq[4], ccq[4], ddq[4];
  CCTK_REAL rho0, rho_uu[4], rho_vv[4], rho_uu2[4], rho_vv2[4], rho_L[4], rho_R[4];
  CCTK_REAL pp0, pp_uu[4], pp_vv[4], pp_uu2[4], pp_vv2[4], pp_L[4], pp_R[4];
  CCTK_REAL vv_x0, vv_x_uu[4], vv_x_vv[4], vv_x_uu2[4], vv_x_vv2[4], vv_x_L[4], vv_x_R[4];
  CCTK_REAL vv_y0, vv_y_uu[4], vv_y_vv[4], vv_y_uu2[4], vv_y_vv2[4], vv_y_L[4], vv_y_R[4];
  CCTK_REAL vv_z0, vv_z_uu[4], vv_z_vv[4], vv_z_uu2[4], vv_z_vv2[4], vv_z_L[4], vv_z_R[4];
  CCTK_REAL bb_x0, bb_x_uu[4], bb_x_vv[4], bb_x_uu2[4], bb_x_vv2[4], bb_x_L[4], bb_x_R[4];
  CCTK_REAL bb_y0, bb_y_uu[4], bb_y_vv[4], bb_y_uu2[4], bb_y_vv2[4], bb_y_L[4], bb_y_R[4];
  CCTK_REAL bb_z0, bb_z_uu[4], bb_z_vv[4], bb_z_uu2[4], bb_z_vv2[4], bb_z_L[4], bb_z_R[4];
  CCTK_REAL uu_L[4][4], uu_R[4][4], vv_L[4][4], vv_R[4][4], bb_L[4][4], bb_R[4][4];
  CCTK_REAL lapse, shift[4], conf, shift2, gg[4][4], dtg, igg3[4][4], beta_L[4][4], beta_R[4][4];
  CCTK_REAL rho_star_L[4], rho_star_R[4], beta2_L[4], beta2_R[4], omega_L, omega_R;
  CCTK_REAL Ttt_mhd_L[4], Ttx_mhd_L[4], Tty_mhd_L[4], Ttz_mhd_L[4], delta0[4], delta1[4], delta2[4];
  CCTK_REAL Txx_mhd_L[4], Txy_mhd_L[4], Txz_mhd_L[4], Tyy_mhd_L[4], Tyz_mhd_L[4], Tzz_mhd_L[4];
  CCTK_REAL Ttt_mhd_R[4], Ttx_mhd_R[4], Tty_mhd_R[4], Ttz_mhd_R[4], vb2_L[4], vb2_R[4];
  CCTK_REAL Txx_mhd_R[4], Txy_mhd_R[4], Txz_mhd_R[4], Tyy_mhd_R[4], Tyz_mhd_R[4], Tzz_mhd_R[4];
  CCTK_REAL eps_L[4], eps_R[4], uu_t_L[4], uu_t_R[4], uu_x_L[4], uu_x_R[4], uu_y_L[4], uu_y_R[4], uu_z_L[4], uu_z_R[4];
  CCTK_REAL Speed_Sound_L[4], Speed_Sound_R[4], Speed_Alfven_L[4], Speed_Alfven_R[4]; 
  CCTK_REAL rho_max, rho_min, pp_max, pp_min, vv_x_max, vv_x_min, vv_y_max, vv_y_min;
  CCTK_REAL vv_z_max, vv_z_min,  bb_x_max, bb_x_min,  bb_y_max, bb_y_min,  bb_z_max, bb_z_min;
  
  if (cctk_nghostzones[0] >= 2) {
  istart = cctk_nghostzones[0];
  iend   = cctk_lsh[ 0 ] - cctk_nghostzones[0];
  }
  if (cctk_nghostzones[1] >= 2) {
  jstart = cctk_nghostzones[1];
  jend   = cctk_lsh[ 1 ] - cctk_nghostzones[1];
  }
  if (cctk_nghostzones[2] >= 2) {
  kstart = cctk_nghostzones[2];
  kend   = cctk_lsh[ 2 ] - cctk_nghostzones[2];
  }
  
  if (cctk_nghostzones[0] < 2) {
  istart = 2;
  iend   = cctk_lsh[ 0 ] - 2;
  }
  if (cctk_nghostzones[1] < 2) {
  jstart = 2;
  jend   = cctk_lsh[ 1 ] - 2;
  }
  if (cctk_nghostzones[2] < 2) {
  kstart = 2;
  kend   = cctk_lsh[ 2 ] - 2;
  }
  
  for(k=0; k < cctk_lsh[2]; k++)
  {
	  for(j=0; j < cctk_lsh[1]; j++)
	  {
		  for(i=0; i < cctk_lsh[0]; i++)
          {
  
  index = CCTK_GFINDEX3D( cctkGH, i, j, k);
  
  F_rho_x[index] = rho_star[index]*vv_x[index];
  F_rho_y[index] = rho_star[index]*vv_y[index];
  F_rho_z[index] = rho_star[index]*vv_z[index];
  
  F_tau_x[index] = alpha[index]*sqrtdetg[index]*Ttx_mhd[index]-rho_star[index]*vv_x[index];
  F_tau_y[index] = alpha[index]*sqrtdetg[index]*Tty_mhd[index]-rho_star[index]*vv_y[index];
  F_tau_z[index] = alpha[index]*sqrtdetg[index]*Ttz_mhd[index]-rho_star[index]*vv_z[index];
  
  F_SS_xx[index] = sqrtdetg[index]*(betax[index]*Ttx_mhd[index]+gxx[index]*Txx_mhd[index]+gxy[index]*Txy_mhd[index]+gxz[index]*Txz_mhd[index]);
  F_SS_yx[index] = sqrtdetg[index]*(betax[index]*Tty_mhd[index]+gxx[index]*Txy_mhd[index]+gxy[index]*Tyy_mhd[index]+gxz[index]*Tyz_mhd[index]);
  F_SS_zx[index] = sqrtdetg[index]*(betax[index]*Ttz_mhd[index]+gxx[index]*Txz_mhd[index]+gxy[index]*Tyz_mhd[index]+gxz[index]*Tzz_mhd[index]);
  F_SS_xy[index] = sqrtdetg[index]*(betay[index]*Ttx_mhd[index]+gxy[index]*Txx_mhd[index]+gyy[index]*Txy_mhd[index]+gyz[index]*Txz_mhd[index]);
  F_SS_yy[index] = sqrtdetg[index]*(betay[index]*Tty_mhd[index]+gxy[index]*Txy_mhd[index]+gyy[index]*Tyy_mhd[index]+gyz[index]*Tyz_mhd[index]);
  F_SS_zy[index] = sqrtdetg[index]*(betay[index]*Ttz_mhd[index]+gxy[index]*Txz_mhd[index]+gyy[index]*Tyz_mhd[index]+gyz[index]*Tzz_mhd[index]);
  F_SS_xz[index] = sqrtdetg[index]*(betaz[index]*Ttx_mhd[index]+gxz[index]*Txx_mhd[index]+gyz[index]*Txy_mhd[index]+gzz[index]*Txz_mhd[index]);
  F_SS_yz[index] = sqrtdetg[index]*(betaz[index]*Tty_mhd[index]+gxz[index]*Txy_mhd[index]+gyz[index]*Tyy_mhd[index]+gzz[index]*Tyz_mhd[index]);
  F_SS_zz[index] = sqrtdetg[index]*(betaz[index]*Ttz_mhd[index]+gxz[index]*Txz_mhd[index]+gyz[index]*Tyz_mhd[index]+gzz[index]*Tzz_mhd[index]);
  
  F_BB_xy[index] = vv_x[index]*bb_y[index] - vv_y[index]*bb_x[index];
  F_BB_yx[index] = vv_y[index]*bb_x[index] - vv_x[index]*bb_y[index];
  F_BB_xz[index] = vv_x[index]*bb_z[index] - vv_z[index]*bb_x[index];
  F_BB_zx[index] = vv_z[index]*bb_x[index] - vv_x[index]*bb_z[index];
  F_BB_yz[index] = vv_y[index]*bb_z[index] - vv_z[index]*bb_y[index];
  F_BB_zy[index] = vv_z[index]*bb_y[index] - vv_y[index]*bb_z[index];
  		}	
  	 }
  }
  
    
  for(k=kstart; k < kend; k++)
  {
	  for(j=jstart; j < jend; j++)
	  {
		  for(i=istart; i < iend; i++)
          {
		index = CCTK_GFINDEX3D( cctkGH, i, j, k );
		uu_index1 = CCTK_GFINDEX3D( cctkGH, i+1, j, k );
        vv_index1 = CCTK_GFINDEX3D( cctkGH, i-1, j, k );
        uu_index2 = CCTK_GFINDEX3D( cctkGH, i, j+1, k );
        vv_index2 = CCTK_GFINDEX3D( cctkGH, i, j-1, k );
		uu_index3 = CCTK_GFINDEX3D( cctkGH, i, j, k+1 );
        vv_index3 = CCTK_GFINDEX3D( cctkGH, i, j, k-1 );
        uu2_index1 = CCTK_GFINDEX3D( cctkGH, i+2, j, k );
        uu2_index2 = CCTK_GFINDEX3D( cctkGH, i, j+2, k );
		uu2_index3 = CCTK_GFINDEX3D( cctkGH, i, j, k+2 );
		vv2_index1 = CCTK_GFINDEX3D( cctkGH, i-2, j, k );
        vv2_index2 = CCTK_GFINDEX3D( cctkGH, i, j-2, k );
		vv2_index3 = CCTK_GFINDEX3D( cctkGH, i, j, k-2 );
		
        
  F_rho_x[index] = 0.0;
  F_rho_y[index] = 0.0;
  F_rho_z[index] = 0.0;
  
  F_tau_x[index] = 0.0;
  F_tau_y[index] = 0.0;
  F_tau_z[index] = 0.0;
  
  F_SS_xx[index] = 0.0;
  F_SS_yx[index] = 0.0;
  F_SS_zx[index] = 0.0;
  F_SS_xy[index] = 0.0;
  F_SS_yy[index] = 0.0;
  F_SS_zy[index] = 0.0;
  F_SS_xz[index] = 0.0;
  F_SS_yz[index] = 0.0;
  F_SS_zz[index] = 0.0;
  
  F_BB_xy[index] = 0.0;
  F_BB_yx[index] = 0.0;
  F_BB_xz[index] = 0.0;
  F_BB_zx[index] = 0.0;
  F_BB_yz[index] = 0.0;
  F_BB_zy[index] = 0.0;
  
  rho0 = rho[index];
  rho_uu[1] = rho[uu_index1];
  rho_vv[1] = rho[vv_index1];
  rho_uu[2] = rho[uu_index2];
  rho_vv[2] = rho[vv_index2];
  rho_uu[3] = rho[uu_index3];
  rho_vv[3] = rho[vv_index3];
  rho_uu2[1] = rho[uu2_index1];
  rho_uu2[2] = rho[uu2_index2];
  rho_uu2[3] = rho[uu2_index3];
  rho_vv2[1] = rho[vv2_index1];
  rho_vv2[2] = rho[vv2_index2];
  rho_vv2[3] = rho[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = rho0-rho_vv[m];
  	bbq[m] = rho_uu[m]-rho0;
  	ccq[m] = rho_uu2[m]-rho_uu[m];
  	ddq[m] = rho_vv[m]-rho_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, rho,index,rho_max);
  qmin(CCTK_PASS_CTOC, rho,index,rho_min);
  
  for (m = 1; m < 4; m++) {
  rho_L[m] = rho0 + 0.5*(rho_uu[m] - rho0) + (delta1[m] - delta2[m])/8.0;
  rho_R[m] = rho_vv[m] + 0.5*(rho0 - rho_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  rho_L[m] = fmax(rho_min,fmin(rho_max,rho_L[m]));
  rho_R[m] = fmin(rho_max,fmax(rho_min,rho_R[m]));
  }
  monotonize(rho0,rho_R,rho_L);
  
  pp0 = pp[index];
  pp_uu[1] = pp[uu_index1];
  pp_vv[1] = pp[vv_index1];
  pp_uu[2] = pp[uu_index2];
  pp_vv[2] = pp[vv_index2];
  pp_uu[3] = pp[uu_index3];
  pp_vv[3] = pp[vv_index3];
  pp_uu2[1] = pp[uu2_index1];
  pp_uu2[2] = pp[uu2_index2];
  pp_uu2[3] = pp[uu2_index3];
  pp_vv2[1] = pp[vv2_index1];
  pp_vv2[2] = pp[vv2_index2];
  pp_vv2[3] = pp[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = pp0-pp_vv[m];
  	bbq[m] = pp_uu[m]-pp0;
  	ccq[m] = pp_uu2[m]-pp_uu[m];
  	ddq[m] = pp_vv[m]-pp_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, pp,index,pp_max);
  qmin(CCTK_PASS_CTOC, pp,index,pp_min);
  
  for (m = 1; m < 4; m++) {
  pp_L[m] = pp0 + 0.5*(pp_uu[m] - pp0) + (delta1[m] - delta2[m])/8.0;
  pp_R[m] = pp_vv[m] + 0.5*(pp0 - pp_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  pp_L[m] = fmax(pp_min,fmin(pp_max,pp_L[m]));
  pp_R[m] = fmin(pp_max,fmax(pp_min,pp_R[m]));
  }
  monotonize(pp0,pp_R,pp_L);
  
  vv_x0 = vv_x[index];
  vv_x_uu[1] = vv_x[uu_index1];
  vv_x_vv[1] = vv_x[vv_index1];
  vv_x_uu[2] = vv_x[uu_index2];
  vv_x_vv[2] = vv_x[vv_index2];
  vv_x_uu[3] = vv_x[uu_index3];
  vv_x_vv[3] = vv_x[vv_index3];
  vv_x_uu2[1] = vv_x[uu2_index1];
  vv_x_uu2[2] = vv_x[uu2_index2];
  vv_x_uu2[3] = vv_x[uu2_index3];
  vv_x_vv2[1] = vv_x[vv2_index1];
  vv_x_vv2[2] = vv_x[vv2_index2];
  vv_x_vv2[3] = vv_x[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vv_x0-vv_x_vv[m];
  	bbq[m] = vv_x_uu[m]-vv_x0;
  	ccq[m] = vv_x_uu2[m]-vv_x_uu[m];
  	ddq[m] = vv_x_vv[m]-vv_x_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vv_x,index,vv_x_max);
  qmin(CCTK_PASS_CTOC, vv_x,index,vv_x_min);
  
  for (m = 1; m < 4; m++) {
  vv_x_L[m] = vv_x0 + 0.5*(vv_x_uu[m] - vv_x0) + (delta1[m] - delta2[m])/8.0;
  vv_x_R[m] = vv_x_vv[m] + 0.5*(vv_x0 - vv_x_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vv_x_L[m] = fmax(vv_x_min,fmin(vv_x_max,vv_x_L[m]));
  vv_x_R[m] = fmin(vv_x_max,fmax(vv_x_min,vv_x_R[m]));
  }
  monotonize(vv_x0,vv_x_R,vv_x_L);
  
  vv_y0 = vv_y[index];
  vv_y_uu[1] = vv_y[uu_index1];
  vv_y_vv[1] = vv_y[vv_index1];
  vv_y_uu[2] = vv_y[uu_index2];
  vv_y_vv[2] = vv_y[vv_index2];
  vv_y_uu[3] = vv_y[uu_index3];
  vv_y_vv[3] = vv_y[vv_index3];
  vv_y_uu2[1] = vv_y[uu2_index1];
  vv_y_uu2[2] = vv_y[uu2_index2];
  vv_y_uu2[3] = vv_y[uu2_index3];
  vv_y_vv2[1] = vv_y[vv2_index1];
  vv_y_vv2[2] = vv_y[vv2_index2];
  vv_y_vv2[3] = vv_y[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vv_y0-vv_y_vv[m];
  	bbq[m] = vv_y_uu[m]-vv_y0;
  	ccq[m] = vv_y_uu2[m]-vv_y_uu[m];
  	ddq[m] = vv_y_vv[m]-vv_y_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vv_y,index,vv_y_max);
  qmin(CCTK_PASS_CTOC, vv_y,index,vv_y_min);
  
  for (m = 1; m < 4; m++) {
  vv_y_L[m] = vv_y0 + 0.5*(vv_y_uu[m] - vv_y0) + (delta1[m] - delta2[m])/8.0;
  vv_y_R[m] = vv_y_vv[m] + 0.5*(vv_y0 - vv_y_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vv_y_L[m] = fmax(vv_y_min,fmin(vv_y_max,vv_y_L[m]));
  vv_y_R[m] = fmin(vv_y_max,fmax(vv_y_min,vv_y_R[m]));
  }
  monotonize(vv_y0,vv_y_R,vv_y_L);
  
  vv_z0 = vv_z[index];
  vv_z_uu[1] = vv_z[uu_index1];
  vv_z_vv[1] = vv_z[vv_index1];
  vv_z_uu[2] = vv_z[uu_index2];
  vv_z_vv[2] = vv_z[vv_index2];
  vv_z_uu[3] = vv_z[uu_index3];
  vv_z_vv[3] = vv_z[vv_index3];
  vv_z_uu2[1] = vv_z[uu2_index1];
  vv_z_uu2[2] = vv_z[uu2_index2];
  vv_z_uu2[3] = vv_z[uu2_index3];
  vv_z_vv2[1] = vv_z[vv2_index1];
  vv_z_vv2[2] = vv_z[vv2_index2];
  vv_z_vv2[3] = vv_z[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vv_z0-vv_z_vv[m];
  	bbq[m] = vv_z_uu[m]-vv_z0;
  	ccq[m] = vv_z_uu2[m]-vv_z_uu[m];
  	ddq[m] = vv_z_vv[m]-vv_z_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vv_z,index,vv_z_max);
  qmin(CCTK_PASS_CTOC, vv_z,index,vv_z_min);
  
  for (m = 1; m < 4; m++) {
  vv_z_L[m] = vv_z0 + 0.5*(vv_z_uu[m] - vv_z0) + (delta1[m] - delta2[m])/8.0;
  vv_z_R[m] = vv_z_vv[m] + 0.5*(vv_z0 - vv_z_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vv_z_L[m] = fmax(vv_z_min,fmin(vv_z_max,vv_z_L[m]));
  vv_z_R[m] = fmin(vv_z_max,fmax(vv_z_min,vv_z_R[m]));
  }
  monotonize(vv_z0,vv_z_R,vv_z_L);
  
  bb_x0 = b_x[index];
  bb_x_uu[1] = b_x[uu_index1];
  bb_x_vv[1] = b_x[vv_index1];
  bb_x_uu[2] = b_x[uu_index2];
  bb_x_vv[2] = b_x[vv_index2];
  bb_x_uu[3] = b_x[uu_index3];
  bb_x_vv[3] = b_x[vv_index3];
  bb_x_uu2[1] = b_x[uu2_index1];
  bb_x_uu2[2] = b_x[uu2_index2];
  bb_x_uu2[3] = b_x[uu2_index3];
  bb_x_vv2[1] = b_x[vv2_index1];
  bb_x_vv2[2] = b_x[vv2_index2];
  bb_x_vv2[3] = b_x[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = bb_x0-bb_x_vv[m];
  	bbq[m] = bb_x_uu[m]-bb_x0;
  	ccq[m] = bb_x_uu2[m]-bb_x_uu[m];
  	ddq[m] = bb_x_vv[m]-bb_x_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, b_x,index,bb_x_max);
  qmin(CCTK_PASS_CTOC, b_x,index,bb_x_min);
  
  for (m = 1; m < 4; m++) {
  bb_x_L[m] = bb_x0 + 0.5*(bb_x_uu[m] - bb_x0) + (delta1[m] - delta2[m])/8.0;
  bb_x_R[m] = bb_x_vv[m] + 0.5*(bb_x0 - bb_x_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  bb_x_L[m] = fmax(bb_x_min,fmin(bb_x_max,bb_x_L[m]));
  bb_x_R[m] = fmin(bb_x_max,fmax(bb_x_min,bb_x_R[m]));
  }
  monotonize(bb_x0,bb_x_R,bb_x_L);
  
  bb_y0 = b_y[index];
  bb_y_uu[1] = b_y[uu_index1];
  bb_y_vv[1] = b_y[vv_index1];
  bb_y_uu[2] = b_y[uu_index2];
  bb_y_vv[2] = b_y[vv_index2];
  bb_y_uu[3] = b_y[uu_index3];
  bb_y_vv[3] = b_y[vv_index3];
  bb_y_uu2[1] = b_y[uu2_index1];
  bb_y_uu2[2] = b_y[uu2_index2];
  bb_y_uu2[3] = b_y[uu2_index3];
  bb_y_vv2[1] = b_y[vv2_index1];
  bb_y_vv2[2] = b_y[vv2_index2];
  bb_y_vv2[3] = b_y[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = bb_y0-bb_y_vv[m];
  	bbq[m] = bb_y_uu[m]-bb_y0;
  	ccq[m] = bb_y_uu2[m]-bb_y_uu[m];
  	ddq[m] = bb_y_vv[m]-bb_y_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, b_y,index,bb_y_max);
  qmin(CCTK_PASS_CTOC, b_y,index,bb_y_min);
  
  for (m = 1; m < 4; m++) {
  bb_y_L[m] = bb_y0 + 0.5*(bb_y_uu[m] - bb_y0) + (delta1[m] - delta2[m])/8.0;
  bb_y_R[m] = bb_y_vv[m] + 0.5*(bb_y0 - bb_y_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  bb_y_L[m] = fmax(bb_y_min,fmin(bb_y_max,bb_y_L[m]));
  bb_y_R[m] = fmin(bb_y_max,fmax(bb_y_min,bb_y_R[m]));
  }
  monotonize(bb_y0,bb_y_R,bb_y_L);
  
  bb_z0 = b_z[index];
  bb_z_uu[1] = b_z[uu_index1];
  bb_z_vv[1] = b_z[vv_index1];
  bb_z_uu[2] = b_z[uu_index2];
  bb_z_vv[2] = b_z[vv_index2];
  bb_z_uu[3] = b_z[uu_index3];
  bb_z_vv[3] = b_z[vv_index3];
  bb_z_uu2[1] = b_z[uu2_index1];
  bb_z_uu2[2] = b_z[uu2_index2];
  bb_z_uu2[3] = b_z[uu2_index3];
  bb_z_vv2[1] = b_z[vv2_index1];
  bb_z_vv2[2] = b_z[vv2_index2];
  bb_z_vv2[3] = b_z[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = bb_z0-bb_z_vv[m];
  	bbq[m] = bb_z_uu[m]-bb_z0;
  	ccq[m] = bb_z_uu2[m]-bb_z_uu[m];
  	ddq[m] = bb_z_vv[m]-bb_z_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, b_z,index,bb_z_max);
  qmin(CCTK_PASS_CTOC, b_z,index,bb_z_min);
  
  for (m = 1; m < 4; m++) {
  bb_z_L[m] = bb_z0 + 0.5*(bb_z_uu[m] - bb_z0) + (delta1[m] - delta2[m])/8.0;
  bb_z_R[m] = bb_z_vv[m] + 0.5*(bb_z0 - bb_z_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  bb_z_L[m] = fmax(bb_z_min,fmin(bb_z_max,bb_z_L[m]));
  bb_z_R[m] = fmin(bb_z_max,fmax(bb_z_min,bb_z_R[m]));
  }
  monotonize(bb_z0,bb_z_R,bb_z_L);
  
  vv_L[1][1] = vv_x_L[1]; 
  vv_L[1][2] = vv_x_L[2]; 
  vv_L[1][3] = vv_x_L[3]; 
  vv_L[2][1] = vv_y_L[1]; 
  vv_L[2][2] = vv_y_L[2]; 
  vv_L[2][3] = vv_y_L[3]; 
  vv_L[3][1] = vv_z_L[1]; 
  vv_L[3][2] = vv_z_L[2]; 
  vv_L[3][3] = vv_z_L[3]; 
  
  vv_R[1][1] = vv_x_R[1]; 
  vv_R[1][2] = vv_x_R[2]; 
  vv_R[1][3] = vv_x_R[3]; 
  vv_R[2][1] = vv_y_R[1]; 
  vv_R[2][2] = vv_y_R[2]; 
  vv_R[2][3] = vv_y_R[3]; 
  vv_R[3][1] = vv_z_R[1]; 
  vv_R[3][2] = vv_z_R[2]; 
  vv_R[3][3] = vv_z_R[3]; 
  
  bb_L[1][1] = bb_x_L[1];
  bb_L[1][2] = bb_x_L[2];
  bb_L[1][3] = bb_x_L[3];
  bb_L[2][1] = bb_y_L[1];
  bb_L[2][2] = bb_y_L[2];
  bb_L[2][3] = bb_y_L[3];
  bb_L[3][1] = bb_z_L[1];
  bb_L[3][2] = bb_z_L[2];
  bb_L[3][3] = bb_z_L[3];
  
  bb_R[1][1] = bb_x_R[1];
  bb_R[1][2] = bb_x_R[2];
  bb_R[1][3] = bb_x_R[3];
  bb_R[2][1] = bb_y_R[1];
  bb_R[2][2] = bb_y_R[2];
  bb_R[2][3] = bb_y_R[3];
  bb_R[3][1] = bb_z_R[1];
  bb_R[3][2] = bb_z_R[2];
  bb_R[3][3] = bb_z_R[3];
  
  lapse = alpha[index];
    	
  shift[1] = betax[index];
  shift[2] = betay[index];
  shift[3] = betaz[index];
			
  conf = exp(4.0*phi[index]);
  	        
  gg[0][0] = gtt[index];
  gg[0][1] = gtx[index];
  gg[0][2] = gty[index];
  gg[0][3] = gtt[index];
  gg[1][0] = gg[0][1];
  gg[2][0] = gg[0][2];
  gg[3][0] = gg[0][3];
  gg[1][1] = gxx[index];
  gg[2][2] = gyy[index];
  gg[3][3] = gzz[index];
  gg[1][2] = gxy[index];
  gg[1][3] = gxz[index];
  gg[2][3] = gyz[index];
  gg[2][1] = gxy[index];
  gg[3][1] = gxz[index];
  gg[3][2] = gyz[index];
			
  cofactor(gg,igg3,4);
  
  beta_L[0][1] = 0.0;
  beta_L[0][2] = 0.0;
  beta_L[0][3] = 0.0;
  beta_R[0][1] = 0.0;
  beta_R[0][2] = 0.0;
  beta_R[0][3] = 0.0;
  
  for (m = 1; m < 4; m++) {
  eps_L[m] = 0.0;
  eps_R[m] = 0.0;
  vb2_L[m] = 0.0;
  vb2_R[m] = 0.0;
  uu_t_L[m] = 1.0/alpha[index];
  uu_t_R[m] = 1.0/alpha[index];
  }
  
  for (m = 1; m < 4; m++) {
  	for (n = 1; n < 4; n++) {
  		for (p = 1; p < 4; p++) {
  		vb2_L[m] += gg[n][p]*(vv_L[n][m] + shift[n])*(vv_L[p][m] + shift[p]);
  		vb2_R[m] += gg[n][p]*(vv_R[n][m] + shift[n])*(vv_R[p][m] + shift[p]);
  		}
	}
  }
    
  for (m = 1; m < 4; m++) {
  if (rho_L[m] != 0.0) {
  eps_L[m] = pp_L[m]/(rho_L[m]*(gam-1.0)); 
  }
  if (rho_R[m] != 0.0) {
  eps_R[m] = pp_R[m]/(rho_R[m]*(gam-1.0));
  }
  if (vb2_L[m]<alpha[index]*alpha[index]) { 
  uu_t_L[m] = 1.0/pow(alpha[index]*alpha[index] - vb2_L[m],0.5);
  }
  if (vb2_R[m]<alpha[index]*alpha[index]) { 
  uu_t_R[m] = 1.0/pow(alpha[index]*alpha[index] - vb2_R[m],0.5);
  }
  uu_x_L[m] = (vv_x_L[m] + betax[index])*alpha[index]*uu_t_L[m];
  uu_x_R[m] = (vv_x_R[m] + betax[index])*alpha[index]*uu_t_R[m];
  uu_y_L[m] = (vv_y_L[m] + betay[index])*alpha[index]*uu_t_L[m];
  uu_y_R[m] = (vv_y_R[m] + betay[index])*alpha[index]*uu_t_R[m];
  uu_z_L[m] = (vv_z_L[m] + betaz[index])*alpha[index]*uu_t_L[m];
  uu_z_R[m] = (vv_z_R[m] + betaz[index])*alpha[index]*uu_t_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  	uu_L[1][m] = uu_x_L[m];
  	uu_L[2][m] = uu_y_L[m];
  	uu_L[3][m] = uu_z_L[m];
  	uu_R[1][m] = uu_x_R[m];
  	uu_R[2][m] = uu_y_R[m];
  	uu_R[3][m] = uu_z_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  rho_star_L[m] = alpha[index]*sqrtdetg[index]*rho_L[m]*uu_t_L[m]; 
  rho_star_R[m] = alpha[index]*sqrtdetg[index]*rho_R[m]*uu_t_R[m];
  	for (n = 1; n < 4; n++) {
  		for (p = 1; p < 4; p++) {
  			beta_L[0][m] += (gg[n][p]*uu_L[n][m]*bb_L[p][m])/(sqrt(4.0*pi)*alpha[index]); 
  			beta_R[0][m] += (gg[n][p]*uu_R[n][m]*bb_R[p][m])/(sqrt(4.0*pi)*alpha[index]); 
  		}
  	}
  	for (n = 1; n < 4; n++) {
  		beta_L[n][m] = bb_L[n][m]/(alpha[index]*sqrt(4.0*pi)*uu_t_L[m]) + beta_L[0][m]*(uu_L[n][m]/uu_t_L[m]);
  		beta_R[n][m] = bb_R[n][m]/(alpha[index]*sqrt(4.0*pi)*uu_t_R[m]) + beta_R[0][m]*(uu_R[n][m]/uu_t_R[m]);
  	}
  }
  
for (m = 1; m < 4; m++) {
	beta2_L[m] = 0.0;
	beta2_R[m] = 0.0;

  	for (n = 0; n < 4; n++) {
  		for (p = 0; p < 4; p++) {		
	    	beta2_L[m] += gg[n][p]*beta_L[n][m]*beta_L[p][m];
	    	beta2_R[m] += gg[n][p]*beta_R[n][m]*beta_R[p][m];
	    }
  	}
}
  
  for (m = 1; m < 4; m++) {
  if (rho_L[m] > 0.0) {
  Speed_Sound_L[m] = pow(gam*pp_L[m]/(rho_L[m] + pp_L[m] + rho_L[m]*eps_L[m] + trQ[index]),0.5);
  }
  if (rho_R[m] > 0.0) {
  Speed_Sound_R[m] = pow(gam*pp_R[m]/(rho_R[m] + pp_R[m] + rho_R[m]*eps_R[m] + trQ[index]),0.5);
  }
  if ((rho_L[m] > 0.0) || (beta2_L[m] != 0.0)) {
  Speed_Alfven_L[m] = pow(beta2_L[m]/(rho_L[m] + pp_L[m] + rho_L[m]*eps_L[m] + trQ[index] + beta2_L[m]),0.5);
  }
  if ((rho_R[m] > 0.0) || (beta2_R[m] != 0.0)) {
  Speed_Alfven_R[m] = pow(beta2_R[m]/(rho_R[m] + pp_R[m] + rho_R[m]*eps_R[m] + trQ[index] + beta2_R[m]),0.5);
  }
  
  Cwave_L[m] = sqrt(pow(Speed_Alfven_L[m],2.0) + pow(Speed_Sound_L[m],2.0)*(1.0-pow(Speed_Alfven_L[m],2.0)));
  Cwave_R[m] = sqrt(pow(Speed_Alfven_R[m],2.0) + pow(Speed_Sound_R[m],2.0)*(1.0-pow(Speed_Alfven_R[m],2.0)));
  
  Cwave_max[m] = fmax(0,fmax(Cwave_R[m],Cwave_L[m]));
  Cwave_min[m] = -fmin(0,fmin(-Cwave_R[m],-Cwave_L[m]));
  }
  
  Cwave_x_min[index] = Cwave_min[1];
  Cwave_y_min[index] = Cwave_min[2];
  Cwave_z_min[index] = Cwave_min[3];
  Cwave_x_max[index] = Cwave_max[1];
  Cwave_y_max[index] = Cwave_max[2];
  Cwave_z_max[index] = Cwave_max[3];
  
  	for (m = 1; m < 4; m++) {
		
/* Calculate the Tmunu for MHD */

		omega_L = rho_L[m] + pp_L[m] + rho_L[m]*eps_L[m] + trQ[index];
  		
  		Ttt_mhd_L[m] = (omega_L+beta2_L[m])*uu_t_L[m]*uu_t_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[0][0] - beta_L[0][m]*beta_L[0][m];
  		
  		Ttx_mhd_L[m] = (omega_L+beta2_L[m])*uu_t_L[m]*uu_x_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[0][1] - beta_L[0][m]*beta_L[1][m];
  		Tty_mhd_L[m] = (omega_L+beta2_L[m])*uu_t_L[m]*uu_y_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[0][2] - beta_L[0][m]*beta_L[2][m];
  		Ttz_mhd_L[m] = (omega_L+beta2_L[m])*uu_t_L[m]*uu_z_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[0][3] - beta_L[0][m]*beta_L[3][m];
  		
        Txx_mhd_L[m] = (omega_L+beta2_L[m])*uu_x_L[m]*uu_x_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[1][1] - beta_L[1][m]*beta_L[1][m] + Qxx[index];
  		Tyy_mhd_L[m] = (omega_L+beta2_L[m])*uu_y_L[m]*uu_y_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[2][2] - beta_L[2][m]*beta_L[2][m] + Qyy[index];
  		Tzz_mhd_L[m] = (omega_L+beta2_L[m])*uu_z_L[m]*uu_z_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[3][3] - beta_L[3][m]*beta_L[3][m] + Qzz[index];
  		Txy_mhd_L[m] = (omega_L+beta2_L[m])*uu_x_L[m]*uu_y_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[1][2] - beta_L[1][m]*beta_L[2][m] + Qxy[index];
  		Txz_mhd_L[m] = (omega_L+beta2_L[m])*uu_x_L[m]*uu_z_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[1][3] - beta_L[1][m]*beta_L[3][m] + Qxz[index];
  		Tyz_mhd_L[m] = (omega_L+beta2_L[m])*uu_y_L[m]*uu_z_L[m] + (pp_L[m]+trQ[index]+0.5*beta2_L[m])*igg3[2][3] - beta_L[2][m]*beta_L[3][m] + Qyz[index];
  		
  		omega_R = rho_R[m] + pp_R[m] + rho_R[m]*eps_R[m] + trQ[index];
  		
  		Ttt_mhd_R[m] = (omega_R+beta2_R[m])*uu_t_R[m]*uu_t_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[0][0] - beta_R[0][m]*beta_R[0][m];
  		
  		Ttx_mhd_R[m] = (omega_R+beta2_R[m])*uu_t_R[m]*uu_x_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[0][1] - beta_R[0][m]*beta_R[1][m];
  		Tty_mhd_R[m] = (omega_R+beta2_R[m])*uu_t_R[m]*uu_y_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[0][2] - beta_R[0][m]*beta_R[2][m];
  		Ttz_mhd_R[m] = (omega_R+beta2_R[m])*uu_t_R[m]*uu_z_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[0][3] - beta_R[0][m]*beta_R[3][m];
  		
        Txx_mhd_R[m] = (omega_R+beta2_R[m])*uu_x_R[m]*uu_x_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[1][1] - beta_R[1][m]*beta_R[1][m] + Qxx[index];
  		Tyy_mhd_R[m] = (omega_R+beta2_R[m])*uu_y_R[m]*uu_y_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[2][2] - beta_R[2][m]*beta_R[2][m] + Qyy[index];
  		Tzz_mhd_R[m] = (omega_R+beta2_R[m])*uu_z_R[m]*uu_z_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[3][3] - beta_R[3][m]*beta_R[3][m] + Qzz[index];
  		Txy_mhd_R[m] = (omega_R+beta2_R[m])*uu_x_R[m]*uu_y_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[1][2] - beta_R[1][m]*beta_R[2][m] + Qxy[index];
  		Txz_mhd_R[m] = (omega_R+beta2_R[m])*uu_x_R[m]*uu_z_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[1][3] - beta_R[1][m]*beta_R[3][m] + Qxz[index];
  		Tyz_mhd_R[m] = (omega_R+beta2_R[m])*uu_y_R[m]*uu_z_R[m] + (pp_R[m]+trQ[index]+0.5*beta2_R[m])*igg3[2][3] - beta_R[2][m]*beta_R[3][m] + Qyz[index];
  		
  		if (rho[index] <= 0.0) {
  	
    	Ttt_mhd_L[m] = beta2_L[m]*uu_t_L[m]*uu_t_L[m] - beta_L[0][m]*beta_L[0][m] + 0.5*beta2_L[m]*igg3[0][0];
    	Ttx_mhd_L[m] = beta2_L[m]*uu_t_L[m]*uu_x_L[m] - beta_L[0][m]*beta_L[1][m] + 0.5*beta2_L[m]*igg3[0][1];
  		Tty_mhd_L[m] = beta2_L[m]*uu_t_L[m]*uu_y_L[m] - beta_L[0][m]*beta_L[2][m] + 0.5*beta2_L[m]*igg3[0][2];
  		Ttz_mhd_L[m] = beta2_L[m]*uu_t_L[m]*uu_z_L[m] - beta_L[0][m]*beta_L[3][m] + 0.5*beta2_L[m]*igg3[0][3];
    	Txx_mhd_L[m] = beta2_L[m]*uu_x_L[m]*uu_x_L[m] - beta_L[1][m]*beta_L[1][m] + 0.5*beta2_L[m]*igg3[1][1];
  		Tyy_mhd_L[m] = beta2_L[m]*uu_y_L[m]*uu_y_L[m] - beta_L[2][m]*beta_L[2][m] + 0.5*beta2_L[m]*igg3[2][2];
  		Tzz_mhd_L[m] = beta2_L[m]*uu_z_L[m]*uu_z_L[m] - beta_L[3][m]*beta_L[3][m] + 0.5*beta2_L[m]*igg3[3][3];
  		Txy_mhd_L[m] = beta2_L[m]*uu_x_L[m]*uu_y_L[m] - beta_L[1][m]*beta_L[2][m] + 0.5*beta2_L[m]*igg3[1][2];
  		Txz_mhd_L[m] = beta2_L[m]*uu_x_L[m]*uu_z_L[m] - beta_L[1][m]*beta_L[3][m] + 0.5*beta2_L[m]*igg3[1][3];
  		Tyz_mhd_L[m] = beta2_L[m]*uu_y_L[m]*uu_z_L[m] - beta_L[2][m]*beta_L[3][m] + 0.5*beta2_L[m]*igg3[2][3];
  		
  		Ttt_mhd_R[m] = beta2_R[m]*uu_t_R[m]*uu_t_R[m] - beta_R[0][m]*beta_R[0][m] + 0.5*beta2_R[m]*igg3[0][0];
    	Ttx_mhd_R[m] = beta2_R[m]*uu_t_R[m]*uu_x_R[m] - beta_R[0][m]*beta_R[1][m] + 0.5*beta2_R[m]*igg3[0][1];
  		Tty_mhd_R[m] = beta2_R[m]*uu_t_R[m]*uu_y_R[m] - beta_R[0][m]*beta_R[2][m] + 0.5*beta2_R[m]*igg3[0][2];
  		Ttz_mhd_R[m] = beta2_R[m]*uu_t_R[m]*uu_z_R[m] - beta_R[0][m]*beta_R[3][m] + 0.5*beta2_R[m]*igg3[0][3];
    	Txx_mhd_R[m] = beta2_R[m]*uu_x_R[m]*uu_x_R[m] - beta_R[1][m]*beta_R[1][m] + 0.5*beta2_R[m]*igg3[1][1];
  		Tyy_mhd_R[m] = beta2_R[m]*uu_y_R[m]*uu_y_R[m] - beta_R[2][m]*beta_R[2][m] + 0.5*beta2_R[m]*igg3[2][2];
  		Tzz_mhd_R[m] = beta2_R[m]*uu_z_R[m]*uu_z_R[m] - beta_R[3][m]*beta_R[3][m] + 0.5*beta2_R[m]*igg3[3][3];
  		Txy_mhd_R[m] = beta2_R[m]*uu_x_R[m]*uu_y_R[m] - beta_R[1][m]*beta_R[2][m] + 0.5*beta2_R[m]*igg3[1][2];
  		Txz_mhd_R[m] = beta2_R[m]*uu_x_R[m]*uu_z_R[m] - beta_R[1][m]*beta_R[3][m] + 0.5*beta2_R[m]*igg3[1][3];
  		Tyz_mhd_R[m] = beta2_R[m]*uu_y_R[m]*uu_z_R[m] - beta_R[2][m]*beta_R[3][m] + 0.5*beta2_R[m]*igg3[2][3];
  		
   		}

     }
     
  //Calculate Conserved Variables
  
  for (m = 1; m < 4; m++) {
  U_rho_L[m] = rho_star_L[m];
  U_rho_R[m] = rho_star_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  U_tau_L[m] = alpha[index]*sqrtdetg[index]*Ttt_mhd_L[m] - rho_star_L[m];
  U_tau_R[m] = alpha[index]*sqrtdetg[index]*Ttt_mhd_R[m] - rho_star_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  U_S_x_L[m] = sqrtdetg[index]*(gg[1][0]*Ttt_mhd_L[m]+gg[1][1]*Ttx_mhd_L[m]+gg[1][2]*Tty_mhd_L[m]+gg[1][3]*Ttz_mhd_L[m]);
  U_S_x_R[m] = sqrtdetg[index]*(gg[1][0]*Ttt_mhd_R[m]+gg[1][1]*Ttx_mhd_R[m]+gg[1][2]*Tty_mhd_R[m]+gg[1][3]*Ttz_mhd_R[m]);
  }
  
  for (m = 1; m < 4; m++) {
  U_S_y_L[m] = sqrtdetg[index]*(gg[2][0]*Ttt_mhd_L[m]+gg[2][1]*Ttx_mhd_L[m]+gg[2][2]*Tty_mhd_L[m]+gg[2][3]*Ttz_mhd_L[m]);
  U_S_y_R[m] = sqrtdetg[index]*(gg[2][0]*Ttt_mhd_R[m]+gg[2][1]*Ttx_mhd_R[m]+gg[2][2]*Tty_mhd_R[m]+gg[2][3]*Ttz_mhd_R[m]);
  }

  for (m = 1; m < 4; m++) {
  U_S_z_L[m] = sqrtdetg[index]*(gg[3][0]*Ttt_mhd_L[m]+gg[3][1]*Ttx_mhd_L[m]+gg[3][2]*Tty_mhd_L[m]+gg[3][3]*Ttz_mhd_L[m]);
  U_S_z_R[m] = sqrtdetg[index]*(gg[3][0]*Ttt_mhd_R[m]+gg[3][1]*Ttx_mhd_R[m]+gg[3][2]*Tty_mhd_R[m]+gg[3][3]*Ttz_mhd_R[m]);
  }
 
  for (m = 1; m < 4; m++) {
  U_bb_x_L[m] = sqrtdetg[index]*bb_x_L[m]/alpha[index];
  U_bb_x_R[m] = sqrtdetg[index]*bb_x_R[m]/alpha[index];
  }
   
  for (m = 1; m < 4; m++) {
  U_bb_y_L[m] = sqrtdetg[index]*bb_y_L[m]/alpha[index];
  U_bb_y_R[m] = sqrtdetg[index]*bb_y_R[m]/alpha[index];
  }
   
  for (m = 1; m < 4; m++) {
  U_bb_z_L[m] = sqrtdetg[index]*bb_z_L[m]/alpha[index];
  U_bb_z_R[m] = sqrtdetg[index]*bb_z_R[m]/alpha[index];
  }
  
  //Calc Fluxes
  
  F_rho_x_L = rho_star_L[1]*vv_x_L[1];
  F_rho_y_L = rho_star_L[2]*vv_y_L[2];
  F_rho_z_L = rho_star_L[3]*vv_z_L[3];
  
  F_tau_x_L = alpha[index]*sqrtdetg[index]*Ttx_mhd_L[1]-rho_star_L[1]*vv_x_L[1];
  F_tau_y_L = alpha[index]*sqrtdetg[index]*Tty_mhd_L[2]-rho_star_L[2]*vv_y_L[2];
  F_tau_z_L = alpha[index]*sqrtdetg[index]*Ttz_mhd_L[3]-rho_star_L[3]*vv_z_L[3];
  
  F_SS_xx_L = sqrtdetg[index]*(gg[1][0]*Ttx_mhd_L[1]+gxx[index]*Txx_mhd_L[1]+gxy[index]*Txy_mhd_L[1]+gxz[index]*Txz_mhd_L[1]);
  F_SS_yx_L = sqrtdetg[index]*(gg[1][0]*Tty_mhd_L[2]+gxx[index]*Txy_mhd_L[2]+gxy[index]*Tyy_mhd_L[2]+gxz[index]*Tyz_mhd_L[2]);
  F_SS_zx_L = sqrtdetg[index]*(gg[1][0]*Ttz_mhd_L[3]+gxx[index]*Txz_mhd_L[3]+gxy[index]*Tyz_mhd_L[3]+gxz[index]*Tzz_mhd_L[3]);
  F_SS_xy_L = sqrtdetg[index]*(gg[2][0]*Tty_mhd_L[1]+gxy[index]*Txx_mhd_L[1]+gyy[index]*Txy_mhd_L[1]+gyz[index]*Txz_mhd_L[1]);
  F_SS_yy_L = sqrtdetg[index]*(gg[2][0]*Tty_mhd_L[2]+gxy[index]*Txy_mhd_L[2]+gyy[index]*Tyy_mhd_L[2]+gyz[index]*Tyz_mhd_L[2]);
  F_SS_zy_L = sqrtdetg[index]*(gg[2][0]*Ttz_mhd_L[3]+gxy[index]*Txz_mhd_L[3]+gyy[index]*Tyz_mhd_L[3]+gyz[index]*Tzz_mhd_L[3]);
  F_SS_xz_L = sqrtdetg[index]*(gg[3][0]*Ttx_mhd_L[1]+gxz[index]*Txx_mhd_L[1]+gyz[index]*Txy_mhd_L[1]+gzz[index]*Txz_mhd_L[1]);
  F_SS_yz_L = sqrtdetg[index]*(gg[3][0]*Tty_mhd_L[2]+gxz[index]*Txy_mhd_L[2]+gyz[index]*Tyy_mhd_L[2]+gzz[index]*Tyz_mhd_L[2]);
  F_SS_zz_L = sqrtdetg[index]*(gg[3][0]*Ttz_mhd_L[3]+gxz[index]*Txz_mhd_L[3]+gyz[index]*Tyz_mhd_L[3]+gzz[index]*Tzz_mhd_L[3]);
  
  F_BB_xy_L = vv_x_L[1]*bb_y_L[1] - vv_y_L[1]*bb_x_L[1];
  F_BB_yx_L = vv_y_L[2]*bb_x_L[2] - vv_x_L[2]*bb_y_L[2];
  F_BB_xz_L = vv_x_L[1]*bb_z_L[1] - vv_z_L[1]*bb_x_L[1];
  F_BB_zx_L = vv_z_L[3]*bb_x_L[3] - vv_x_L[3]*bb_z_L[3];
  F_BB_yz_L = vv_y_L[2]*bb_z_L[2] - vv_z_L[2]*bb_y_L[2];
  F_BB_zy_L = vv_z_L[3]*bb_y_L[3] - vv_y_L[3]*bb_z_L[3];
  
  F_rho_x_R = rho_star_R[1]*vv_x_R[1];
  F_rho_y_R = rho_star_R[2]*vv_y_R[2];
  F_rho_z_R = rho_star_R[3]*vv_z_R[3];
  
  F_tau_x_R = alpha[index]*sqrtdetg[index]*Ttx_mhd_R[1]-rho_star_R[1]*vv_x_R[1];
  F_tau_y_R = alpha[index]*sqrtdetg[index]*Tty_mhd_R[2]-rho_star_R[2]*vv_y_R[2];
  F_tau_z_R = alpha[index]*sqrtdetg[index]*Ttz_mhd_R[3]-rho_star_R[3]*vv_z_R[3];
  
  F_SS_xx_R = sqrtdetg[index]*(gg[1][0]*Ttx_mhd_R[1]+gxx[index]*Txx_mhd_R[1]+gxy[index]*Txy_mhd_R[1]+gxz[index]*Txz_mhd_R[1]);
  F_SS_yx_R = sqrtdetg[index]*(gg[1][0]*Tty_mhd_R[2]+gxx[index]*Txy_mhd_R[2]+gxy[index]*Tyy_mhd_R[2]+gxz[index]*Tyz_mhd_R[2]);
  F_SS_zx_R = sqrtdetg[index]*(gg[1][0]*Ttz_mhd_R[3]+gxx[index]*Txz_mhd_R[3]+gxy[index]*Tyz_mhd_R[3]+gxz[index]*Tzz_mhd_R[3]);
  F_SS_xy_R = sqrtdetg[index]*(gg[2][0]*Ttx_mhd_R[1]+gxy[index]*Txx_mhd_R[1]+gyy[index]*Txy_mhd_R[1]+gyz[index]*Txz_mhd_R[1]);
  F_SS_yy_R = sqrtdetg[index]*(gg[2][0]*Tty_mhd_R[2]+gxy[index]*Txy_mhd_R[2]+gyy[index]*Tyy_mhd_R[2]+gyz[index]*Tyz_mhd_R[2]);
  F_SS_zy_R = sqrtdetg[index]*(gg[2][0]*Ttz_mhd_R[3]+gxy[index]*Txz_mhd_R[3]+gyy[index]*Tyz_mhd_R[3]+gyz[index]*Tzz_mhd_R[3]);
  F_SS_xz_R = sqrtdetg[index]*(gg[3][0]*Ttx_mhd_R[1]+gxz[index]*Txx_mhd_R[1]+gyz[index]*Txy_mhd_R[1]+gzz[index]*Txz_mhd_R[1]);
  F_SS_yz_R = sqrtdetg[index]*(gg[3][0]*Tty_mhd_R[2]+gxz[index]*Txy_mhd_R[2]+gyz[index]*Tyy_mhd_R[2]+gzz[index]*Tyz_mhd_R[2]);
  F_SS_zz_R = sqrtdetg[index]*(gg[3][0]*Ttz_mhd_R[3]+gxz[index]*Txz_mhd_R[3]+gyz[index]*Tyz_mhd_R[3]+gzz[index]*Tzz_mhd_R[3]);
  
  F_BB_xy_R = vv_x_R[1]*bb_y_R[1] - vv_y_R[1]*bb_x_R[1];
  F_BB_yx_R = vv_y_R[2]*bb_x_R[2] - vv_x_R[2]*bb_y_R[2];
  F_BB_xz_R = vv_x_R[1]*bb_z_R[1] - vv_z_R[1]*bb_x_R[1];
  F_BB_zx_R = vv_z_R[3]*bb_x_R[3] - vv_x_R[3]*bb_z_R[3];
  F_BB_yz_R = vv_y_R[2]*bb_z_R[2] - vv_z_R[2]*bb_y_R[2];
  F_BB_zy_R = vv_z_R[3]*bb_y_R[3] - vv_y_R[3]*bb_z_R[3];
  
  
  if ((Cwave_min[1] + Cwave_max[1] != 0.0) && (Cwave_min[2] + Cwave_max[2] != 0.0) && (Cwave_min[3] + Cwave_max[3] != 0.0)) {
  F_rho_x[index] = (Cwave_min[1]*F_rho_x_R+Cwave_max[1]*F_rho_x_L - Cwave_min[1]*Cwave_max[1]*(U_rho_R[1]-U_rho_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_rho_y[index] = (Cwave_min[2]*F_rho_y_R+Cwave_max[2]*F_rho_y_L - Cwave_min[2]*Cwave_max[2]*(U_rho_R[2]-U_rho_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_rho_z[index] = (Cwave_min[3]*F_rho_z_R+Cwave_max[3]*F_rho_z_L - Cwave_min[3]*Cwave_max[3]*(U_rho_R[3]-U_rho_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  
  F_tau_x[index] = (Cwave_min[1]*F_tau_x_R+Cwave_max[1]*F_tau_x_L - Cwave_min[1]*Cwave_max[1]*(U_tau_R[1]-U_tau_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_tau_y[index] = (Cwave_min[2]*F_tau_y_R+Cwave_max[2]*F_tau_y_L - Cwave_min[2]*Cwave_max[2]*(U_tau_R[2]-U_tau_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_tau_z[index] = (Cwave_min[3]*F_tau_z_R+Cwave_max[3]*F_tau_z_L - Cwave_min[3]*Cwave_max[3]*(U_tau_R[3]-U_tau_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  
  F_SS_xx[index] = (Cwave_min[1]*F_SS_xx_R+Cwave_max[1]*F_SS_xx_L - Cwave_min[1]*Cwave_max[1]*(U_S_x_R[1]-U_S_x_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SS_yx[index] = (Cwave_min[2]*F_SS_yx_R+Cwave_max[2]*F_SS_yx_L - Cwave_min[2]*Cwave_max[2]*(U_S_x_R[2]-U_S_x_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SS_zx[index] = (Cwave_min[3]*F_SS_zx_R+Cwave_max[3]*F_SS_zx_L - Cwave_min[3]*Cwave_max[3]*(U_S_x_R[3]-U_S_x_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_SS_xy[index] = (Cwave_min[1]*F_SS_xy_R+Cwave_max[1]*F_SS_xy_L - Cwave_min[1]*Cwave_max[1]*(U_S_y_R[1]-U_S_y_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SS_yy[index] = (Cwave_min[2]*F_SS_yy_R+Cwave_max[2]*F_SS_yy_L - Cwave_min[2]*Cwave_max[2]*(U_S_y_R[2]-U_S_y_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SS_zy[index] = (Cwave_min[3]*F_SS_zy_R+Cwave_max[3]*F_SS_zy_L - Cwave_min[3]*Cwave_max[3]*(U_S_y_R[3]-U_S_y_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_SS_xz[index] = (Cwave_min[1]*F_SS_xz_R+Cwave_max[1]*F_SS_xz_L - Cwave_min[1]*Cwave_max[1]*(U_S_z_R[1]-U_S_z_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SS_yz[index] = (Cwave_min[2]*F_SS_yz_R+Cwave_max[2]*F_SS_yz_L - Cwave_min[2]*Cwave_max[2]*(U_S_z_R[2]-U_S_z_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SS_zz[index] = (Cwave_min[3]*F_SS_zz_R+Cwave_max[3]*F_SS_zz_L - Cwave_min[3]*Cwave_max[3]*(U_S_z_R[3]-U_S_z_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  
  F_BB_yx[index] = (Cwave_min[2]*F_BB_yx_R+Cwave_max[2]*F_BB_yx_L - Cwave_min[2]*Cwave_max[2]*(U_bb_x_R[2]-U_bb_x_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_BB_zx[index] = (Cwave_min[3]*F_BB_zx_R+Cwave_max[3]*F_BB_zx_L - Cwave_min[3]*Cwave_max[3]*(U_bb_x_R[3]-U_bb_x_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_BB_xy[index] = (Cwave_min[1]*F_BB_xy_R+Cwave_max[1]*F_BB_xy_L - Cwave_min[1]*Cwave_max[1]*(U_bb_y_R[1]-U_bb_y_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_BB_zy[index] = (Cwave_min[3]*F_BB_zy_R+Cwave_max[3]*F_BB_zy_L - Cwave_min[3]*Cwave_max[3]*(U_bb_y_R[3]-U_bb_y_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_BB_xz[index] = (Cwave_min[1]*F_BB_xz_R+Cwave_max[1]*F_BB_xz_L - Cwave_min[1]*Cwave_max[1]*(U_bb_z_R[1]-U_bb_z_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_BB_yz[index] = (Cwave_min[2]*F_BB_yz_R+Cwave_max[2]*F_BB_yz_L - Cwave_min[2]*Cwave_max[2]*(U_bb_z_R[2]-U_bb_z_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  }
  
           }
	  }
  }        

}

void No_Reconstruction(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  // Declare and initialize variables
  CCTK_INT i, j, k, index;
  
  for(k=0; k < cctk_lsh[2]; k++)
  {
	  for(j=0; j < cctk_lsh[1]; j++)
	  {
		  for(i=0; i < cctk_lsh[0]; i++)
          {
  
  index = CCTK_GFINDEX3D( cctkGH, i, j, k);
  
  F_rho_x[index] = rho_star[index]*vv_x[index];
  F_rho_y[index] = rho_star[index]*vv_y[index];
  F_rho_z[index] = rho_star[index]*vv_z[index];
  
  F_tau_x[index] = alpha[index]*sqrtdetg[index]*Ttx_mhd[index]-rho_star[index]*vv_x[index];
  F_tau_y[index] = alpha[index]*sqrtdetg[index]*Tty_mhd[index]-rho_star[index]*vv_y[index];
  F_tau_z[index] = alpha[index]*sqrtdetg[index]*Ttz_mhd[index]-rho_star[index]*vv_z[index];
  
  F_SS_xx[index] = sqrtdetg[index]*(gtx[index]*Ttx_mhd[index]+gxx[index]*Txx_mhd[index]+gxy[index]*Txy_mhd[index]+gxz[index]*Txz_mhd[index]);
  F_SS_yx[index] = sqrtdetg[index]*(gtx[index]*Tty_mhd[index]+gxx[index]*Txy_mhd[index]+gxy[index]*Tyy_mhd[index]+gxz[index]*Tyz_mhd[index]);
  F_SS_zx[index] = sqrtdetg[index]*(gtx[index]*Ttz_mhd[index]+gxx[index]*Txz_mhd[index]+gxy[index]*Tyz_mhd[index]+gxz[index]*Tzz_mhd[index]);
  F_SS_xy[index] = sqrtdetg[index]*(gty[index]*Ttx_mhd[index]+gxy[index]*Txx_mhd[index]+gyy[index]*Txy_mhd[index]+gyz[index]*Txz_mhd[index]);
  F_SS_yy[index] = sqrtdetg[index]*(gty[index]*Tty_mhd[index]+gxy[index]*Txy_mhd[index]+gyy[index]*Tyy_mhd[index]+gyz[index]*Tyz_mhd[index]);
  F_SS_zy[index] = sqrtdetg[index]*(gty[index]*Ttz_mhd[index]+gxy[index]*Txz_mhd[index]+gyy[index]*Tyz_mhd[index]+gyz[index]*Tzz_mhd[index]);
  F_SS_xz[index] = sqrtdetg[index]*(gtz[index]*Ttx_mhd[index]+gxz[index]*Txx_mhd[index]+gyz[index]*Txy_mhd[index]+gzz[index]*Txz_mhd[index]);
  F_SS_yz[index] = sqrtdetg[index]*(gtz[index]*Tty_mhd[index]+gxz[index]*Txy_mhd[index]+gyz[index]*Tyy_mhd[index]+gzz[index]*Tyz_mhd[index]);
  F_SS_zz[index] = sqrtdetg[index]*(gtz[index]*Ttz_mhd[index]+gxz[index]*Txz_mhd[index]+gyz[index]*Tyz_mhd[index]+gzz[index]*Tzz_mhd[index]);
  
  F_BB_xy[index] = vv_x[index]*bb_y[index] - vv_y[index]*bb_x[index];
  F_BB_yx[index] = vv_y[index]*bb_x[index] - vv_x[index]*bb_y[index];
  F_BB_xz[index] = vv_x[index]*bb_z[index] - vv_z[index]*bb_x[index];
  F_BB_zx[index] = vv_z[index]*bb_x[index] - vv_x[index]*bb_z[index];
  F_BB_yz[index] = vv_y[index]*bb_z[index] - vv_z[index]*bb_y[index];
  F_BB_zy[index] = vv_z[index]*bb_y[index] - vv_y[index]*bb_z[index];

  		}	
  	 }
  }      

}

void ludcmp(CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *dd) 
{  
  CCTK_INT i,imax,j,k;
  CCTK_REAL big, dum, sum, temp, vv[4];
  
  for (i = 0; i < 4; i++) {
    vv[i] = 0.0;
  }
  
  imax = -1;

  *dd = 1.0; 
  for (i = 0; i < 4; i++) {
    big = 0.0; 
    for (j = 0; j < 4; j++) {
      if ((temp = fabs(aa[i][j])) > big) {
          big = temp; 
      }
    }
    if (big != 0.0) {
          vv[i] = 1.0/big; 
    } else {
      //fprintf(stderr,"Singular Matrix in Routine LUDCMP\n");
      //exit(1);
    }    
  }
  for (j = 0; j < 4; j++) {
    for (i = 0; i < j; i++) {
      sum = aa[i][j]; 
      for (k = 0; k <i; k++) {
         sum -=  aa[i][k]*aa[k][j]; 
      }
      aa[i][j] = sum; 
    }
    big= 0.0; 
    for (i = j; i < 4; i++) {
      sum = aa[i][j]; 
      for (k = 0; k < j; k++) {
         sum -=  aa[i][k]*aa[k][j]; 
      }
      aa[i][j] = sum; 
      if ( (dum= vv[i]*fabs(sum)) >=  big){
	    big = dum; 
	    imax = i; 
	    //if(imax < 0) {imax = 0;}
	    //if(imax > 3) {imax = 3;}
      }
    }
    if (j != imax) {
      for (k = 0; k < 4; k++) {
	     dum = aa[imax][k]; 
	     aa[imax][k] = aa[j][k]; 
	     aa[j][k] = dum; 
      }
      *dd = -(*dd); 
      vv[imax] = vv[j]; 
    }
    indx[j] = imax; 
    if (aa[j][j] == 0.0) {
       aa[j][j] = TINY; 
    }
    if (j != 3) {
      dum = 1.0/(aa[j][j]); 
      for (i = j+1; i < 4; i++) {
         aa[i][j] *=  dum; 
      }
    }
  }
}

void lubksb(CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *bb)
{
	CCTK_INT i,ii=-1,ip,j;
	CCTK_REAL sum;

	for (i=0;i<4;i++) {
		ip=indx[i];
		sum=bb[ip];
		bb[ip]=bb[i];
		if (ii >= 0) {
			for (j=ii;j<i;j++) {
			   sum -= aa[i][j]*bb[j]; 
			}
		} else if (sum != 0.0) {
		  ii=i; 
		}
		bb[i]=sum;
	}
	for (i=3;i>=0;i--) {
	   sum=bb[i];
		if (i < 4) {
		 for (j=i+1;j<4;j++) {
		   sum -= aa[i][j]*bb[j];
		 }
		}
	  bb[i]=sum/aa[i][i];
	}
}

#undef TINY
#undef LIGHT
#undef SIE
#undef Rho_Min
