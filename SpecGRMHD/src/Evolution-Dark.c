/*@@ Calculates the Tmunu and the RHS of Dark Matter @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

#define TINY 1.0e-30
#define Density 6.6742867e-11  
#define LIGHT 2.99792458e8
#define Rad_Const 7.565767e-16 

void Dark_RHS(CCTK_ARGUMENTS);
void Dark_Tmunu(CCTK_ARGUMENTS);
void Dark_Boundaries(CCTK_ARGUMENTS);
void Reconstruction_Dark(CCTK_ARGUMENTS);
void No_Reconstruction_Dark(CCTK_ARGUMENTS);
void SpecDeriv_Vector_Derivative( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void ludcmp( CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *dd );
void lubksb( CCTK_REAL **aa, CCTK_INT *indx, CCTK_REAL *bb );
void Clean_Read_Scalar(CCTK_ARGUMENTS, CCTK_REAL *scalar); 
int sgn(CCTK_REAL v);
void MC(CCTK_REAL *aaq,CCTK_REAL *bbq,CCTK_REAL *delta);
void monotonize(CCTK_REAL U, CCTK_REAL *Ur, CCTK_REAL *Ul);
CCTK_REAL determinant(CCTK_REAL a[4][4], CCTK_REAL k);
void cofactor(CCTK_REAL num[4][4], CCTK_REAL inverse[4][4], CCTK_REAL f);
void transpose(CCTK_REAL num[4][4], CCTK_REAL fac[4][4], CCTK_REAL inverse[4][4], CCTK_REAL r);
void qmax(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_max);
void qmin(CCTK_ARGUMENTS, CCTK_REAL *qq, CCTK_INT index, CCTK_REAL qq_min); 

/* Begin Dark_RHS */

void Dark_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
/*  Declare local variables */

  CCTK_INT i,j,k,m,n,p,four,handle;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez, index, ierr;
  CCTK_REAL ss,rho_star_dark_dot,tau_dark_dot,shift2,aa_avg,aa_ratio, EERad0;
  CCTK_REAL SD_x_dot,SD_y_dot,SD_z_dot,conf,sqrtdtg,Kavg,lapse_avg,Rad_avg;
  
/* Declare Arrays */

  CCTK_REAL ***Flux_rho_dark,***Flux_tau_dark,****Flux_SSD;
  CCTK_REAL kk3[4][4], igg3[4][4], detg[4], bssn_dk[4], T_dark[4][4];
  CCTK_REAL bssn_dkij[4][4][4],gg4[4][4][4],lapse[4],bssn_beta[4][4];
  

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
    
    Flux_rho_dark = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    Flux_tau_dark = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    Flux_SSD = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
                  
    for(m=0; m < 4; m++) {
       Flux_rho_dark[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       Flux_tau_dark[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       Flux_SSD[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
      
         for(n=0; n < 4; n++) {
            Flux_rho_dark[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
            Flux_tau_dark[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
            Flux_SSD[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
      
               for(p=0; p < 4; p++) {
                  Flux_SSD[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
               }
         }
    }
    
    if (reconstruction) {
		Reconstruction_Dark(CCTK_PASS_CTOC); 
	}
	
	if (!reconstruction) {
    	No_Reconstruction_Dark(CCTK_PASS_CTOC);
    }
	
	if (clean_scalars) {
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_dark_x); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_dark_y); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_rho_dark_z); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_dark_x); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_dark_y); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_tau_dark_z);
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_xx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_yx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_xz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_xy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_yy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_zy); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_zx); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_yz); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, F_SSD_zz); 
	}

for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
/* Keep it Physical */

		if ((rho_star_dark[index] < 0.0) || (isnan(rho_star_dark[index]))) {
		    rho_star_dark[index] = 0.0;
		}
		
		if(isnan(F_rho_dark_x[index])) {F_rho_dark_x[index] = 0.0;} 
		if(isnan(F_rho_dark_y[index])) {F_rho_dark_y[index] = 0.0;} 
		if(isnan(F_rho_dark_z[index])) {F_rho_dark_z[index] = 0.0;}
		if(isnan(F_tau_dark_x[index])) {F_tau_dark_x[index] = 0.0;} 
		if(isnan(F_tau_dark_y[index])) {F_tau_dark_y[index] = 0.0;} 
		if(isnan(F_tau_dark_z[index])) {F_tau_dark_z[index] = 0.0;} 
		if(isnan(F_SSD_xx[index])) {F_SSD_xx[index] = 0.0;} 
		if(isnan(F_SSD_yy[index])) {F_SSD_yy[index] = 0.0;} 
		if(isnan(F_SSD_zz[index])) {F_SSD_zz[index] = 0.0;}
		if(isnan(F_SSD_xy[index])) {F_SSD_xy[index] = 0.0;} 
		if(isnan(F_SSD_xz[index])) {F_SSD_xz[index] = 0.0;} 
		if(isnan(F_SSD_yz[index])) {F_SSD_yz[index] = 0.0;}
		if(isnan(F_SSD_yx[index])) {F_SSD_yx[index] = 0.0;} 
		if(isnan(F_SSD_zx[index])) {F_SSD_zx[index] = 0.0;} 
		if(isnan(F_SSD_zy[index])) {F_SSD_zy[index] = 0.0;}
		    
/* Define Arrays */
        
        Flux_rho_dark[0][1][index] = F_rho_dark_x[index];
        Flux_rho_dark[0][2][index] = F_rho_dark_y[index];
        Flux_rho_dark[0][3][index] = F_rho_dark_z[index];
        
        Flux_tau_dark[0][1][index] = F_tau_dark_x[index];
        Flux_tau_dark[0][2][index] = F_tau_dark_y[index];
        Flux_tau_dark[0][3][index] = F_tau_dark_z[index];

        Flux_SSD[0][1][1][index] = F_SSD_xx[index];
        Flux_SSD[0][2][1][index] = F_SSD_yx[index];
        Flux_SSD[0][1][3][index] = F_SSD_xz[index];
        Flux_SSD[0][1][2][index] = F_SSD_xy[index];
        Flux_SSD[0][2][2][index] = F_SSD_yy[index];
        Flux_SSD[0][3][2][index] = F_SSD_zy[index];
        Flux_SSD[0][3][1][index] = F_SSD_zx[index];
        Flux_SSD[0][2][3][index] = F_SSD_yz[index];
        Flux_SSD[0][3][3][index] = F_SSD_zz[index];
        
        }
	}
}

/* Calc derivatives */
 
  if (CCTK_Equals(diff,"Spectral")) {
 
  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, Flux_rho_dark); 
  
  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, Flux_tau_dark); 
  
  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, Flux_SSD); 
  
  }
  
  if (CCTK_Equals(diff,"Spectral_1d")) {
 
  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, Flux_rho_dark); 
  
  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, Flux_tau_dark); 
  
  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, Flux_SSD);
  
  }
  
  if (CCTK_Equals(diff,"Finite")) {
 
  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, Flux_rho_dark); 
  
  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, Flux_tau_dark); 
  
  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, Flux_SSD);
  
  }
  
  if (CCTK_Equals(diff,"Finite4")) {
 
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, Flux_rho_dark); 
  
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, Flux_tau_dark); 
  
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, Flux_SSD);
  
  }
  
  if (CCTK_Equals(diff,"SCR3")) {
 
  SpecDeriv_Vector_Derivative_SCR3(CCTK_PASS_CTOC, Flux_rho_dark); 
  
  SpecDeriv_Vector_Derivative_SCR3(CCTK_PASS_CTOC, Flux_tau_dark); 
  
  SpecDeriv_Tensor_Derivative_SCR3(CCTK_PASS_CTOC, Flux_SSD); 
  
  }
  
	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&Kavg, 1, CCTK_VarIndex("MHD_Analysis::TrK_avg"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&lapse_avg, 1, CCTK_VarIndex("MHD_Analysis::alp_avg"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&Rad_avg, 1, CCTK_VarIndex("SpecGRMHD::EErad"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&aa_avg, 1, CCTK_VarIndex("MHD_Analysis::aa_out_avg"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&aa_ratio, 1, CCTK_VarIndex("MHD_Analysis::aa_ratio_avg"));
	
	if (aa_ratio == 0.0) {aa_ratio = 1.0;}
	if (aa_avg == 0.0) {aa_avg = 1.0;}
	if (lapse_avg == 0.0) {lapse_avg = 1.0;}
	if (Rad_avg == 0.0) {Rad_avg = 1.0;}
	if (Kavg == 0.0) {Kavg = 1.0;}
  
/* Do the */

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
        
        detg[0] =      gxx[index] * gyy[index] * gzz[index] 
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
        
        T_dark[0][0] = Ttt_dark[index];
        T_dark[0][1] = Ttx_dark[index];
        T_dark[0][2] = Tty_dark[index];
        T_dark[0][3] = Ttz_dark[index];
        T_dark[1][0] = Ttx_dark[index];
        T_dark[2][0] = Tty_dark[index];
        T_dark[3][0] = Ttz_dark[index];
        
        T_dark[1][1] = Txx_dark[index];
        T_dark[2][2] = Tyy_dark[index];
        T_dark[3][3] = Tzz_dark[index];
        T_dark[1][2] = Txy_dark[index];
        T_dark[1][3] = Txz_dark[index];
        T_dark[2][3] = Tyz_dark[index];
        T_dark[2][1] = Txy_dark[index];
        T_dark[3][1] = Txz_dark[index];
        T_dark[3][2] = Tyz_dark[index];
        
/* Initialize temp variables */
  
 	ss = 0.0;
	
	for(m=1; m < 4; m++) {
      for(n=1; n < 4; n++) {
	  ss += sqrtdtg*(T_dark[0][0]*bssn_beta[0][m]*bssn_beta[0][n] + 
	        2.0*T_dark[0][m]*bssn_beta[0][n] + T_dark[m][n])*kk3[m][n]; 
	  }
	 }
	
	for(m=1; m < 4; m++) {
	 ss += -sqrtdtg*(T_dark[0][0]*bssn_beta[0][m] + T_dark[0][m])*lapse[m];
    }
   	    
/* Calculate RHS */
        
        rho_star_dark_dot = 0.0;
        tau_dark_dot = 0.0;
        SD_x_dot = 0.0;
        SD_y_dot = 0.0;
        SD_z_dot = 0.0;
        
        for(m=1; m < 4; m++) {
         rho_star_dark_dot += -Flux_rho_dark[m][m][index];
         
         tau_dark_dot += -Flux_tau_dark[m][m][index];
        }
        
        rho_star_dark_dt[index] = rho_star_dark_dot;
        
        tau_dark_dt[index] = tau_dark_dot + ss; //*pow(aa_ratio,-1.5); 
        
        for(m=0; m < 4; m++) {
         for(n=0; n < 4; n++) {        
             SD_x_dot += 0.5*sqrtdtg*T_dark[m][n]*gg4[1][m][n];
             SD_y_dot += 0.5*sqrtdtg*T_dark[m][n]*gg4[2][m][n];
             SD_z_dot += 0.5*sqrtdtg*T_dark[m][n]*gg4[3][m][n];
         }
        }
        
        for(m=1; m < 4; m++) { 
          SD_x_dot += -Flux_SSD[m][m][1][index];
          SD_y_dot += -Flux_SSD[m][m][2][index];
          SD_z_dot += -Flux_SSD[m][m][3][index];
        }
        
        SD_x_dt[index] = SD_x_dot;
        SD_y_dt[index] = SD_y_dot;
        SD_z_dt[index] = SD_z_dot;
        
        EERad0 = irad*Density*pow(aa_avg/(aa_ratio*aa_ratio),-4.0);
        
        if ((EErad[index] < 0.99*EERad0) && (EErad[index] > 1.01*EERad0)) {
        	EErad[index] = EERad0;
        }
        
        if (EErad[index] > 1.0e-80) {
        	EErad_dt[index] = (4.0/3.0)*(aa_avg*pow(aa_ratio,-1.5))*Kavg*EErad[index]; 
        }
        
        if ((EErad[index] <= 1.0e-80) || (Kavg > 0.0)) {
        	EErad_dt[index] = 0.0;
        	EErad[index] = irad*Density*pow(aa_avg,-4.0);
        }
        
        if (TrK[index] >= TrK_max) {
		
		rho_star_dark_dt[index] = 0.0;
        tau_dark_dt[index] = 0.0;
        SD_x_dt[index] = 0.0;
        SD_y_dt[index] = 0.0;
        SD_z_dt[index] = 0.0;
        
		}
		
        //EErad_dt[index] = 0.0;
        //EErad[index] = irad*Density*pow(aa_avg,-4.0);
        
        }
      }
   }
    
   for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           free(Flux_SSD[m][n][p]);
        }
      }
    }
    
    for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
        free(Flux_SSD[m][n]);
        free(Flux_rho_dark[m][n]);
        free(Flux_tau_dark[m][n]);
      }
    }
    
    for(m=0; m < 4; m++) {
      free(Flux_SSD[m]);
      free(Flux_rho_dark[m]);
      free(Flux_tau_dark[m]);
    }
   
   free(Flux_rho_dark);
   free(Flux_tau_dark);
   free(Flux_SSD);
}

void Dark_Tmunu(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
/* Declare local variables */
  CCTK_INT i,j,k,m,n,q,index,ierr,handle;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez;
  CCTK_REAL SD[4],T_dark[4][4],shift[4],lapse;
  CCTK_REAL shift2, sqrtdtg, pi=4.0*atan(1.0);
  CCTK_REAL u[4], u_l[4], uu2;
  CCTK_REAL vv3[4], gg[4][4], igg3[4][4], conf, dtg;
  
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

 /* Do the */

	for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
		if ((rho_star_dark[index] < 0.0) || (isnan(rho_star_dark[index]))) {
		    rho_star_dark[index] = 0.0;
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

/* Begin Tmunu Calculation */

  		SD[1] = SD_x[index];
  		SD[2] = SD_y[index];
  		SD[3] = SD_z[index];
  	
  	for(m=0; m < 4; m++) {
  	    u[m] = 0.0; 
  	}

  	for(m=1; m < 4; m++) {
       for(n=1; n < 4; n++) {
       	if (rho_star_dark[index] > 0.0) {
  	      u[m] += igg3[m][n]*SD[n]/rho_star_dark[index]; 
  	  	}
  	   }
  	}
  	
  	for(m=1; m < 4; m++) {
	   if (fabs(u[m]) >= 1.0) {
	   		u[m] = sgn(u[m])*Speed_Sound[index]; // or 0.0?
	   }
	}
	
	for(m=0; m < 4; m++) {
		if (isnan(u[m])) {u[m] = 0.0;}
	}

	uu2 = 0.0;	
	
 	for(m=1; m < 4; m++) {
     for(n=1; n < 4; n++) {
 	  uu2 += gg[m][n]*u[m]*u[n];
 	 }
 	}
 	
  	u[0] = sqrt(1.0 + fabs(uu2))/alpha[index];
  	
  	rho_dark[index] = rho_star_dark[index]/(sqrtdtg*u[0]);
  	
  	if ((rho_dark[index] <= 0.0) || (isnan(rho_dark[index]))) {
		rho_dark[index] = 0.0;
		u[0] = 1.0/alpha[index];
		u[1] = 0.0;
		u[2] = 0.0;
		u[3] = 0.0;
		uu2 = 0.0;
	}

	for(m=1; m < 4; m++) {
  	 vv3[m] = u[m]/(sqrt(1.0 + fabs(uu2))) - shift[m];
  	 if (isnan(vv3[m])) {vv3[m] = 0.0;}
	}  
	
	vvd_x[index] = vv3[1];
	vvd_y[index] = vv3[2];
	vvd_z[index] = vv3[3];
	
	uud_t[index] = u[0];
	uud_x[index] = u[1];
	uud_y[index] = u[2];
	uud_z[index] = u[3];
		
/* Calculate Spacial part of Tmunu for Dark Matter and Dark Energy*/

  		Ttt_dark[index] = rho_dark[index] * uud_t[index] * uud_t[index];
  		Ttx_dark[index] = rho_dark[index] * uud_t[index] * uud_x[index];
  		Tty_dark[index] = rho_dark[index] * uud_t[index] * uud_y[index];
  		Ttz_dark[index] = rho_dark[index] * uud_t[index] * uud_z[index];
  		Txx_dark[index] = rho_dark[index] * uud_x[index] * uud_x[index];
  		Tyy_dark[index] = rho_dark[index] * uud_y[index] * uud_y[index];
  		Tzz_dark[index] = rho_dark[index] * uud_z[index] * uud_z[index];
  		Txy_dark[index] = rho_dark[index] * uud_x[index] * uud_y[index];
  		Txz_dark[index] = rho_dark[index] * uud_x[index] * uud_z[index];
  		Tyz_dark[index] = rho_dark[index] * uud_y[index] * uud_z[index];
  		
  		if (isnan(Ttt_dark[index])) {Ttt_dark[index] = 0.0;}
   		if (isnan(Ttx_dark[index])) {Ttx_dark[index] = 0.0;}
   		if (isnan(Tty_dark[index])) {Tty_dark[index] = 0.0;}
   		if (isnan(Ttz_dark[index])) {Ttz_dark[index] = 0.0;}
   		if (isnan(Txx_dark[index])) {Txx_dark[index] = 0.0;}
   		if (isnan(Tyy_dark[index])) {Tyy_dark[index] = 0.0;}
   		if (isnan(Tzz_dark[index])) {Tzz_dark[index] = 0.0;}
   		if (isnan(Txy_dark[index])) {Txy_dark[index] = 0.0;}
   		if (isnan(Txz_dark[index])) {Txz_dark[index] = 0.0;}
   		if (isnan(Tyz_dark[index])) {Tyz_dark[index] = 0.0;}
  		
  		Ttt_vac[index] = -cosmo_constant*igg3[0][0]/(8.0*pi) + (4.0/3.0)*EErad[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[0][0]/3.0;
  		Ttx_vac[index] = -cosmo_constant*igg3[0][1]/(8.0*pi) - (4.0/3.0)*EErad[index]*betax[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[0][1]/3.0;
  		Tty_vac[index] = -cosmo_constant*igg3[0][2]/(8.0*pi) - (4.0/3.0)*EErad[index]*betay[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[0][2]/3.0;
  		Ttz_vac[index] = -cosmo_constant*igg3[0][3]/(8.0*pi) - (4.0/3.0)*EErad[index]*betaz[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[0][3]/3.0;
  		Txx_vac[index] = -cosmo_constant*igg3[1][1]/(8.0*pi) + (4.0/3.0)*EErad[index]*pow(betax[index]/alpha[index],2.0) + EErad[index]*igg3[1][1]/3.0;
  		Tyy_vac[index] = -cosmo_constant*igg3[2][2]/(8.0*pi) + (4.0/3.0)*EErad[index]*pow(betay[index]/alpha[index],2.0) + EErad[index]*igg3[2][2]/3.0;
  		Tzz_vac[index] = -cosmo_constant*igg3[3][3]/(8.0*pi) + (4.0/3.0)*EErad[index]*pow(betaz[index]/alpha[index],2.0) + EErad[index]*igg3[3][3]/3.0;
  		Txy_vac[index] = -cosmo_constant*igg3[1][2]/(8.0*pi) + (4.0/3.0)*EErad[index]*betax[index]*betay[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[1][2]/3.0;
  		Txz_vac[index] = -cosmo_constant*igg3[1][3]/(8.0*pi) + (4.0/3.0)*EErad[index]*betax[index]*betaz[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[1][3]/3.0;
  		Tyz_vac[index] = -cosmo_constant*igg3[2][3]/(8.0*pi) + (4.0/3.0)*EErad[index]*betay[index]*betaz[index]*pow(alpha[index],-2.0) + EErad[index]*igg3[2][3]/3.0;
  		      
  		if (isnan(Ttt_vac[index])) {Ttt_vac[index] = 0.0;}
   		if (isnan(Ttx_vac[index])) {Ttx_vac[index] = 0.0;}
   		if (isnan(Tty_vac[index])) {Tty_vac[index] = 0.0;}
   		if (isnan(Ttz_vac[index])) {Ttz_vac[index] = 0.0;}
   		if (isnan(Txx_vac[index])) {Txx_vac[index] = 0.0;}
   		if (isnan(Tyy_vac[index])) {Tyy_vac[index] = 0.0;}
   		if (isnan(Tzz_vac[index])) {Tzz_vac[index] = 0.0;}
   		if (isnan(Txy_vac[index])) {Txy_vac[index] = 0.0;}
   		if (isnan(Txz_vac[index])) {Txz_vac[index] = 0.0;}
   		if (isnan(Tyz_vac[index])) {Tyz_vac[index] = 0.0;}
   		
		}
     }
  } 
    
}

void Reconstruction_Dark(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  // Declare and initialize variables
  CCTK_INT i, j, k, m, n, p, ierr;
  CCTK_INT index,istart,jstart,kstart,iend,jend,kend;
  CCTK_INT uu_index1, vv_index1, uu_index2, vv_index2, uu_index3, vv_index3;
  CCTK_INT uu2_index1, uu2_index2, uu2_index3, vv2_index1, vv2_index2, vv2_index3;
  
  CCTK_REAL pi = 4.0*atan(1.0), Cwave_max[4], Cwave_min[4];
  CCTK_REAL Urho_dark, Urho_dark_uu[4], Urho_dark_vv[4], Urho_dark_uu2[4], Urho_dark_vv2[4], U_rho_dark_L[4], U_rho_dark_R[4];
  CCTK_REAL Utau_dark, Utau_dark_uu[4], Utau_dark_vv[4], Utau_dark_uu2[4], Utau_dark_vv2[4], U_tau_dark_L[4], U_tau_dark_R[4];
  CCTK_REAL USD_x, USD_x_uu[4], USD_x_vv[4], USD_x_uu2[4], USD_x_vv2[4], U_SD_x_L[4], U_SD_x_R[4];
  CCTK_REAL USD_y, USD_y_uu[4], USD_y_vv[4], USD_y_uu2[4], USD_y_vv2[4], U_SD_y_L[4], U_SD_y_R[4];
  CCTK_REAL USD_z, USD_z_uu[4], USD_z_vv[4], USD_z_uu2[4], USD_z_vv2[4], U_SD_z_L[4], U_SD_z_R[4];
  CCTK_REAL F_rho_dark_x_L, F_rho_dark_x_R, F_rho_dark_y_L, F_rho_dark_y_R, F_rho_dark_z_L, F_rho_dark_z_R;
  CCTK_REAL F_tau_dark_x_L, F_tau_dark_x_R, F_tau_dark_y_L, F_tau_dark_y_R, F_tau_dark_z_L, F_tau_dark_z_R;
  CCTK_REAL F_SSD_xx_L, F_SSD_xx_R, F_SSD_xy_L, F_SSD_xy_R, F_SSD_yx_L, F_SSD_yx_R;
  CCTK_REAL F_SSD_xz_L, F_SSD_xz_R, F_SSD_zx_L, F_SSD_zx_R, F_SSD_yy_L, F_SSD_yy_R;
  CCTK_REAL F_SSD_yz_L, F_SSD_yz_R, F_SSD_zy_L, F_SSD_zy_R, F_SSD_zz_L, F_SSD_zz_R;
  CCTK_REAL rho_dark0, rho_dark_uu[4], rho_dark_vv[4], rho_dark_uu2[4], rho_dark_vv2[4], rho_dark_L[4], rho_dark_R[4];
  CCTK_REAL vvd_x0, vvd_x_uu[4], vvd_x_vv[4], vvd_x_uu2[4], vvd_x_vv2[4], vvd_x_L[4], vvd_x_R[4];
  CCTK_REAL vvd_y0, vvd_y_uu[4], vvd_y_vv[4], vvd_y_uu2[4], vvd_y_vv2[4], vvd_y_L[4], vvd_y_R[4];
  CCTK_REAL vvd_z0, vvd_z_uu[4], vvd_z_vv[4], vvd_z_uu2[4], vvd_z_vv2[4], vvd_z_L[4], vvd_z_R[4];
  CCTK_REAL uud_L[4][4], uud_R[4][4], vvd_L[4][4], vvd_R[4][4], aaq[4], bbq[4], ccq[4], ddq[4];
  CCTK_REAL lapse, shift[4], conf, shift2, gg[4][4], rho_dark_star_L[4], rho_dark_star_R[4];
  CCTK_REAL Ttt_dark_L[4], Ttx_dark_L[4], Tty_dark_L[4], Ttz_dark_L[4], delta0[4], delta1[4], delta2[4];
  CCTK_REAL Txx_dark_L[4], Txy_dark_L[4], Txz_dark_L[4], Tyy_dark_L[4], Tyz_dark_L[4], Tzz_dark_L[4];
  CCTK_REAL Ttt_dark_R[4], Ttx_dark_R[4], Tty_dark_R[4], Ttz_dark_R[4], vb2_L[4], vb2_R[4];
  CCTK_REAL Txx_dark_R[4], Txy_dark_R[4], Txz_dark_R[4], Tyy_dark_R[4], Tyz_dark_R[4], Tzz_dark_R[4];
  CCTK_REAL uud_t_L[4], uud_t_R[4], uud_x_L[4], uud_x_R[4], uud_y_L[4], uud_y_R[4], uud_z_L[4], uud_z_R[4];
  CCTK_REAL rho_dark_max, rho_dark_min, vvd_x_max, vvd_x_min, vvd_y_max, vvd_y_min, vvd_z_max, vvd_z_min;
  
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
  
  F_rho_dark_x[index] = rho_star_dark[index]*vvd_x[index];
  F_rho_dark_y[index] = rho_star_dark[index]*vvd_y[index];
  F_rho_dark_z[index] = rho_star_dark[index]*vvd_z[index];
  
  F_tau_dark_x[index] = alpha[index]*sqrtdetg[index]*Ttx_dark[index]-rho_star_dark[index]*vvd_x[index];
  F_tau_dark_y[index] = alpha[index]*sqrtdetg[index]*Tty_dark[index]-rho_star_dark[index]*vvd_y[index];
  F_tau_dark_z[index] = alpha[index]*sqrtdetg[index]*Ttz_dark[index]-rho_star_dark[index]*vvd_z[index];
  
  F_SSD_xx[index] = sqrtdetg[index]*(gtx[index]*Ttx_dark[index]+gxx[index]*Txx_dark[index]+gxy[index]*Txy_dark[index]+gxz[index]*Txz_dark[index]);
  F_SSD_yx[index] = sqrtdetg[index]*(gtx[index]*Tty_dark[index]+gxx[index]*Txy_dark[index]+gxy[index]*Tyy_dark[index]+gxz[index]*Tyz_dark[index]);
  F_SSD_zx[index] = sqrtdetg[index]*(gtx[index]*Ttz_dark[index]+gxx[index]*Txz_dark[index]+gxy[index]*Tyz_dark[index]+gxz[index]*Tzz_dark[index]);
  F_SSD_xy[index] = sqrtdetg[index]*(gty[index]*Ttx_dark[index]+gxy[index]*Txx_dark[index]+gyy[index]*Txy_dark[index]+gyz[index]*Txz_dark[index]);
  F_SSD_yy[index] = sqrtdetg[index]*(gty[index]*Tty_dark[index]+gxy[index]*Txy_dark[index]+gyy[index]*Tyy_dark[index]+gyz[index]*Tyz_dark[index]);
  F_SSD_zy[index] = sqrtdetg[index]*(gty[index]*Ttz_dark[index]+gxy[index]*Txz_dark[index]+gyy[index]*Tyz_dark[index]+gyz[index]*Tzz_dark[index]);
  F_SSD_xz[index] = sqrtdetg[index]*(gtz[index]*Ttx_dark[index]+gxz[index]*Txx_dark[index]+gyz[index]*Txy_dark[index]+gzz[index]*Txz_dark[index]);
  F_SSD_yz[index] = sqrtdetg[index]*(gtz[index]*Tty_dark[index]+gxz[index]*Txy_dark[index]+gyz[index]*Tyy_dark[index]+gzz[index]*Tyz_dark[index]);
  F_SSD_zz[index] = sqrtdetg[index]*(gtz[index]*Ttz_dark[index]+gxz[index]*Txz_dark[index]+gyz[index]*Tyz_dark[index]+gzz[index]*Tzz_dark[index]);

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
		
        
  F_rho_dark_x[index] = 0.0;
  F_rho_dark_y[index] = 0.0;
  F_rho_dark_z[index] = 0.0;
  
  F_tau_dark_x[index] = 0.0;
  F_tau_dark_y[index] = 0.0;
  F_tau_dark_z[index] = 0.0;
  
  F_SSD_xx[index] = 0.0;
  F_SSD_yx[index] = 0.0;
  F_SSD_zx[index] = 0.0;
  F_SSD_xy[index] = 0.0;
  F_SSD_yy[index] = 0.0;
  F_SSD_zy[index] = 0.0;
  F_SSD_xz[index] = 0.0;
  F_SSD_yz[index] = 0.0;
  F_SSD_zz[index] = 0.0;
  
  rho_dark0 = rho_dark[index];
  rho_dark_uu[1] = rho_dark[uu_index1];
  rho_dark_vv[1] = rho_dark[vv_index1];
  rho_dark_uu[2] = rho_dark[uu_index2];
  rho_dark_vv[2] = rho_dark[vv_index2];
  rho_dark_uu[3] = rho_dark[uu_index3];
  rho_dark_vv[3] = rho_dark[vv_index3];
  rho_dark_uu2[1] = rho_dark[uu2_index1];
  rho_dark_uu2[2] = rho_dark[uu2_index2];
  rho_dark_uu2[3] = rho_dark[uu2_index3];
  rho_dark_vv2[1] = rho_dark[vv2_index1];
  rho_dark_vv2[2] = rho_dark[vv2_index2];
  rho_dark_vv2[3] = rho_dark[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = rho_dark0-rho_dark_vv[m];
  	bbq[m] = rho_dark_uu[m]-rho_dark0;
  	ccq[m] = rho_dark_uu2[m]-rho_dark_uu[m];
  	ddq[m] = rho_dark_vv[m]-rho_dark_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, rho_dark,index,rho_dark_max);
  qmin(CCTK_PASS_CTOC, rho_dark,index,rho_dark_min);
  
  for (m = 1; m < 4; m++) {
  rho_dark_L[m] = rho_dark0 + 0.5*(rho_dark_uu[m] - rho_dark0) + (delta1[m] - delta2[m])/8.0;
  rho_dark_R[m] = rho_dark_vv[m] + 0.5*(rho_dark0 - rho_dark_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  rho_dark_L[m] = fmax(rho_dark_min,fmin(rho_dark_max,rho_dark_L[m]));
  rho_dark_R[m] = fmin(rho_dark_max,fmax(rho_dark_min,rho_dark_R[m]));
  }
  monotonize(rho_dark0,rho_dark_R,rho_dark_L);
  
  vvd_x0 = vvd_x[index];
  vvd_x_uu[1] = vvd_x[uu_index1];
  vvd_x_vv[1] = vvd_x[vv_index1];
  vvd_x_uu[2] = vvd_x[uu_index2];
  vvd_x_vv[2] = vvd_x[vv_index2];
  vvd_x_uu[3] = vvd_x[uu_index3];
  vvd_x_vv[3] = vvd_x[vv_index3];
  vvd_x_uu2[1] = vvd_x[uu2_index1];
  vvd_x_uu2[2] = vvd_x[uu2_index2];
  vvd_x_uu2[3] = vvd_x[uu2_index3];
  vvd_x_vv2[1] = vvd_x[vv2_index1];
  vvd_x_vv2[2] = vvd_x[vv2_index2];
  vvd_x_vv2[3] = vvd_x[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vvd_x0-vvd_x_vv[m];
  	bbq[m] = vvd_x_uu[m]-vvd_x0;
  	ccq[m] = vvd_x_uu2[m]-vvd_x_uu[m];
  	ddq[m] = vvd_x_vv[m]-vvd_x_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vvd_x,index,vvd_x_max);
  qmin(CCTK_PASS_CTOC, vvd_x,index,vvd_x_min);
  
  for (m = 1; m < 4; m++) {
  vvd_x_L[m] = vvd_x0 + 0.5*(vvd_x_uu[m] - vvd_x0) + (delta1[m] - delta2[m])/8.0;
  vvd_x_R[m] = vvd_x_vv[m] + 0.5*(vvd_x0 - vvd_x_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vvd_x_L[m] = fmax(vvd_x_min,fmin(vvd_x_max,vvd_x_L[m]));
  vvd_x_R[m] = fmin(vvd_x_max,fmax(vvd_x_min,vvd_x_R[m]));
  }
  monotonize(vvd_x0,vvd_x_R,vvd_x_L);
  
  vvd_y0 = vvd_y[index];
  vvd_y_uu[1] = vvd_y[uu_index1];
  vvd_y_vv[1] = vvd_y[vv_index1];
  vvd_y_uu[2] = vvd_y[uu_index2];
  vvd_y_vv[2] = vvd_y[vv_index2];
  vvd_y_uu[3] = vvd_y[uu_index3];
  vvd_y_vv[3] = vvd_y[vv_index3];
  vvd_y_uu2[1] = vvd_y[uu2_index1];
  vvd_y_uu2[2] = vvd_y[uu2_index2];
  vvd_y_uu2[3] = vvd_y[uu2_index3];
  vvd_y_vv2[1] = vvd_y[vv2_index1];
  vvd_y_vv2[2] = vvd_y[vv2_index2];
  vvd_y_vv2[3] = vvd_y[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vvd_y0-vvd_y_vv[m];
  	bbq[m] = vvd_y_uu[m]-vvd_y0;
  	ccq[m] = vvd_y_uu2[m]-vvd_y_uu[m];
  	ddq[m] = vvd_y_vv[m]-vvd_y_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vvd_y,index,vvd_y_max);
  qmin(CCTK_PASS_CTOC, vvd_y,index,vvd_y_min);
  
  for (m = 1; m < 4; m++) {
  vvd_y_L[m] = vvd_y0 + 0.5*(vvd_y_uu[m] - vvd_y0) + (delta1[m] - delta2[m])/8.0;
  vvd_y_R[m] = vvd_y_vv[m] + 0.5*(vvd_y0 - vvd_y_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vvd_y_L[m] = fmax(vvd_y_min,fmin(vvd_y_max,vvd_y_L[m]));
  vvd_y_R[m] = fmin(vvd_y_max,fmax(vvd_y_min,vvd_y_R[m]));
  }
  monotonize(vvd_y0,vvd_y_R,vvd_y_L);
  
  vvd_z0 = vvd_z[index];
  vvd_z_uu[1] = vvd_z[uu_index1];
  vvd_z_vv[1] = vvd_z[vv_index1];
  vvd_z_uu[2] = vvd_z[uu_index2];
  vvd_z_vv[2] = vvd_z[vv_index2];
  vvd_z_uu[3] = vvd_z[uu_index3];
  vvd_z_vv[3] = vvd_z[vv_index3];
  vvd_z_uu2[1] = vvd_z[uu2_index1];
  vvd_z_uu2[2] = vvd_z[uu2_index2];
  vvd_z_uu2[3] = vvd_z[uu2_index3];
  vvd_z_vv2[1] = vvd_z[vv2_index1];
  vvd_z_vv2[2] = vvd_z[vv2_index2];
  vvd_z_vv2[3] = vvd_z[vv2_index3];
  
  for (m = 1; m < 4; m++) {
  	aaq[m] = vvd_z0-vvd_z_vv[m];
  	bbq[m] = vvd_z_uu[m]-vvd_z0;
  	ccq[m] = vvd_z_uu2[m]-vvd_z_uu[m];
  	ddq[m] = vvd_z_vv[m]-vvd_z_vv2[m];
  }
  
  MC(aaq, bbq, delta1);
  MC(bbq, ccq, delta2);
  MC(ddq, aaq, delta0);
  
  qmax(CCTK_PASS_CTOC, vvd_z,index,vvd_z_max);
  qmin(CCTK_PASS_CTOC, vvd_z,index,vvd_z_min);
  
  for (m = 1; m < 4; m++) {
  vvd_z_L[m] = vvd_z0 + 0.5*(vvd_z_uu[m] - vvd_z0) + (delta1[m] - delta2[m])/8.0;
  vvd_z_R[m] = vvd_z_vv[m] + 0.5*(vvd_z0 - vvd_z_vv[m]) + (delta0[m] - delta1[m])/8.0;
  
  vvd_z_L[m] = fmax(vvd_z_min,fmin(vvd_z_max,vvd_z_L[m]));
  vvd_z_R[m] = fmin(vvd_z_max,fmax(vvd_z_min,vvd_z_R[m]));
  }
  monotonize(vvd_z0,vvd_z_R,vvd_z_L);
  
  vvd_L[1][1] = vvd_x_L[1]; 
  vvd_L[1][2] = vvd_x_L[2]; 
  vvd_L[1][3] = vvd_x_L[3]; 
  vvd_L[2][1] = vvd_y_L[1]; 
  vvd_L[2][2] = vvd_y_L[2]; 
  vvd_L[2][3] = vvd_y_L[3]; 
  vvd_L[3][1] = vvd_z_L[1]; 
  vvd_L[3][2] = vvd_z_L[2]; 
  vvd_L[3][3] = vvd_z_L[3]; 
  
  vvd_R[1][1] = vvd_x_R[1]; 
  vvd_R[1][2] = vvd_x_R[2]; 
  vvd_R[1][3] = vvd_x_R[3]; 
  vvd_R[2][1] = vvd_y_R[1]; 
  vvd_R[2][2] = vvd_y_R[2]; 
  vvd_R[2][3] = vvd_y_R[3]; 
  vvd_R[3][1] = vvd_z_R[1]; 
  vvd_R[3][2] = vvd_z_R[2]; 
  vvd_R[3][3] = vvd_z_R[3]; 
  
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
  
  for (m = 1; m < 4; m++) {
  vb2_L[m] = 0.0;
  vb2_R[m] = 0.0;
  uud_t_L[m] = 1.0/alpha[index];
  uud_t_R[m] = 1.0/alpha[index];
  }
  
  for (m = 1; m < 4; m++) {
  	for (n = 1; n < 4; n++) {
  		for (p = 1; p < 4; p++) {
  		vb2_L[m] += gg[n][p]*(vvd_L[n][m] + shift[n])*(vvd_L[p][m] + shift[p]);
  		vb2_R[m] += gg[n][p]*(vvd_R[n][m] + shift[n])*(vvd_R[p][m] + shift[p]);
  		}
	}
  }
    
  for (m = 1; m < 4; m++) {
  if (vb2_L[m]<alpha[index]*alpha[index]) { 
  uud_t_L[m] = 1.0/pow(alpha[index]*alpha[index] - vb2_L[m],0.5);
  }
  if (vb2_R[m]<alpha[index]*alpha[index]) { 
  uud_t_R[m] = 1.0/pow(alpha[index]*alpha[index] - vb2_R[m],0.5);
  }
  uud_x_L[m] = (vvd_x_L[m] + betax[index])*alpha[index]*uud_t_L[m];
  uud_x_R[m] = (vvd_x_R[m] + betax[index])*alpha[index]*uud_t_R[m];
  uud_y_L[m] = (vvd_y_L[m] + betay[index])*alpha[index]*uud_t_L[m];
  uud_y_R[m] = (vvd_y_R[m] + betay[index])*alpha[index]*uud_t_R[m];
  uud_z_L[m] = (vvd_z_L[m] + betaz[index])*alpha[index]*uud_t_L[m];
  uud_z_R[m] = (vvd_z_R[m] + betaz[index])*alpha[index]*uud_t_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  	uud_L[1][m] = uud_x_L[m];
  	uud_L[2][m] = uud_y_L[m];
  	uud_L[3][m] = uud_z_L[m];
  	uud_R[1][m] = uud_x_R[m];
  	uud_R[2][m] = uud_y_R[m];
  	uud_R[3][m] = uud_z_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  rho_dark_star_L[m] = alpha[index]*sqrtdetg[index]*rho_dark_L[m]*uud_t_L[m]; 
  rho_dark_star_R[m] = alpha[index]*sqrtdetg[index]*rho_dark_R[m]*uud_t_R[m];
  }
  
  Cwave_min[1] = Cwave_x_min[index];
  Cwave_min[2] = Cwave_y_min[index];
  Cwave_min[3] = Cwave_z_min[index];
  Cwave_max[1] = Cwave_x_max[index];
  Cwave_max[2] = Cwave_y_max[index];
  Cwave_max[3] = Cwave_z_max[index];
  
  	for (m = 1; m < 4; m++) {
		
/* Calculate the Tmunu for MHD */

		Ttt_dark_L[m] = rho_dark_L[m] * uud_t_L[m] * uud_t_L[m];
  		Ttx_dark_L[m] = rho_dark_L[m] * uud_t_L[m] * uud_x_L[m];
  		Tty_dark_L[m] = rho_dark_L[m] * uud_t_L[m] * uud_y_L[m];
  		Ttz_dark_L[m] = rho_dark_L[m] * uud_t_L[m] * uud_z_L[m];
  		Txx_dark_L[m] = rho_dark_L[m] * uud_x_L[m] * uud_x_L[m];
  		Tyy_dark_L[m] = rho_dark_L[m] * uud_y_L[m] * uud_y_L[m];
  		Tzz_dark_L[m] = rho_dark_L[m] * uud_z_L[m] * uud_z_L[m];
  		Txy_dark_L[m] = rho_dark_L[m] * uud_x_L[m] * uud_y_L[m];
  		Txz_dark_L[m] = rho_dark_L[m] * uud_x_L[m] * uud_z_L[m];
  		Tyz_dark_L[m] = rho_dark_L[m] * uud_y_L[m] * uud_z_L[m];
  		
  		Ttt_dark_R[m] = rho_dark_R[m] * uud_t_R[m] * uud_t_R[m];
  		Ttx_dark_R[m] = rho_dark_R[m] * uud_t_R[m] * uud_x_R[m];
  		Tty_dark_R[m] = rho_dark_R[m] * uud_t_R[m] * uud_y_R[m];
  		Ttz_dark_R[m] = rho_dark_R[m] * uud_t_R[m] * uud_z_R[m];
  		Txx_dark_R[m] = rho_dark_R[m] * uud_x_R[m] * uud_x_R[m];
  		Tyy_dark_R[m] = rho_dark_R[m] * uud_y_R[m] * uud_y_R[m];
  		Tzz_dark_R[m] = rho_dark_R[m] * uud_z_R[m] * uud_z_R[m];
  		Txy_dark_R[m] = rho_dark_R[m] * uud_x_R[m] * uud_y_R[m];
  		Txz_dark_R[m] = rho_dark_R[m] * uud_x_R[m] * uud_z_R[m];
  		Tyz_dark_R[m] = rho_dark_R[m] * uud_y_R[m] * uud_z_R[m];

     }
     
  //Calculate Conserved Variables
  
  for (m = 1; m < 4; m++) {
  U_rho_dark_L[m] = rho_dark_star_L[m];
  U_rho_dark_R[m] = rho_dark_star_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  U_tau_dark_L[m] = alpha[index]*sqrtdetg[index]*Ttt_dark_L[m] - rho_dark_star_L[m];
  U_tau_dark_R[m] = alpha[index]*sqrtdetg[index]*Ttt_dark_R[m] - rho_dark_star_R[m];
  }
  
  for (m = 1; m < 4; m++) {
  U_SD_x_L[m] = sqrtdetg[index]*(gg[1][0]*Ttt_dark_L[m]+gg[1][1]*Ttx_dark_L[m]+gg[1][2]*Tty_dark_L[m]+gg[1][3]*Ttz_dark_L[m]);
  U_SD_x_R[m] = sqrtdetg[index]*(gg[1][0]*Ttt_dark_R[m]+gg[1][1]*Ttx_dark_R[m]+gg[1][2]*Tty_dark_R[m]+gg[1][3]*Ttz_dark_R[m]);
  }
  
  for (m = 1; m < 4; m++) {
  U_SD_y_L[m] = sqrtdetg[index]*(gg[2][0]*Ttt_dark_L[m]+gg[2][1]*Ttx_dark_L[m]+gg[2][2]*Tty_dark_L[m]+gg[2][3]*Ttz_dark_L[m]);
  U_SD_y_R[m] = sqrtdetg[index]*(gg[2][0]*Ttt_dark_R[m]+gg[2][1]*Ttx_dark_R[m]+gg[2][2]*Tty_dark_R[m]+gg[2][3]*Ttz_dark_R[m]);
  }

  for (m = 1; m < 4; m++) {
  U_SD_z_L[m] = sqrtdetg[index]*(gg[3][0]*Ttt_dark_L[m]+gg[3][1]*Ttx_dark_L[m]+gg[3][2]*Tty_dark_L[m]+gg[3][3]*Ttz_dark_L[m]);
  U_SD_z_R[m] = sqrtdetg[index]*(gg[3][0]*Ttt_dark_R[m]+gg[3][1]*Ttx_dark_R[m]+gg[3][2]*Tty_dark_R[m]+gg[3][3]*Ttz_dark_R[m]);
  }
  
  //Calc Fluxes
  
  F_rho_dark_x_L = rho_dark_star_L[1]*vvd_x_L[1];
  F_rho_dark_y_L = rho_dark_star_L[2]*vvd_y_L[2];
  F_rho_dark_z_L = rho_dark_star_L[3]*vvd_z_L[3];
  
  F_tau_dark_x_L = alpha[index]*sqrtdetg[index]*Ttx_dark_L[1]-rho_dark_star_L[1]*vvd_x_L[1];
  F_tau_dark_y_L = alpha[index]*sqrtdetg[index]*Tty_dark_L[2]-rho_dark_star_L[2]*vvd_y_L[2];
  F_tau_dark_z_L = alpha[index]*sqrtdetg[index]*Ttz_dark_L[3]-rho_dark_star_L[3]*vvd_z_L[3];
  
  F_SSD_xx_L = sqrtdetg[index]*(gg[1][0]*Ttx_dark_L[1]+gxx[index]*Txx_dark_L[1]+gxy[index]*Txy_dark_L[1]+gxz[index]*Txz_dark_L[1]);
  F_SSD_yx_L = sqrtdetg[index]*(gg[1][0]*Tty_dark_L[2]+gxx[index]*Txy_dark_L[2]+gxy[index]*Tyy_dark_L[2]+gxz[index]*Tyz_dark_L[2]);
  F_SSD_zx_L = sqrtdetg[index]*(gg[1][0]*Ttz_dark_L[3]+gxx[index]*Txz_dark_L[3]+gxy[index]*Tyz_dark_L[3]+gxz[index]*Tzz_dark_L[3]);
  F_SSD_xy_L = sqrtdetg[index]*(gg[2][0]*Tty_dark_L[1]+gxy[index]*Txx_dark_L[1]+gyy[index]*Txy_dark_L[1]+gyz[index]*Txz_dark_L[1]);
  F_SSD_yy_L = sqrtdetg[index]*(gg[2][0]*Tty_dark_L[2]+gxy[index]*Txy_dark_L[2]+gyy[index]*Tyy_dark_L[2]+gyz[index]*Tyz_dark_L[2]);
  F_SSD_zy_L = sqrtdetg[index]*(gg[2][0]*Ttz_dark_L[3]+gxy[index]*Txz_dark_L[3]+gyy[index]*Tyz_dark_L[3]+gyz[index]*Tzz_dark_L[3]);
  F_SSD_xz_L = sqrtdetg[index]*(gg[3][0]*Ttx_dark_L[1]+gxz[index]*Txx_dark_L[1]+gyz[index]*Txy_dark_L[1]+gzz[index]*Txz_dark_L[1]);
  F_SSD_yz_L = sqrtdetg[index]*(gg[3][0]*Tty_dark_L[2]+gxz[index]*Txy_dark_L[2]+gyz[index]*Tyy_dark_L[2]+gzz[index]*Tyz_dark_L[2]);
  F_SSD_zz_L = sqrtdetg[index]*(gg[3][0]*Ttz_dark_L[3]+gxz[index]*Txz_dark_L[3]+gyz[index]*Tyz_dark_L[3]+gzz[index]*Tzz_dark_L[3]);
  
  F_rho_dark_x_R = rho_dark_star_R[1]*vvd_x_R[1];
  F_rho_dark_y_R = rho_dark_star_R[2]*vvd_y_R[2];
  F_rho_dark_z_R = rho_dark_star_R[3]*vvd_z_R[3];
  
  F_tau_dark_x_R = alpha[index]*sqrtdetg[index]*Ttx_dark_R[1]-rho_dark_star_R[1]*vvd_x_R[1];
  F_tau_dark_y_R = alpha[index]*sqrtdetg[index]*Tty_dark_R[2]-rho_dark_star_R[2]*vvd_y_R[2];
  F_tau_dark_z_R = alpha[index]*sqrtdetg[index]*Ttz_dark_R[3]-rho_dark_star_R[3]*vvd_z_R[3];
  
  F_SSD_xx_R = sqrtdetg[index]*(gg[1][0]*Ttx_dark_R[1]+gxx[index]*Txx_dark_R[1]+gxy[index]*Txy_dark_R[1]+gxz[index]*Txz_dark_R[1]);
  F_SSD_yx_R = sqrtdetg[index]*(gg[1][0]*Tty_dark_R[2]+gxx[index]*Txy_dark_R[2]+gxy[index]*Tyy_dark_R[2]+gxz[index]*Tyz_dark_R[2]);
  F_SSD_zx_R = sqrtdetg[index]*(gg[1][0]*Ttz_dark_R[3]+gxx[index]*Txz_dark_R[3]+gxy[index]*Tyz_dark_R[3]+gxz[index]*Tzz_dark_R[3]);
  F_SSD_xy_R = sqrtdetg[index]*(gg[2][0]*Ttx_dark_R[1]+gxy[index]*Txx_dark_R[1]+gyy[index]*Txy_dark_R[1]+gyz[index]*Txz_dark_R[1]);
  F_SSD_yy_R = sqrtdetg[index]*(gg[2][0]*Tty_dark_R[2]+gxy[index]*Txy_dark_R[2]+gyy[index]*Tyy_dark_R[2]+gyz[index]*Tyz_dark_R[2]);
  F_SSD_zy_R = sqrtdetg[index]*(gg[2][0]*Ttz_dark_R[3]+gxy[index]*Txz_dark_R[3]+gyy[index]*Tyz_dark_R[3]+gyz[index]*Tzz_dark_R[3]);
  F_SSD_xz_R = sqrtdetg[index]*(gg[3][0]*Ttx_dark_R[1]+gxz[index]*Txx_dark_R[1]+gyz[index]*Txy_dark_R[1]+gzz[index]*Txz_dark_R[1]);
  F_SSD_yz_R = sqrtdetg[index]*(gg[3][0]*Tty_dark_R[2]+gxz[index]*Txy_dark_R[2]+gyz[index]*Tyy_dark_R[2]+gzz[index]*Tyz_dark_R[2]);
  F_SSD_zz_R = sqrtdetg[index]*(gg[3][0]*Ttz_dark_R[3]+gxz[index]*Txz_dark_R[3]+gyz[index]*Tyz_dark_R[3]+gzz[index]*Tzz_dark_R[3]);
  
  
  if ((Cwave_min[1] + Cwave_max[1] != 0.0) && (Cwave_min[2] + Cwave_max[2] != 0.0) && (Cwave_min[3] + Cwave_max[3] != 0.0)) {
  F_rho_dark_x[index] = (Cwave_min[1]*F_rho_dark_x_R+Cwave_max[1]*F_rho_dark_x_L - Cwave_min[1]*Cwave_max[1]*(U_rho_dark_R[1]-U_rho_dark_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_rho_dark_y[index] = (Cwave_min[2]*F_rho_dark_y_R+Cwave_max[2]*F_rho_dark_y_L - Cwave_min[2]*Cwave_max[2]*(U_rho_dark_R[2]-U_rho_dark_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_rho_dark_z[index] = (Cwave_min[3]*F_rho_dark_z_R+Cwave_max[3]*F_rho_dark_z_L - Cwave_min[3]*Cwave_max[3]*(U_rho_dark_R[3]-U_rho_dark_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  
  F_tau_dark_x[index] = (Cwave_min[1]*F_tau_dark_x_R+Cwave_max[1]*F_tau_dark_x_L - Cwave_min[1]*Cwave_max[1]*(U_tau_dark_R[1]-U_tau_dark_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_tau_dark_y[index] = (Cwave_min[2]*F_tau_dark_y_R+Cwave_max[2]*F_tau_dark_y_L - Cwave_min[2]*Cwave_max[2]*(U_tau_dark_R[2]-U_tau_dark_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_tau_dark_z[index] = (Cwave_min[3]*F_tau_dark_z_R+Cwave_max[3]*F_tau_dark_z_L - Cwave_min[3]*Cwave_max[3]*(U_tau_dark_R[3]-U_tau_dark_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  
  F_SSD_xx[index] = (Cwave_min[1]*F_SSD_xx_R+Cwave_max[1]*F_SSD_xx_L - Cwave_min[1]*Cwave_max[1]*(U_SD_x_R[1]-U_SD_x_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SSD_yx[index] = (Cwave_min[2]*F_SSD_yx_R+Cwave_max[2]*F_SSD_yx_L - Cwave_min[2]*Cwave_max[2]*(U_SD_x_R[2]-U_SD_x_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SSD_zx[index] = (Cwave_min[3]*F_SSD_zx_R+Cwave_max[3]*F_SSD_zx_L - Cwave_min[3]*Cwave_max[3]*(U_SD_x_R[3]-U_SD_x_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_SSD_xy[index] = (Cwave_min[1]*F_SSD_xy_R+Cwave_max[1]*F_SSD_xy_L - Cwave_min[1]*Cwave_max[1]*(U_SD_y_R[1]-U_SD_y_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SSD_yy[index] = (Cwave_min[2]*F_SSD_yy_R+Cwave_max[2]*F_SSD_yy_L - Cwave_min[2]*Cwave_max[2]*(U_SD_y_R[2]-U_SD_y_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SSD_zy[index] = (Cwave_min[3]*F_SSD_zy_R+Cwave_max[3]*F_SSD_zy_L - Cwave_min[3]*Cwave_max[3]*(U_SD_y_R[3]-U_SD_y_L[3]))/(Cwave_min[3]+Cwave_max[3]);
  F_SSD_xz[index] = (Cwave_min[1]*F_SSD_xz_R+Cwave_max[1]*F_SSD_xz_L - Cwave_min[1]*Cwave_max[1]*(U_SD_z_R[1]-U_SD_z_L[1]))/(Cwave_min[1]+Cwave_max[1]);
  F_SSD_yz[index] = (Cwave_min[2]*F_SSD_yz_R+Cwave_max[2]*F_SSD_yz_L - Cwave_min[2]*Cwave_max[2]*(U_SD_z_R[2]-U_SD_z_L[2]))/(Cwave_min[2]+Cwave_max[2]);
  F_SSD_zz[index] = (Cwave_min[3]*F_SSD_zz_R+Cwave_max[3]*F_SSD_zz_L - Cwave_min[3]*Cwave_max[3]*(U_SD_z_R[3]-U_SD_z_L[3]))/(Cwave_min[3]+Cwave_max[3]);

  }
  
           }
	  }
  }        

}


void No_Reconstruction_Dark(CCTK_ARGUMENTS) 
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
  
  F_rho_dark_x[index] = rho_star_dark[index]*vvd_x[index];
  F_rho_dark_y[index] = rho_star_dark[index]*vvd_y[index];
  F_rho_dark_z[index] = rho_star_dark[index]*vvd_z[index];
  
  F_tau_dark_x[index] = alpha[index]*sqrtdetg[index]*Ttx_dark[index]-rho_star_dark[index]*vvd_x[index];
  F_tau_dark_y[index] = alpha[index]*sqrtdetg[index]*Tty_dark[index]-rho_star_dark[index]*vvd_y[index];
  F_tau_dark_z[index] = alpha[index]*sqrtdetg[index]*Ttz_dark[index]-rho_star_dark[index]*vvd_z[index];
  
  F_SSD_xx[index] = sqrtdetg[index]*(gtx[index]*Ttx_dark[index]+gxx[index]*Txx_dark[index]+gxy[index]*Txy_dark[index]+gxz[index]*Txz_dark[index]);
  F_SSD_yx[index] = sqrtdetg[index]*(gtx[index]*Tty_dark[index]+gxx[index]*Txy_dark[index]+gxy[index]*Tyy_dark[index]+gxz[index]*Tyz_dark[index]);
  F_SSD_zx[index] = sqrtdetg[index]*(gtx[index]*Ttz_dark[index]+gxx[index]*Txz_dark[index]+gxy[index]*Tyz_dark[index]+gxz[index]*Tzz_dark[index]);
  F_SSD_xy[index] = sqrtdetg[index]*(gty[index]*Ttx_dark[index]+gxy[index]*Txx_dark[index]+gyy[index]*Txy_dark[index]+gyz[index]*Txz_dark[index]);
  F_SSD_yy[index] = sqrtdetg[index]*(gty[index]*Tty_dark[index]+gxy[index]*Txy_dark[index]+gyy[index]*Tyy_dark[index]+gyz[index]*Tyz_dark[index]);
  F_SSD_zy[index] = sqrtdetg[index]*(gty[index]*Ttz_dark[index]+gxy[index]*Txz_dark[index]+gyy[index]*Tyz_dark[index]+gyz[index]*Tzz_dark[index]);
  F_SSD_xz[index] = sqrtdetg[index]*(gtz[index]*Ttx_dark[index]+gxz[index]*Txx_dark[index]+gyz[index]*Txy_dark[index]+gzz[index]*Txz_dark[index]);
  F_SSD_yz[index] = sqrtdetg[index]*(gtz[index]*Tty_dark[index]+gxz[index]*Txy_dark[index]+gyz[index]*Tyy_dark[index]+gzz[index]*Tyz_dark[index]);
  F_SSD_zz[index] = sqrtdetg[index]*(gtz[index]*Ttz_dark[index]+gxz[index]*Txz_dark[index]+gyz[index]*Tyz_dark[index]+gzz[index]*Tzz_dark[index]);
  		}	
  	 }
  }      

}

#undef TINY
#undef Density
#undef LIGHT
#undef Rad_Const
