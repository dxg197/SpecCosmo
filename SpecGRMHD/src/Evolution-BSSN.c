/*@@ Calculates RHS of GRMHD equations @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define SIE 1.11265006e-17
#define Speed_Light 2.99792458e8
#define Small 1.0e-10

/* Calc RHS of bssn Eqns */

void bssn_RHS(CCTK_ARGUMENTS);
void bssn_Boundaries(CCTK_ARGUMENTS);
void SpecDeriv_Scalar_Derivative( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Vector_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Scalar_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative2( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );


void bssn_RHS(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
/* Declare local variables */

  CCTK_INT a,b,c,d,i,j,k,l,m,n,p,q,four,ierr;
  CCTK_INT sizex, sizey, sizez, index,handle;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_REAL shift2,Ricci_Scalar,pi,conf,conf_inv,zmin,zmax,lz,wc,t;
  CCTK_REAL bj0,bj1,by0,by1,lam,aa_avg,K_avg,aa_ratio;
  CCTK_REAL T_total[4][4],igg3[4][4],Si[4],Ricci[4][4],Sij[4][4];
  CCTK_REAL bssn_Chris1[4][4][4],bssn_Chris2[4][4][4];
  CCTK_REAL dtg,Lie_gxx,Lie_gyy,Lie_gzz,Aij_Scalar,ddphi,dtg_inv;
  CCTK_REAL Lie_phi,Lie_gxy,Lie_gxz,Lie_gyz,alphai[4][4],bssn_beta[4][4][4];
  CCTK_REAL AilAljxx,AilAljyy,AilAljzz,AilAljxy;
  CCTK_REAL STFxx,STFyy,STFzz,STFxy,STFxz,STFyz;
  CCTK_REAL AilAljxz,AilAljyz,Lie_Axx,Lie_Ayy;
  CCTK_REAL Lie_Azz,Lie_Axy,Lie_Axz,Lie_Ayz,AijAij;
  CCTK_REAL Lie_K,inv_Aij[4][4],S_Scalar,Ttotal[4][4],gg4[4][4];
  CCTK_REAL Aij[4][4],AijTF[4][4],Chris[4][4][4],Ricci_phi[4][4]; 
  CCTK_REAL Lie_gamma_x,Lie_gamma_y,Lie_gamma_z;
  CCTK_REAL bssn_gamma_x_RHS1,bssn_gamma_x_RHS2;
  CCTK_REAL bssn_gamma_y_RHS1,bssn_gamma_y_RHS2;
  CCTK_REAL bssn_gamma_z_RHS1,bssn_gamma_z_RHS2; 
  CCTK_REAL gamma_x, gamma_y, gamma_z,lapse_avg;
  CCTK_REAL rho_star_avg, rho_star_dark_avg;
  
/* Declare Arrays */
  
  CCTK_REAL *****bssn_g,****bssn_ginv,***lapse,***bssn_phi,***bssn_bb;
  CCTK_REAL ****bssn_A,****bssn_Ainv,**bssn_K,****shift,***gamma;
  
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
    pi = 4.0*atan(1.0);
    
    bssn_phi = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    bssn_K = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    gamma = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    lapse = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    bssn_A = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    bssn_Ainv = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    shift = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    bssn_bb = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    bssn_g = (CCTK_REAL *****)malloc(sizeof(CCTK_REAL ****)*four);
    bssn_ginv = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
                   
    for(m=0; m < 4; m++) {
       bssn_phi[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       lapse[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       bssn_K[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       gamma[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       bssn_A[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       bssn_Ainv[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       shift[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       bssn_bb[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       bssn_g[m] = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL **)*four);
       bssn_ginv[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
        
        for(n=0; n < 4; n++) {
          bssn_phi[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
          lapse[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
          gamma[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
          bssn_bb[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
          shift[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_A[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_Ainv[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_ginv[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_g[m][n] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL *)*four);
        
            for(p=0; p < 4; p++) {
              shift[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_A[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_Ainv[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_ginv[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_g[m][n][p] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
        
                for(q=0; q < 4; q++) {
                   bssn_g[m][n][p][q] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
                }
            }
        }
    }
    
    zmin = 0.0;
    zmax = 0.0;
    
    handle = CCTK_ReductionHandle("maximum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmax, 1, CCTK_VarIndex("grid::z"));
	
	handle = CCTK_ReductionHandle("minimum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmin, 1, CCTK_VarIndex("grid::z"));
 
  	zmin += CCTK_DELTA_SPACE(2)*(cctk_nghostzones[2]-0.5);
  	zmax += -CCTK_DELTA_SPACE(2)*(cctk_nghostzones[2]-0.5);
  	lz = zmax - zmin;
    
    handle = CCTK_ReductionHandle("average");
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&aa_avg, 1, CCTK_VarIndex("MHD_Analysis::aa_out_avg"));
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&lapse_avg, 1, CCTK_VarIndex("MHD_Analysis::alp_avg"));
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&K_avg, 1, CCTK_VarIndex("MHD_Analysis::TrK_avg"));
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&aa_ratio, 1, CCTK_VarIndex("MHD_Analysis::aa_ratio_avg"));
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&rho_star_avg, 1, CCTK_VarIndex("SpecGRMHD::rho_star"));
    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, (void *)&rho_star_dark_avg, 1, CCTK_VarIndex("SpecGRMHD::rho_star_dark"));
    
    if ((aa_ratio == 0.0) || (isnan(aa_ratio))) {aa_ratio = 1.0;}
    if ((aa_avg == 0.0) || (isnan(aa_avg))) {aa_avg = 1.0;}
    if ((K_avg == 0.0) || (isnan(K_avg))) {K_avg = 1.0;}
    if ((lapse_avg == 0.0) || (isnan(lapse_avg))) {lapse_avg = 1.0;}

/* Do the */

for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
    
/* Calc bssn metric */
  
    bssn_g[0][0][1][1][index] = bssn_gxx[index];
    bssn_g[0][0][2][2][index] = bssn_gyy[index];
    bssn_g[0][0][3][3][index] = bssn_gzz[index];
    bssn_g[0][0][1][2][index] = bssn_gxy[index];
    bssn_g[0][0][1][3][index] = bssn_gxz[index];
    bssn_g[0][0][2][3][index] = bssn_gyz[index];
    bssn_g[0][0][2][1][index] = bssn_gxy[index];
    bssn_g[0][0][3][1][index] = bssn_gxz[index];
    bssn_g[0][0][3][2][index] = bssn_gyz[index];
    
    dtg =      bssn_gxx[index] * bssn_gyy[index] * bssn_gzz[index] 
           + 2.0*bssn_gxy[index] * bssn_gxz[index] * bssn_gyz[index] 
           -   bssn_gxx[index] * pow(bssn_gyz[index],2.0)           
           -   bssn_gyy[index] * pow(bssn_gxz[index],2.0)           
           -   bssn_gzz[index] * pow(bssn_gxy[index],2.0);
    
    detg1[index] = dtg - 1.0;
    
    bssn_ginv[0][1][1][index] = (bssn_gyy[index] * bssn_gzz[index] - bssn_gyz[index] * bssn_gyz[index]) / dtg;
    bssn_ginv[0][2][2][index] = (bssn_gxx[index] * bssn_gzz[index] - bssn_gxz[index] * bssn_gxz[index]) / dtg;
    bssn_ginv[0][3][3][index] = (bssn_gxx[index] * bssn_gyy[index] - bssn_gxy[index] * bssn_gxy[index]) / dtg;                         
    bssn_ginv[0][1][2][index] = (bssn_gxz[index] * bssn_gyz[index] - bssn_gxy[index] * bssn_gzz[index]) / dtg;
    bssn_ginv[0][1][3][index] = (bssn_gxy[index] * bssn_gyz[index] - bssn_gxz[index] * bssn_gyy[index]) / dtg;
    bssn_ginv[0][2][3][index] = (bssn_gxz[index] * bssn_gxy[index] - bssn_gyz[index] * bssn_gxx[index]) / dtg; 
    bssn_ginv[0][2][1][index] = bssn_ginv[0][1][2][index];
    bssn_ginv[0][3][1][index] = bssn_ginv[0][1][3][index];
    bssn_ginv[0][3][2][index] = bssn_ginv[0][2][3][index];
    
    dtg_inv =      bssn_ginv[0][1][1][index] * bssn_ginv[0][2][2][index] * bssn_ginv[0][3][3][index] 
           	 + 2.0*bssn_ginv[0][1][2][index] * bssn_ginv[0][1][3][index] * bssn_ginv[0][2][3][index] 
               -   bssn_ginv[0][1][1][index] * pow(bssn_ginv[0][2][3][index],2.0)           
               -   bssn_ginv[0][2][2][index] * pow(bssn_ginv[0][1][3][index],2.0)           
               -   bssn_ginv[0][3][3][index] * pow(bssn_ginv[0][1][2][index],2.0);
    
/* Calc Trace K */
	
	conf = exp(4.0*phi[index]); 
	conf_inv = exp(-4.0*phi[index]); 
    
    bssn_K[0][index] = TrK[index];
	                  
/* Calc Trace-Free Extrinsic Curvature */
  
    bssn_A[0][1][1][index] = bssn_Axx[index];
    bssn_A[0][2][2][index] = bssn_Ayy[index];
    bssn_A[0][3][3][index] = bssn_Azz[index];
    bssn_A[0][1][2][index] = bssn_Axy[index];
    bssn_A[0][1][3][index] = bssn_Axz[index];
    bssn_A[0][2][3][index] = bssn_Ayz[index];
    bssn_A[0][2][1][index] = bssn_Axy[index];
    bssn_A[0][3][1][index] = bssn_Axz[index];
    bssn_A[0][3][2][index] = bssn_Ayz[index];
    
    for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
		        bssn_Ainv[0][m][n][index] = 0.0;
		    }
		}
		
		for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
               for(p=1; p < 4; p++) {
                  for(q=1; q < 4; q++) {
		             bssn_Ainv[0][m][n][index] += bssn_ginv[0][m][p][index]*bssn_ginv[0][n][q][index]*
		                                          bssn_A[0][p][q][index];
		          }
		       }
            }
        }
    
    AATF[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) { 
          AATF[index] += bssn_ginv[0][m][n][index]*bssn_A[0][m][n][index];
        }
    }
    
/* Calc gamma and shift */
    
    gamma[0][1][index] = bssn_gamma_x[index];
    gamma[0][2][index] = bssn_gamma_y[index];
    gamma[0][3][index] = bssn_gamma_z[index];
    
    bssn_phi[0][0][index] = phi[index];
    
    if (CCTK_Equals(gauge_condition,"geodesic")) {
    	alpha[index] = 1.0;
    }
    
    if (CCTK_Equals(gauge_condition,"scale")) {
    	alpha[index] = aa_avg;
    }
    	
    if (CCTK_Equals(gauge_condition,"geod")) {
    	alpha[index] = exp(12.0*phi[index]*sigma);
    	//alpha_x[index] = 12.0*bssn_dk_x[index]*sigma*exp(12.0*phi[index]*sigma);
    	//alpha_y[index] = 12.0*bssn_dk_y[index]*sigma*exp(12.0*phi[index]*sigma);
    	//alpha_z[index] = 12.0*bssn_dk_z[index]*sigma*exp(12.0*phi[index]*sigma);
    }
    	
    if (CCTK_Equals(gauge_condition,"gwave")) {
    	wc = nw*2.0*pi/lz;
    
    	alpha[index] = exp(0.5*Amp*sin(wc*(cctk_time + z[index])));
    	//alpha_z[index] = 0.5*Amp*wc*cos(wc*(cctk_time + z[index]))*alpha[index];
    }
    	
    if (CCTK_Equals(gauge_condition,"gowdy")) {
    	
    	if (CCTK_EQUALS(gowdytyp,"expand")) {
           	t = cctk_time;
        } else {
           	t = exp(-cctk_time);
        }
          	
        bj0 = j0(t);
   	 	bj1 = j1(t);
   	    by0 = y0(2.0*t);
   	    by1 = y1(2.0*t);
   	    
        lam = pow(Amp,2.0)*(-t*bj0*bj1*pow(cos(z[index]),2.0) + 0.5*pow(t*bj0,2.0) + 0.5*pow(t*bj1,2.0));   
	    
    	alpha[index] = exp(0.25*lam)*pow(t,-0.25);
    	//alpha_z[index] = 2.0*pow(Amp,2.0)*t*bj0*bj1*cos(z[index])*sin(z[index])*alpha[index];
    }
		
    lapse[0][0][index] = alpha[index];
    
    shift[0][0][1][index] = betax[index];
    shift[0][0][2][index] = betay[index];
    shift[0][0][3][index] = betaz[index];
    
    bssn_bb[0][1][index] = bssn_bb_x[index];
    bssn_bb[0][2][index] = bssn_bb_y[index];
    bssn_bb[0][3][index] = bssn_bb_z[index];
    
         }
     }
 }
 
/* Calc derivatives */

  if (CCTK_Equals(diff,"Spectral")) {
  
/* lapse derivatives */

  SpecDeriv_Scalar_Derivative2(CCTK_PASS_CTOC, lapse); 
  
/* gamma derivatives */

  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, gamma); 
  
/* bssn inverse metric derivatives */

  SpecDeriv_Scalar_Derivative2(CCTK_PASS_CTOC, bssn_phi);
  SpecDeriv_Tensor_Derivative2(CCTK_PASS_CTOC, bssn_g);
  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, bssn_ginv);

/* shift derivatives */
 
  SpecDeriv_Vector_Derivative2(CCTK_PASS_CTOC, shift);
  SpecDeriv_Vector_Derivative(CCTK_PASS_CTOC, bssn_bb);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, bssn_A);
  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Spectral_1d")) {
  
/* lapse derivatives */

  SpecDeriv_Scalar_Derivative2_1d(CCTK_PASS_CTOC, lapse); 
  
/* gamma derivatives */

  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, gamma); 
  
/* bssn inverse metric derivatives */

  SpecDeriv_Scalar_Derivative2_1d(CCTK_PASS_CTOC, bssn_phi);
  SpecDeriv_Tensor_Derivative2_1d(CCTK_PASS_CTOC, bssn_g);
  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, bssn_ginv);

/* shift derivatives */
 
  SpecDeriv_Vector_Derivative2_1d(CCTK_PASS_CTOC, shift);
  SpecDeriv_Vector_Derivative_1d(CCTK_PASS_CTOC, bssn_bb);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, bssn_A);
  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Finite")) {
  
/* lapse derivatives */

  SpecDeriv_Scalar_Derivative2_FD(CCTK_PASS_CTOC, lapse);  
  
/* gamma derivatives */

  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, gamma); 
  
/* bssn inverse metric derivatives */

  SpecDeriv_Scalar_Derivative2_FD(CCTK_PASS_CTOC, bssn_phi);
  SpecDeriv_Tensor_Derivative2_FD(CCTK_PASS_CTOC, bssn_g); 
  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, bssn_ginv); 

/* shift derivatives */
 
  SpecDeriv_Vector_Derivative2_FD(CCTK_PASS_CTOC, shift);
  SpecDeriv_Vector_Derivative_FD(CCTK_PASS_CTOC, bssn_bb);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, bssn_A);
  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Finite4")) {
  
/* lapse derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, lapse); 
  
/* gamma derivatives */

  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, gamma); 
  
/* bssn inverse metric derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, bssn_phi);
  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_g); 
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_ginv); 

/* shift derivatives */
 
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, shift);
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, bssn_bb);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_A);
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"SCR3")) {
  
/* lapse derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, lapse); 
  
/* gamma derivatives */

  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, gamma); 
  
/* bssn inverse metric derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, bssn_phi);
  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_g); 
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_ginv); 

/* shift derivatives */
 
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, shift);
  SpecDeriv_Vector_Derivative_FD4(CCTK_PASS_CTOC, bssn_bb);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_A);
  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, bssn_K); 
  
  }

for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
	
	if (isnan(lapse[0][1][index])) {lapse[0][1][index] = 0.0;}
	if (isnan(lapse[0][2][index])) {lapse[0][2][index] = 0.0;}
	if (isnan(lapse[0][3][index])) {lapse[0][3][index] = 0.0;}
			
	alpha_x[index] = lapse[0][1][index];
    alpha_y[index] = lapse[0][2][index];
    alpha_z[index] = lapse[0][3][index];
    
    if (isnan(shift[0][1][1][index])) {shift[0][1][1][index] = 0.0;}
    if (isnan(shift[0][2][1][index])) {shift[0][2][1][index] = 0.0;}
    if (isnan(shift[0][3][1][index])) {shift[0][3][1][index] = 0.0;}
    if (isnan(shift[0][1][2][index])) {shift[0][1][2][index] = 0.0;}
    if (isnan(shift[0][2][2][index])) {shift[0][2][2][index] = 0.0;}
    if (isnan(shift[0][3][2][index])) {shift[0][3][2][index] = 0.0;}
    if (isnan(shift[0][1][3][index])) {shift[0][1][3][index] = 0.0;}
    if (isnan(shift[0][2][3][index])) {shift[0][2][3][index] = 0.0;}
    if (isnan(shift[0][3][3][index])) {shift[0][3][3][index] = 0.0;}
    
    betax1[index] = shift[0][1][1][index];
    betax2[index] = shift[0][2][1][index];
    betax3[index] = shift[0][3][1][index];
    betay1[index] = shift[0][1][2][index];
    betay2[index] = shift[0][2][2][index];
    betay3[index] = shift[0][3][2][index];
    betaz1[index] = shift[0][1][3][index];
    betaz2[index] = shift[0][2][3][index];
    betaz3[index] = shift[0][3][3][index];

	if (isnan(bssn_phi[0][1][index])) {bssn_phi[0][1][index] = 0.0;}
	if (isnan(bssn_phi[0][2][index])) {bssn_phi[0][2][index] = 0.0;}
	if (isnan(bssn_phi[0][3][index])) {bssn_phi[0][3][index] = 0.0;}
	
    bssn_dk_x[index] = bssn_phi[0][1][index];
    bssn_dk_y[index] = bssn_phi[0][2][index];
    bssn_dk_z[index] = bssn_phi[0][3][index];
    
    if (isnan(bssn_g[0][1][1][1][index])) {bssn_g[0][1][1][1][index] = 0.0;}
    if (isnan(bssn_g[0][2][1][1][index])) {bssn_g[0][2][1][1][index] = 0.0;}
    if (isnan(bssn_g[0][3][1][1][index])) {bssn_g[0][3][1][1][index] = 0.0;}
    if (isnan(bssn_g[0][1][2][2][index])) {bssn_g[0][1][2][2][index] = 0.0;}
    if (isnan(bssn_g[0][2][2][2][index])) {bssn_g[0][2][2][2][index] = 0.0;}
    if (isnan(bssn_g[0][3][2][2][index])) {bssn_g[0][3][2][2][index] = 0.0;}
    if (isnan(bssn_g[0][1][3][3][index])) {bssn_g[0][1][3][3][index] = 0.0;}
    if (isnan(bssn_g[0][2][3][3][index])) {bssn_g[0][2][3][3][index] = 0.0;}
    if (isnan(bssn_g[0][3][3][3][index])) {bssn_g[0][3][3][3][index] = 0.0;}
    if (isnan(bssn_g[0][1][1][2][index])) {bssn_g[0][1][1][2][index] = 0.0;}
    if (isnan(bssn_g[0][2][1][2][index])) {bssn_g[0][2][1][2][index] = 0.0;}
    if (isnan(bssn_g[0][3][1][2][index])) {bssn_g[0][3][1][2][index] = 0.0;}
    if (isnan(bssn_g[0][1][1][3][index])) {bssn_g[0][1][1][3][index] = 0.0;}
    if (isnan(bssn_g[0][2][1][3][index])) {bssn_g[0][2][1][3][index] = 0.0;}
    if (isnan(bssn_g[0][3][1][3][index])) {bssn_g[0][3][1][3][index] = 0.0;}
    if (isnan(bssn_g[0][1][2][3][index])) {bssn_g[0][1][2][3][index] = 0.0;}
    if (isnan(bssn_g[0][2][2][3][index])) {bssn_g[0][2][2][3][index] = 0.0;}
    if (isnan(bssn_g[0][3][2][3][index])) {bssn_g[0][3][2][3][index] = 0.0;}
    
    bssn_dkij_xxx[index] = bssn_g[0][1][1][1][index];
    bssn_dkij_yxx[index] = bssn_g[0][2][1][1][index];
    bssn_dkij_zxx[index] = bssn_g[0][3][1][1][index];
    bssn_dkij_xyy[index] = bssn_g[0][1][2][2][index];
    bssn_dkij_yyy[index] = bssn_g[0][2][2][2][index];
    bssn_dkij_zyy[index] = bssn_g[0][3][2][2][index];
    bssn_dkij_xzz[index] = bssn_g[0][1][3][3][index];
    bssn_dkij_yzz[index] = bssn_g[0][2][3][3][index];
    bssn_dkij_zzz[index] = bssn_g[0][3][3][3][index];
    bssn_dkij_xxy[index] = bssn_g[0][1][1][2][index];
    bssn_dkij_yxy[index] = bssn_g[0][2][1][2][index];
    bssn_dkij_zxy[index] = bssn_g[0][3][1][2][index];
    bssn_dkij_xxz[index] = bssn_g[0][1][1][3][index];
    bssn_dkij_yxz[index] = bssn_g[0][2][1][3][index];
    bssn_dkij_zxz[index] = bssn_g[0][3][1][3][index];
    bssn_dkij_xyz[index] = bssn_g[0][1][2][3][index];
    bssn_dkij_yyz[index] = bssn_g[0][2][2][3][index];
    bssn_dkij_zyz[index] = bssn_g[0][3][2][3][index];
			
	conf = exp(4.0*phi[index]); 
	conf_inv = exp(-4.0*phi[index]); 

/* Calc metric */
    
    gg4[0][0] = gtt[index];
    gg4[0][1] = gtx[index];
    gg4[0][2] = gty[index]; 
    gg4[0][3] = gtz[index];
    gg4[1][0] = gg4[0][1];
    gg4[2][0] = gg4[0][2];
    gg4[3][0] = gg4[0][3];
    gg4[1][1] = gxx[index];
    gg4[2][2] = gyy[index];
    gg4[3][3] = gzz[index];
    gg4[1][2] = gxy[index];
    gg4[1][3] = gxz[index];
    gg4[2][3] = gyz[index];
    gg4[2][1] = gxy[index];
    gg4[3][1] = gxz[index];
    gg4[3][2] = gyz[index];
    			
/* Calc curvature */
     
     for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           Ricci[m][n] = 0.0;
           Ricci_phi[m][n] = 0.0;
           for(p=1; p < 4; p++) {
              bssn_Chris1[m][n][p] = 0.0;
              bssn_Chris2[m][n][p] = 0.0;
           }
        }
     }
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           for(p=1; p < 4; p++) {
              for(q=1; q < 4; q++) {
     				bssn_Chris1[m][n][p] += 0.5*bssn_ginv[0][m][q][index]*(bssn_g[0][n][p][q][index] + 
     										bssn_g[0][p][n][q][index] - bssn_g[0][q][n][p][index]);
               }
           }
        }
     }
     
     for(a=1; a < 4; a++) {
        for(b=1; b < 4; b++) {
           for(c=1; c < 4; c++) {
           	   for(d=1; d < 4; d++) {
              		bssn_Chris2[a][b][c] += bssn_g[0][0][a][d][index]*bssn_Chris1[d][b][c];  
     		   }
           }
        }
     }

     Chris_xxx[index] = bssn_Chris1[1][1][1];
     Chris_xxy[index] = bssn_Chris1[1][1][2];
     Chris_xxz[index] = bssn_Chris1[1][1][3];
     Chris_xyx[index] = bssn_Chris1[1][2][1];
     Chris_xyy[index] = bssn_Chris1[1][2][2];
     Chris_xyz[index] = bssn_Chris1[1][2][3];
     Chris_xzx[index] = bssn_Chris1[1][3][1];
     Chris_xzy[index] = bssn_Chris1[1][3][2];
     Chris_xzz[index] = bssn_Chris1[1][3][3];
     Chris_yxx[index] = bssn_Chris1[2][1][1];
     Chris_yxy[index] = bssn_Chris1[2][1][2];
     Chris_yxz[index] = bssn_Chris1[2][1][3];
     Chris_yyx[index] = bssn_Chris1[2][2][1];
     Chris_yyy[index] = bssn_Chris1[2][2][2];
     Chris_yyz[index] = bssn_Chris1[2][2][3];
     Chris_yzx[index] = bssn_Chris1[2][3][1];
     Chris_yzy[index] = bssn_Chris1[2][3][2];
     Chris_yzz[index] = bssn_Chris1[2][3][3];
     Chris_zxx[index] = bssn_Chris1[3][1][1];
     Chris_zxy[index] = bssn_Chris1[3][1][2];
     Chris_zxz[index] = bssn_Chris1[3][1][3];
     Chris_zyx[index] = bssn_Chris1[3][2][1];
     Chris_zyy[index] = bssn_Chris1[3][2][2];
     Chris_zyz[index] = bssn_Chris1[3][2][3];
     Chris_zzx[index] = bssn_Chris1[3][3][1];
     Chris_zzy[index] = bssn_Chris1[3][3][2];
     Chris_zzz[index] = bssn_Chris1[3][3][3];
     
/* Calc Ricci Tensor */

    for(a=1; a < 4; a++) {
        for(b=1; b < 4; b++) {
           for(c=1; c < 4; c++) {
              for(l=1; l < 4; l++) {
                 for(m=1; m < 4; m++) {
     
     			   Ricci[a][b] += bssn_ginv[0][c][l][index]*(bssn_Chris1[m][a][c]*bssn_Chris2[m][l][b] +
     			                  bssn_Chris1[m][c][b]*bssn_Chris2[a][m][l] +
     			                  bssn_Chris1[m][c][a]*bssn_Chris2[b][m][l]);
                   }
                   
                   Ricci[a][b] += -0.5*bssn_ginv[0][c][l][index]*bssn_g[c][l][a][b][index];
               }
               
               Ricci[a][b] += 0.5*gamma[0][c][index]*bssn_Chris2[a][b][c] +
     			              0.5*gamma[0][c][index]*bssn_Chris2[b][a][c] +
     			              0.5*bssn_g[0][0][c][a][index]*gamma[b][c][index] + 
     				          0.5*bssn_g[0][0][c][b][index]*gamma[a][c][index];
           }
        }
     }
   
   Ricci_xx[index] = Ricci[1][1];
   Ricci_yy[index] = Ricci[2][2];
   Ricci_zz[index] = Ricci[3][3];
   Ricci_xy[index] = Ricci[1][2];
   Ricci_xz[index] = Ricci[1][3];
   Ricci_yz[index] = Ricci[2][3];
   
   
/* Calc Ricci Phi */

		ddphi = 0.0;
   
   		for(m=1; m < 4; m++) {
     		for(n=1; n < 4; n++) {
   				for(l=1; l < 4; l++) {
   					ddphi += -bssn_ginv[0][m][n][index]*bssn_Chris1[l][m][n]*bssn_phi[0][l][index];
   				}
   			ddphi += bssn_ginv[0][m][n][index]*bssn_phi[m][n][index];
   	 		}
   		}
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
            	for(l=1; l < 4; l++) {
            		Ricci_phi[m][n] += 2.0*bssn_Chris1[l][m][n]*bssn_phi[0][l][index];
            			for(a=1; a < 4; a++) {
            				Ricci_phi[m][n] += -4.0*bssn_g[0][0][m][n][index]*bssn_ginv[0][l][a][index]*bssn_phi[0][a][index]*bssn_phi[0][l][index];
            			}
            	}
        		Ricci_phi[m][n] += -2.0*bssn_phi[m][n][index] - 2.0*bssn_g[0][0][m][n][index]*ddphi +
        		                   4.0*bssn_phi[0][m][index]*bssn_phi[0][n][index];
            }
        }
	
   Ricci_phi_xx[index] = Ricci_phi[1][1];
   Ricci_phi_yy[index] = Ricci_phi[2][2];
   Ricci_phi_zz[index] = Ricci_phi[3][3];
   Ricci_phi_xy[index] = Ricci_phi[1][2];
   Ricci_phi_xz[index] = Ricci_phi[1][3];
   Ricci_phi_yz[index] = Ricci_phi[2][3];
   
   Ricci_Scalar = 0.0;
   
   for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
             Ricci_Scalar += conf_inv*bssn_ginv[0][m][n][index]*(Ricci[m][n]+Ricci_phi[m][n]);
        }
   }
   
/* Stress Tensor Terms */

   if (add_DM) {
   
   Ttotal[0][0] = Ttt_mhd[index] + Ttt_dark[index] + Ttt_vac[index];
   Ttotal[0][1] = Ttx_mhd[index] + Ttx_dark[index] + Ttx_vac[index];
   Ttotal[0][2] = Tty_mhd[index] + Tty_dark[index] + Tty_vac[index];
   Ttotal[0][3] = Ttz_mhd[index] + Ttz_dark[index] + Ttz_vac[index];
   Ttotal[1][1] = Txx_mhd[index] + Txx_dark[index] + Txx_vac[index];
   Ttotal[1][2] = Txy_mhd[index] + Txy_dark[index] + Txy_vac[index];
   Ttotal[1][3] = Txz_mhd[index] + Txz_dark[index] + Txz_vac[index];
   Ttotal[2][2] = Tyy_mhd[index] + Tyy_dark[index] + Tyy_vac[index];
   Ttotal[2][3] = Tyz_mhd[index] + Tyz_dark[index] + Tyz_vac[index];
   Ttotal[3][3] = Tzz_mhd[index] + Tzz_dark[index] + Tzz_vac[index];
   
   }
   
   if (!add_DM) {
   
   Ttotal[0][0] = Ttt_mhd[index];
   Ttotal[0][1] = Ttx_mhd[index];
   Ttotal[0][2] = Tty_mhd[index];
   Ttotal[0][3] = Ttz_mhd[index];
   Ttotal[1][1] = Txx_mhd[index];
   Ttotal[1][2] = Txy_mhd[index];
   Ttotal[1][3] = Txz_mhd[index];
   Ttotal[2][2] = Tyy_mhd[index];
   Ttotal[2][3] = Tyz_mhd[index];
   Ttotal[3][3] = Tzz_mhd[index];
   
   }
   
   Ttotal[1][0] = Ttotal[0][1];
   Ttotal[2][0] = Ttotal[0][2];
   Ttotal[3][0] = Ttotal[0][3];
   
   Ttotal[2][1] = Ttotal[1][2];
   Ttotal[3][1] = Ttotal[1][3];
   Ttotal[3][2] = Ttotal[2][3];
   
   for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) { 
           T_total[m][n] = 0.0;
        }
   }
  
   for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) { 
           T_total[0][0] += Ttotal[m][n]*gg4[0][m]*gg4[0][n];
           T_total[0][1] += Ttotal[m][n]*gg4[0][m]*gg4[1][n];
           T_total[0][2] += Ttotal[m][n]*gg4[0][m]*gg4[2][n];
           T_total[0][3] += Ttotal[m][n]*gg4[0][m]*gg4[3][n];
           T_total[1][1] += Ttotal[m][n]*gg4[1][m]*gg4[1][n];
           T_total[2][2] += Ttotal[m][n]*gg4[2][m]*gg4[2][n];
           T_total[3][3] += Ttotal[m][n]*gg4[3][m]*gg4[3][n];
           T_total[1][2] += Ttotal[m][n]*gg4[1][m]*gg4[2][n];
           T_total[1][3] += Ttotal[m][n]*gg4[1][m]*gg4[3][n];
           T_total[2][3] += Ttotal[m][n]*gg4[2][m]*gg4[3][n];
        }
   }
           
   T_total[1][0] = T_total[0][1];
   T_total[2][0] = T_total[0][2];
   T_total[3][0] = T_total[0][3];
   T_total[2][1] = T_total[1][2];
   T_total[3][1] = T_total[1][3];
   T_total[3][2] = T_total[2][3];
   
   rho_total[index] = 8.0*pi*pow(alpha[index],2.0)*Ttotal[0][0];
     
   for(m=0; m < 4; m++) {
   	Si[m] = 0.0;
   }
     
   for(n=1; n < 4; n++) {
   	for(m=0; m < 4; m++) {
   		Si[n] += 8.0*pi*alpha[index]*gg4[n][m]*Ttotal[0][m];
   	}
   }
   
   Si_x[index] = Si[1];
   Si_y[index] = Si[2];
   Si_z[index] = Si[3];
   
   Sij_xx[index] = 8.0*pi*T_total[1][1];
   Sij_yy[index] = 8.0*pi*T_total[2][2];
   Sij_zz[index] = 8.0*pi*T_total[3][3];
   Sij_xy[index] = 8.0*pi*T_total[1][2];
   Sij_xz[index] = 8.0*pi*T_total[1][3];
   Sij_yz[index] = 8.0*pi*T_total[2][3];
   
/* Calculate Constraints */

    AijAij = 0.0;
        
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
		    AijAij += bssn_Ainv[0][m][n][index]*bssn_A[0][m][n][index];
    	}
    }
    
    ham[index] = Ricci_Scalar + 2.0*pow(TrK[index],2.0)/3.0 - AijAij - 2.0*rho_total[index];
    
    if (rho_total[index] != 0.0) {
    	ham[index] = ham[index]/(2.0*rho_total[index]);
    }
    
    momx[index] = 0.0;
    momy[index] = 0.0;
    momz[index] = 0.0;
    
    for(m=1; m < 4; m++) {
	    momx[index] += bssn_Ainv[m][1][m][index] + 6.0*bssn_Ainv[0][1][m][index]*bssn_phi[0][m][index] - 
	                   bssn_ginv[0][1][m][index]*(2.0*bssn_K[m][index]/3.0 + Si[m]);
	    momy[index] += bssn_Ainv[m][2][m][index] + 6.0*bssn_Ainv[0][2][m][index]*bssn_phi[0][m][index] - 
	                   bssn_ginv[0][2][m][index]*(2.0*bssn_K[m][index]/3.0 + Si[m]);
	    momz[index] += bssn_Ainv[m][3][m][index] + 6.0*bssn_Ainv[0][3][m][index]*bssn_phi[0][m][index] - 
	                   bssn_ginv[0][3][m][index]*(2.0*bssn_K[m][index]/3.0 + Si[m]);
	}
	
	for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           momx[index] += bssn_Ainv[0][m][n][index]*bssn_Chris1[1][m][n];
           momy[index] += bssn_Ainv[0][m][n][index]*bssn_Chris1[2][m][n];
           momz[index] += bssn_Ainv[0][m][n][index]*bssn_Chris1[3][m][n];
        }
   }
	
	
	gamma_x = 0.0;
	gamma_y = 0.0;
	gamma_z = 0.0;
	
	for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           gamma_x += bssn_ginv[0][m][n][index]*bssn_Chris1[1][m][n];
           gamma_y += bssn_ginv[0][m][n][index]*bssn_Chris1[2][m][n];
           gamma_z += bssn_ginv[0][m][n][index]*bssn_Chris1[3][m][n];
        }
   }
   
   //bssn_gamma_x[index] = gamma_x;
   //bssn_gamma_y[index] = gamma_y;
   //bssn_gamma_z[index] = gamma_z;
    
    dGx[index] = 0.0; 
    dGy[index] = 0.0;
    dGz[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
        
           dGx[index] += bssn_ginv[0][m][n][index]*bssn_Chris1[1][m][n];
           dGy[index] += bssn_ginv[0][m][n][index]*bssn_Chris1[2][m][n];
           dGz[index] += bssn_ginv[0][m][n][index]*bssn_Chris1[3][m][n];
           
        }
   }
    
	dGx[index] += -bssn_gamma_x[index];
    dGy[index] += -bssn_gamma_y[index];
    dGz[index] += -bssn_gamma_z[index];
    
    dHx[index] = 0.0; 
    dHy[index] = 0.0;
    dHz[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        dHx[index] += -bssn_ginv[m][1][m][index];
        dHy[index] += -bssn_ginv[m][2][m][index];
        dHz[index] += -bssn_ginv[m][3][m][index];
    }
    
    dHx[index] += -bssn_gamma_x[index];
    dHy[index] += -bssn_gamma_y[index];
    dHz[index] += -bssn_gamma_z[index];
    
    if (CCTK_Equals(slicing,"scale") && (! fix_lapse) && (TrK[index] < 0.0)) {
    	Ek[index] = (lapse_avg - aa_avg)/aa_avg;
    } else {
    	Ek[index] = 0.0;
    }
    
    
    AATF[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) { 
          AATF[index] += bssn_ginv[0][m][n][index]*bssn_A[0][m][n][index];
        }
    }
    
        }
	}
}

/* Begin Evolution */
   
/* Set up shorthands */

    istart = cctk_nghostzones[0];
    jstart = cctk_nghostzones[1];
    kstart = cctk_nghostzones[2];
  
    iend = cctk_lsh[0] - cctk_nghostzones[0];
    jend = cctk_lsh[1] - cctk_nghostzones[1];
    kend = cctk_lsh[2] - cctk_nghostzones[2];

/* Do the */

for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );

/* Calc shift squared */
  
  	 conf = exp(4.0*phi[index]);
  	 conf_inv = exp(-4.0*phi[index]); 
  	 
     bssn_Chris1[1][1][1] = Chris_xxx[index];
     bssn_Chris1[1][1][2] = Chris_xxy[index];
     bssn_Chris1[1][1][3] = Chris_xxz[index];
     bssn_Chris1[1][2][1] = Chris_xyx[index];
     bssn_Chris1[1][2][2] = Chris_xyy[index];
     bssn_Chris1[1][2][3] = Chris_xyz[index];
     bssn_Chris1[1][3][1] = Chris_xzx[index];
     bssn_Chris1[1][3][2] = Chris_xzy[index];
     bssn_Chris1[1][3][3] = Chris_xzz[index];
     bssn_Chris1[2][1][1] = Chris_yxx[index];
     bssn_Chris1[2][1][2] = Chris_yxy[index];
     bssn_Chris1[2][1][3] = Chris_yxz[index];
     bssn_Chris1[2][2][1] = Chris_yyx[index];
     bssn_Chris1[2][2][2] = Chris_yyy[index];
     bssn_Chris1[2][2][3] = Chris_yyz[index];
     bssn_Chris1[2][3][1] = Chris_yzx[index];
     bssn_Chris1[2][3][2] = Chris_yzy[index];
     bssn_Chris1[2][3][3] = Chris_yzz[index];
     bssn_Chris1[3][1][1] = Chris_zxx[index];
     bssn_Chris1[3][1][2] = Chris_zxy[index];
     bssn_Chris1[3][1][3] = Chris_zxz[index];
     bssn_Chris1[3][2][1] = Chris_zyx[index];
     bssn_Chris1[3][2][2] = Chris_zyy[index];
     bssn_Chris1[3][2][3] = Chris_zyz[index];
     bssn_Chris1[3][3][1] = Chris_zzx[index];
     bssn_Chris1[3][3][2] = Chris_zzy[index];
     bssn_Chris1[3][3][3] = Chris_zzz[index];
     
     for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           for(p=1; p < 4; p++) {
     			Chris[m][n][p] = bssn_Chris1[m][n][p];
           }
        }
     }
     
     for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           for(p=1; p < 4; p++) {
              for(q=1; q < 4; q++) {
     				Chris[m][n][p] += -2.0*bssn_g[0][0][n][p][index]*bssn_ginv[0][m][q][index]*bssn_phi[0][q][index];
               }
               if (m == n) {
               		Chris[m][n][p] += 2.0*bssn_phi[0][p][index];
               }
               if (m == p) {
               		Chris[m][n][p] += 2.0*bssn_phi[0][n][index];
               }
           }
        }
     }
     
     Ricci[1][1] = Ricci_xx[index];
   	 Ricci[2][2] = Ricci_yy[index];
     Ricci[3][3] = Ricci_zz[index];
     Ricci[1][2] = Ricci_xy[index];
     Ricci[1][3] = Ricci_xz[index];
     Ricci[2][3] = Ricci_yz[index];
     Ricci[2][1] = Ricci_xy[index];
     Ricci[3][1] = Ricci_xz[index];
     Ricci[3][2] = Ricci_yz[index];
     
     Ricci_phi[1][1] = Ricci_phi_xx[index];
   	 Ricci_phi[2][2] = Ricci_phi_yy[index];
     Ricci_phi[3][3] = Ricci_phi_zz[index];
     Ricci_phi[1][2] = Ricci_phi_xy[index];
     Ricci_phi[1][3] = Ricci_phi_xz[index];
     Ricci_phi[2][3] = Ricci_phi_yz[index];
     Ricci_phi[2][1] = Ricci_phi_xy[index];
     Ricci_phi[3][1] = Ricci_phi_xz[index];
     Ricci_phi[3][2] = Ricci_phi_yz[index];
    
     Sij[1][1] = Sij_xx[index];
     Sij[2][2] = Sij_yy[index];
     Sij[3][3] = Sij_zz[index];
     Sij[1][2] = Sij_xy[index];
     Sij[1][3] = Sij_xz[index];
     Sij[2][3] = Sij_yz[index];
     Sij[2][1] = Sij_xy[index];
     Sij[3][1] = Sij_xz[index];
     Sij[3][2] = Sij_yz[index];
    
     Si[1] = Si_x[index];
     Si[2] = Si_y[index];
     Si[3] = Si_z[index];
     
     /* Calc metric */
    
    gg4[0][0] = gtt[index];
    gg4[0][1] = gtx[index];
    gg4[0][2] = gty[index]; 
    gg4[0][3] = gtz[index];
    gg4[1][0] = gg4[0][1];
    gg4[2][0] = gg4[0][2];
    gg4[3][0] = gg4[0][3];
    gg4[1][1] = gxx[index];
    gg4[2][2] = gyy[index];
    gg4[3][3] = gzz[index];
    gg4[1][2] = gxy[index];
    gg4[1][3] = gxz[index];
    gg4[2][3] = gyz[index];
    gg4[2][1] = gxy[index];
    gg4[3][1] = gxz[index];
    gg4[3][2] = gyz[index];
    
    igg3[0][0] = igtt[index];
    igg3[0][1] = igtx[index];
    igg3[0][2] = igty[index];
    igg3[0][3] = igtz[index];
    igg3[1][1] = igxx[index];
    igg3[2][2] = igyy[index];
    igg3[3][3] = igzz[index];
    igg3[1][2] = igxy[index];
    igg3[1][3] = igxz[index];
    igg3[2][3] = igyz[index];
    igg3[0][1] = igtx[index];
    igg3[0][2] = igty[index];
    igg3[0][3] = igtz[index];
    igg3[2][1] = igxy[index];
    igg3[3][1] = igxz[index];
    igg3[3][2] = igyz[index];
     
     for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) {
        	for (p=1; p < 4; p++) {
        		bssn_beta[m][n][p] = shift[m][n][p][index];
        	}
        }
     }
     
     for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) {
        	alphai[m][n] = lapse[m][n][index];
        }
     }
     
    S_Scalar = 0.0;
		
	for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
		    S_Scalar += igg3[m][n]*Sij[m][n];
        }
    }
    
    STFxx = Sij[1][1] - gg4[1][1]*S_Scalar/3.0;
	STFyy = Sij[2][2] - gg4[2][2]*S_Scalar/3.0;
	STFzz = Sij[3][3] - gg4[3][3]*S_Scalar/3.0;
	STFxy = Sij[1][2] - gg4[1][2]*S_Scalar/3.0;
	STFxz = Sij[1][3] - gg4[1][3]*S_Scalar/3.0;
	STFyz = Sij[2][3] - gg4[2][3]*S_Scalar/3.0;
        
    AijAij = 0.0;
        
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
		    AijAij += bssn_Ainv[0][m][n][index]*bssn_A[0][m][n][index];
        }
    }
    
/* Calc RHS of K */
		
		Lie_K = 0.0;
        
        for(m=1; m < 4; m++) {
            Lie_K += shift[0][0][m][index]*bssn_K[m][index];
        }
        
        bssn_K_dt[index] = 0.0;
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
            	for(p=1; p < 4; p++) {
            		bssn_K_dt[index] += igg3[m][n]*Chris[p][m][n]*alphai[0][p];
            	}
            }
        }
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
				bssn_K_dt[index] += -igg3[m][n]*alphai[m][n]; 
		    }
		}
		            
		bssn_K_dt[index] += alphai[0][0]*(AijAij + pow(TrK[index],2.0)/3.0) + 0.5*alphai[0][0]*
							(rho_total[index] + S_Scalar) + Lie_K;
		
if (CCTK_Equals(slicing,"scale") && (! fix_lapse)) {

/* Calc RHS of lapse */

		alpha_dt[index] = 0.0;
		
		if (K_avg < 0) { 
		for(m=1; m < 4; m++) {
			alpha_dt[index] += bssn_beta[0][0][m]*alphai[0][m] + (1.0/3.0)*alphai[0][0]*bssn_beta[0][m][m];
		}
		
		alpha_dt[index] += -(1.0/3.0)*pow(aa_avg,2.0)*K_avg;
		} 
		
		//alpha_dt[index] += -(1.0/3.0)*pow(alphai[0][0],2.0)*TrK[index];
	
}

if (CCTK_Equals(slicing,"1+log") && (! fix_lapse)) {

/* Calc RHS of lapse */

		alpha_dt[index] = 0.0;
		
		for(m=1; m < 4; m++) {
			alpha_dt[index] += bssn_beta[0][0][m]*alphai[0][m];
		}
		
		alpha_dt[index] += -2.0*alphai[0][0]*TrK[index];
		
}

if (CCTK_Equals(slicing,"harmonic") && (! fix_lapse)) {

/* Calc RHS of lapse */

		alpha_dt[index] = 0.0;
		
		for(m=1; m < 4; m++) {
			alpha_dt[index] += bssn_beta[0][0][m]*alphai[0][m];
		}
		
		alpha_dt[index] += -pow(alphai[0][0],2.0)*TrK[index];
		
}

		if (fix_lapse) {
			alpha_dt[index] = 0.0;
		}

/* Calc RHS of phi */
		
		Lie_phi = 0.0;
        
        for(m=1; m < 4; m++) {
            Lie_phi += shift[0][0][m][index]*bssn_phi[0][m][index] + bssn_beta[0][m][m]/6.0;
        }
		
		phi_dt[index] = -alphai[0][0]*TrK[index]/6.0 + Lie_phi;
                               
/* Calc RHS of metric */
		
		Lie_gxx = 0.0;
		Lie_gyy = 0.0;
		Lie_gzz = 0.0;
		Lie_gxy = 0.0;
		Lie_gxz = 0.0;
		Lie_gyz = 0.0;
        
        for(m=1; m < 4; m++) {
            Lie_gxx += bssn_beta[0][0][m]*bssn_g[0][m][1][1][index] - 2.0*bssn_g[0][0][1][1][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][1][index]*bssn_beta[1][0][m] + 
                      bssn_g[0][0][m][1][index]*bssn_beta[1][0][m];
            Lie_gyy += bssn_beta[0][0][m]*bssn_g[0][m][2][2][index] - 2.0*bssn_g[0][0][2][2][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][2][index]*bssn_beta[2][0][m] + 
                      bssn_g[0][0][m][2][index]*bssn_beta[2][0][m];
            Lie_gzz += bssn_beta[0][0][m]*bssn_g[0][m][3][3][index] - 2.0*bssn_g[0][0][3][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][3][index]*bssn_beta[3][0][m] + 
                      bssn_g[0][0][m][3][index]*bssn_beta[3][0][m];
            Lie_gxy += bssn_beta[0][0][m]*bssn_g[0][m][1][2][index] - 2.0*bssn_g[0][0][1][2][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][1][index]*bssn_beta[2][0][m] + 
                      bssn_g[0][0][m][2][index]*bssn_beta[1][0][m];
            Lie_gxz += bssn_beta[0][0][m]*bssn_g[0][m][1][3][index] - 2.0*bssn_g[0][0][1][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][1][index]*bssn_beta[3][0][m] + 
                      bssn_g[0][0][m][3][index]*bssn_beta[1][0][m];
            Lie_gyz += bssn_beta[0][0][m]*bssn_g[0][m][2][3][index] - 2.0*bssn_g[0][0][2][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_g[0][0][m][2][index]*bssn_beta[3][0][m] + 
                      bssn_g[0][0][m][3][index]*bssn_beta[2][0][m];   
        }
        
		bssn_gxx_dt[index] = -2.0*alphai[0][0]*bssn_Axx[index] + Lie_gxx;
		bssn_gyy_dt[index] = -2.0*alphai[0][0]*bssn_Ayy[index] + Lie_gyy;
		bssn_gzz_dt[index] = -2.0*alphai[0][0]*bssn_Azz[index] + Lie_gzz;
		bssn_gxy_dt[index] = -2.0*alphai[0][0]*bssn_Axy[index] + Lie_gxy;
		bssn_gxz_dt[index] = -2.0*alphai[0][0]*bssn_Axz[index] + Lie_gxz;
		bssn_gyz_dt[index] = -2.0*alphai[0][0]*bssn_Ayz[index] + Lie_gyz;

/* Calc AilAlj */
		
		AilAljxx = 0.0; 
		AilAljyy = 0.0; 
		AilAljzz = 0.0; 
		AilAljxy = 0.0; 
		AilAljxz = 0.0; 
		AilAljyz = 0.0; 
		
		for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
		           AilAljxx += bssn_ginv[0][n][m][index]*bssn_A[0][1][n][index]*bssn_A[0][m][1][index];
		           AilAljyy += bssn_ginv[0][n][m][index]*bssn_A[0][2][n][index]*bssn_A[0][m][2][index];
		           AilAljzz += bssn_ginv[0][n][m][index]*bssn_A[0][3][n][index]*bssn_A[0][m][3][index];
		           AilAljxy += bssn_ginv[0][n][m][index]*bssn_A[0][1][n][index]*bssn_A[0][m][2][index];
		           AilAljxz += bssn_ginv[0][n][m][index]*bssn_A[0][1][n][index]*bssn_A[0][m][3][index];
		           AilAljyz += bssn_ginv[0][n][m][index]*bssn_A[0][2][n][index]*bssn_A[0][m][3][index];
            }
        }
        
/* Calc RHS of Aij */
        
		Lie_Axx = 0.0;
		Lie_Ayy = 0.0;
		Lie_Azz = 0.0;
		Lie_Axy = 0.0;
		Lie_Axz = 0.0;
		Lie_Ayz = 0.0;
        
        for(m=1; m < 4; m++) {
            Lie_Axx += bssn_beta[0][0][m]*bssn_A[m][1][1][index] - 2.0*bssn_A[0][1][1][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][1][index]*bssn_beta[1][0][m] + 
                      bssn_A[0][m][1][index]*bssn_beta[1][0][m];
            Lie_Ayy += bssn_beta[0][0][m]*bssn_A[m][2][2][index] - 2.0*bssn_A[0][2][2][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][2][index]*bssn_beta[2][0][m] + 
                      bssn_A[0][m][2][index]*bssn_beta[2][0][m];
            Lie_Azz += bssn_beta[0][0][m]*bssn_A[m][3][3][index] - 2.0*bssn_A[0][3][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][3][index]*bssn_beta[3][0][m] + 
                      bssn_A[0][m][3][index]*bssn_beta[3][0][m];
            Lie_Axy += bssn_beta[0][0][m]*bssn_A[m][1][2][index] - 2.0*bssn_A[0][1][2][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][1][index]*bssn_beta[2][0][m] + 
                      bssn_A[0][m][2][index]*bssn_beta[1][0][m];
            Lie_Axz += bssn_beta[0][0][m]*bssn_A[m][1][3][index] - 2.0*bssn_A[0][1][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][1][index]*bssn_beta[3][0][m] + 
                      bssn_A[0][m][3][index]*bssn_beta[1][0][m];
            Lie_Ayz += bssn_beta[0][0][m]*bssn_A[m][2][3][index] - 2.0*bssn_A[0][2][3][index] * 
                      bssn_beta[m][0][m]/3.0 + bssn_A[0][m][2][index]*bssn_beta[3][0][m] + 
                      bssn_A[0][m][3][index]*bssn_beta[2][0][m];          
        }
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
        		Aij[m][n] = 0.0;
        		AijTF[m][n] = 0.0;
        	}
        }
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
        		Aij[m][n] += alphai[0][0]*(Ricci[m][n]+Ricci_phi[m][n]) - alphai[m][n];
        		for(p=1; p < 4; p++) {
        			Aij[m][n] += Chris[p][m][n]*alphai[0][p];
        		}
            }
        }
        
        Aij_Scalar = 0.0;
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
        		Aij_Scalar += igg3[m][n]*Aij[m][n]; 
        	}
        }
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
        		AijTF[m][n] = Aij[m][n] - gg4[m][n]*Aij_Scalar/3.0;
        	}
        }
        
		bssn_Axx_dt[index] = conf_inv*(AijTF[1][1] - alphai[0][0]*STFxx) +
		                     alphai[0][0]*(TrK[index]*bssn_A[0][1][1][index] - 2.0*AilAljxx) + Lie_Axx;
	    bssn_Ayy_dt[index] = conf_inv*(AijTF[2][2] - alphai[0][0]*STFyy) +  
	                         alphai[0][0]*(TrK[index]*bssn_A[0][2][2][index] - 2.0*AilAljyy) + Lie_Ayy;
	    bssn_Azz_dt[index] = conf_inv*(AijTF[3][3] - alphai[0][0]*STFzz) + 
	                         alphai[0][0]*(TrK[index]*bssn_A[0][3][3][index] - 2.0*AilAljzz) + Lie_Azz;
		bssn_Axy_dt[index] = conf_inv*(AijTF[1][2] - alphai[0][0]*STFxy) + 
		                     alphai[0][0]*(TrK[index]*bssn_A[0][1][2][index] - 2.0*AilAljxy) + Lie_Axy;
	    bssn_Axz_dt[index] = conf_inv*(AijTF[1][3] - alphai[0][0]*STFxz) + 
	                         alphai[0][0]*(TrK[index]*bssn_A[0][1][3][index] - 2.0*AilAljxz) + Lie_Axz;
	    bssn_Ayz_dt[index] = conf_inv*(AijTF[2][3] - alphai[0][0]*STFyz) + 
	                         alphai[0][0]*(TrK[index]*bssn_A[0][2][3][index] - 2.0*AilAljyz) + Lie_Ayz;
                                   
/* Calc RHS of Gamma */
        
        bssn_gamma_x_RHS1  = 0.0;
        bssn_gamma_y_RHS1  = 0.0;
        bssn_gamma_z_RHS1  = 0.0;
        bssn_gamma_x_RHS2  = 0.0;
        bssn_gamma_y_RHS2  = 0.0;
        bssn_gamma_z_RHS2  = 0.0;
        
		for(m=1; m < 4; m++) { 
		    bssn_gamma_x_RHS1  += -2.0*(alphai[0][0]*bssn_Ainv[m][1][m][index] + bssn_Ainv[0][1][m][index]*alphai[0][m]);
        	bssn_gamma_y_RHS1  += -2.0*(alphai[0][0]*bssn_Ainv[m][2][m][index] + bssn_Ainv[0][2][m][index]*alphai[0][m]);
        	bssn_gamma_z_RHS1  += -2.0*(alphai[0][0]*bssn_Ainv[m][3][m][index] + bssn_Ainv[0][3][m][index]*alphai[0][m]);
		    
		    /*bssn_gamma_x_RHS1 += -2.0*bssn_Ainv[0][1][m][index]*alphai[0][m] + 2.0*alphai[0][0] *
		                          (-2.0*bssn_ginv[0][1][m][index]*bssn_K[m][index]/3.0 +  
		                          6.0*bssn_Ainv[0][1][m][index]*bssn_phi[0][m][index] -
		                          bssn_ginv[0][1][m][index]*Si[m]);
		    bssn_gamma_y_RHS1 += -2.0*bssn_Ainv[0][2][m][index]*alphai[0][m] + 2.0*alphai[0][0] *
		                          (-2.0*bssn_ginv[0][2][m][index]*bssn_K[m][index]/3.0 +  
		                          6.0*bssn_Ainv[0][2][m][index]*bssn_phi[0][m][index] -
		                          bssn_ginv[0][2][m][index]*Si[m]);
		    bssn_gamma_z_RHS1 += -2.0*bssn_Ainv[0][3][m][index]*alphai[0][m] + 2.0*alphai[0][0] *
		                          (-2.0*bssn_ginv[0][3][m][index]*bssn_K[m][index]/3.0 +  
		                          6.0*bssn_Ainv[0][3][m][index]*bssn_phi[0][m][index] -
		                          bssn_ginv[0][3][m][index]*Si[m]);
		   
            for(n=1; n < 4; n++) {
		         bssn_gamma_x_RHS2 += 2.0*alphai[0][0]*(bssn_Chris1[1][m][n] * bssn_Ainv[0][n][m][index]);
		         bssn_gamma_y_RHS2 += 2.0*alphai[0][0]*(bssn_Chris1[2][m][n] * bssn_Ainv[0][n][m][index]);
		         bssn_gamma_z_RHS2 += 2.0*alphai[0][0]*(bssn_Chris1[3][m][n] * bssn_Ainv[0][n][m][index]);
		    }*/
        }
        
        Lie_gamma_x = 0.0;
        Lie_gamma_y = 0.0;
        Lie_gamma_z = 0.0;  
        
        for(m=1; m < 4; m++) {
            Lie_gamma_x += bssn_beta[0][0][m]*gamma[m][1][index] - gamma[0][m][index]*bssn_beta[0][m][1] +
                            2.0*gamma[0][1][index]*bssn_beta[0][m][m]/3.0;
            Lie_gamma_y += bssn_beta[0][0][m]*gamma[m][2][index] - gamma[0][m][index]*bssn_beta[0][m][2] +
                            2.0*gamma[0][2][index]*bssn_beta[0][m][m]/3.0;
            Lie_gamma_z += bssn_beta[0][0][m]*gamma[m][3][index] - gamma[0][m][index]*bssn_beta[0][m][3] +
                            2.0*gamma[0][3][index]*bssn_beta[0][m][m]/3.0;
            
            for(n=1; n < 4; n++) {
                Lie_gamma_x += bssn_ginv[0][m][n][index]*bssn_beta[m][n][1] + bssn_ginv[0][1][m][index]*bssn_beta[m][n][n]/3.0;
                Lie_gamma_y += bssn_ginv[0][m][n][index]*bssn_beta[m][n][2] + bssn_ginv[0][2][m][index]*bssn_beta[m][n][n]/3.0;
                Lie_gamma_z += bssn_ginv[0][m][n][index]*bssn_beta[m][n][3] + bssn_ginv[0][3][m][index]*bssn_beta[m][n][n]/3.0;
                            
            }
        }
		bssn_gamma_x_dt[index] = bssn_gamma_x_RHS1 + bssn_gamma_x_RHS2 + Lie_gamma_x;
		bssn_gamma_y_dt[index] = bssn_gamma_y_RHS1 + bssn_gamma_y_RHS2 + Lie_gamma_y;
		bssn_gamma_z_dt[index] = bssn_gamma_z_RHS1 + bssn_gamma_z_RHS2 + Lie_gamma_z;
		
		if (! fix_shift) {
		                    
/* Calc RHS of shift */

		betax_dt[index] = 0.0;
		betay_dt[index] = 0.0;
		betaz_dt[index] = 0.0;
		
		for(m=1; m < 4; m++) {
			betax_dt[index] += bssn_beta[0][0][m]*bssn_beta[0][m][1];
			betay_dt[index] += bssn_beta[0][0][m]*bssn_beta[0][m][2];
			betaz_dt[index] += bssn_beta[0][0][m]*bssn_beta[0][m][3];
		}
		
		betax_dt[index] += eta*bssn_bb_x[index];
		betay_dt[index] += eta*bssn_bb_y[index];
		betaz_dt[index] += eta*bssn_bb_z[index];
		
/* Calc RHS of B */

		bssn_bb_x_dt[index] = 0.0;
		bssn_bb_y_dt[index] = 0.0;
		bssn_bb_z_dt[index] = 0.0;
		
		for(m=1; m < 4; m++) {
			bssn_bb_x_dt[index] += bssn_beta[0][0][m]*bssn_bb[m][1][index];
			bssn_bb_y_dt[index] += bssn_beta[0][0][m]*bssn_bb[m][2][index];
			bssn_bb_z_dt[index] += bssn_beta[0][0][m]*bssn_bb[m][3][index];
		}
		
		bssn_bb_x_dt[index] += H_parm*bssn_gamma_x_dt[index] - bssn_bb_x[index];
		bssn_bb_y_dt[index] += H_parm*bssn_gamma_x_dt[index] - bssn_bb_y[index];
		bssn_bb_z_dt[index] += H_parm*bssn_gamma_x_dt[index] - bssn_bb_z[index];
		
		}
		
		if (TrK[index] >= TrK_max) {
		
			phi_dt[index] = 0.0;
			bssn_gxx_dt[index] = 0.0;
			bssn_gyy_dt[index] = 0.0;
			bssn_gzz_dt[index] = 0.0;
			bssn_gxy_dt[index] = 0.0;
			bssn_gxz_dt[index] = 0.0;
			bssn_gyz_dt[index] = 0.0;
			
			bssn_K_dt[index] = 0.0;
			bssn_Axx_dt[index] = 0.0;
			bssn_Ayy_dt[index] = 0.0;
			bssn_Azz_dt[index] = 0.0;
			bssn_Axy_dt[index] = 0.0;
			bssn_Axz_dt[index] = 0.0;
			bssn_Ayz_dt[index] = 0.0;
			
			bssn_gamma_x_dt[index] = 0.0;
			bssn_gamma_y_dt[index] = 0.0;
			bssn_gamma_z_dt[index] = 0.0;
			
		}
		
		if (fix_shift) {
			bssn_bb_x_dt[index] = 0.0; 
			bssn_bb_y_dt[index] = 0.0; 
			bssn_bb_z_dt[index] = 0.0; 
			betax_dt[index] = 0.0;
			betay_dt[index] = 0.0;
			betaz_dt[index] = 0.0; 
		}
                    
        }
     }
  }
	
   for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           for(q=0; q < 4; q++) {
               free(bssn_g[m][n][p][q]);
           }
        }
      }
    }
    
   for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           free(bssn_g[m][n][p]);
           free(bssn_ginv[m][n][p]);
           free(bssn_A[m][n][p]);
           free(bssn_Ainv[m][n][p]);
           free(shift[m][n][p]);
        }
      }
   }
    
   for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        free(bssn_g[m][n]);
        free(bssn_ginv[m][n]);
        free(bssn_A[m][n]);
        free(bssn_Ainv[m][n]);
        free(shift[m][n]);
        free(bssn_bb[m][n]);
        free(gamma[m][n]);
        free(lapse[m][n]);
        free(bssn_phi[m][n]);
     }
   }
    
  for(m=0; m < 4; m++) {
     free(bssn_g[m]);
     free(bssn_phi[m]);
     free(bssn_ginv[m]);
     free(bssn_bb[m]);
     free(bssn_A[m]);
     free(bssn_Ainv[m]);
     free(shift[m]);
     free(gamma[m]);
     free(lapse[m]);
     free(bssn_K[m]);
  }

  free(bssn_g);
  free(bssn_phi);
  free(bssn_ginv);
  free(bssn_bb);
  free(bssn_A);
  free(bssn_Ainv);
  free(shift);
  free(gamma);
  free(lapse);
  free(bssn_K);
}

/* Calculate Boundary Conditions */

void bssn_Boundaries(CCTK_ARGUMENTS)
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
    (cctkGH, CCTK_ALL_FACES, 1, -1, "SpecGRMHD::bssn_variables", bound);
  
  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }
  
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "SpecGRMHD::bssn_variables_Misc", bound);
  
  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }
  
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "SpecGRMHD::adm_bssn_cons", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }

  return;
}

