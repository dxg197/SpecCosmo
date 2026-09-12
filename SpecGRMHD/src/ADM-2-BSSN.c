/*@@ ADM to bssn conversion functions for thorn SpecGRMHD @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

#define SIE 1.11265006e-17
#define Speed_Light 2.99792458e8
 
/* Calc bssn Variables */

void calcbssn(CCTK_ARGUMENTS);
void SpecDeriv_Scalar_Derivative( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Scalar_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative2( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Scalar_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative2_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Scalar_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_FD( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );
void SpecDeriv_Scalar_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void SpecDeriv_Scalar_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Vector_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Vector_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor );
void SpecDeriv_Tensor_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );

  
void calcbssn(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
/* Declare local variables */

  CCTK_INT a,b,c,d,i,j,k,l,m,n,p,q,four,index;
  CCTK_INT sizex, sizey, sizez, ierr, handle;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_REAL dtg,shift2,Ricci_Scalar,con_fac,pi,ddphi;
  CCTK_REAL Sijxx,Sijyy,Sijzz,Sijxy,Sijxz,Sijyz;
  CCTK_REAL bj0,bj1,by0,by1,lam,zmin,zmax,lz,wc,t,scale_factor;
  CCTK_REAL T_total[4][4],igg3[4][4],Si[4],Ricci[4][4],gg4[4][4];
  CCTK_REAL bssn_Chris1[4][4][4],bssn_Chris2[4][4][4],kk3[4][4];
  CCTK_REAL bssn_gamma[4][4],Ttotal[4][4],bssn_A[4][4],AijAij,Ricci_phi[4][4];
   
/* Declare Arrays */

  CCTK_REAL *****bssn_g,*****bssn_ginv,**lapse,****bssn_beta;
  CCTK_REAL ****bssn_Ainv,**bssn_K,***bssn_phi;
  
  
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
    
    lapse = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    bssn_K = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    bssn_phi = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
    bssn_ginv = (CCTK_REAL *****)malloc(sizeof(CCTK_REAL ****)*four);
    bssn_Ainv = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    bssn_beta = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
    bssn_g = (CCTK_REAL *****)malloc(sizeof(CCTK_REAL ****)*four);
                   
    for(m=0; m < 4; m++) {
       lapse[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       bssn_K[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       bssn_phi[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       bssn_ginv[m] = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
       bssn_Ainv[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       bssn_beta[m] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
       bssn_g[m] = (CCTK_REAL ****)malloc(sizeof(CCTK_REAL ***)*four);
        
        for(n=0; n < 4; n++) {
          bssn_phi[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
          bssn_ginv[m][n] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
          bssn_Ainv[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_beta[m][n] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
          bssn_g[m][n] = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
        
            for(p=0; p < 4; p++) {
              bssn_Ainv[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_beta[m][n][p] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
              bssn_ginv[m][n][p] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
              bssn_g[m][n][p] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
              
              for(q=0; q < 4; q++) {
                   bssn_ginv[m][n][p][q] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
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

  	zmax +=  CCTK_DELTA_SPACE(2)*(cctk_nghostzones[2]-0.5);
  	zmin += -CCTK_DELTA_SPACE(2)*(cctk_nghostzones[2]-0.5);
  	lz = zmax - zmin;
  	
  	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &scale_factor, 1, CCTK_VarIndex("MHD_Analysis::aa_out_avg"));

/* Do the */

for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
        
/* Calc det of ADM metric */
    
    dtg =      gxx[index] * gyy[index] * gzz[index] 
           + 2.0*gxy[index] * gxz[index] * gyz[index] 
           -   gxx[index] * pow(gyz[index],2.0)           
           -   gyy[index] * pow(gxz[index],2.0)           
           -   gzz[index] * pow(gxy[index],2.0);
         
/* Calc inverse of ADM metric */

    igxx[index] = (gyy[index] * gzz[index] - pow(gyz[index],2.0) ) / dtg;
    igyy[index] = (gxx[index] * gzz[index] - pow(gxz[index],2.0) ) / dtg;
    igzz[index] = (gxx[index] * gyy[index] - pow(gxy[index],2.0) ) / dtg;
    igxy[index] = (gxz[index] * gyz[index] - gxy[index] * gzz[index]) / dtg;
    igxz[index] = (gxy[index] * gyz[index] - gxz[index] * gyy[index]) / dtg;
    igyz[index] = (gxz[index] * gxy[index] - gyz[index] * gxx[index]) / dtg;
		
    igg3[1][1] = igxx[index];
    igg3[2][2] = igyy[index];
    igg3[3][3] = igzz[index];
    igg3[1][2] = igxy[index];
    igg3[1][3] = igxz[index];
    igg3[2][3] = igyz[index];
    igg3[2][1] = igxy[index];
    igg3[3][1] = igxz[index];
    igg3[3][2] = igyz[index];
    
/* Calc phi */
    
    phi[index] = log(dtg)/12.0;
    
    bssn_phi[0][0][index] = phi[index];
    
    con_fac = exp(-4.0*phi[index]);
    
/* Calc bssn metric */
  
    bssn_g[0][0][1][1][index] = con_fac * gxx[index];
    bssn_g[0][0][2][2][index] = con_fac * gyy[index];
    bssn_g[0][0][3][3][index] = con_fac * gzz[index];
    bssn_g[0][0][1][2][index] = con_fac * gxy[index];
    bssn_g[0][0][1][3][index] = con_fac * gxz[index];
    bssn_g[0][0][2][3][index] = con_fac * gyz[index];
    bssn_g[0][0][2][1][index] = con_fac * gxy[index];
    bssn_g[0][0][3][1][index] = con_fac * gxz[index];
    bssn_g[0][0][3][2][index] = con_fac * gyz[index];
    
    bssn_gxx[index] = bssn_g[0][0][1][1][index];
    bssn_gyy[index] = bssn_g[0][0][2][2][index];
    bssn_gzz[index] = bssn_g[0][0][3][3][index];
    bssn_gxy[index] = bssn_g[0][0][1][2][index];
    bssn_gxz[index] = bssn_g[0][0][1][3][index];
    bssn_gyz[index] = bssn_g[0][0][2][3][index];
    
    dtg =      bssn_gxx[index] * bssn_gyy[index] * bssn_gzz[index] 
           + 2.0*bssn_gxy[index] * bssn_gxz[index] * bssn_gyz[index] 
           -   bssn_gxx[index] * pow(bssn_gyz[index],2.0)           
           -   bssn_gyy[index] * pow(bssn_gxz[index],2.0)           
           -   bssn_gzz[index] * pow(bssn_gxy[index],2.0);
    
    detg1[index] = dtg - 1.0;
    
    bssn_ginv[0][0][1][1][index] = (bssn_gyy[index] * bssn_gzz[index] - bssn_gyz[index] * bssn_gyz[index]) / dtg;
    bssn_ginv[0][0][2][2][index] = (bssn_gxx[index] * bssn_gzz[index] - bssn_gxz[index] * bssn_gxz[index]) / dtg;
    bssn_ginv[0][0][3][3][index] = (bssn_gxx[index] * bssn_gyy[index] - bssn_gxy[index] * bssn_gxy[index]) / dtg;                         
    bssn_ginv[0][0][1][2][index] = (bssn_gxz[index] * bssn_gyz[index] - bssn_gxy[index] * bssn_gzz[index]) / dtg;
    bssn_ginv[0][0][1][3][index] = (bssn_gxy[index] * bssn_gyz[index] - bssn_gxz[index] * bssn_gyy[index]) / dtg;
    bssn_ginv[0][0][2][3][index] = (bssn_gxz[index] * bssn_gxy[index] - bssn_gyz[index] * bssn_gxx[index]) / dtg; 
    bssn_ginv[0][0][2][1][index] = (bssn_gxz[index] * bssn_gyz[index] - bssn_gxy[index] * bssn_gzz[index]) / dtg;
    bssn_ginv[0][0][3][1][index] = (bssn_gxy[index] * bssn_gyz[index] - bssn_gxz[index] * bssn_gyy[index]) / dtg;
    bssn_ginv[0][0][3][2][index] = (bssn_gxz[index] * bssn_gxy[index] - bssn_gyz[index] * bssn_gxx[index]) / dtg;   
    
/* Calc Trace K */

    kk3[1][1] = kxx[index];
    kk3[2][2] = kyy[index];
    kk3[3][3] = kzz[index];
    kk3[1][2] = kxy[index];
    kk3[1][3] = kxz[index];
    kk3[2][3] = kyz[index];
    kk3[2][1] = kxy[index];
    kk3[3][1] = kxz[index];
    kk3[3][2] = kyz[index];
     
	bssn_K[0][index] = 0.0;
	
	for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
	        bssn_K[0][index] += igg3[m][n]*kk3[m][n];
	    }
	}
	
	TrK[index] = bssn_K[0][index];
	                  
/* Calc Trace-Free Extrinsic Curvature */
  
    bssn_Axx[index] = con_fac * (kxx[index] - gxx[index]*TrK[index]/3.0);
    bssn_Ayy[index] = con_fac * (kyy[index] - gyy[index]*TrK[index]/3.0);
    bssn_Azz[index] = con_fac * (kzz[index] - gzz[index]*TrK[index]/3.0);
    bssn_Axy[index] = con_fac * (kxy[index] - gxy[index]*TrK[index]/3.0);
    bssn_Axz[index] = con_fac * (kxz[index] - gxz[index]*TrK[index]/3.0);
    bssn_Ayz[index] = con_fac * (kyz[index] - gyz[index]*TrK[index]/3.0);
    
    bssn_A[1][1] = bssn_Axx[index];
    bssn_A[2][2] = bssn_Ayy[index];
    bssn_A[3][3] = bssn_Azz[index];
    bssn_A[1][2] = bssn_Axy[index];
    bssn_A[1][3] = bssn_Axz[index];
    bssn_A[2][3] = bssn_Ayz[index];
    bssn_A[2][1] = bssn_Axy[index];
    bssn_A[3][1] = bssn_Axz[index];
    bssn_A[3][2] = bssn_Ayz[index];
       
       for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
		       bssn_Ainv[0][m][n][index] = 0.0;
		    }
		}
    
       for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
               for(p=1; p < 4; p++) {
                  for(q=1; q < 4; q++) {
		             bssn_Ainv[0][m][n][index] += bssn_ginv[0][0][m][p][index]*bssn_ginv[0][0][n][q][index]*
		                                          bssn_A[p][q];
		          }
		       }
            }
        }
    
/* Calc lapse and shift */
    
    if (CCTK_Equals(gauge_condition,"geodesic")) {
    	alpha[index] = 1.0;
    }
    
    if (CCTK_Equals(gauge_condition,"scale")) {
    	alpha[index] = aa0;
    }
    	
    if (CCTK_Equals(gauge_condition,"geod")) {
        alpha[index] = exp(12.0*phi[index]*sigma);
    }
    	
    if (CCTK_Equals(gauge_condition,"gwave")) {
    
    	wc = nw*2.0*pi/lz;
    
    	alpha[index] = exp(0.5*Amp*sin(wc*(cctk_time + z[index])));
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
    }
    
    sqrtdetg[index] = alpha[index]*exp(6.0*phi[index]);
    
    lapse[0][index] = alpha[index];
    
    bssn_beta[0][0][1][index] = betax[index];
    bssn_beta[0][0][2][index] = betay[index];
    bssn_beta[0][0][3][index] = betaz[index];
    
         }
     }
 }
 
 printf("Calculate Derivatives \n");
 
/* Calc derivatives */

  if (CCTK_Equals(diff,"Spectral")) {
  
  printf("Use Sepectral Methods (3 Dimensions) \n");
  
/* bssn phi derivatives */

  SpecDeriv_Scalar_Derivative2(CCTK_PASS_CTOC, bssn_phi); 
  
/* bssn metric derivatives */

  SpecDeriv_Tensor_Derivative2(CCTK_PASS_CTOC, bssn_g); 
  SpecDeriv_Tensor_Derivative2(CCTK_PASS_CTOC, bssn_ginv); 
  
/* gauge derivatives */
 
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, lapse);
  SpecDeriv_Vector_Derivative2(CCTK_PASS_CTOC, bssn_beta);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Spectral_1d")) {
  
  printf("Use Sepectral Methods (One Dimension - z only) \n");
  
/* bssn phi derivatives */

  SpecDeriv_Scalar_Derivative2_1d(CCTK_PASS_CTOC, bssn_phi); 
  
/* bssn metric derivatives */

  SpecDeriv_Tensor_Derivative2_1d(CCTK_PASS_CTOC, bssn_g);
  SpecDeriv_Tensor_Derivative2_1d(CCTK_PASS_CTOC, bssn_ginv);
  
/* gauge derivatives */
 
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, lapse);
  SpecDeriv_Vector_Derivative2_1d(CCTK_PASS_CTOC, bssn_beta);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_1d(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_1d(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Finite")) {
  
  printf("Use Finite Differencing (Second Order) \n");
  
/* bssn phi derivatives */

  SpecDeriv_Scalar_Derivative2_FD(CCTK_PASS_CTOC, bssn_phi); 
  
/* bssn metric derivatives */

  SpecDeriv_Tensor_Derivative2_FD(CCTK_PASS_CTOC, bssn_g); 
  SpecDeriv_Tensor_Derivative2_FD(CCTK_PASS_CTOC, bssn_ginv); 
  
/* gauge derivatives */
 
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, lapse);
  SpecDeriv_Vector_Derivative2_FD(CCTK_PASS_CTOC, bssn_beta);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_FD(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_FD(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"Finite4")) {
  
  printf("Use Finite Differencing (Fourth Order) \n");
  
/* bssn phi derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, bssn_phi); 
  
/* bssn metric derivatives */

  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_g); \
  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_ginv); 
  
/* gauge derivatives */
 
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, lapse);
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, bssn_beta);
  
/* extrinsic curvature derivative */  

  SpecDeriv_Tensor_Derivative_FD4(CCTK_PASS_CTOC, bssn_Ainv);
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, bssn_K); 
  
  }
  
  if (CCTK_Equals(diff,"SCR3")) {
  
  printf("Use SCR3 Finite Differencing (Fourth Order)\n");
  
/* bssn phi derivatives */

  SpecDeriv_Scalar_Derivative2_FD4(CCTK_PASS_CTOC, bssn_phi); 
  
/* bssn metric derivatives */

  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_g); \
  SpecDeriv_Tensor_Derivative2_FD4(CCTK_PASS_CTOC, bssn_ginv); 
  
/* gauge derivatives */
 
  SpecDeriv_Scalar_Derivative_FD4(CCTK_PASS_CTOC, lapse);
  SpecDeriv_Vector_Derivative2_FD4(CCTK_PASS_CTOC, bssn_beta);
  
/* extrinsic curvature derivative */  

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
 
/* Calc shift squared */
  
  	shift2 = gxx[index]*betax[index]*betax[index] + 2.0*gxy[index]*betax[index]*betay[index] + 
  	         gyy[index]*betay[index]*betay[index] + 2.0*gxz[index]*betax[index]*betaz[index] + 
  	         gzz[index]*betaz[index]*betaz[index] + 2.0*gyz[index]*betay[index]*betaz[index];
    
/* Calc metric */
    
    gg4[0][0] = -pow(alpha[index],2.0) + shift2;
    gg4[0][1] = gxx[index]*betax[index]+gxy[index]*betay[index]+gxz[index]*betaz[index];
    gg4[0][2] = gxy[index]*betax[index]+gyy[index]*betay[index]+gyz[index]*betaz[index];
    gg4[0][3] = gxz[index]*betax[index]+gyz[index]*betay[index]+gzz[index]*betaz[index];
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
        
/* Calc Connections */

     for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           for(p=1; p < 4; p++) {
              bssn_Chris1[m][n][p] = 0.0;
              bssn_Chris2[m][n][p] = 0.0;
           }
        }
     }
     
     
     for(a=1; a < 4; a++) {
        for(b=1; b < 4; b++) {
           for(c=1; c < 4; c++) {
              for(d=1; d < 4; d++) {
                 bssn_Chris1[a][b][c] += 0.5*bssn_ginv[0][0][a][d][index]*(bssn_g[0][b][c][d][index] + 
                                         bssn_g[0][c][b][d][index] - bssn_g[0][d][b][c][index]);
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
        
/* Calc Gamma */
   
   for(m=0; m < 4; m++) {
        for(n=0; n < 4; n++) {
   			bssn_gamma[m][n] = 0.0;
   			Ricci[m][n] = 0.0;
   		}
   }
   
   for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           bssn_gamma[0][1] += bssn_ginv[0][0][m][n][index]*bssn_Chris1[1][m][n];
           bssn_gamma[0][2] += bssn_ginv[0][0][m][n][index]*bssn_Chris1[2][m][n];
           bssn_gamma[0][3] += bssn_ginv[0][0][m][n][index]*bssn_Chris1[3][m][n];
        }
   }
   
   for(m=1; m < 4; m++) {
     for(n=1; n < 4; n++) {
        for(p=1; p < 4; p++) {
           bssn_gamma[m][n] += -bssn_ginv[m][p][n][p][index];
        }
     }
   }
   
   bssn_gamma_x[index] = bssn_gamma[0][1];
   bssn_gamma_y[index] = bssn_gamma[0][2];
   bssn_gamma_z[index] = bssn_gamma[0][3]; 
   
   bssn_bb_x[index] = 0.0;
   bssn_bb_y[index] = 0.0;
   bssn_bb_z[index] = 0.0;
        
/* Calc Ricci Tensor */

    for(a=1; a < 4; a++) {
        for(b=1; b < 4; b++) {
           for(c=1; c < 4; c++) {
              for(l=1; l < 4; l++) {
                 for(m=1; m < 4; m++) {
     
     			   Ricci[a][b] += bssn_ginv[0][0][c][l][index]*(bssn_Chris1[m][c][a]*bssn_Chris2[b][m][l] +
     			                  bssn_Chris1[m][c][b]*bssn_Chris2[a][m][l] +
     			                  bssn_Chris1[m][c][a]*bssn_Chris2[b][m][l]);
                   }
                   
                   Ricci[a][b] += -0.5*bssn_ginv[0][0][c][l][index]*bssn_g[c][l][a][b][index];
               }
               
               Ricci[a][b] += 0.5*bssn_gamma[0][c]*bssn_Chris2[a][b][c] +
     			              0.5*bssn_gamma[0][c]*bssn_Chris2[b][a][c] +
     			              0.5*bssn_g[0][0][c][a][index]*bssn_gamma[b][c] + 
     				          0.5*bssn_g[0][0][c][b][index]*bssn_gamma[a][c];
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
   					ddphi += -bssn_ginv[0][0][m][n][index]*bssn_Chris1[l][m][n]*bssn_phi[0][l][index];
   				}
   			ddphi += bssn_ginv[0][0][m][n][index]*bssn_phi[m][n][index];
   	 		}
   		}
        
        for(m=1; m < 4; m++) {
            for(n=1; n < 4; n++) {
            	for(l=1; l < 4; l++) {
            		Ricci_phi[m][n] += 2.0*bssn_Chris1[l][m][n]*bssn_phi[0][l][index];
            			for(a=1; a < 4; a++) {
            				Ricci_phi[m][n] += -4.0*bssn_g[0][0][m][n][index]*bssn_ginv[0][0][l][a][index]*bssn_phi[0][a][index]*bssn_phi[0][l][index];
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
             Ricci_Scalar += con_fac*bssn_ginv[0][0][m][n][index]*(Ricci[m][n]+Ricci_phi[m][n]);
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
  
   for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) { 
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
   
   rho_total[index] = 8.0*pi*pow(lapse[0][index],2.0)*Ttotal[0][0];
     
   for(m=0; m < 4; m++) {
   	Si[m] = 0.0;
   }
     
   for(n=1; n < 4; n++) {
   	for(m=1; m < 4; m++) {
   		Si[n] += 8.0*pi*alpha[index]*(gg4[n][m]*Ttotal[0][m]);
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
            for(p=1; p < 4; p++) {
                for(q=1; q < 4; q++) {
		            AijAij += bssn_g[0][0][m][n][index]*bssn_g[0][0][p][q][index]*
		                      bssn_Ainv[0][m][p][index]*bssn_Ainv[0][n][q][index];
		        }
		    }
    	}
    }
    
   ddphi = 0.0;
   
   for(m=1; m < 4; m++) {
     for(n=1; n < 4; n++) {
   		for(l=1; l < 4; l++) {
   			ddphi += -bssn_ginv[0][0][m][n][index]*bssn_Chris1[l][m][n]*bssn_phi[0][l][index];
   		}
   		ddphi += bssn_ginv[0][0][m][n][index]*bssn_phi[m][n][index];
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
	                   bssn_ginv[0][0][1][m][index]*(2.0*bssn_K[m][index]/3.0 - exp(-4.0*phi[index])*Si[m]);
	    momy[index] += bssn_Ainv[m][2][m][index] + 6.0*bssn_Ainv[0][2][m][index]*bssn_phi[0][m][index] - 
	                   bssn_ginv[0][0][2][m][index]*(2.0*bssn_K[m][index]/3.0 - exp(-4.0*phi[index])*Si[m]);
	    momz[index] += bssn_Ainv[m][3][m][index] + 6.0*bssn_Ainv[0][3][m][index]*bssn_phi[0][m][index] - 
	                   bssn_ginv[0][0][3][m][index]*(2.0*bssn_K[m][index]/3.0 - exp(-4.0*phi[index])*Si[m]);
	}
	
	for(m=1; m < 4; m++) {
		for(n=1; n < 4; n++) {
	        momx[index] += bssn_Chris1[1][m][n]*bssn_Ainv[0][m][n][index];
	        momy[index] += bssn_Chris1[2][m][n]*bssn_Ainv[0][m][n][index];
	        momz[index] += bssn_Chris1[3][m][n]*bssn_Ainv[0][m][n][index];
	    }
	}
    
    dGx[index] = 0.0; 
    dGy[index] = 0.0;
    dGz[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) {
           dGx[index] += -bssn_ginv[0][0][m][n][index]*bssn_Chris1[1][m][n];
           dGy[index] += -bssn_ginv[0][0][m][n][index]*bssn_Chris1[2][m][n];
           dGz[index] += -bssn_ginv[0][0][m][n][index]*bssn_Chris1[3][m][n];
        }
   }
    
    dGx[index] += bssn_gamma_x[index];
    dGy[index] += bssn_gamma_y[index];
    dGz[index] += bssn_gamma_z[index];
    
    dHx[index] = 0.0; 
    dHy[index] = 0.0;
    dHz[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        dHx[index] += bssn_ginv[0][m][1][m][index];
        dHy[index] += bssn_ginv[0][m][2][m][index];
        dHz[index] += bssn_ginv[0][m][3][m][index];
    }
    
    dHx[index] += bssn_gamma_x[index];
    dHy[index] += bssn_gamma_y[index];
    dHz[index] += bssn_gamma_z[index];
    
    AATF[index] = 0.0;
    
    for(m=1; m < 4; m++) {
        for(n=1; n < 4; n++) { 
          AATF[index] += bssn_g[0][0][m][n][index]*bssn_Ainv[0][m][n][index];
        }
    }
    
        }
	}
}
  
  for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           for(q=0; q < 4; q++) {
               free(bssn_g[m][n][p][q]);
               free(bssn_ginv[m][n][p][q]);
           }
        }
     }
  }
  
  for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        for(p=0; p < 4; p++) {
           free(bssn_g[m][n][p]);
           free(bssn_ginv[m][n][p]);
           free(bssn_Ainv[m][n][p]);
           free(bssn_beta[m][n][p]);
        }
     }
  }
  
  for(m=0; m < 4; m++) {
     for(n=0; n < 4; n++) {
        free(bssn_g[m][n]);
        free(bssn_ginv[m][n]);
        free(bssn_Ainv[m][n]);
        free(bssn_beta[m][n]);
        free(bssn_phi[m][n]);
     }
  }   
   
  for(m=0; m < 4; m++) {  
    free(bssn_g[m]);
    free(bssn_ginv[m]);
    free(lapse[m]);
    free(bssn_Ainv[m]);
    free(bssn_beta[m]);
    free(bssn_phi[m]);
    free(bssn_K[m]);
  }

  free(bssn_g);
  free(bssn_ginv);
  free(lapse);
  free(bssn_Ainv);
  free(bssn_K);
  free(bssn_phi);
  free(bssn_beta);
}

