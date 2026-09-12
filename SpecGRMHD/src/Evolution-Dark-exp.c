/*@@ Calculates the Tmunu for Dark Matter/Energy using the Theory by Ma and Wang physics.gen-ph/1206.5078 @@*/
/* Assumes Geodesic Slicing */

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

void Dark_RHS(CCTK_ARGUMENTS);
void Dark_Tmunu(CCTK_ARGUMENTS);
void Dark_Boundaries(CCTK_ARGUMENTS);

/* Begin Dark_Tmunu */

void Dark_Tmunu(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
/* Declare local variables */
  CCTK_INT i,j,k,m,n,q,a,b,c,d,four,index;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez;
  CCTK_REAL detg,gg4[4][4],igg4[4][4],conf,dtg,pi;
  CCTK_REAL beta2,Tdark[4][4],T_dark[4][4],T_mhd[4][4];
  
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
        
   beta2 = conf*(bssn_gxx[index]*betax[index]*betax[index] + 2.0*bssn_gxy[index]*betax[index]*betay[index] + 
  	             bssn_gyy[index]*betay[index]*betay[index] + 2.0*bssn_gxz[index]*betax[index]*betaz[index] + 
  	             bssn_gzz[index]*betaz[index]*betaz[index] + 2.0*bssn_gyz[index]*betay[index]*betaz[index]);
    
/* Calc metric */
    
    gg4[0][0] = -pow(alpha[index],2) + beta2;
    gg4[0][1] = betax[index];
    gg4[0][2] = betay[index];
    gg4[0][3] = betaz[index];
    gg4[1][0] = betax[index];
    gg4[2][0] = betay[index];
    gg4[3][0] = betaz[index];
    gg4[1][1] = conf*bssn_gxx[index];
    gg4[2][2] = conf*bssn_gyy[index];
    gg4[3][3] = conf*bssn_gzz[index];
    gg4[1][2] = conf*bssn_gxy[index];
    gg4[1][3] = conf*bssn_gxz[index];
    gg4[2][3] = conf*bssn_gyz[index];
    gg4[2][1] = conf*bssn_gxy[index];
    gg4[3][1] = conf*bssn_gxz[index];
    gg4[3][2] = conf*bssn_gyz[index];
    
/* Calc det of ADM metric */
    
    dtg =      gxx[index] * gyy[index] * gzz[index] 
           + 2*gxy[index] * gxz[index] * gyz[index] 
           -   gxx[index] * pow(gyz[index],2)           
           -   gyy[index] * pow(gxz[index],2)           
           -   gzz[index] * pow(gxy[index],2);
         
/* Calc inverse of ADM metric */

        igg4[1][1] = -(gg4[0][0] * gyy[index] * gzz[index] - gg4[0][0] * pow(gyz[index],2) + 
		               gg4[0][2] * gg4[0][3] * gyz[index] + gg4[0][2] * gg4[0][3] * gyz[index] -
		               gg4[0][2] * gg4[0][2] * gzz[index] - gg4[0][3] * gg4[0][3] * gyy[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[2][2] = -(gg4[0][0] * gxx[index] * gzz[index] - gg4[0][0] * pow(gxz[index],2) + 
		               gg4[0][1] * gg4[0][3] * gxz[index] + gg4[0][1] * gg4[0][3] * gxz[index] -
		               gg4[0][1] * gg4[0][1] * gzz[index] - gg4[0][3] * gg4[0][3] * gxx[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[3][3] = -(gg4[0][0] * gxx[index] * gyy[index] - gg4[0][0] * pow(gxy[index],2) + 
		               gg4[0][1] * gg4[0][2] * gxy[index] + gg4[0][1] * gg4[0][2] * gxy[index] -
		               gg4[0][1] * gg4[0][1] * gyy[index] - gg4[0][2] * gg4[0][2] * gxx[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[1][2] = -(gg4[0][0] * gxz[index] * gyz[index] - gg4[0][0] * gxy[index] * gzz[index] + 
		               gg4[0][1] * gg4[0][2] * gzz[index] + gg4[0][3] * gg4[0][3] * gxy[index] -
		               gg4[0][2] * gg4[0][3] * gxz[index] - gg4[0][3] * gg4[0][1] * gyz[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[1][3] = -(gg4[0][0] * gxy[index] * gyz[index] - gg4[0][0] * gxz[index] * gyy[index] + 
		               gg4[0][2] * gg4[0][2] * gxz[index] + gg4[0][3] * gg4[0][1] * gyy[index] -
		               gg4[0][2] * gg4[0][1] * gyz[index] - gg4[0][3] * gg4[0][2] * gxy[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[2][3] = -(gg4[0][0] * gxz[index] * gxy[index] - gg4[0][0] * gyz[index] * gxx[index] + 
		               gg4[0][1] * gg4[0][1] * gyz[index] + gg4[0][3] * gg4[0][2] * gxx[index] -
		               gg4[0][1] * gg4[0][2] * gxz[index] - gg4[0][3] * gg4[0][1] * gxy[index]) / 
		               (alpha[index]*alpha[index]*dtg);
		               
        igg4[2][1] = igg4[1][2];
        
        igg4[3][1] = igg4[1][3];
        
        igg4[3][2] = igg4[2][3];
        
        igg4[0][0] = -(gxx[index]*gyy[index]*gzz[index] + gxy[index]*gyz[index]*gxz[index] + 
        			   gxz[index]*gxy[index]*gyz[index] - gxx[index]*gyz[index]*gyz[index] - 
        			   gxy[index]*gxy[index]*gzz[index] - gxz[index]*gyy[index]*gxz[index]) /
                      (alpha[index]*alpha[index]*dtg);
                      
        igg4[0][1] = -(gg4[0][1]*gyz[index]*gyz[index] + gg4[0][2]*gxy[index]*gzz[index] + 
        			   gg4[0][3]*gyy[index]*gxz[index] - gg4[0][1]*gyy[index]*gzz[index] - 
        			   gg4[0][2]*gyz[index]*gxz[index] - gg4[0][3]*gxy[index]*gyz[index]) /
                      (alpha[index]*alpha[index]*dtg);
                      
        igg4[0][2] = -(gg4[0][1]*gxy[index]*gzz[index] + gg4[0][2]*gxz[index]*gxz[index] + 
        			   gg4[0][3]*gxx[index]*gyz[index] - gg4[0][1]*gxz[index]*gyz[index] - 
        			   gg4[0][2]*gxx[index]*gzz[index] - gg4[0][3]*gxy[index]*gxz[index]) /
                      (alpha[index]*alpha[index]*dtg);
                      
        igg4[0][3] = -(gg4[0][1]*gxz[index]*gyy[index] + gg4[0][2]*gxx[index]*gyz[index] + 
        			   gg4[0][3]*gxy[index]*gxy[index] - gg4[0][1]*gxy[index]*gyz[index] - 
        			   gg4[0][2]*gxz[index]*gxy[index] - gg4[0][3]*gxx[index]*gyy[index]) /
                      (alpha[index]*alpha[index]*dtg);
                      
        igg4[1][0] = igg4[0][1];
        
        igg4[2][0] = igg4[0][2];
        
        igg4[3][0] = igg4[0][3];
        
/* Calculate scalar stress-tensor */

		T_mhd[0][0] = Ttt_mhd[index];
		T_mhd[0][1] = Ttx_mhd[index];
		T_mhd[0][2] = Tty_mhd[index];
		T_mhd[0][3] = Ttz_mhd[index];	
		T_mhd[1][1] = Txx_mhd[index];
		T_mhd[2][2] = Tyy_mhd[index];
		T_mhd[3][3] = Tzz_mhd[index];
		T_mhd[1][2] = Txy_mhd[index];
		T_mhd[1][3] = Txz_mhd[index];
		T_mhd[2][3] = Tyz_mhd[index];
				
		T_mhd[1][0] = Ttx_mhd[index];
		T_mhd[2][0] = Tty_mhd[index];
		T_mhd[3][0] = Ttz_mhd[index];
		
		T_mhd[2][1] = Txy_mhd[index];
		T_mhd[3][1] = Txz_mhd[index];
		T_mhd[3][2] = Tyz_mhd[index];
		
		TT[index] = 0.0;

		for(m=0; m < 4; m++) {
          for(n=0; n < 4; n++) { 
             TT[index] += T_mhd[m][n]*gg4[m][n];
             T_dark[m][n] = 0.0;
          }
        }
		
/* Calculate Spacial part of Tmunu for Dark Matter */

  		Ttt_dark[index] = -igg4[0][0]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Ttx_dark[index] = -igg4[0][1]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Tty_dark[index] = -igg4[0][2]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Ttz_dark[index] = -igg4[0][3]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Txx_dark[index] = -igg4[1][1]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Tyy_dark[index] = -igg4[2][2]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Tzz_dark[index] = -igg4[3][3]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Txy_dark[index] = -igg4[1][2]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Txz_dark[index] = -igg4[1][3]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		Tyz_dark[index] = -igg4[2][3]*(RR[index]/(8.0*pi) + TT[index])/4.0;
  		
  		Tdark[0][0] = Ttt_dark[index];
        Tdark[0][1] = Ttx_dark[index];
        Tdark[0][2] = Tty_dark[index];
        Tdark[0][3] = Ttz_dark[index];
        Tdark[1][1] = Txx_dark[index];
        Tdark[1][2] = Txy_dark[index];
        Tdark[1][3] = Txz_dark[index];
        Tdark[2][2] = Tyy_dark[index];
        Tdark[2][3] = Tyz_dark[index];
        Tdark[3][3] = Tzz_dark[index];
   
        Tdark[1][0] = Tdark[0][1];
        Tdark[2][0] = Tdark[0][2];
        Tdark[3][0] = Tdark[0][3];
   
        Tdark[2][1] = Tdark[1][2];
        Tdark[3][1] = Tdark[1][3];
        Tdark[3][2] = Tdark[2][3];
   
  		for(m=0; m < 4; m++) {
          for(n=0; n < 4; n++) { 
           T_dark[0][0] += Tdark[m][n]*gg4[0][m]*gg4[0][n];
           T_dark[0][1] += Tdark[m][n]*gg4[0][m]*gg4[1][n];
           T_dark[0][2] += Tdark[m][n]*gg4[0][m]*gg4[2][n];
           T_dark[0][3] += Tdark[m][n]*gg4[0][m]*gg4[3][n];
           T_dark[1][1] += Tdark[m][n]*gg4[1][m]*gg4[1][n];
           T_dark[2][2] += Tdark[m][n]*gg4[2][m]*gg4[2][n];
           T_dark[3][3] += Tdark[m][n]*gg4[3][m]*gg4[3][n];
           T_dark[1][2] += Tdark[m][n]*gg4[1][m]*gg4[2][n];
           T_dark[1][3] += Tdark[m][n]*gg4[1][m]*gg4[3][n];
           T_dark[2][3] += Tdark[m][n]*gg4[2][m]*gg4[3][n];
          }
        }
   
   T_dark[1][0] = T_dark[0][1];
   T_dark[2][0] = T_dark[0][2];
   T_dark[3][0] = T_dark[0][3];
   
   T_dark[2][1] = T_dark[1][2];
   T_dark[3][1] = T_dark[1][3];
   T_dark[3][2] = T_dark[2][3];
  		
  	rho_dark[index] = (T_dark[0][0] - 2.0*betax[index]*T_dark[0][1] - 2.0*betay[index]*T_dark[0][2] - 
                      2.0*betaz[index]*T_dark[0][3] + betax[index]*betax[index]*T_dark[1][1] + 
                      betay[index]*betay[index]*T_dark[2][2] + betaz[index]*betaz[index]*T_dark[3][3] + 
                      2.0*(betax[index]*betay[index]*T_dark[1][2] + betax[index]*betaz[index]*T_dark[1][3] +
                      betay[index]*betaz[index]*T_dark[2][3]))/pow(alpha[index],2);
  		      
		}
     }
  } 
  
}

