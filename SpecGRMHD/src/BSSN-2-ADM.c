/*@@ bssn to ADM conversion functions for thorn SpecGRMHD @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Calc ADM Variables */
  
void calcadm(CCTK_ARGUMENTS);
int sgn(CCTK_REAL v);
CCTK_REAL determinant(CCTK_REAL a[4][4], CCTK_REAL k);
void cofactor(CCTK_REAL num[4][4], CCTK_REAL inverse[4][4], CCTK_REAL f);
void transpose(CCTK_REAL num[4][4], CCTK_REAL fac[4][4], CCTK_REAL inverse[4][4], CCTK_REAL r);

void calcadm(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
/* Declare local variables */

  CCTK_INT i,j,k;
  CCTK_INT sizex, sizey, sizez, index;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_REAL dtg, conf, shift2, gg[4][4], igg[4][4],four;
  
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
    
    four = 4.0;

/* Do the */

for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
	conf = exp(4.0*phi[index]);
        
/* Calculate ADM metric */

	gxx[index] = conf*bssn_gxx[index];
    gyy[index] = conf*bssn_gyy[index];
    gzz[index] = conf*bssn_gzz[index];
    gxy[index] = conf*bssn_gxy[index];
    gxz[index] = conf*bssn_gxz[index];
    gyz[index] = conf*bssn_gyz[index];

  	shift2 = conf*(bssn_gxx[index]*betax[index]*betax[index] + 2.0*bssn_gxy[index]*betax[index]*betay[index] + 
  	       	bssn_gyy[index]*betay[index]*betay[index] + 2.0*bssn_gxz[index]*betax[index]*betaz[index] + 
  	       	bssn_gzz[index]*betaz[index]*betaz[index] + 2.0*bssn_gyz[index]*betay[index]*betaz[index]);
  	        
  	gtt[index] = -pow(alpha[index],2.0) + shift2;
  	gtx[index] = gxx[index]*betax[index]+gxy[index]*betay[index]+gxz[index]*betaz[index];
  	gty[index] = gxy[index]*betax[index]+gyy[index]*betay[index]+gyz[index]*betaz[index]; 
  	gtz[index] = gxz[index]*betax[index]+gyz[index]*betay[index]+gzz[index]*betaz[index];
  	
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
    
/* Calculate ADM Extrinsic Curvature */
      
    kxx[index] = conf*(bssn_Axx[index] + bssn_gxx[index]*TrK[index]/3.0);
    kyy[index] = conf*(bssn_Ayy[index] + bssn_gyy[index]*TrK[index]/3.0);
    kzz[index] = conf*(bssn_Azz[index] + bssn_gzz[index]*TrK[index]/3.0);
    kxy[index] = conf*(bssn_Axy[index] + bssn_gxy[index]*TrK[index]/3.0);
    kxz[index] = conf*(bssn_Axz[index] + bssn_gxz[index]*TrK[index]/3.0);
    kyz[index] = conf*(bssn_Ayz[index] + bssn_gyz[index]*TrK[index]/3.0);
    
/* Calc inverse of ADM metric */

    dtg =      gxx[index] * gyy[index] * gzz[index] 
           + 2.0*gxy[index] * gxz[index] * gyz[index] 
           -   gxx[index] * pow(gyz[index],2.0)           
           -   gyy[index] * pow(gxz[index],2.0)           
           -   gzz[index] * pow(gxy[index],2.0);
           
    cofactor(gg,igg,4);

    igtt[index] = igg[0][0];
    igtx[index] = igg[0][1];
    igty[index] = igg[0][2];
    igtz[index] = igg[0][3];
    igxx[index] = igg[1][1];
    igyy[index] = igg[2][2];
    igzz[index] = igg[3][3];
    igxy[index] = igg[1][2];
    igxz[index] = igg[1][3];
    igyz[index] = igg[2][3];
    
    sqrtdetg[index] = alpha[index] * sqrt(fabs(dtg));
    
        }
     }
  }
}

