#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

void SpecDeriv_Tensor_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative2_FD( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, m, n, p, q, istart, jstart, kstart, iend, jend, kend;
    CCTK_INT uu_index[4][4], vv_index[4][4], index;
    CCTK_INT uv_index[4][4], vu_index[4][4];
    CCTK_INT uuindex, vvindex, uvindex, vuindex;
	istart = cctk_nghostzones[0];
	jstart = cctk_nghostzones[1];
	kstart = cctk_nghostzones[2];
	iend   = cctk_lsh[ 0 ] - cctk_nghostzones[0];
	jend   = cctk_lsh[ 1 ] - cctk_nghostzones[1];
	kend   = cctk_lsh[ 2 ] - cctk_nghostzones[2];
    
    for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
            {
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
				uu_index[1][1] = CCTK_GFINDEX3D( cctkGH, i+1, j, k );
                vv_index[1][1] = CCTK_GFINDEX3D( cctkGH, i-1, j, k );
                uu_index[2][2] = CCTK_GFINDEX3D( cctkGH, i, j+1, k );
                vv_index[2][2] = CCTK_GFINDEX3D( cctkGH, i, j-1, k );
				uu_index[3][3] = CCTK_GFINDEX3D( cctkGH, i, j, k+1 );
                vv_index[3][3] = CCTK_GFINDEX3D( cctkGH, i, j, k-1 );
                uu_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+1, j+1, k );
                vv_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-1, j-1, k );
                uv_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+1, j-1, k );
                vu_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-1, j+1, k );
                uu_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+1, j, k+1 );
                vv_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-1, j, k-1 );
                uv_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+1, j, k-1 );
                vu_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-1, j, k+1 );
                uu_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+1, k+1 );
                vv_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-1, k-1 );
                uv_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+1, k-1 );
                vu_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-1, k+1 );
                
                uu_index[2][1] = uu_index[1][2];
                vv_index[2][1] = vv_index[1][2];
                uv_index[2][1] = uv_index[1][2];
                vu_index[2][1] = vu_index[1][2];
                uu_index[3][1] = uu_index[1][3];
                vv_index[3][1] = vv_index[1][3];
                uv_index[3][1] = uv_index[1][3];
                vu_index[3][1] = vu_index[1][3];
                uu_index[3][2] = uu_index[2][3];
                vv_index[3][2] = vv_index[2][3];
                uv_index[3][2] = uv_index[2][3];
                vu_index[3][2] = vu_index[2][3];
                
                
                for (m = 1; m < 4; m++) {
                   for (n = 1; n < 4; n++) {
                      for (p = 1; p < 4; p++) {
                   uuindex = uu_index[m][m];
                   vvindex = vv_index[m][m];
                   tensor2[0][m][n][p][index] = (tensor2[0][0][n][p][uuindex] - tensor2[0][0][n][p][vvindex]) / 
                                                                                  (2.0*CCTK_DELTA_SPACE(m-1));
                   tensor2[m][0][n][p][index] = tensor2[0][m][n][p][index];   
                                          
                   for (q = 1; q < 4; q++) {
                      if (q != m) {
                      uuindex = uu_index[m][q];
                      vvindex = vv_index[m][q];
                      uvindex = uv_index[m][q];
                      vuindex = vu_index[m][q];
                      tensor2[m][q][n][p][index] = (tensor2[0][0][n][p][uuindex] - tensor2[0][0][n][p][uvindex] -
                                                    tensor2[0][0][n][p][vuindex] + tensor2[0][0][n][p][vvindex]) / 
                                                               (4.0*CCTK_DELTA_SPACE(m-1)*CCTK_DELTA_SPACE(q-1));
                      }
                   } 
                   uuindex = uu_index[m][m];
                   vvindex = vv_index[m][m];
                   tensor2[m][m][n][p][index] = (tensor2[0][0][n][p][uuindex] - 2.0*tensor2[0][0][n][p][index] +
                                                 tensor2[0][0][n][p][vvindex]) / (CCTK_DELTA_SPACE(m-1)*CCTK_DELTA_SPACE(m-1));
                      }
                   }
                }                                                      
			}
		}
	}
}