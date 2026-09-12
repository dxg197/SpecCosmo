#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

void SpecDeriv_Tensor_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative2_FD4( CCTK_ARGUMENTS, CCTK_REAL *****tensor2 )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, m, n, p, q, istart, jstart, kstart, iend, jend, kend;
    CCTK_INT uu_index[4][4], vv_index[4][4], index;
    CCTK_INT uv_index[4][4], vu_index[4][4];
    CCTK_INT uu_2_index[4][4], vv_2_index[4][4];
    CCTK_INT u2v_index[4][4], v2u_index[4][4];
    CCTK_INT uv2_index[4][4], vu2_index[4][4];
    CCTK_INT u2v2_index[4][4], v2u2_index[4][4];
    CCTK_INT u2u_index[4][4], v2v_index[4][4];
    CCTK_INT uu2_index[4][4], vv2_index[4][4];
    CCTK_INT u2u2_index[4][4], v2v2_index[4][4];
    CCTK_INT uuindex, vvindex, uvindex, vuindex;
    CCTK_INT uu_2index, vv_2index;
    CCTK_INT u2vindex, v2uindex, uv2index, vu2index;
    CCTK_INT u2v2index, v2u2index, u2uindex, v2vindex;
    CCTK_INT uu2index, vv2index, u2u2index, v2v2index;

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
        
        uu_2_index[1][1] = CCTK_GFINDEX3D( cctkGH, i+2, j, k );
        vv_2_index[1][1] = CCTK_GFINDEX3D( cctkGH, i-2, j, k );
        uu_2_index[2][2] = CCTK_GFINDEX3D( cctkGH, i, j+2, k );
        vv_2_index[2][2] = CCTK_GFINDEX3D( cctkGH, i, j-2, k );
		uu_2_index[3][3] = CCTK_GFINDEX3D( cctkGH, i, j, k+2 );
        vv_2_index[3][3] = CCTK_GFINDEX3D( cctkGH, i, j, k-2 );
        uu_2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+2, j+2, k );
        vv_2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-2, j-2, k );
        uu_2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+2, j, k+2 );
        vv_2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-2, j, k-2 );
		uu_2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+2, k+2 );
        vv_2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-2, k-2 );
        
        uu_2_index[2][1] = uu_2_index[1][2];
        vv_2_index[2][1] = vv_2_index[1][2];
        uu_2_index[3][1] = uu_2_index[1][3];
        vv_2_index[3][1] = vv_2_index[1][3];
        uu_2_index[3][2] = uu_2_index[2][3];
        vv_2_index[3][2] = vv_2_index[2][3];
		
        u2v_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+2, j-1, k );
        v2u_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-2, j+1, k );
        uv2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+1, j-2, k );
        vu2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-1, j+2, k );
        u2v2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+2, j-2, k );
        v2u2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-2, j+2, k );
        u2u_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+2, j+1, k );
        v2v_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-2, j-1, k );
        uu2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+1, j+2, k );
        vv2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-1, j-2, k );
        u2u2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i+2, j+2, k );
        v2v2_index[1][2] = CCTK_GFINDEX3D( cctkGH, i-2, j-2, k );
        
        u2v_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+2, j, k-1 );
        v2u_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-2, j, k+1 );
        uv2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+1, j, k-2 );
        vu2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-1, j, k+2 );
        u2v2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+2, j, k-2 );
        v2u2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-2, j, k+2 );
        u2u_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+2, j, k+1 );
        v2v_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-2, j, k-1 );
        uu2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+1, j, k+2 );
        vv2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-1, j, k-2 );
        u2u2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i+2, j, k+2 );
        v2v2_index[1][3] = CCTK_GFINDEX3D( cctkGH, i-2, j, k-2 );
        
        u2v_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+2, k-1 );
        v2u_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-2, k+1 );
        uv2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+1, k-2 );
        vu2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-1, k+2 );
        u2v2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+2, k-2 );
        v2u2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-2, k+2 );
        u2u_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+2, k+1 );
        v2v_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-2, k-1 );
        uu2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+1, k+2 );
        vv2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-1, k-2 );
        u2u2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j+2, k+2 );
        v2v2_index[2][3] = CCTK_GFINDEX3D( cctkGH, i, j-2, k-2 );
                
        u2v_index[2][1] = u2v_index[1][2];
        v2u_index[2][1] = v2u_index[1][2];
        uv2_index[2][1] = uv2_index[1][2];
        vu2_index[2][1] = vu2_index[1][2];
        u2v2_index[2][1] = u2v2_index[1][2];
        v2u2_index[2][1] = v2u2_index[1][2];
        u2u_index[2][1] = u2u_index[1][2];
        v2v_index[2][1] = v2v_index[1][2];
        uu2_index[2][1] = uu2_index[1][2];
        vv2_index[2][1] = vv2_index[1][2];
        u2u2_index[2][1] = u2u2_index[1][2];
        v2v2_index[2][1] = v2v2_index[1][2];
                
        u2v_index[3][1] = u2v_index[1][3];
        v2u_index[3][1] = v2u_index[1][3];
        uv2_index[3][1] = uv2_index[1][3];
        vu2_index[3][1] = vu2_index[1][3];
        u2v2_index[3][1] = u2v2_index[1][3];
        v2u2_index[3][1] = v2u2_index[1][3];
        u2u_index[3][1] = u2u_index[1][3];
        v2v_index[3][1] = v2v_index[1][3];
        uu2_index[3][1] = uu2_index[1][3];
        vv2_index[3][1] = vv2_index[1][3];
        u2u2_index[3][1] = u2u2_index[1][3];
        v2v2_index[3][1] = v2v2_index[1][3];
                
        vu_index[3][2] = vu_index[2][3];
        u2v_index[3][2] = u2v_index[2][3];
        v2u_index[3][2] = v2u_index[2][3];
        uv2_index[3][2] = uv2_index[2][3];
        vu2_index[3][2] = vu2_index[2][3];
        u2v2_index[3][2] = u2v2_index[2][3];
        v2u2_index[3][2] = v2u2_index[2][3];
        u2u_index[3][2] = u2u_index[2][3];
        v2v_index[3][2] = v2v_index[2][3];
        uu2_index[3][2] = uu2_index[2][3];
        vv2_index[3][2] = vv2_index[2][3];
        u2u2_index[3][2] = u2u2_index[2][3];
        v2v2_index[3][2] = v2v2_index[2][3];
                
                
        for (m = 1; m < 4; m++) {
            for (n = 1; n < 4; n++) {
                for (p = 1; p < 4; p++) {
                   uuindex = uu_index[m][m];
                   vvindex = vv_index[m][m];
                   vv_2index = vv_2_index[m][m];
                   uu_2index = uu_2_index[m][m];

                   tensor2[0][m][n][p][index] = (tensor2[0][0][n][p][vv_2index] - 
                                               8.0*tensor2[0][0][n][p][vvindex] + 
                                               8.0*tensor2[0][0][n][p][uuindex] - 
                                               tensor2[0][0][n][p][uu_2index]) / 
                                               (12.0*CCTK_DELTA_SPACE(m-1));

                   tensor2[m][0][n][p][index] = tensor2[0][m][n][p][index];   
                                          
                   for (q = 1; q < 4; q++) {
                      if (q != m) {
                      uuindex = uu_index[m][q];
                      vvindex = vv_index[m][q];
                      uvindex = uv_index[m][q];
                      vuindex = vu_index[m][q];
                      u2vindex = u2v_index[m][q];
                      v2uindex = v2u_index[m][q];
                      uv2index = uv2_index[m][q];
                      vu2index = vu2_index[m][q];
                      u2v2index = u2v2_index[m][q];
                      v2u2index = v2u2_index[m][q];
                      u2uindex = u2u_index[m][q];
                      v2vindex = v2v_index[m][q];
                      uu2index = uu2_index[m][q];
                      vv2index = vv2_index[m][q];
                      u2u2index = u2u2_index[m][q];
                      v2v2index = v2v2_index[m][q];
                      
                      tensor2[m][q][n][p][index] = ( tensor2[0][0][n][p][v2v2index] - 
                                                    8.0 * tensor2[0][0][n][p][v2vindex] + 
                                                    8.0 * tensor2[0][0][n][p][v2uindex] - 
                                                    tensor2[0][0][n][p][v2u2index] - 
                                                    8.0 * tensor2[0][0][n][p][vv2index] + 
                                                    64.0 * tensor2[0][0][n][p][vvindex] - 
                                                    64.0 * tensor2[0][0][n][p][vuindex] + 
                                                    8.0 * tensor2[0][0][n][p][vu2index] + 
                                                    8.0 * tensor2[0][0][n][p][uv2index] - 
                                                    64.0 * tensor2[0][0][n][p][uvindex] + 
                                                    64.0 * tensor2[0][0][n][p][uuindex] - 
                                                    8.0 * tensor2[0][0][n][p][uu2index] - 
                                                    tensor2[0][0][n][p][u2v2index] + 
                                                    8.0 * tensor2[0][0][n][p][u2vindex] - 
                                                    8.0 * tensor2[0][0][n][p][u2uindex] + 
                                                    tensor2[0][0][n][p][u2u2index] ) * 
                                                    ( 1.0 / (144.0 * CCTK_DELTA_SPACE(m-1) * 
                                                    CCTK_DELTA_SPACE(q-1) ) );

                      }
                   } 

                   uuindex = uu_index[m][m];
                   vvindex = vv_index[m][m];
                   vv_2index = vv_2_index[m][m];
                   uu_2index = uu_2_index[m][m];
                   
                   tensor2[m][m][n][p][index] = ( -1.0 * tensor2[0][0][n][p][vv_2index] + 
                                                 16.0 * tensor2[0][0][n][p][vvindex] - 
                                                 30.0 * tensor2[0][0][n][p][index] + 
                                                 16.0 * tensor2[0][0][n][p][uuindex] - 
                                                 tensor2[0][0][n][p][uu_2index] ) / 
                                                 ( 12.0 * CCTK_DELTA_SPACE(m-1) * CCTK_DELTA_SPACE(m-1) );


                      }
                   }
                }                                                      
			}
		}
	}
}
