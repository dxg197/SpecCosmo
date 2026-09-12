#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative_FD4( CCTK_ARGUMENTS, CCTK_REAL ****tensor )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, m, n, p, istart, jstart, kstart, iend, jend, kend;
    CCTK_INT uu_index[4], vv_index[4],index;
    CCTK_INT uu2_index[4], vv2_index[4];
    CCTK_INT uuindex, vvindex;
    CCTK_INT uu2index, vv2index;

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
		uu_index[1] = CCTK_GFINDEX3D( cctkGH, i+1, j, k );
        vv_index[1] = CCTK_GFINDEX3D( cctkGH, i-1, j, k );
        uu_index[2] = CCTK_GFINDEX3D( cctkGH, i, j+1, k );
        vv_index[2] = CCTK_GFINDEX3D( cctkGH, i, j-1, k );
		uu_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k+1 );
        vv_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k-1 );
        uu2_index[1] = CCTK_GFINDEX3D( cctkGH, i+2, j, k );
        vv2_index[1] = CCTK_GFINDEX3D( cctkGH, i-2, j, k );
        uu2_index[2] = CCTK_GFINDEX3D( cctkGH, i, j+2, k );
        vv2_index[2] = CCTK_GFINDEX3D( cctkGH, i, j-2, k );
		uu2_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k+2 );
        vv2_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k-2 );

        for (m = 1; m < 4; m++) {
            for (n = 1; n < 4; n++) {
                for (p = 1; p < 4; p++) {
                    uuindex = uu_index[m];
                    vvindex = vv_index[m];
                    uu2index = uu2_index[m];
                    vv2index = vv2_index[m];

                    tensor[m][n][p][index] = (tensor[0][n][p][vv2index] - 8.0*tensor[0][n][p][vvindex] + 
                                             8.0*tensor[0][n][p][uuindex] - tensor[0][n][p][uu2index]) / 
                                              (12.0*CCTK_DELTA_SPACE(m-1)); 

                      } 
                   }
                }
			}
		}
	}
}

