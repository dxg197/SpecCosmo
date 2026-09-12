/*@@ Calculates New Domain Size @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "CoordBase.h"
#include "util_Table.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

void Expansion(CCTK_ARGUMENTS);

void Expansion(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
/*  Declare local variables */

	char coord_name[16];
	CCTK_INT ierr, i, j, k, idx, handle, coord_handle;
	CCTK_REAL upper[3], lower[3];
	CCTK_REAL xmin, xmax, ymin, ymax, zmin, zmax;
	CCTK_REAL scale_factor;
	
	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &scale_factor, 1, CCTK_VarIndex("MHD_Analysis::aa_ratio_avg"));

	handle = CCTK_ReductionHandle("maximum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmax, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymax, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmax, 1, CCTK_VarIndex("grid::z"));
	
	handle = CCTK_ReductionHandle("minimum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmin, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymin, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmin, 1, CCTK_VarIndex("grid::z"));

	if (!isnan(scale_factor)) {
	
	upper[0] = scale_factor*xmax;
	upper[1] = scale_factor*ymax;
	upper[2] = scale_factor*zmax;
	
	lower[0] = scale_factor*xmin;
	lower[1] = scale_factor*ymin;
	lower[2] = scale_factor*zmin;
	
	cctkGH->cctk_delta_space[0] = scale_factor*CCTK_DELTA_SPACE(0);
	cctkGH->cctk_delta_space[1] = scale_factor*CCTK_DELTA_SPACE(1);
	cctkGH->cctk_delta_space[2] = scale_factor*CCTK_DELTA_SPACE(2);
	
	cctkGH->cctk_origin_space[0] = lower[0];
	cctkGH->cctk_origin_space[1] = lower[1];
	cctkGH->cctk_origin_space[2] = lower[2];
	
	}
	
	/* cart3d */
  	for (i = 0; i < 3; i++) {
    sprintf(coord_name, "%c", 'x' + i);

    coord_handle = Coord_CoordHandle(cctkGH, coord_name, "cart3d");
    if (coord_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Error retreiving coordinate handle for '%s' of cart3d",
                 coord_name);
    }
    sprintf(coord_name, "grid::%c", 'x' + i);
    ierr += Util_TableSetReal(coord_handle, lower[i], "COMPMIN");
    ierr += Util_TableSetReal(coord_handle, upper[i], "COMPMAX");
    ierr += Util_TableSetReal(coord_handle, cctk_delta_space[i], "DELTA");
    
  }
 
  for (int k = 0; k < cctk_lsh[2]; k++) {
    for (int j = 0; j < cctk_lsh[1]; j++) {
      for (int i = 0; i < cctk_lsh[0]; i++) {
        int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        x[idx] = CCTK_DELTA_SPACE(0) * (i + cctk_lbnd[0]) + CCTK_ORIGIN_SPACE(0);
        y[idx] = CCTK_DELTA_SPACE(1) * (j + cctk_lbnd[1]) + CCTK_ORIGIN_SPACE(1);
        z[idx] = CCTK_DELTA_SPACE(2) * (k + cctk_lbnd[2]) + CCTK_ORIGIN_SPACE(2);
        r[idx] = sqrt((x[idx]*x[idx]) + (y[idx]*y[idx]) + (z[idx]*z[idx]));
      }
    }
  }
      
}
    

