/*@@ Calculates the output data for Plasma Field @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine analysis2(CCTK_ARGUMENTS)
    use constants_analysis
    implicit none
	
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    CCTK_REAL results
    integer ierr, handle, vindex

    ! -----------------
    ! Do the 
    ! Calculate avg values

    call CCTK_ReductionHandle(handle, "average")
    results = 0.0

    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of rho_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rad_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rad_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of rad_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_dark_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_dark_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of rho_dark_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out_total")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_out_total_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of rho_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of temp_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_rel")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_out_rel_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of temp_out_rel_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_nonrel")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_out_nonrel_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of temp_out_nonrel_avg failed!");
    endif
   
    call CCTK_VarIndex(vindex, "mhd_analysis::pp_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    pp_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of pp_out_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "mhd_analysis::b_x_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_x_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of b_x_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::b_y_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_y_out_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of b_y_out_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::b_z_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_z_out_avg = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of b_z_out_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "SpecGRMHD::vv_x")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    vv_x_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of vv_x_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::vv_y")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    vv_y_avg = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of vv_y_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "SpecGRMHD::vv_z")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    vv_z_avg = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of vv_z_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "SpecGRMHD::gxx")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gxx_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gxx_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::gyy")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gyy_avg = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gyy_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "SpecGRMHD::gzz")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gzz_avg = results
 
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gzz_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::gxy")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gxy_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gxy_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::gxz")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gxz_avg = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gxz_avg failed!");
    endif

    call CCTK_VarIndex(vindex, "SpecGRMHD::gyz")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gyz_avg = results
 
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gyz_avg failed!");
    endif
   
    call CCTK_VarIndex(vindex, "SpecGRMHD::sqrtdetg")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    sqrtdetg_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of sqrtdetg_avg failed!");
    endif
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::alpha")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    alp_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of alp_avg failed!");
    endif
    
end subroutine analysis2
