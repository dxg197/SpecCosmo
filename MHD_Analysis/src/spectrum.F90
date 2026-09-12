/*@@ Calculates the output data for GWs @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine analysis4(CCTK_ARGUMENTS)
    use constants_analysis
    implicit none
	
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    integer ierr, handle, vindex
    integer i,j,k
    integer istart,jstart,kstart,iend,jend,kend
    CCTK_REAL results, npstrain, npgwspec

    ! -----------------
    ! Do the 
    ! Calculate avg values
    
    istart = 1 
    jstart = 1 
    kstart = 1 

    iend = cctk_lsh(1) 
    jend = cctk_lsh(2) 
    kend = cctk_lsh(3) 
    
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
             
    		if (gwspec_PSD(i,j,k).ne.0.0) then
        		npgwspec_temp(i,j,k) = 1.0
        	else
        		npgwspec_temp(i,j,k) = 0.0
    		endif
        
    		if (geostrain(i,j,k).ne.0.0) then
        		npstrain_temp(i,j,k) = 1.0
        	else
        		npstrain_temp(i,j,k) = 0.0
    		endif
    
    	  end do
       end do
    end do

    call CCTK_ReductionHandle(handle, "sum")
    results = 0.0
    
    call CCTK_VarIndex(vindex, "mhd_analysis::npstrain_temp")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    npstrain = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of npstrain_temp failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::npgwspec_temp")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    npgwspec = results
    
    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of npgwspec_temp failed!");
    endif

    call CCTK_VarIndex(vindex, "mhd_analysis::geostrain")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    strain_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of geostrain failed!");
    endif
    
    call CCTK_VarIndex(vindex, "mhd_analysis::gwspec_PSD")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    gwspec_avg = results

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of gwspec_PSD failed!");
    endif
        
    if ((npstrain > 0.0).and.(npgwspec > 0.0)) then
      strain_avg = strain_avg/npstrain
	  gwspec_avg = gwspec_avg/npgwspec
	endif
 
    
end subroutine analysis4
