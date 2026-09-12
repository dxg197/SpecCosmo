! MHD_Analysis

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module temp_analysis
  implicit none
  private
  public calc_temp
  !
contains
  !
  subroutine calc_temp (aa1, bb1, cc1, dd1, ee1, ff1, temp_o, temp_rel, temp_nonrel)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(out) :: temp_o
    CCTK_REAL, intent(in)  :: aa1, bb1, cc1, dd1, ee1, ff1
    CCTK_REAL xx1, xx2, xx3, temp_rel, temp_nonrel
    
    !
    temp_rel = ((dd1 - aa1)/cc1)**0.25
    
    temp_nonrel = ee1/bb1
    
    IF ((dd1-aa1) > aa1) THEN
    temp_o = temp_rel
    ELSE
    temp_o = temp_nonrel
    END IF
    
    !print '("temp_rel = "E50.3)',temp_rel
    
    !print '("temp_nonrel = "E50.3)',temp_nonrel
    
    
    !xx1 = 5.0396842*(dd1-aa1)
    !xx2 = (SQRT(729.0*(bb1**4.0)*(cc1**2.0)-6912.0*(cc1**3.0)*((dd1-aa1)**3.0)) + 27.0*(bb1**2.0)*cc1)**(0.33333333)
    !xx3 = 3.77976315*cc1
    
    !temp_o = 0.5*SQRT(2.0*bb1/(cc1*SQRT(-xx1/xx2+xx2/xx3)) + xx1/xx2 - xx2/xx3)) - 0.5*SQRT(-xx1/xx2 + xx2/xx3)
    
    !
  end subroutine calc_temp
  !
end module temp_analysis

