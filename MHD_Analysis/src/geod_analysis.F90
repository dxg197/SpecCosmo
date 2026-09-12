! MHD_Anlaysis

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



module geod_analysis
  implicit none
  private
  public makegeod
  !
contains
  !
  subroutine makegeod(phi, alph, beta)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in) :: phi
    CCTK_REAL, intent(out) :: alph, beta(3)
    !
    alph = exp(12.d0*phi*sigma)
    !
    ! betai = 0
    beta(:) = 0
    !
  end subroutine makegeod
  !
end module geod_analysis