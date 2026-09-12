! MHD_Anlaysis

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



module scale_analysis
  use coords_analysis
  use constants_analysis
  implicit none
  private
  public makescale
  !
contains
  !
  subroutine makescale(alph, beta, scale)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in) :: scale
    CCTK_REAL, intent(out) :: alph, beta(3)
    !
    ! alpha = scale factor squared
    alph = scale
    !
    ! betai = 0
    beta(:) = 0
    !
  end subroutine makescale
  !
end module scale_analysis
