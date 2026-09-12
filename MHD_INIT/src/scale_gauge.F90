! MHD_Init
! scale_gauge: Gauge depends on Scale Factor:

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



module scale_gauge
  use constants_init
  implicit none
  private
  public makescale
  !
contains
  !
  subroutine makescale(alph, beta)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(out) :: alph, beta(3)
    !
    ! alpha = inital scale factor
    alph = aa0
    !
    ! betai = 0
    beta(:) = 0
    !
  end subroutine makescale
  !
end module scale_gauge
