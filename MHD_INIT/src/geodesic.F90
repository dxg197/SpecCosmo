! MHD_Init
! geodesic: Geodesic gauge condition:
! Provide initial data and boundary conditions for a geodesic gauge.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



module geodesic
  implicit none
  private
  public makegeodesic
  !
contains
  !
  subroutine makegeodesic(alph, beta)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(out) :: alph, beta(3)
    !
    ! alpha = 1.0
    alph = 1.0
    !
    ! betai = 0
    beta(:) = 0
    !
  end subroutine makegeodesic
  !
end module geodesic
