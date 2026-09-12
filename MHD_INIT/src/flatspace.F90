! MHD_Init
! flatspace: Analytic solution for flat space:
! Provide initial data and boundary conditions for a flat space.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module flatspace
  implicit none
  private
  public makeflatspace
  !
contains
  !
  subroutine makeflatspace(g, kk)
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer i
    !
    ! gij = deltaij
    g(:,:) = 0.0
    forall (i=1:3) g(i,i) = 1.0
    !
    ! kij = 0.0
    kk(:,:) = 0.0
    !
  end subroutine makeflatspace
  !
end module flatspace
