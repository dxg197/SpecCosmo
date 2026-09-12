! MHD_Init
! Bondi: Initial Data for Linear Plane Waves:
! Written by David Garrison

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module bondi 
  use constants_init
  implicit none
  private
  public make_bondi 
  !
contains
  !
  subroutine make_bondi(xx, time, lz, g, kk)
    DECLARE_CCTK_PARAMETERS
    !
    CCTK_REAL plus, plusdot
    CCTK_REAL wc, v 
    CCTK_REAL, intent(in) :: xx(3), time, lz 
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer i
    !
    wc = nw*2.0*pi/lz
    v = time + xx(3) 
    !
    g(:,:) = 0.0
    forall (i=1:3) g(i,i) = 1.0 
    !
    kk(:,:) = 0.0
    !
    plus = iplusA*SIN(wc*v)
    plusdot = wc*iplusA*COS(wc*v)
    !
    !
    g(1,1) = g(1,1) + plus
    !
    g(2,2) = g(2,2) - plus
    !
    g(3,3) = 1.0
    !
    kk(1,1) = -0.5*plusdot
    !
    kk(2,2) =  0.5*plusdot
    !
    kk(3,3) = 0.0
    !
    !
  end subroutine make_bondi 
  !
end module bondi  
