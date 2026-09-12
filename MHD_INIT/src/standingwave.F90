! MHD_Init
! Standing Wave: Initial Data for Linear Standing Waves:
! Written by David Garrison

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module standingwave 
  use constants_init
  implicit none
  private
  public make_swave 
  !
contains
  !
  subroutine make_swave(xx, time, lz, g, kk)
    DECLARE_CCTK_PARAMETERS
    !
    CCTK_REAL plus, plusdot, cross, crossdot
    CCTK_REAL wc, HH_geo
    CCTK_REAL, intent(in) :: xx(3), time, lz 
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer i
    !
    wc = nw*2.0*pi/lz
    HH_geo = HH/(3.085678D19*Speed_Light)
    !
    g(:,:) = 0.0
    forall (i=1:3) g(i,i) = aa0**2 
    !
    kk(:,:) = 0.0
    forall (i=1:3) kk(i,i) = -HH_geo*aa0**2
    !
    plus = iplusA*SIN(wc*xx(3))*COS(wc*time)
    plusdot = -wc*iplusA*SIN(wc*xx(3))*SIN(wc*time)
    !
    !
    g(1,1) = g(1,1) + plus*aa0**2
    !
    g(2,2) = g(2,2) - plus*aa0**2
    !
    kk(1,1) = kk(1,1) - 0.5*plusdot*aa0**2 - HH_geo*plus*aa0**2
    !
    kk(2,2) = kk(2,2) + 0.5*plusdot*aa0**2 + HH_geo*plus*aa0**2
    !
    cross = icrossB*SIN(wc*xx(3))*COS(wc*time)
    crossdot = -wc*icrossB*SIN(wc*xx(3))*SIN(wc*time) 
    !
    !
    g(1,2) = cross*aa0**2
    !
    g(2,1) = cross*aa0**2
    !
    kk(1,2) = -0.5*crossdot*aa0**2 - HH_geo*cross*aa0**2
    !
    kk(2,1) = -0.5*crossdot*aa0**2 - HH_geo*cross*aa0**2
    !
    !
  end subroutine make_swave 
  !
end module standingwave  
