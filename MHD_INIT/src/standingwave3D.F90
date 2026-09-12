! MHD_Init
! Standing Wave 3D: Initial Data for Linear Standing Waves:
! Written by David Garrison

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module standingwave3D 
  use constants_init
  implicit none
  private
  public make_swave3D
  !
contains
  !
  subroutine make_swave3D(CCTK_ARGUMENTS, xx, time, lx, ly, lz, g, kk)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    !
    CCTK_REAL plusx, plusdotx, crossx, crossdotx
    CCTK_REAL plusy, plusdoty, crossy, crossdoty
    CCTK_REAL plusz, plusdotz, crossz, crossdotz
    CCTK_REAL wcx, wcy, wcz 
    CCTK_REAL, intent(in) :: xx(3), time, lx, ly, lz 
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer i,n,nn
    !
    !
    g(:,:) = 0.0
    forall (i=1:3) g(i,i) = 1.0 
    !
    kk(:,:) = 0.0
    !
    nn = cctk_gsh(3)
    !
    if (cctk_gsh(2) < cctk_gsh(3)) then
       nn = cctk_gsh(2)
       if (cctk_gsh(1) < cctk_gsh(2)) then
          nn = cctk_gsh(1)
       end if
    end if
    if (cctk_gsh(1) < cctk_gsh(3)) then
       nn = cctk_gsh(1)
       if (cctk_gsh(2) < cctk_gsh(1)) then
          nn = cctk_gsh(2)
       end if
    end if
    !
    do n = 2,nn/4,2
    !
    wcx = n*pi/lx
    wcy = n*pi/ly
    wcz = n*pi/lz
    !
    plusx = iplusA*SIN(wcx*xx(1))*COS(wcx*time)
    plusdotx = -wcx*iplusA*SIN(wcx*xx(1))*SIN(wcx*time)
    plusy = iplusA*SIN(wcy*xx(2))*COS(wcy*time)
    plusdoty = -wcy*iplusA*SIN(wcy*xx(2))*SIN(wcy*time)
    plusz = iplusA*SIN(wcz*xx(3))*COS(wcz*time)
    plusdotz = -wcz*iplusA*SIN(wcz*xx(3))*SIN(wcz*time)
    !
    !
    g(1,1) = g(1,1) + plusz - plusy
    !
    g(2,2) = g(2,2) - plusz + plusx
    !
    g(3,3) = g(3,3) + plusy - plusx
    !
    kk(1,1) = kk(1,1) - 0.5*plusdotz + 0.5*plusdoty
    !
    kk(2,2) = kk(2,2) + 0.5*plusdotz - 0.5*plusdotx
    !
    kk(3,3) = kk(3,3) - 0.5*plusdoty + 0.5*plusdotx
    !
    !
    crossx = icrossB*SIN(wcx*xx(1))*COS(wcx*time)
    crossdotx = -wcx*icrossB*SIN(wcx*xx(1))*SIN(wcx*time) 
    crossy = icrossB*SIN(wcy*xx(2))*COS(wcy*time)
    crossdoty = -wcy*icrossB*SIN(wcy*xx(2))*SIN(wcy*time) 
    crossz = icrossB*SIN(wcz*xx(3))*COS(wcz*time)
    crossdotz = -wcz*icrossB*SIN(wcz*xx(3))*SIN(wcz*time) 
    !
    !
    g(1,2) = g(1,2) + crossz
    g(2,1) = g(2,1) + crossz
    !
    g(1,3) = g(1,3) + crossy
    g(3,1) = g(3,1) + crossy
    !
    g(2,3) = g(2,3) + crossx
    g(3,2) = g(3,2) + crossx
    !
    kk(1,2) = kk(1,2) - 0.5*crossdotz
    kk(2,1) = kk(2,1) - 0.5*crossdotz
    !
    kk(1,3) = kk(1,3) - 0.5*crossdoty
    kk(3,1) = kk(3,1) - 0.5*crossdoty
    !
    kk(2,3) = kk(2,3) - 0.5*crossdotx
    kk(3,2) = kk(3,2) - 0.5*crossdotx
    !
    end do
    !
  end subroutine make_swave3D 
  !
end module standingwave3D  
