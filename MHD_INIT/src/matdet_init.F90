#include "cctk.h"

module matdet_init
  implicit none
  private
  
  public calc_symdet4
  
contains
  
  subroutine calc_symdet4 (g4, dtg4)
    CCTK_REAL, intent(in)  :: g4(4,4)
    CCTK_REAL, intent(out) :: dtg4
    CCTK_REAL :: tmp(4,4)
    integer   :: ipiv(4)
    integer   :: info
    integer   :: nperms
    integer   :: i
    character :: msg*1000
    
    tmp = g4
    !call sytrf (4, 4, tmp, 4, ipiv, info)
    
    if (info < 0) then
       write (msg, '("Error in call to SYTRF, info=",i4)') info
       call CCTK_WARN (1, trim(msg))
    end if
    
    if (info > 0) then
       dtg4 = 0
       return
    end if
    
    nperms = 0
    do i=1,4
       if (mod(ipiv(i),2) /= mod(i,2)) nperms = nperms+1
    end do
    if (mod(nperms,2) /=0) call CCTK_WARN (0, "internal error")
    nperms = nperms / 2
    
    dtg4 = 1
    if (mod(nperms,2)/=0) dtg4 = -1
    do i=1,4
       dtg4 = dtg4 * tmp(i,i)
    end do
  end subroutine calc_symdet4
  
end module matdet_init
