! MHD_Init
! Kasner: Analytic Kasner Solution:
! Written by David Garrison

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module kasner 
  implicit none
  private
  public makekasner 
  !
contains
  !
  subroutine makekasner (time, g, kk)
    DECLARE_CCTK_PARAMETERS 
    CCTK_REAL, intent(in) :: time
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3) 
    CCTK_REAL p1, p2, p3
    integer i
    !
    g(:,:) = 0.0
    forall (i=1:3) g(i,i) = 1.0 
    !
    ! kij = 0
    kk(:,:) = 0.0
    !
    p1 = (kasner_param**2-1.0)/(kasner_param**2+3.0)
    p2 = 2.0*(1.0+kasner_param)/(kasner_param**2+3.0)
    p3 = 2.0*(1.0-kasner_param)/(kasner_param**2+3.0)
    !
    g(1,1) = (time+1.0)**(2.0*p1) 
     
    g(2,2) = (time+1.0)**(2.0*p2) 
    
    g(3,3) = (time+1.0)**(2.0*p3) 
    
    kk(1,1) = -p1*(time+1.0)**(2.0*p1-1.0) 
    
    kk(2,2) = -p2*(time+1.0)**(2.0*p2-1.0)  
    
    kk(3,3) = -p3*(time+1.0)**(2.0*p3-1.0)   
    ! 
  end subroutine makekasner  
  !
end module kasner  
