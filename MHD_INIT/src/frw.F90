! MHD_Init
! FRW: Initial Data for FRW spacetime with linear GW:
! Written by David Garrison
 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module frw
  use cactus
  use constants_init
  implicit none
  private
  public make_frw
  
contains
  
  subroutine make_frw(CCTK_ARGUMENTS, g, kk, i,j,k)
  	DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    INTEGER, intent(in) :: i,j,k
    CCTK_REAL plus, plusdot, cross, crossdot, HH_geo
    CCTK_REAL pert(3,3), hhdot(3,3), rand, gwamp, rho_0, BB_0, pp_0
    CCTK_REAL bb3(3), eps, beta(3), shift2, rho_c, kw,rho_d0,rad0 
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer l,m,n, zero
    
    zero = 0
    
    pert(:,:) = 0.0
    
    hhdot(:,:) = 0.0
    
    bb3(1) = iBx*Gauss
    bb3(2) = iBy*Gauss
    bb3(3) = iBz*Gauss
    
    beta(1) = bb3(1)/(alpha(i,j,k)*sqrt(4.0*pi))
    beta(2) = bb3(2)/(alpha(i,j,k)*sqrt(4.0*pi))
    beta(3) = bb3(3)/(alpha(i,j,k)*sqrt(4.0*pi))
 
	BB_0 = beta(1)*beta(1) + beta(2)*beta(2) + beta(3)*beta(3)
	
	rho_0 = irho*Density
	rho_d0 = irho_dark*Density
	rad0 = irad*Density
	
	 if (rho_0.ne.zero) then
       !eps = (eps_cof*Rad_Const*SIE*(itemp**4))/irho
       eps = eps_cof*Rad_Const*SIE*(itemp**4)/irho - irho_dark/irho - irad/irho - 1.0
    else
       eps = 0.0
    end if
	
	pp_0 = (gam - 1.0)*eps*rho_0
	
	rho_c = rho_0*(1.0 + eps) + BB_0*(1.0-0.5/alpha(i,j,k)**2) + rho_d0 + rad0
	
	HH_geo = sqrt(8.0*pi*Flatness_Ratio*rho_c/3.0)*aa0
    
    !HH_geo = HH/(3.085678D19)
    !HH_geo = HH/(9.25062992217D27)
    
    !kw = HH*aa/100.0
    
    !gwamp = pAmp /(aa0*sqrt(2.0*kw))
    
    gwamp = Amp
    
    call RANDOM_SEED()
    
    ! for standing wave in the xy plane

	call RANDOM_NUMBER(rand)
    		plus = gwamp*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    		plusdot = gwamp*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    		cross = gwamp*2.0*(rand - 0.5)
    		
    call RANDOM_NUMBER(rand)
    		crossdot = gwamp*2.0*(rand - 0.5)
             
    		pert(1,1) = pert(1,1) + plus
    		pert(2,2) = pert(2,2) - plus
    
    		hhdot(1,1) = hhdot(1,1) + plusdot
    		hhdot(2,2) = hhdot(2,2) - plusdot
    
    		pert(1,2) = pert(1,2) + cross
    		pert(2,1) = pert(2,1) + cross
    
    		hhdot(1,2) = hhdot(1,2) + crossdot
    		hhdot(2,1) = hhdot(2,1) + crossdot
    
    ! for standing wave in the xz plane
     
    call RANDOM_NUMBER(rand)
    		plus = gwamp*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    		plusdot = gwamp*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    		cross = gwamp*2.0*(rand - 0.5)
    		
    call RANDOM_NUMBER(rand)
    		crossdot = gwamp*2.0*(rand - 0.5)
    		          
    		pert(1,1) = pert(1,1) - plus
    		pert(3,3) = pert(3,3) + plus
    
    		hhdot(1,1) = hhdot(1,1) - plusdot
    		hhdot(3,3) = hhdot(3,3) + plusdot
    
    		pert(1,3) = pert(1,3) + cross
    		pert(3,1) = pert(3,1) + cross
    
    		hhdot(1,3) = hhdot(1,3) + crossdot
    		hhdot(3,1) = hhdot(3,1) + crossdot
    
    ! for standing wave in the yz plane
    
    call RANDOM_NUMBER(rand)
    		plus = gwamp*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    		plusdot = gwamp*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    		cross = gwamp*2.0*(rand - 0.5)
    		
    call RANDOM_NUMBER(rand)
    		crossdot = gwamp*2.0*(rand - 0.5)
    		          
    		pert(2,2) = pert(2,2) + plus
    		pert(3,3) = pert(3,3) - plus
    
    		hhdot(2,2) = hhdot(2,2) + plusdot
    		hhdot(3,3) = hhdot(3,3) - plusdot
    
    		pert(2,3) = pert(2,3) + cross
    		pert(3,2) = pert(3,2) + cross
    
    		hhdot(2,3) = hhdot(2,3) + crossdot
    		hhdot(3,2) = hhdot(3,2) + crossdot
    
    g(:,:) = 0.0
    forall (l=1:3) g(l,l) = aa0**2 
    
    kk(:,:) = 0.0
    forall (l=1:3) kk(l,l) = -HH_geo*aa0**2/alpha(i,j,k)
    
    g(1,1) = g(1,1) + pert(1,1)*aa0**2
    g(2,2) = g(2,2) + pert(2,2)*aa0**2
    g(3,3) = g(3,3) + pert(3,3)*aa0**2
    g(1,2) = g(1,2) + pert(1,2)*aa0**2
    g(1,3) = g(1,3) + pert(1,3)*aa0**2
    g(2,3) = g(2,3) + pert(2,3)*aa0**2
    g(2,1) = g(1,2)
    g(3,1) = g(1,3)
    g(3,2) = g(2,3)
    
    kk(1,1) = kk(1,1) - 0.5*hhdot(1,1)*aa0**2/alpha(i,j,k) - HH_geo*pert(1,1)*aa0**2/alpha(i,j,k)
    kk(2,2) = kk(2,2) - 0.5*hhdot(2,2)*aa0**2/alpha(i,j,k) - HH_geo*pert(2,2)*aa0**2/alpha(i,j,k)
    kk(3,3) = kk(3,3) - 0.5*hhdot(3,3)*aa0**2/alpha(i,j,k) - HH_geo*pert(3,3)*aa0**2/alpha(i,j,k)
    kk(1,2) = kk(1,2) - 0.5*hhdot(1,2)*aa0**2/alpha(i,j,k) - HH_geo*pert(1,2)*aa0**2/alpha(i,j,k)
    kk(1,3) = kk(1,3) - 0.5*hhdot(1,3)*aa0**2/alpha(i,j,k) - HH_geo*pert(1,3)*aa0**2/alpha(i,j,k)
    kk(2,3) = kk(2,3) - 0.5*hhdot(2,3)*aa0**2/alpha(i,j,k) - HH_geo*pert(2,3)*aa0**2/alpha(i,j,k)
    kk(2,1) = kk(1,2)
    kk(3,1) = kk(1,3)
    kk(3,2) = kk(2,3)
    
  end subroutine make_frw
  
end module frw 
