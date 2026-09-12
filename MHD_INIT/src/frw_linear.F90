! MHD_Init
! FRW: Initial Data for FRW spacetime with linear GW:
! Written by David Garrison
 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module frw_linear
  use cactus
  use constants_init
  use cctk_Parameter
  implicit none
  private
  public make_frw_lin 
  
contains
  
  subroutine make_frw_lin(CCTK_ARGUMENTS, xx, time, lx, ly, lz, g, kk, phase, freq_max, i,j,k)
  	DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    INTEGER, intent(in) :: i,j,k, freq_max
    CCTK_REAL plus, plusdot, cross, crossdot, HH_geo, phi0(48), conv, fct, rand, rand1, rand2
    CCTK_REAL wcx, wcy, wcz, w0, pert(3,3), hhdot(3,3), Amp0,rho_d0,rad0, temp0 
    CCTK_REAL bb3(3), eps, beta(3), shift2, rho_c, rho_0, BB_0, pp_0, rho_tot, golden
    CCTK_REAL, intent(in) :: xx(3), time, lx, ly, lz, phase(48)
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer l,m,n,q, ierr, zero, freq_total, seed(8),seed1,seed2
    
    zero = 0
    
    wcx = 2.0*pi/(lx*Meter)
    
    wcy = 2.0*pi/(ly*Meter)
    
    wcz = 2.0*pi/(lz*Meter)
    
    conv = 180/pi
    
    golden = 2.39996
    
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
	temp0 = itemp
	
	 if (rho_0.ne.zero) then
       !eps = (eps_cof*Rad_Const*SIE*(itemp**4))/irho
       eps = eps_cof*Rad_Const*SIE*(itemp**4)/irho - irho_dark/irho - irad/irho - 1.0
    else
       eps = 0.0
    end if
	
	pp_0 = (gam - 1.0)*eps*rho_0
	
	rho_tot = rho_0 + rho_d0
	
	rho_c = rho_0*(1.0 + eps) + BB_0*(1.0-0.5/alpha(i,j,k)**2) + rho_d0 + rad0
	
	HH_geo = sqrt(8.0*pi*Flatness_Ratio*rho_c/3.0)*aa0
	
	Hubble0 = HH_geo
    
    !HH_geo = HH/(3.085678D19)
    !HH_geo = HH/(9.25062992217D27)
    
    freq_total = 0
    
    do n = 1, freq_max
    
    	freq_total = freq_total + n
    
    end do
    
    if (random.ne.0) then
    
    call RANDOM_SEED()
    
    Amp0 = rho_pert*ts_ratio
    
    call RANDOM_NUMBER(rand)
    
    plus = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    plusdot = Amp0*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    
    cross = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    crossdot = Amp0*2.0*(rand - 0.5)
             
    pert(1,1) = pert(1,1) + plus
    pert(2,2) = pert(2,2) - plus
    
    hhdot(1,1) = hhdot(1,1) + plusdot
    hhdot(2,2) = hhdot(2,2) - plusdot
    
    pert(1,2) = pert(1,2) + cross
    pert(2,1) = pert(2,1) + cross
    
    hhdot(1,2) = hhdot(1,2) + crossdot
    hhdot(2,1) = hhdot(2,1) + crossdot
    
    call RANDOM_NUMBER(rand)
    
    plus = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    plusdot = Amp0*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    
    cross = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    crossdot = Amp0*2.0*(rand - 0.5)
    		          
    pert(1,1) = pert(1,1) - plus
    pert(3,3) = pert(3,3) + plus
    
    hhdot(1,1) = hhdot(1,1) - plusdot
    hhdot(3,3) = hhdot(3,3) + plusdot
    
    pert(1,3) = pert(1,3) + cross
    pert(3,1) = pert(3,1) + cross
    
    hhdot(1,3) = hhdot(1,3) + crossdot
    hhdot(3,1) = hhdot(3,1) + crossdot
    		
    call RANDOM_NUMBER(rand)
    
    plus = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    plusdot = Amp0*2.0*(rand - 0.5)
        
    call RANDOM_NUMBER(rand)
    
    cross = Amp0*2.0*(rand - 0.5)
    
    call RANDOM_NUMBER(rand)
    
    crossdot = Amp0*2.0*(rand - 0.5)
    		          
    pert(2,2) = pert(2,2) + plus
    pert(3,3) = pert(3,3) - plus
    
    hhdot(2,2) = hhdot(2,2) + plusdot
    hhdot(3,3) = hhdot(3,3) - plusdot
    
    pert(2,3) = pert(2,3) + cross
    pert(3,2) = pert(3,2) + cross
    
    hhdot(2,3) = hhdot(2,3) + crossdot
    hhdot(3,2) = hhdot(3,2) + crossdot
    
    end if
    
    if (random.eq.0) then
    q = 2
    
    ! for standing wave in the xy plane
    do n = 1, freq_max
    
    	w0 = n*wcz
    	seed1 = INT(xx(1))
    	seed2 = INT(xx(2))
        call RANDOM_SEED(seed1)
        call RANDOM_NUMBER(rand1)
        call RANDOM_SEED(seed2)
        call RANDOM_NUMBER(rand2)
   
        Amp0 = rho_pert*ts_ratio*(w0/HH_geo)**tensor_spec
        
        phi0(43) = (w0*(xx(3)+lz*(phase(43)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(44) = (w0*(xx(3)+lz*(phase(44)+rand1+rand2-1.5)))*Meter + golden*(n-1)

    		plus = Amp0*COS(phi0(43))
    		plusdot = Amp0*w0*SIN(phi0(43))
        
    		cross = Amp0*COS(phi0(44))
    		crossdot = Amp0*w0*SIN(phi0(44))
             
    		pert(1,1) = pert(1,1) + plus
    		pert(2,2) = pert(2,2) - plus
    
    		hhdot(1,1) = hhdot(1,1) + plusdot
    		hhdot(2,2) = hhdot(2,2) - plusdot
    
    		pert(1,2) = pert(1,2) + cross
    		pert(2,1) = pert(2,1) + cross
    
    		hhdot(1,2) = hhdot(1,2) + crossdot
    		hhdot(2,1) = hhdot(2,1) + crossdot
    
    end do
    
    ! for standing wave in the xz plane
    do n = 1, freq_max
    
        w0 = n*wcy
        seed1 = INT(xx(1))
    	seed2 = INT(xx(3))
        call RANDOM_SEED(seed1)
        call RANDOM_NUMBER(rand1)
        call RANDOM_SEED(seed2)
        call RANDOM_NUMBER(rand2)
    
        Amp0 = rho_pert*ts_ratio*(w0/HH_geo)**tensor_spec

        phi0(45) = (w0*(xx(2)+ly*(phase(45)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(46) = (w0*(xx(2)+ly*(phase(46)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    		
    		plus = Amp0*COS(phi0(45))
    		plusdot = Amp0*w0*SIN(phi0(45))
        
    		cross = Amp0*COS(phi0(46))
    		crossdot = Amp0*w0*SIN(phi0(46))
    		          
    		pert(1,1) = pert(1,1) - plus
    		pert(3,3) = pert(3,3) + plus
    
    		hhdot(1,1) = hhdot(1,1) - plusdot
    		hhdot(3,3) = hhdot(3,3) + plusdot
    
    		pert(1,3) = pert(1,3) + cross
    		pert(3,1) = pert(3,1) + cross
    
    		hhdot(1,3) = hhdot(1,3) + crossdot
    		hhdot(3,1) = hhdot(3,1) + crossdot
    
    end do
    
    ! for standing wave in the yz plane
    do n = 1, freq_max
    
        w0 = n*wcx
        seed1 = INT(xx(2))
    	seed2 = INT(xx(3))
        call RANDOM_SEED(seed1)
        call RANDOM_NUMBER(rand1)
        call RANDOM_SEED(seed2)
        call RANDOM_NUMBER(rand2)
    
        Amp0 = rho_pert*ts_ratio*(w0/HH_geo)**tensor_spec

        phi0(47) = (w0*(xx(1)+lx*(phase(47)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(48) = (w0*(xx(1)+lx*(phase(48)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    		
    		plus = Amp0*COS(phi0(47))
    		plusdot = Amp0*w0*SIN(phi0(47))
        
    		cross = Amp0*COS(phi0(48))
    		crossdot = Amp0*w0*SIN(phi0(48))
    		          
    		pert(2,2) = pert(2,2) + plus
    		pert(3,3) = pert(3,3) - plus
    
    		hhdot(2,2) = hhdot(2,2) + plusdot
    		hhdot(3,3) = hhdot(3,3) - plusdot
    
    		pert(2,3) = pert(2,3) + cross
    		pert(3,2) = pert(3,2) + cross
    
    		hhdot(2,3) = hhdot(2,3) + crossdot
    		hhdot(3,2) = hhdot(3,2) + crossdot
        
    end do
    
    end if
    
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
    
  end subroutine make_frw_lin 
  
end module frw_linear 
