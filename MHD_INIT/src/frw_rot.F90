! MHD_Init
! FRW: Initial Data for FRW spacetime with rotating GW:
! Written by David Garrison
 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

module frw_rot
  use cactus
  use constants_init
  implicit none
  private
  public make_frw_rot
  
contains
  
  subroutine make_frw_rot(CCTK_ARGUMENTS, xx, time, lx, ly, lz, g, kk, phase, freq_max, i,j,k)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    INTEGER, intent(in) :: i,j,k, freq_max
    CCTK_REAL plus, plusdot, cross, crossdot, HH_geo, phi0(48), conv, temp0, rand1, rand2
    CCTK_REAL wcx, wcy, wcz, w0, pert(3,3), hhdot(3,3), fct,rho_d0,rad0, rand 
    CCTK_REAL hhl(3,3), Ql(3,3), hhr(3,3), Qr(3,3), time0, Amp0, rho_tot
    CCTK_REAL bb3(3), eps, beta(3), shift2, rho_c, rho_0, BB_0, pp_0, golden
    CCTK_REAL, intent(in) :: xx(3), time, lx, ly, lz, phase(48) 
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    integer l,m,n,p, zero, seed(8),seed1,seed2
    
    zero = 0
    
    wcx = 2.0*pi/(lx*Meter)
    
    wcy = 2.0*pi/(ly*Meter)
    
    wcz = 2.0*pi/(lz*Meter)
    
    conv = 180/pi
    
    hhl(:,:) = 0.0
    
    Ql(:,:) = 0.0
    
    hhr(:,:) = 0.0
    
    Qr(:,:) = 0.0
    
    golden = 2.39996
    
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
    
    !HH_geo = HH/(3.085678D19)
    !HH_geo = HH/(9.25062992217D27)
    
    time0 = -1e-20 
    p = 2
    
    ! for standing wave in the xy plane
    do n = 1, freq_max
    
    	w0 = n*wcz
    	seed1 = INT(xx(1))
    	seed2 = INT(xx(2))
        call RANDOM_SEED(seed1)
        call RANDOM_NUMBER(rand1)
        call RANDOM_SEED(seed2)
        call RANDOM_NUMBER(rand2)
    
        Amp0 = pAmp*sqrt(3.0*(rho_tot*rho_pert)**2*ts_ratio/(4.0*pi**2*(rho_tot+3.0*pp_0)**2*(w0/HH_geo)**3))
        	
    	phi0(43) = (w0*(xx(3)+lz*(phase(43)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(44) = (w0*(xx(3)+lz*(phase(44)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	    
    		plus = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(43))
    		plusdot = Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(43))
    
    		cross = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(44))
    		crossdot = -Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(44))
    		          
    		hhl(1,1) = hhl(1,1) + plus
    		hhl(2,2) = hhl(2,2) - plus
    
    		Ql(1,1) = Ql(1,1) + plusdot
    		Ql(2,2) = Ql(2,2) - plusdot
    
    		hhl(1,2) = hhl(1,2) + cross
    		hhl(2,1) = hhl(2,1) + cross
    
    		Ql(1,2) = Ql(1,2) + crossdot
    		Ql(2,1) = Ql(2,1) + crossdot
    
    	
    		plus = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(43))
    		plusdot = Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(43))
    
    		cross = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(44))
    		crossdot = -Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(44))
    		
    		hhr(1,1) = hhr(1,1) + plus
    		hhr(2,2) = hhr(2,2) - plus
    
    		Qr(1,1) = Qr(1,1) + plusdot
    		Qr(2,2) = Qr(2,2) - plusdot
    
    		hhr(1,2) = hhr(1,2) - cross
    		hhr(2,1) = hhr(2,1) - cross
    
    		Qr(1,2) = Qr(1,2) - crossdot
    		Qr(2,1) = Qr(2,1) - crossdot
    
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
    
       	Amp0 = pAmp*sqrt(3.0*(rho_tot*rho_pert)**2*ts_ratio/(4.0*pi**2*(rho_tot+3.0*pp_0)**2*(w0/HH_geo)**3))
        
        phi0(45) = (w0*(xx(2)+ly*(phase(45)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(46) = (w0*(xx(2)+ly*(phase(46)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    		plus = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(45))
    		plusdot = Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(45))
    
    		cross = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(46))
    		crossdot = -Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(46))
    		          
    		hhl(1,1) = hhl(1,1) - plus
    		hhl(3,3) = hhl(3,3) + plus
    
    		Ql(1,1) = Ql(1,1) - plusdot
    		Ql(3,3) = Ql(3,3) + plusdot
    
    		hhl(1,3) = hhl(1,3) + cross
    		hhl(3,1) = hhl(3,1) + cross
    
    		Ql(1,3) = Ql(1,3) + crossdot
    		Ql(3,1) = Ql(3,1) + crossdot
    		
    
    		plus = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(45))
    		plusdot = Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(45))
    
    		cross = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(46))
    		crossdot = -Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(46))
    		          
    		hhr(1,1) = hhr(1,1) - plus
    		hhr(3,3) = hhr(3,3) + plus
    
    		Qr(1,1) = Qr(1,1) - plusdot
    		Qr(3,3) = Qr(3,3) + plusdot
    
    		hhr(1,3) = hhr(1,3) - cross
    		hhr(3,1) = hhr(3,1) - cross
    
    		Qr(1,3) = Qr(1,3) - crossdot
    		Qr(3,1) = Qr(3,1) - crossdot
    
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
    
       	Amp0 = pAmp*sqrt(3.0*(rho_tot*rho_pert)**2*ts_ratio/(4.0*pi**2*(rho_tot+3.0*pp_0)**2*(w0/HH_geo)**3))
        
        phi0(47) = (w0*(xx(1)+lx*(phase(47)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(48) = (w0*(xx(1)+lx*(phase(48)+rand1+rand2-1.5)))*Meter + golden*(n-1)

    		plus = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(47))
    		plusdot = Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(47))
    
    		cross = Amp0*EXP(w0*TH*time0*Second)*COS(phi0(48))
    		crossdot = -Amp0*w0*EXP(w0*TH*time0*Second)*SIN(phi0(48))
    		          
    		hhl(2,2) = hhl(2,2) + plus
    		hhl(3,3) = hhl(3,3) - plus
    
    		Ql(2,2) = Ql(2,2) + plusdot
    		Ql(3,3) = Ql(3,3) - plusdot
    
    		hhl(2,3) = hhl(2,3) + cross
    		hhl(3,2) = hhl(3,2) + cross
    
    		Ql(2,3) = Ql(2,3) + crossdot
    		Ql(3,2) = Ql(3,2) + crossdot
    
    	
    		plus = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(47))
    		plusdot = Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(47))
    
    		cross = Amp0*EXP(-w0*TH*time0*Second)*COS(phi0(48))
    		crossdot = -Amp0*w0*EXP(-w0*TH*time0*Second)*SIN(phi0(48))
    		          
    		hhr(2,2) = hhr(2,2) + plus
    		hhr(3,3) = hhr(3,3) - plus
    
    		Qr(2,2) = Qr(2,2) + plusdot
    		Qr(3,3) = Qr(3,3) - plusdot
    
    		hhr(2,3) = hhr(2,3) - cross
    		hhr(3,2) = hhr(3,2) - cross
    
    		Qr(2,3) = Qr(2,3) - crossdot
    		Qr(3,2) = Qr(3,2) - crossdot
    
    end do

    pert(1,1) = (hhl(1,1) + hhr(1,1))/sqrt(2.0)
    pert(2,2) = (hhl(2,2) + hhr(2,2))/sqrt(2.0)
    pert(3,3) = (hhl(3,3) + hhr(3,3))/sqrt(2.0)
    pert(1,2) = (hhl(1,2) + hhr(1,2))/sqrt(2.0)
    pert(1,3) = (hhl(1,3) + hhr(1,3))/sqrt(2.0)
    pert(2,3) = (hhl(2,3) + hhr(2,3))/sqrt(2.0)
    
    hhdot(1,1) = (Ql(1,1) + Qr(1,1))/sqrt(2.0)
    hhdot(2,2) = (Ql(2,2) + Qr(2,2))/sqrt(2.0)
    hhdot(3,3) = (Ql(3,3) + Qr(3,3))/sqrt(2.0)
    hhdot(1,2) = (Ql(1,2) + Qr(1,2))/sqrt(2.0)
    hhdot(1,3) = (Ql(1,3) + Qr(1,3))/sqrt(2.0)
    hhdot(2,3) = (Ql(2,3) + Qr(2,3))/sqrt(2.0)

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
    
  end subroutine make_frw_rot 
  
end module frw_rot 
