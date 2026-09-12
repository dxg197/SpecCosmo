! MHD_Init
! Spectrum of Initial Data for Early Universe
! Written by David Garrison
 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module spec_init
  use constants_init
  use pointwise_init
  use matdet_init
  use matinv_init
  use tensor_init
  use coords
  implicit none
  private
  public make_spec_init 
  !
contains
  !
  subroutine make_spec_init(CCTK_ARGUMENTS, xx, time, phase, freq_max, lx, ly, lz, i,j,k)
  	DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    INTEGER, intent(in) :: i,j,k, freq_max
    INTEGER :: shape(3), pos, m, n, p, fct, kx, ky, kz, seed(14),seed1,seed2
    CCTK_REAL, intent(in) :: xx(3), time, lx, ly, lz, phase(48)
    CCTK_REAL :: rhostar0, rhodstar0, tau0, s0(3), sd0(3), bb0(3), vv0(3), Tmunu(0:3,0:3)
    CCTK_REAL :: sqrtg, detg, g(3,3), g_inv(0:3,0:3), gg(0:3,0:3), zero, rand, rand1, rand2
    CCTK_REAL :: uu4(0:3), bb3(3), rho0, rho_d0, eps, pp0, tau_dark0, Tdark(0:3,0:3)
    CCTK_REAL :: shift(3), hp, beta(0:3), uu2, kappa, rho_c, Tvac(0:3,0:3)
    CCTK_REAL :: mass, vel, bfield, mass_d, maxk1, maxk2, maxk3, temp0, temp_p, mass0
    CCTK_REAL :: wcx, wcy, wcz, w0, conv, phi0(48), one, rad0, dx, dy, dz, HH_geo, golden
 	
 	! Calculate GR variables using TGRTensor
 	
 	shape(:) = cctk_lsh(:)
 	call calc_position(shape,i,j,k,pos)
 
 	call get_vector(betax,betay,betaz,shift,pos)
 	call get_tensor(gxx,gxy,gxz,gyy,gyz,gzz,g,pos)
 	
 	! Initialize temp variables
    ! ------------------------
 	  gg(:,:)  = 0.0
 	  uu4(:) = 0.0
 	  bb3(:) = 0.0
 	  beta(:) = 0.0
 	  beta2(i,j,k) = 0.0
 	  Tmunu(:,:) = 0.0 
 	  Tdark(:,:) = 0.0 
 	  s0(:) = 0.0
 	  sd0(:) = 0.0
 	  zero = 0.0
 	  uu2 = 0.0
 	  one = 1.0
 	  golden = 2.39996
 		                
        do m = 1,3
         do n = 1,3
          gg(m,n)  = g(m,n)
         end do
        end do
        
        call calc_det(g,detg)
        sqrtg = alpha(i,j,k)*sqrt(detg)
        
        gg(1,0) = shift(1)
        gg(2,0) = shift(2)
        gg(3,0) = shift(3)
        gg(0,1) = gg(1,0)
        gg(0,2) = gg(2,0)
        gg(0,3) = gg(3,0)
        do m = 1,3
         do n = 1,3
          gg(0,0) = gg(0,0) + g(m,n)*shift(m)*shift(n)
         end do
        end do
        gg(0,0) = gg(0,0) - alpha(i,j,k)**2
        
        call calc_syminv4(gg,g_inv)   
        
    wcx = 2.0*pi/(lx*Meter)
    
    wcy = 2.0*pi/(ly*Meter)
    
    wcz = 2.0*pi/(lz*Meter)
    
    dx = cctk_delta_space(1)
    
    dy = cctk_delta_space(2)
    
    dz = cctk_delta_space(3)
    
    conv = 180/pi
    
    bb3(1) = iBx*Gauss
    bb3(2) = iBy*Gauss
    bb3(3) = iBz*Gauss
    rho0 = irho*Density
    rho_d0 = irho_dark*Density
    temp0 = itemp
    rad0 = irad*Density
    
    if (random.ne.0) then
    	call RANDOM_SEED()
    	call RANDOM_NUMBER(rand)
		rho0 = (1.0 + rho_pert*2.0*(rand - 0.5))*rho0
		call RANDOM_NUMBER(rand)
    	rho_d0 = (1.0 + rho_dark_pert*2.0*(rand - 0.5))*rho_d0
    	call RANDOM_NUMBER(rand)
    	temp0 = (1.0 + temp_pert*2.0*(rand - 0.5))*temp0
    	call RANDOM_NUMBER(rand)
    	uu4(1) = maxvel*2.0*(rand - 0.5)
    	call RANDOM_NUMBER(rand)
    	uu4(2) = maxvel*2.0*(rand - 0.5)
    	call RANDOM_NUMBER(rand)
    	uu4(3) = maxvel*2.0*(rand - 0.5)
    	call RANDOM_NUMBER(rand)
    	bb3(1) = maxbb*Gauss*2.0*(rand - 0.5)
    	call RANDOM_NUMBER(rand)
    	bb3(2) = maxbb*Gauss*2.0*(rand - 0.5)
    	call RANDOM_NUMBER(rand)
    	bb3(3) = maxbb*Gauss*2.0*(rand - 0.5)
    end if
    
    if (irho.ne.zero) then
       eps = eps_cof*Rad_Const*SIE*(itemp**4)/irho - irho_dark/irho - irad/irho - 1.0
    else
       eps = 0.0
    end if
    
    if (eps < 0.0) then
    	eps = 0.0
    end if
    
	rho_c = rho0*(1.0 + eps) + rho_d0 + rad0
	
	HH_geo = sqrt(8.0*pi*rho_c/3.0)*aa0
	
	maxk1 = 0.0
    maxk2 = 0.0
    maxk3 = 0.0
    
    do kx = 0, freq_max
    	
    	maxk1 = maxk1 + SQRT(2.0)*(kx*(wcx+wcy+wcz)/(3.0*HH_geo))**scalar_spec
    		
    	maxk2 = maxk2 + SQRT(2.0)*(kx*(wcx+wcy+wcz)/(3.0*HH_geo))**vector_spec
    		
    	maxk3 = maxk3 + SQRT(2.0)*(kx*(wcx+wcy+wcz)/(3.0*HH_geo))**velocity_spec
    
    end do
    
    mass = rho_pert*irho*Density/maxk1
    mass_d = rho_dark_pert*irho_dark*Density/maxk1
    temp_p = temp_pert*temp0/maxk1
        
    bfield = maxbb*Gauss/maxk2
    
    vel = maxvel/maxk3
    
    do m = 1, 42
    	phi0(m) = 0
    end do
    p = 14
    
    ! for standing wave in the xy plane
    do n = 1, freq_max
    
    	w0 = n*wcz
    	seed1 = INT(xx(1))
    	seed2 = INT(xx(2))
        call RANDOM_SEED(seed1)
        call RANDOM_NUMBER(rand1)
        call RANDOM_SEED(seed2)
        call RANDOM_NUMBER(rand2)
        
    	phi0(1) = (w0*(xx(3)+lz*(phase(1)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(2) = (w0*(xx(3)+lz*(phase(2)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(3) = (w0*(xx(3)+lz*(phase(3)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(4) = (w0*(xx(3)+lz*(phase(4)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(5) = (w0*(xx(3)+lz*(phase(5)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(6) = (w0*(xx(3)+lz*(phase(6)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(7) = (w0*(xx(3)+lz*(phase(7)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    	phi0(22) = (w0*(xx(3)+lz*(phase(22)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(23) = (w0*(xx(3)+lz*(phase(23)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(24) = (w0*(xx(3)+lz*(phase(24)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(25) = (w0*(xx(3)+lz*(phase(25)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(26) = (w0*(xx(3)+lz*(phase(26)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(27) = (w0*(xx(3)+lz*(phase(27)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(28) = (w0*(xx(3)+lz*(phase(28)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    	!print '(/,3x,"seed1 = ",i10, 3x, "seed1 = ",i10, 3x, "rand1 = ",f7.4, 3x, "rand2 = ",f7.4, 3x, "phi0(1) = ",f7.4)',seed1, seed2, rand1, rand2, phi0(1)
    	
    	if (random.eq.0) then
			rho0 = rho0 + mass*COS(phi0(1))*(w0/HH_geo)**scalar_spec
    		rho_d0 = rho_d0 + mass_d*COS(phi0(2))*(w0/HH_geo)**scalar_spec
    		temp0 = temp0 + temp_p*COS(phi0(3))*(w0/HH_geo)**scalar_spec
    		uu4(1) = uu4(1) + vel*COS(phi0(4))*(w0/HH_geo)**velocity_spec
    		uu4(2) = uu4(2) + vel*COS(phi0(5))*(w0/HH_geo)**velocity_spec 
    		uu4(3) = uu4(3) + vel*COS(phi0(6))*(w0/HH_geo)**velocity_spec
    		bb3(1) = bb3(1) + bfield*COS(phi0(7))*(w0/HH_geo)**vector_spec
    		bb3(2) = bb3(2) + bfield*COS(phi0(22))*(w0/HH_geo)**vector_spec
    		bb3(3) = bb3(3) + bfield*COS(phi0(23))*(w0/HH_geo)**vector_spec
    	end if
    
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
        phi0(8) = (w0*(xx(2)+ly*(phase(8)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(9) = (w0*(xx(2)+ly*(phase(9)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(10) = (w0*(xx(2)+ly*(phase(10)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(11) = (w0*(xx(2)+ly*(phase(11)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(12) = (w0*(xx(2)+ly*(phase(12)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(13) = (w0*(xx(2)+ly*(phase(13)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(14) = (w0*(xx(2)+ly*(phase(14)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    	phi0(29) = (w0*(xx(2)+ly*(phase(29)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(30) = (w0*(xx(2)+ly*(phase(30)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(31) = (w0*(xx(2)+ly*(phase(31)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(32) = (w0*(xx(2)+ly*(phase(32)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(33) = (w0*(xx(2)+ly*(phase(33)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(34) = (w0*(xx(2)+ly*(phase(34)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(35) = (w0*(xx(2)+ly*(phase(35)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    	if (random.eq.0) then
			rho0 = rho0 + mass*COS(phi0(8))*(w0/HH_geo)**scalar_spec
    		rho_d0 = rho_d0 + mass_d*COS(phi0(9))*(w0/HH_geo)**scalar_spec
    		temp0 = temp0 + temp_p*COS(phi0(10))*(w0/HH_geo)**scalar_spec
    		uu4(1) = uu4(1) + vel*COS(phi0(11))*(w0/HH_geo)**velocity_spec
    		uu4(2) = uu4(2) + vel*COS(phi0(12))*(w0/HH_geo)**velocity_spec
    		uu4(3) = uu4(3) + vel*COS(phi0(13))*(w0/HH_geo)**velocity_spec
    		bb3(1) = bb3(1) + bfield*COS(phi0(14))*(w0/HH_geo)**vector_spec
    		bb3(2) = bb3(2) + bfield*COS(phi0(29))*(w0/HH_geo)**vector_spec
    		bb3(3) = bb3(3) + bfield*COS(phi0(30))*(w0/HH_geo)**vector_spec
    	end if
    		
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
        phi0(15) = (w0*(xx(1)+lx*(phase(15)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(16) = (w0*(xx(1)+lx*(phase(16)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(17) = (w0*(xx(1)+lx*(phase(17)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(18) = (w0*(xx(1)+lx*(phase(18)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(19) = (w0*(xx(1)+lx*(phase(19)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(20) = (w0*(xx(1)+lx*(phase(20)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(21) = (w0*(xx(1)+lx*(phase(21)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	
    	phi0(36) = (w0*(xx(1)+lx*(phase(36)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(37) = (w0*(xx(1)+lx*(phase(37)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(38) = (w0*(xx(1)+lx*(phase(38)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(39) = (w0*(xx(1)+lx*(phase(39)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(40) = (w0*(xx(1)+lx*(phase(40)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(41) = (w0*(xx(1)+lx*(phase(41)+rand1+rand2-1.5)))*Meter + golden*(n-1)
    	phi0(42) = (w0*(xx(1)+lx*(phase(42)+rand1+rand2-1.5)))*Meter + golden*(n-1)
        
        if (random.eq.0) then
			rho0 = rho0 + mass*COS(phi0(15))*(w0/HH_geo)**scalar_spec
    		rho_d0 = rho_d0 + mass_d*COS(phi0(16))*(w0/HH_geo)**scalar_spec
    		temp0 = temp0 + temp_p*COS(phi0(17))*(w0/HH_geo)**scalar_spec
    		uu4(1) = uu4(1) + vel*COS(phi0(18))*(w0/HH_geo)**velocity_spec
    		uu4(2) = uu4(2) + vel*COS(phi0(19))*(w0/HH_geo)**velocity_spec
    		uu4(3) = uu4(3) + vel*COS(phi0(20))*(w0/HH_geo)**velocity_spec
    		bb3(1) = bb3(1) + bfield*COS(phi0(21))*(w0/HH_geo)**vector_spec
    		bb3(2) = bb3(2) + bfield*COS(phi0(36))*(w0/HH_geo)**vector_spec
    		bb3(3) = bb3(3) + bfield*COS(phi0(37))*(w0/HH_geo)**vector_spec
    	end if
    		
    end do
    
    do m = 1,3
      do n = 1,3
        uu2  = uu2 + gg(m,n)*uu4(m)*uu4(n)
      end do
    end do
    
    uu4(0) = sqrt(1.0 + uu2)/alpha(i,j,k)
    
    do m = 1,3
      do n = 1,3
  	    beta(0) = beta(0) + gg(m,n)*uu4(m)*bb3(n)/(sqrt(4.0*pi)*alpha(i,j,k))
      end do
    end do
    
	do m = 1,3	
  	 beta(m) = (bb3(m)/alpha(i,j,k) + sqrt(4.0*pi)*beta(0)*uu4(m))/(sqrt(4.0*pi)*uu4(0))
	end do
 
	do m = 0,3
	 do n = 0,3  		
	  beta2(i,j,k) = beta2(i,j,k) + gg(m,n)*beta(m)*beta(n)
	 end do
	end do
    
    if (rho0.ne.0.0) then
       eps = (eps_cof*Rad_Const*SIE*Density*(temp0**4.0) - rho_d0 - rad0)/rho0 - 1.0
    else
       eps = 0.0
    end if
    
    if (eps < 0.0) then
    	eps = 0.0
    end if
    
    pp0 = (gam-1.0)*eps*rho0
    
    if (rho0.ne.zero) then
    	hp = 1.0 + eps + pp0/rho0
    else
    	hp = 1.0 + gam*eps
    end if

    ! Calc Initial gridfunctions
        
    rhostar0 = sqrtg*rho0*uu4(0)
    rhodstar0 = sqrtg*rho_d0
       
    do m = 1,3   
     bb0(m) = sqrt(detg)*bb3(m)
     vv0(m) = uu4(m)/(alpha(i,j,k)*uu4(0))-shift(m)
    end do
    
    do m = 0,3
     do n = 0,3
      Tmunu(m,n) = (rho0*hp+beta2(i,j,k))*uu4(m)*uu4(n)+(pp0+0.5*beta2(i,j,k))*g_inv(m,n)-beta(m)*beta(n)
     end do
    end do
    
    do m = 0,3
     do n = 0,3
    	Tvac(m,n) = cosmo_constant*g_inv(m,n)/(8.0*pi)
   	 end do
    end do
    
    Tvac(0,0) = Tvac(0,0) + rad0;
    Tvac(1,1) = Tvac(1,1) + rad0/3.0;
    Tvac(2,2) = Tvac(2,2) + rad0/3.0;
    Tvac(3,3) = Tvac(3,3) + rad0/3.0;
    
    Tdark(0,0) = rho_d0
    
    do m = 1,3
     do n = 0,3
      s0(m) = s0(m) + sqrtg*gg(m,n)*Tmunu(0,n)
      sd0(m) = sd0(m) + sqrtg*gg(m,n)*Tdark(0,n)
     end do
    end do
    
    tau0 = alpha(i,j,k)*sqrtg*Tmunu(0,0)-rhostar0
    tau_dark0 = alpha(i,j,k)*sqrtg*Tdark(0,0)-rhodstar0
    
    ! Calc other gridfunctions
    
      call set_scalar(rhostar0, rho_star, pos)
      call set_scalar(eps, epsilon, pos)
      call set_scalar(uu4(0), uu_t, pos)
      call set_vector(vv0, vv_x, vv_y, vv_z, pos)
      call set_scalar(rho0, rho, pos)
      call set_scalar(pp0, pp,pos)
      call set_vector(bb0, bb_x, bb_y, bb_z, pos)
      call set_vector(bb3, b_x, b_y, b_z, pos)
      call set_scalar(tau0, tau, pos)
      call set_vector(s0, S_x, S_y, S_z, pos)
      call set_scalar(Tmunu(0,0), Ttt_mhd, pos)
      call set_scalar(Tmunu(0,1), Ttx_mhd, pos)
      call set_scalar(Tmunu(0,2), Tty_mhd, pos)
      call set_scalar(Tmunu(0,3), Ttz_mhd, pos)
      call set_scalar(Tmunu(1,1), Txx_mhd, pos)
      call set_scalar(Tmunu(1,2), Txy_mhd, pos)
      call set_scalar(Tmunu(1,3), Txz_mhd, pos)
      call set_scalar(Tmunu(2,2), Tyy_mhd, pos)
      call set_scalar(Tmunu(2,3), Tyz_mhd, pos)
      call set_scalar(Tmunu(3,3), Tzz_mhd, pos)
      
      call set_scalar(rhodstar0, rho_star_dark, pos)
      call set_scalar(one, uud_t, pos)
      call set_vector(-shift, vvd_x, vvd_y, vvd_z, pos)
      call set_scalar(rho_d0, rho_dark, pos)
      call set_scalar(tau_dark0, tau_dark, pos)
      call set_vector(sd0, SD_x, SD_y, SD_z, pos)
      call set_scalar(Tdark(0,0), Ttt_dark, pos)
      call set_scalar(Tdark(0,1), Ttx_dark, pos)
      call set_scalar(Tdark(0,2), Tty_dark, pos)
      call set_scalar(Tdark(0,3), Ttz_dark, pos)
      call set_scalar(Tdark(1,1), Txx_dark, pos)
      call set_scalar(Tdark(1,2), Txy_dark, pos)
      call set_scalar(Tdark(1,3), Txz_dark, pos)
      call set_scalar(Tdark(2,2), Tyy_dark, pos)
      call set_scalar(Tdark(2,3), Tyz_dark, pos)
      call set_scalar(Tdark(3,3), Tzz_dark, pos)
      
      call set_scalar(Tvac(0,0), Ttt_vac, pos)
      call set_scalar(Tvac(0,1), Ttx_vac, pos)
      call set_scalar(Tvac(0,2), Tty_vac, pos)
      call set_scalar(Tvac(0,3), Ttz_vac, pos)
      call set_scalar(Tvac(1,1), Txx_vac, pos)
      call set_scalar(Tvac(1,2), Txy_vac, pos)
      call set_scalar(Tvac(1,3), Txz_vac, pos)
      call set_scalar(Tvac(2,2), Tyy_vac, pos)
      call set_scalar(Tvac(2,3), Tyz_vac, pos)
      call set_scalar(Tvac(3,3), Tzz_vac, pos)
      call set_scalar(rad0, EErad, pos)
    
  end subroutine make_spec_init 
  
end module spec_init