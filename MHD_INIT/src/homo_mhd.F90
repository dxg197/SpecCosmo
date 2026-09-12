! MHD_Init
! Homogenous MHD Field

! $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module homo_mhd
  use constants_init
  use pointwise_init
  use matdet_init
  use matinv_init
  use tensor_init
  use coords
  implicit none
  private
  public makehomomhd
  !
contains
  !
  subroutine makehomomhd(CCTK_ARGUMENTS,i,j,k,volume)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS 
    !    
    INTEGER, intent(in) :: i,j,k
    INTEGER :: shape(3), pos, m, n
    CCTK_REAL, intent(in) :: volume
    CCTK_REAL :: rhostar0, rhodstar0, tau0, s0(3), sd0(3), bb0(3), vv0(3), Tmunu(0:3,0:3)
    CCTK_REAL :: sqrtg, detg, g(3,3), g_inv(0:3,0:3), gg(0:3,0:3), zero, one
    CCTK_REAL :: uu4(0:3), bb3(3), rho0, rho_d0, eps, pp0, tau_dark0, Tdark(0:3,0:3)
    CCTK_REAL :: shift(3), hp, beta(0:3), uu2, rho_c, Tvac(0:3,0:3), rad0
    !
    ! Get Position using TGRTensor
 
 	shape(:) = cctk_lsh(:)
 	call calc_position(shape,i,j,k,pos)
    
    !
    ! Calculate GR variables using TGRTensor
 
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
        
    ! Calc variables in geometerized units
    
    uu4(0) = 1.0/alpha(i,j,k)
    bb3(1) = iBx*Gauss
    bb3(2) = iBy*Gauss
    bb3(3) = iBz*Gauss
    rho0   = irho*Density
    rho_d0 = irho_dark*Density
    rad0 = irad*Density
    
    beta(0) = 0.0

	do m = 1,3	
  	 beta(m) = bb3(m)/(alpha(i,j,k)*sqrt(4.0*pi)*uu4(0))
	end do
 
	do m = 0,3
	 do n = 0,3  		
	  beta2(i,j,k) = beta2(i,j,k) + gg(m,n)*beta(m)*beta(n)
	 end do
	end do
    
    if (rho0.ne.zero) then
       eps = (eps_cof*Rad_Const*SIE*(itemp**4) - irho_dark)/irho - 1.0
    else
       eps = 0.0
    end if
    
    if ((nonrel.ne.zero).and.(rho0.ne.zero)) then
       	eps = Boltz_Const*SIE*itemp/((gam-1.0)*irho*volume)
    end if
    
    pp0    = (gam-1.0)*eps*rho0
    hp     = 1.0 + gam*eps 

    ! Calc Initial gridfunctions
        
    rhostar0 = sqrtg*rho0*uu4(0)
    rhodstar0 = sqrtg*rho_d0
       
    do m = 1,3   
     bb0(m) = sqrt(detg)*bb3(m)
     vv0(m) = uu4(m)/uu4(0)-shift(m)
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
    
  end subroutine makehomomhd
  !
end module homo_mhd
