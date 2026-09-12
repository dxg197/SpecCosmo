! MHD_Init
! Switch-on Slow Rarefaction

! $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module slow_rarefaction
  use constants_init
  use pointwise_init
  use matdet_init
  use matinv_init
  use tensor_init
  implicit none
  private
  public makeslowrarefaction
  !
contains
  !
  subroutine makeslowrarefaction(CCTK_ARGUMENTS,i,j,k)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS 
    !    
    INTEGER, intent(in) :: i,j,k
    INTEGER :: shape(3), pos, m, n
    CCTK_REAL :: rhostar0, rhodstar0, tau0, s0(3), sd0(3), bb0(3), vv0(3), Tmunu(0:3,0:3)
    CCTK_REAL :: sqrtg, detg, g(3,3), g_inv(0:3,0:3), gg(0:3,0:3), zero
    CCTK_REAL :: uu4(0:3), bb3(3), rho0, rho_d0, eps, pp0, tau_dark0, Tdark(0:3,0:3)
    CCTK_REAL :: shift(3), hp, beta(0:3), uu2
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
 	  zero = 0.0
 	  uu2 = 0.0
 	  sd0(:) = 0.0
 		                
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
    
    if (z(i,j,k).eq.0) then
    
    ! Calc variables in geometerized units
    
    uu4(3) = maxvel*(-0.3825)
    uu4(2) = maxvel*(-0.693)
    uu4(1) = 0.0
    
    do m = 1,3
      do n = 1,3
        uu2  = uu2 + g_inv(m,n)*uu4(m)*uu4(n)
      end do
    end do
    
    uu4(0) = sqrt(1.0 + uu2)/alpha(i,j,k)
    bb3(3) = 1.0*SIE
    bb3(2) = 0.511*SIE
    bb3(1) = 0.0
    rho0   = 5.89e-3*SIE
    pp0 = 0.55*SIE
    
    end if
    
    if (z(i,j,k).lt.0) then
    
    ! Calc variables in geometerized units
    
    uu4(3) = maxvel*(-0.765)
    uu4(2) = maxvel*(-1.386)
    uu4(1) = 0.0
    
    do m = 1,3
      do n = 1,3
        uu2  = uu2 + g_inv(m,n)*uu4(m)*uu4(n)
      end do
    end do
    
    uu4(0) = sqrt(1.0 + uu2)/alpha(i,j,k)
    bb3(3) = 1.0*SIE
    bb3(2) = 1.022*SIE
    bb3(1) = 0.0
    rho0   = 1.78e-3*SIE
    pp0 = 0.1*SIE
    
    end if
    
    if (z(i,j,k).gt.0) then
    
    ! Calc variables in geometerized units
    
    uu4(3) = 0.0
    uu4(2) = 0.0
    uu4(1) = 0.0
    
    do m = 1,3
      do n = 1,3
        uu2  = uu2 + g_inv(m,n)*uu4(m)*uu4(n)
      end do
    end do
    
    uu4(0) = sqrt(1.0 + uu2)/alpha(i,j,k)
    bb3(3) = 1.0*SIE
    bb3(2) = 0.0
    bb3(1) = 0.0
    rho0   = 0.01*SIE
    pp0 = 1.0*SIE
    
    end if
    
    rho_d0 = irho_dark*Density
    eps = pp0/((gam-1.0)*rho0)
    
    hp = 1.0 + gam*eps 

    do m = 1,3
      do n = 1,3
  	    beta(0) = beta(0) + gg(m,n)*uu4(m)*bb3(n)/alpha(i,j,k)
      end do
    end do
    
	do m = 1,3	
  	 beta(m) = (bb3(m)/alpha(i,j,k) + beta(0)*uu4(m))/uu4(0)
	end do
 
	do m = 0,3
	 do n = 0,3  		
	  beta2(i,j,k) = beta2(i,j,k) + gg(m,n)*beta(m)*beta(n)
	 end do
	end do

    ! Calc Initial gridfunctions
        
    rhostar0 = sqrtg*rho0*uu4(0)
    rhodstar0 = sqrtg*rho_d0*uu4(0)
       
    do m = 1,3   
     bb0(m) = sqrt(detg)*bb3(m)
     vv0(m) = uu4(m)/uu4(0)-shift(m)
    end do
    
    do m = 0,3
     do n = 0,3
      Tmunu(m,n) = (rho0*hp+beta2(i,j,k))*uu4(m)*uu4(n)+(pp0+0.5*beta2(i,j,k))*g_inv(m,n)-beta(m)*beta(n)
      Tdark(m,n) = rho_d0*uu4(m)*uu4(n)
     end do
    end do
    
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
      call set_scalar(uu4(1), uu_x, pos)
      call set_scalar(uu4(2), uu_y, pos)
      call set_scalar(uu4(3), uu_z, pos)
      call set_vector(vv0, vv_x, vv_y, vv_z, pos)
      call set_scalar(rho0, rho, pos)
      call set_scalar(pp0, pp,pos)
      call set_vector(bb0*sqrt(4.0*pi), bb_x, bb_y, bb_z, pos)
      call set_vector(bb3*sqrt(4.0*pi), b_x, b_y, b_z, pos)
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
      call set_scalar(beta(0), beta_t, pos)
      call set_scalar(beta(1), beta_x, pos)
      call set_scalar(beta(2), beta_y, pos)
      call set_scalar(beta(3), beta_z, pos)
      
      call set_scalar(rhodstar0, rho_star_dark, pos)
      call set_scalar(uu4(0), uud_t, pos)
      call set_vector(vv0, vvd_x, vvd_y, vvd_z, pos)
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
    
  end subroutine makeslowrarefaction
  !
end module slow_rarefaction
