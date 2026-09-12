! MHD_Init
! Empty Space           

! $Header:$

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module empty    
  use constants_init
  use pointwise_init
  use tensor_init
  implicit none
  private
  public makeempty
  !
contains
  !
  subroutine makeempty(CCTK_ARGUMENTS,i,j,k)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS 
    !    
    INTEGER, intent(in) :: i,j,k
    INTEGER :: shape(3), pos, m, n
    CCTK_REAL :: rhostar0, rhodstar0, tau0, s0(3), sd0(3), bb0(3), vv0(3), Tmunu(0:3,0:3)
    CCTK_REAL :: uu4(0:3), bb3(3), rho0, eps, pp0, tau_dark0, Tdark(0:3,0:3)
    CCTK_REAL :: hp, beta(0:3)
    !
    ! Get Position using TGRTensor
 
    shape(:) = cctk_lsh(:)
    call calc_position(shape,i,j,k,pos)
    
    ! Initialize temp variables
    ! ------------------------
      uu4(:) = 0.0
      bb3(:) = 0.0
      beta(:) = 0.0
      beta2(i,j,k) = 0.0
                
    ! Calc variables in geometerized units
    
      uu4(0) = 1.0
      rho0  = 0.0           
      eps   = 0.0 
      pp0   = 0.0
      hp    = 1.0

    ! Calc Initial gridfunctions
        
      rhostar0 = 0.0
      rhodstar0 = 0.0
      bb0(:) = 0.0    
      vv0(:) = 0.0
      Tmunu(:,:) = 0.0
      Tdark(:,:) = 0.0
      tau0 = 0.0 
      s0(:) = 0.0  
      
      tau_dark0 = 0.0
      sd0(:) = 0.0
    
    ! Calc other gridfunctions
    
      call set_scalar(rhostar0, rho_star, pos)
      call set_scalar(eps, epsilon, pos)
      call set_scalar(uu4(0), uu_t, pos)
      call set_vector(vv0, vv_x, vv_y, vv_z, pos)
      call set_scalar(rho0, rho, pos)
      call set_scalar(pp0, pp, pos)
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
      call set_scalar(uu4(0), uud_t, pos)
      call set_vector(vv0, vvd_x, vvd_y, vvd_z, pos)
      call set_scalar(rho0, rho_dark, pos)
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
    
  end subroutine makeempty   
  !
end module empty    
