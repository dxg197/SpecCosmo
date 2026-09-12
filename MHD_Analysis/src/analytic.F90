/*@@ Calculates the Analytic Solution for Plasma Field @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

module analytic
  implicit none
  private
  public analytic_sol
  
contains

subroutine analytic_sol(CCTK_ARGUMENTS)
    use cactus_analysis
    use constants_analysis
    use pointwise_analysis
    implicit none
	
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    !  Declare local variables
    INTEGER   :: i,j,k, ierr, handle
    INTEGER   :: shape(3), pos, vindex
    INTEGER   :: istart,jstart,kstart,iend,jend,kend
    CCTK_REAL :: zmin, zmax, lz, b2, Sx, Sy, Sz, EE, Cm, wc
    CCTK_REAL :: va_x, va_y, va_z, ua_x, ua_y, ua_z, Cs2, vp_x, vp_y, vp_z
    CCTK_REAL :: omegam1, omegam2, um1_x, um1_y, um1_z, um2_x, um2_y, um2_z
    CCTK_REAL :: Ca, Ck, aaa, a11, a12, a21, a22, c1, c2, am1, am2
    CCTK_REAL :: bb3(3), rho0, eps0, pp0, time, vv_x_a_d, vv_y_a_d, vv_z_a_d
    CCTK_REAL :: rho_a_d, rho_e, b_x_a_d, b_x_e, b_y_a_d, b_y_e, b_z_a_d, b_z_e
    CCTK_REAL :: vv_x_e, vv_y_e, vv_z_e, vv_x_ep, vv_y_ep, vv_z_ep
    CCTK_REAL :: b_x_ep, b_y_ep, b_z_ep, rho_ep, zero

    ! Set up shorthands
    ! -----------------

      istart = 1
      jstart = 1
      kstart = 1
  
      iend = cctk_lsh(1)
      jend = cctk_lsh(2)
      kend = cctk_lsh(3)  
      
      shape(:) = cctk_lsh(:)
      
      zero = 0.0
      
      call CCTK_ReductionHandle(handle, "maximum")
      call CCTK_VarIndex(vindex, "grid::z")
      call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmax, 1, vindex)
    
      call CCTK_ReductionHandle(handle, "minimum")
      call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmin, 1, vindex)
      
      zmin = zmin + (cctk_nghostzones(3)-1)*cctk_delta_space(3)+0.5d0*cctk_delta_space(3)
      zmax = zmax - (cctk_nghostzones(3)-1)*cctk_delta_space(3)-0.5d0*cctk_delta_space(3)
      lz = zmax - zmin 
      
      ! Calc variables in geometerized units
    
        bb3(1) = iBx*Gauss
        bb3(2) = iBy*Gauss
        bb3(3) = iBz*Gauss
        rho0  = irho*Density
        eps0  = eps_cof*Rad_Const*(itemp**4)*SIE/irho                        
        pp0   = (gam-1.0)*eps0*rho0
        
          wc = nw*2.0*pi/lz
          
          b2 = bb3(1)**2 + bb3(2)**2 + bb3(3)**2
 
          Sx = (bb3(1)*iplusA + bb3(2)*icrossB)*bb3(3)/(4.0*pi)
          
          Sy = (bb3(1)*icrossB - bb3(2)*iplusA)*bb3(3)/(4.0*pi)
          
          Sz = (-(bb3(1)**2 - bb3(2)**2)*iplusA - 2.0*bb3(1)*bb3(2)*icrossB)/(4.0*pi)
          
          EE = rho0*(1.0+eps0) + pp0 + b2/(4.0*pi)
          
          Cm = sqrt((gam*pp0 + b2/(4.0*pi))/EE)
          
          va_x = bb3(1)/sqrt(4.0*pi*EE)
          
          va_y = bb3(2)/sqrt(4.0*pi*EE)
          
          va_z = bb3(3)/sqrt(4.0*pi*EE)
          
          vp_x = -va_z*Sz*va_x/(EE*(1.0 - (va_x**2+va_y**2+va_z**2))*(1.0 - va_z**2)) + Sx/(EE*(1.0 - va_z**2))
          
          vp_y = -va_z*Sz*va_y/(EE*(1.0 - (va_x**2+va_y**2+va_z**2))*(1.0 - va_z**2)) + Sy/(EE*(1.0 - va_z**2))
          
          vp_z = Sz/(EE*(1.0 - Cm**2))
          
          ua_x = -wc*va_y
          
          ua_y = wc*va_x
          
          ua_z = 0.0
          
          Cs2 = gam*pp0/(rho0 + pp0 + rho0*eps0)
          
          omegam1 = sqrt(0.5*((wc*Cm)**2 + Cs2*(wc*va_z)**2 + sqrt(((wc*Cm)**2 + Cs2*(wc*va_z)**2)**2 - 4.0*Cs2*(wc**2*va_z)**2)))
          
          omegam2 = sqrt(0.5*((wc*Cm)**2 + Cs2*(wc*va_z)**2 - sqrt(((wc*Cm)**2 + Cs2*(wc*va_z)**2)**2 - 4.0*Cs2*(wc**2*va_z)**2)))
          
          um1_x = va_x
          
          um1_y = va_y
          
          um1_z = va_z + omegam1**2*(1.0 - (va_x**2 + va_y**2 + va_z**2))*wc/((omegam1**2-wc**2)*(wc*va_z))
          
          um2_x = va_x
          
          um2_y = va_y
          
          um2_z = va_z + omegam2**2*(1.0 - (va_x**2 + va_y**2 + va_z**2))*wc/((omegam2**2-wc**2)*(wc*va_z))
          
          Ca = wc*(ua_x*Sx + ua_y*Sy + ua_z*Sz)/EE
          
          Ck = wc*(wc*Sz)/EE
          
          aaa = (Ca - wc*(ua_x*vp_x + ua_y*vp_y + ua_z*vp_z))/((ua_x**2 + ua_y**2 + ua_y**2)*(wc*va_z))
          
          a11 = omegam1*(va_x*um1_x + va_y*um1_y + va_z*um1_z)
          
          a12 = omegam2*(va_x*um2_x + va_y*um2_y + va_z*um2_z)
          
          a21 = omegam1*(wc*um1_z)
          
          a22 = omegam2*(wc*um2_z)
          
          c1 = -wc*(va_x*vp_x + va_y*vp_y + va_z*vp_z)
          
          c2 = Ck - wc*(wc*vp_z)
          
          am1 = (c1*a22 - c2*a12)/(a11*a22 - a12*a21)
          
          am2 = (c2*a11 - c1*a21)/(a11*a22 - a12*a21)

    ! -----------------
    ! Do the 
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
          
    ! Get Position and Offsets using TGRTensor
          call calc_position(shape,i,j,k,pos)

 	  time = cctk_time   
          
          if (iBz.eq.zero) then
          
            rho_a_d = rho0*Sz*(cos(wc*Cm*time) - cos(wc*time))*sin(wc*z(i,j,k))/(EE*(1.0 - Cm**2))
 
            b_x_a_d = bb3(1)*Sz*(cos(wc*Cm*time) - cos(wc*time))*sin(wc*z(i,j,k))/(EE*(1.0 - Cm**2))
 
            b_y_a_d = bb3(2)*Sz*(cos(wc*Cm*time) - cos(wc*time))*sin(wc*z(i,j,k))/(EE*(1.0 - Cm**2))

            b_z_a_d = bb3(3)*Sz*(cos(wc*Cm*time) - cos(wc*time))*sin(wc*z(i,j,k))/(EE*(1.0 - Cm**2))

            vv_x_a_d = zero 

            vv_y_a_d = zero

            vv_z_a_d = -Sz*(Cm*sin(wc*Cm*time) - sin(wc*time))*cos(wc*z(i,j,k))/(EE*(1.0 - Cm**2))
          
          else
          
            rho_a_d = wc*rho0*(am1*um1_z*(1.0 - cos(omegam1*time))/omegam1 + &
                      am2*um2_z*(1.0 - cos(omegam2*time))/omegam2 + &
                      vp_z*(1.0 - cos(wc*time))/wc)*sin(wc*z(i,j,k))
 
            b_x_a_d = wc*(-bb3(3)*aaa*ua_x*(1.0 - cos((wc*va_z)*time))/(wc*va_z) + &
                      am1*(um1_z*bb3(1) - bb3(3)*um1_x)*(1.0 - cos(omegam1*time))/omegam1 + &
                      am2*(um2_z*bb3(1) - bb3(3)*um2_x)*(1.0 - cos(omegam2*time))/omegam2 + &
                      (vp_z*bb3(1) - bb3(3)*vp_x)*(1.0 - cos(wc*time))/wc)*sin(wc*z(i,j,k))
 
            b_y_a_d = wc*(-bb3(3)*aaa*ua_y*(1.0 - cos((wc*va_z)*time))/(wc*va_z) + &
                      am1*(um1_z*bb3(2) - bb3(3)*um1_y)*(1.0 - cos(omegam1*time))/omegam1 + &
                      am2*(um2_z*bb3(2) - bb3(3)*um2_y)*(1.0 - cos(omegam2*time))/omegam2 + &
                      (vp_z*bb3(2) - bb3(3)*vp_y)*(1.0 - cos(wc*time))/wc)*sin(wc*z(i,j,k))

            b_z_a_d = wc*(-bb3(3)*aaa*ua_z*(1.0 - cos((wc*va_z)*time))/(wc*va_z) + &
                      am1*(um1_z*bb3(3) - bb3(3)*um1_z)*(1.0 - cos(omegam1*time))/omegam1 + &
                      am2*(um2_z*bb3(3) - bb3(3)*um2_z)*(1.0 - cos(omegam2*time))/omegam2 + &
                      (vp_z*bb3(3) - bb3(3)*vp_z)*(1.0 - cos(wc*time))/wc)*sin(wc*z(i,j,k))

            vv_x_a_d = vp_x*cos(wc*z(i,j,k))*sin(wc*time) + aaa*ua_x*cos(wc*z(i,j,k))*sin((wc*va_z)*time) + &
                       am1*um1_x*cos(wc*z(i,j,k))*sin(omegam1*time) + am2*um2_x*cos(wc*z(i,j,k))*sin(omegam2*time)

            vv_y_a_d = vp_y*cos(wc*z(i,j,k))*sin(wc*time) + aaa*ua_y*cos(wc*z(i,j,k))*sin((wc*va_z)*time) + &
                       am1*um1_y*cos(wc*z(i,j,k))*sin(omegam1*time) + am2*um2_y*cos(wc*z(i,j,k))*sin(omegam2*time)

            vv_z_a_d = vp_z*cos(wc*z(i,j,k))*sin(wc*time) + aaa*ua_z*cos(wc*z(i,j,k))*sin((wc*va_z)*time) + &
                       am1*um1_z*cos(wc*z(i,j,k))*sin(omegam1*time) + am2*um2_z*cos(wc*z(i,j,k))*sin(omegam2*time)
          
          end if
             
          call set_scalar(rho_a_d/Density, rho_a_diff, pos)
          call set_scalar(b_x_a_d/Gauss, b_x_a_diff, pos)
          call set_scalar(b_y_a_d/Gauss, b_y_a_diff, pos)
          call set_scalar(b_z_a_d/Gauss, b_z_a_diff, pos)
          call set_scalar(vv_x_a_d*Speed_Light, vv_x_a_diff, pos)
          call set_scalar(vv_y_a_d*Speed_Light, vv_y_a_diff, pos)
          call set_scalar(vv_z_a_d*Speed_Light, vv_z_a_diff, pos)
          
          end do
       end do
    end do
    
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
          
          call calc_position(shape,i,j,k,pos)
          
          rho_e = rho_out_diff(i,j,k) - rho_a_diff(i,j,k)
          if (ABS(rho_a_diff(i,j,k)) > 1.0e-14) then 
            rho_ep = rho_e/rho_a_diff(i,j,k)
          else
            rho_ep = zero
          end if
          
          b_x_e = b_x_out_diff(i,j,k) - b_x_a_diff(i,j,k)
          if (ABS(b_x_a_diff(i,j,k)) > 1.0e-14) then 
            b_x_ep = b_x_e/b_x_a_diff(i,j,k)
          else 
            b_x_ep = zero
          end if
          
          b_y_e = b_y_out_diff(i,j,k) - b_y_a_diff(i,j,k)
          if (ABS(b_y_a_diff(i,j,k)) > 1.0e-14) then
            b_y_ep = b_y_e/b_y_a_diff(i,j,k)
          else
            b_y_ep = zero
          end if
          
          b_z_e = b_z_out_diff(i,j,k) - b_z_a_diff(i,j,k)
          if (ABS(b_z_a_diff(i,j,k)) > 1.0e-14) then
            b_z_ep = b_z_e/b_z_a_diff(i,j,k)
          else
            b_z_ep = zero
          end if
          
          vv_x_e = vv_x_out_diff(i,j,k) - vv_x_a_diff(i,j,k)
          if (ABS(vv_x_a_diff(i,j,k)) > 1.0e-14) then
            vv_x_ep = vv_x_e/vv_x_a_diff(i,j,k)
          else
            vv_x_ep = zero
          end if
          
          vv_y_e = vv_y_out_diff(i,j,k) - vv_y_a_diff(i,j,k)
          if (ABS(vv_y_a_diff(i,j,k)) > 1.0e-14) then
            vv_y_ep = vv_y_e/vv_y_a_diff(i,j,k)
          else
            vv_y_ep = zero
          end if
          
          vv_z_e = vv_z_out_diff(i,j,k) - vv_z_a_diff(i,j,k)
          if (ABS(vv_z_a_diff(i,j,k)) > 1.0e-14) then
            vv_z_ep = vv_z_e/vv_z_a_diff(i,j,k)
          else
            vv_z_ep = zero
          end if
          
          call set_scalar(rho_e, rho_error, pos)
          call set_scalar(b_x_e, b_x_error, pos)
          call set_scalar(b_y_e, b_y_error, pos)
          call set_scalar(b_z_e, b_z_error, pos)
          call set_scalar(vv_x_e, vv_x_error, pos)
          call set_scalar(vv_y_e, vv_y_error, pos)
          call set_scalar(vv_z_e, vv_z_error, pos)

          call set_scalar(rho_ep, rho_error_percent, pos)
          call set_scalar(b_x_ep, b_x_error_percent, pos)
          call set_scalar(b_y_ep, b_y_error_percent, pos)
          call set_scalar(b_z_ep, b_z_error_percent, pos)
          call set_scalar(vv_x_ep, vv_x_error_percent, pos)
          call set_scalar(vv_y_ep, vv_y_error_percent, pos)
          call set_scalar(vv_z_ep, vv_z_error_percent, pos)


          end do
       end do
    end do

end subroutine analytic_sol

end module analytic
