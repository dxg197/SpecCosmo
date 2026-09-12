/*@@ Calculates the output data for Plasma Field @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine analysis1(CCTK_ARGUMENTS)
    use cactus_analysis
    use coords_analysis
    use constants_analysis
    use pointwise_analysis
    use temp_analysis
    implicit none
	
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    CCTK_REAL magT_o, divb_o, divb_f, rho_dark_o, volume, EErad_o
    CCTK_REAL rho_o, temp_o, pp_o, b_x_o, b_y_o, b_z_o, bb2_o
    CCTK_REAL xmax, xmin, ymax, ymin, zmax, zmin, int_energy
    CCTK_REAL gg(3,3), b(3), b_l(3), beta(0:3), bu, u(3)
    CCTK_REAL mdensity, edensity, boltz, radc, mass, results
    CCTK_REAL press_mass, density_boltz, mass_p, mass_avg, temp_rel, temp_nonrel
    
    integer i,j,k,m,n
    integer shape(3), pos
    integer istart,jstart,kstart,iend,jend,kend
    integer ierr, handle, vindex

    ! Set up shorthands
    ! -----------------

      istart = 1
      jstart = 1 
      kstart = 1 

      iend = cctk_lsh(1) 
      jend = cctk_lsh(2)
      kend = cctk_lsh(3) 

      shape(:) = cctk_lsh(:)
      
    call CCTK_VarIndex(vindex, "grid::x")
    call CCTK_ReductionHandle(handle, "maximum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xmax, 1, vindex)
    call CCTK_ReductionHandle(handle, "minimum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xmin, 1, vindex)
    call CCTK_VarIndex(vindex, "grid::y")
    call CCTK_ReductionHandle(handle, "maximum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, ymax, 1, vindex)
    call CCTK_ReductionHandle(handle, "minimum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, ymin, 1, vindex)
    call CCTK_VarIndex(vindex, "grid::z")
    call CCTK_ReductionHandle(handle, "maximum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmax, 1, vindex)
    call CCTK_ReductionHandle(handle, "minimum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmin, 1, vindex)
  
      volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
      
    call CCTK_ReductionHandle(handle, "average")
    results = 0.0

    call CCTK_VarIndex(vindex, "SpecGRMHD::rho")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    mass_out_mean = results*volume/Density

    if (ierr.ne.0) then
       call CCTK_WARN(1,"Reduction of rho_out_avg failed!");
    endif

    ! -----------------
    ! Do the 
    do k = kstart, kend
       do j = jstart, jend
          do i = istart, iend
             
          ! Get Position and Offsets using TGRTensor
 	        call calc_position(shape,i,j,k,pos)
 	        
 	   	if (rho_star(i,j,k) < 0.0) then
		    rho_star(i,j,k) = 0.0
		end if
		
		if (rho_star_dark(i,j,k) < 0.0) then
		    rho_star_dark(i,j,k) = 0.0
		end if
 	        
 	      gg(1,1) = gxx(i,j,k)
    	  gg(2,2) = gyy(i,j,k)
          gg(3,3) = gzz(i,j,k)
          gg(1,2) = gxy(i,j,k)
          gg(1,3) = gxz(i,j,k)
          gg(2,3) = gyz(i,j,k)
          gg(2,1) = gxy(i,j,k)
          gg(3,1) = gxz(i,j,k)
          gg(3,2) = gyz(i,j,k)
          
          b(1) = b_x(i,j,k)
          b(2) = b_y(i,j,k)
          b(3) = b_z(i,j,k)
          
          u(1) = uu_x(i,j,k)
          u(2) = uu_y(i,j,k)
          u(3) = uu_z(i,j,k)
          
          b_l(:) = 0.0
          beta(:) = 0.0
          bu = 0.0
 	        
 	do m = 1, 3
       do n = 1, 3	
 	      b_l(m) = b_l(m) + gg(m,n)*b(n)	
 	   end do
 	end do
    
    do m = 1, 3
	 	bu = bu + b_l(m)*u(m)
	end do

  	beta(0) = bu/alpha(i,j,k)

	do m = 1, 3
  	 	beta(m) = (b(m) + alpha(i,j,k)*beta(0)*u(m))/uu_t(i,j,k)
	end do
          
          rho_o = rho(i,j,k)/Density
          
          rho_dark_o = rho_dark(i,j,k)/Density
          
          EErad_o = EErad(i,j,k)/Density
          
          pp_o = pp(i,j,k)/Pressure

 	      b_x_o = beta(1)/Gauss
 	      
 	      b_y_o = beta(2)/Gauss
 	      
 	      b_z_o = beta(3)/Gauss
 	      
 	      bb2_o = (beta(1)**2.0 + beta(2)**2.0 + beta(3)**2.0)/Density
 	      
 	      mass = 9.10938e-31
 	      
 	      mass_p = 1.67262192e-27
 	      
 	      mass_avg = 6.3063784271e-30
 	      
 	      mdensity = (rho_o + rho_dark_o)/SIE
 	      
 	      edensity = (epsilon(i,j,k)*rho_o + 0.5*bb2_o + EErad_o + rho_o + rho_dark_o)/SIE 
 	      
 	      boltz = 1.5*rho_o*Boltz_Const/(eps_cof*mass)
 	      
 	      radc = eps_cof*Rad_Const
 	      
 	      int_energy = epsilon(i,j,k)*rho_o*Speed_Light**2.0 
 	      
 	      press_mass = pp_o*mass_avg
 	      
 	      density_boltz = rho_o*Boltz_Const
 	      
 	      call calc_temp(mdensity, density_boltz, radc, edensity, press_mass, int_energy, temp_o, temp_rel, temp_nonrel)
 	      
 	      
 	      magT_o = Ttt_mhd(i,j,k) + Txx_mhd(i,j,k) + Tyy_mhd(i,j,k) + Tzz_mhd(i,j,k) + &
 	      		   Ttt_dark(i,j,k) + Txx_dark(i,j,k) + Tyy_dark(i,j,k) + Tzz_dark(i,j,k) + &
 	      		   Ttt_vac(i,j,k) + Txx_vac(i,j,k) + Tyy_vac(i,j,k) + Tzz_vac(i,j,k)
 	               
 	      divb_o = divb(i,j,k)/Gauss
 	      
 	      if ((bb_x(i,j,k)**2.0 + bb_y(i,j,k)**2.0 + bb_z(i,j,k)**2.0).ne.0.0) then
 	          divb_f = divb(i,j,k)/(bb_x(i,j,k)**2.0 + bb_y(i,j,k)**2.0 + bb_z(i,j,k)**2.0)**(0.5)
 	      end if
              
              call set_scalar(rho_o, rho_out, pos)
              call set_scalar(rho_dark_o, rho_dark_out, pos)
              call set_scalar((rho_o+rho_dark_o+EErad_o), rho_out_total, pos)
              call set_scalar(EErad_o, rad_out, pos)
              call set_scalar(temp_o, temp_out, pos)
              call set_scalar(temp_rel, temp_out_rel, pos)
              call set_scalar(temp_nonrel, temp_out_nonrel, pos)
              call set_scalar(pp_o, pp_out, pos)
              call set_scalar(b_x_o, b_x_out, pos)
              call set_scalar(b_y_o, b_y_out, pos)
              call set_scalar(b_z_o, b_z_out, pos)
              call set_scalar(divb_o, divb_out, pos) 
              call set_scalar(divb_f, divb_frac, pos) 
              call set_scalar(magT_o, magT, pos) 
              call set_scalar(int_energy, ienergy, pos)
 	      
          end do
       end do
    end do
    
end subroutine analysis1
