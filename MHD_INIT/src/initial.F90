/*@@ Calculates the Initial Conditions for the GRMHD equations @@*/
 
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine grmhd_initial(CCTK_ARGUMENTS)

  ! Create the initial data, and the initial gauge.
  ! Initial data and initial gauge are created everywhere, even
  ! on the boundary.  The normal boundary conditions and gauge
  ! conditions are not applied afterwards.  The result is stored
  ! on the current time level.
  
  use bondi
  use cactus
  use constants_init
  use coords
  use empty
  use flatspace
  use geodesic
  use gowdy 
  use gowdylapse
  use gaugewave
  use kasner
  use homo_mhd
  use turb_mhd
  use pointwise_init
  use standingwave
  use standingwave3D
  use frw_linear
  use frw_rot
  use noise
  use alfven
  use collision
  use compound
  use fast_rarefaction
  use fast_shock
  use shock_tube1
  use shock_tube2
  use slow_rarefaction
  use slow_shock
  use periodic_shock
  use spec_init
  use frw
  use scale_gauge
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  !
  integer, parameter :: init_space_none       = 1
  integer, parameter :: init_space_flatspace  = 2
  integer, parameter :: init_space_bondi      = 3
  integer, parameter :: init_space_swave      = 4
  integer, parameter :: init_space_frw_linear = 5
  integer, parameter :: init_space_frw_rot    = 6
  integer, parameter :: init_space_kasner     = 7
  integer, parameter :: init_space_gowdy      = 8
  integer, parameter :: init_space_gaugewave  = 9
  integer, parameter :: init_space_swave3D    = 10
  integer, parameter :: init_space_frw        = 11
  !
  integer, parameter :: init_gauge_none       = 1
  integer, parameter :: init_gauge_geodesic   = 2
  integer, parameter :: init_gauge_gowdy      = 3 
  integer, parameter :: init_gauge_gwave      = 4
  integer, parameter :: init_gauge_scale      = 5
  !
  !
  integer, parameter :: init_mhd_none         = 1
  integer, parameter :: init_mhd_homo         = 2
  integer, parameter :: init_mhd_turb         = 3
  integer, parameter :: init_mhd_empty        = 4
  integer, parameter :: init_mhd_alfven       = 5
  integer, parameter :: init_mhd_collision    = 6
  integer, parameter :: init_mhd_compound     = 7
  integer, parameter :: init_mhd_fast_rare    = 8
  integer, parameter :: init_mhd_fast_shock   = 9
  integer, parameter :: init_mhd_shock1       = 10
  integer, parameter :: init_mhd_shock2       = 11
  integer, parameter :: init_mhd_slow_rare    = 12
  integer, parameter :: init_mhd_slow_shock   = 13
  integer, parameter :: init_mhd_periodic     = 14
  integer, parameter :: init_mhd_spec_init    = 15
  !
  CCTK_REAL xx(3)
  CCTK_REAL offset(3), radius_max
  CCTK_REAL xmax, xmin, lx, ymax, ymin, ly, zmax, zmin, lz, phase(48)
  CCTK_REAL g(3,3), kk(3,3),alph, beta(3), zero, Tmunudark(0:3,0:3), volume
  CCTK_REAL rhostar0, rhodstar0, tau0, s0(3), bb0(3), vv0(3), Tmunumhd(0:3,0:3)
  CCTK_REAL sizex,sizey,sizez,ghostx,ghosty,ghostz
  !
  integer init_space, init_gauge, init_mhd, seed(48)
  integer i,j,k,q, istart, jstart, kstart, vindex, handle
  integer shape(3), pos, ierr, npoints, nghost, freq_max
  !
  !
  istart = 1 
  jstart = 1 
  kstart = 1 
  q = 48
  univ_age = itime
  univ_age_p = itime
  !	
  
  !call CCTK_VarIndex(vindex, "grid::x")
  !call CCTK_ReductionHandle(handle, "maximum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xmax, 1, vindex)
  !call CCTK_ReductionHandle(handle, "minimum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, xmin, 1, vindex)
  !call CCTK_VarIndex(vindex, "grid::y")
  !call CCTK_ReductionHandle(handle, "maximum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, ymax, 1, vindex)
  !call CCTK_ReductionHandle(handle, "minimum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, ymin, 1, vindex)
  !call CCTK_VarIndex(vindex, "grid::z")
  !call CCTK_ReductionHandle(handle, "maximum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmax, 1, vindex)
  !call CCTK_ReductionHandle(handle, "minimum")
  !call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmin, 1, vindex)
  
  call cctk_coordrange(ierr, cctkgh, xmin, xmax, -1, "x", "cart3d")
  xmin = xmin0 + 0.5*cctk_delta_space(1)*(2*cctk_nghostzones(1)-1)
  xmax = xmax0 - 0.5*cctk_delta_space(1)*(2*cctk_nghostzones(1)-1)
  lx = xmax0 - xmin0 
  call cctk_coordrange(ierr, cctkgh, ymin, ymax, -1, "y", "cart3d")
  ymin = ymin0 + 0.5*cctk_delta_space(2)*(2*cctk_nghostzones(2)-1)
  ymax = ymax0 - 0.5*cctk_delta_space(2)*(2*cctk_nghostzones(2)-1)
  ly = ymax0 - ymin0 
  call cctk_coordrange(ierr, cctkgh, zmin, zmax, -1, "z", "cart3d")
  zmin = zmin0 + 0.5*cctk_delta_space(3)*(2*cctk_nghostzones(3)-1)
  zmax = zmax0 - 0.5*cctk_delta_space(3)*(2*cctk_nghostzones(3)-1)
  lz = zmax0 - zmin0 
  !
  shape(:) = cctk_lsh(:)
  zero = 0.0
  
  IF (ABS(xmin)==xmax .AND. ABS(ymin)==ymax .AND. ABS(zmin)==zmax) THEN
  radius_max = SQRT(xmax*xmax+ymax*ymax+zmax*zmax)
  ELSE
  radius_max = 0.0
  END IF
  
  volume = lx*ly*lz
  
  sizex = cctk_gsh(1)
  sizey = cctk_gsh(2)
  sizez = cctk_gsh(3)
  
  ghostx = cctk_nghostzones(1)
  ghosty = cctk_nghostzones(2)
  ghostz = cctk_nghostzones(3)
  
  npoints = MAX(sizex, sizey, sizez)
  
  nghost = MIN(ghostx, ghosty, ghostz)
  
  freq_max = (npoints - 2*nghost)/4
  !
  if (CCTK_EQUALS(initial_gauge, "none")) then
     init_gauge = init_gauge_none
  else if (CCTK_EQUALS(initial_gauge, "geodesic")) then
     init_gauge = init_gauge_geodesic
  else if (CCTK_EQUALS(initial_gauge, "gowdy")) then
     init_gauge = init_gauge_gowdy
  else if (CCTK_EQUALS(initial_gauge, "gwave")) then
     init_gauge = init_gauge_gwave
  else if (CCTK_EQUALS(initial_gauge, "scale")) then
     init_gauge = init_gauge_scale
  else
     call CCTK_WARN(0, "Illegal initial gauge selected")
  end if
  !
  !
  if (init_gauge /= init_gauge_none) then
     do k = kstart, cctk_lsh(3)
        do j = jstart, cctk_lsh(2)
           do i = istart, cctk_lsh(1)
               !
               ! get coordinates
               call getcoords(cctk_time, x,y,z, xx, i,j,k)
               call calc_position(shape,i,j,k,pos)
               !
               ! create initial gauge
               select case (init_gauge)
               case (init_gauge_geodesic)
                  call makegeodesic(alph, beta)
               case (init_gauge_gowdy)
                  call glapse(cctk_time, xx(3), alph, beta, lz)
               case (init_gauge_gwave)
                  call makegwavegauge(xx, cctk_time, lz, alph, beta)
               case (init_gauge_scale)
                  call makescale(alph, beta)
               case default
                  call CCTK_WARN(0, "internal error")
               end select
               !
               call getoffset(x,y,z, offset, i,j,k)
               beta = beta + offset
               !
               call makegaugenoise(alph, beta)
               !
               call set_scalar(alph, alpha, pos)
               call set_vector(beta, betax, betay, betaz, pos)
               !
            end do
         end do
      end do
  end if
  !
  if (CCTK_EQUALS(initial_data, "none")) then
     init_space = init_space_none
  else if (CCTK_EQUALS(initial_data, "flatspace")) then
     init_space = init_space_flatspace
  else if (CCTK_EQUALS(initial_data, "bondi")) then
     init_space = init_space_bondi
  else if (CCTK_EQUALS(initial_data, "swave")) then
     init_space = init_space_swave
  else if (CCTK_EQUALS(initial_data, "swave3D")) then
     init_space = init_space_swave3D
  else if (CCTK_EQUALS(initial_data, "frw_linear")) then
     init_space = init_space_frw_linear
  else if (CCTK_EQUALS(initial_data, "frw_rot")) then
     init_space = init_space_frw_rot
  else if (CCTK_EQUALS(initial_data, "gowdy")) then
     init_space = init_space_gowdy
  else if (CCTK_EQUALS(initial_data, "kasner")) then
     init_space = init_space_kasner  
  else if (CCTK_EQUALS(initial_data, "gaugewave")) then
     init_space = init_space_gaugewave 
  else if (CCTK_EQUALS(initial_data, "frw")) then
     init_space = init_space_frw
  else
     call CCTK_WARN(0, "Illegal initial spacetime selected")
  end if
  !
  do i = 1, 48
  	seed(i) = i**3
  end do
  
  call RANDOM_SEED(SIZE=q)
  call RANDOM_SEED(PUT=seed)
  call RANDOM_NUMBER(HARVEST=phase)
  
  !do i = 1, 48
    !print *, 'phase is', phase(i)
  !end do
  !
  if (init_space /= init_space_none) then
     do k = kstart, cctk_lsh(3)
        do j = jstart, cctk_lsh(2)
           do i = istart, cctk_lsh(1)
              !
              ! get coordinates
              call getcoords(cctk_time, x,y,z, xx, i,j,k)
              call calc_position(shape,i,j,k,pos)
              !
              ! create initial space-time
              select case (init_space)
              case (init_space_flatspace)
                 call makeflatspace(g, kk)
              case (init_space_bondi)
                 call make_bondi(xx, cctk_time, lz, g, kk)
              case (init_space_gowdy)
                 call makegowdy(cctk_time, xx(3), g, kk, lz)
              case (init_space_kasner)
                 call makekasner(cctk_time, g, kk)
              case (init_space_gaugewave)
                 call makegaugewave(xx, cctk_time, lz, g, kk)
              case (init_space_swave)
                 call make_swave(xx, cctk_time, lz, g, kk)
              case (init_space_swave3D)
                 call make_swave3D(CCTK_PASS_FTOF, xx, itime, lx, ly, lz, g, kk)
              case (init_space_frw_linear)
                 call make_frw_lin(CCTK_PASS_FTOF, xx, itime, lx, ly, lz, g, kk, phase, freq_max, i,j,k)
              case (init_space_frw_rot)
                 call make_frw_rot(CCTK_PASS_FTOF, xx, itime, lx, ly, lz, g, kk, phase, freq_max, i,j,k)
              case (init_space_frw)
                 call make_frw(CCTK_PASS_FTOF, g, kk, i,j,k)
              case default
                 call CCTK_WARN(0, "internal error")
              end select
              !
              call makenoise(g, kk)
              !
              call set_tensor(g, gxx, gxy, gxz, gyy, gyz, gzz, pos)
              call set_tensor(kk/alpha(i,j,k), kxx, kxy, kxz, kyy, kyz, kzz, pos)
              !
           end do
        end do
     end do
  end if
  !
  if (CCTK_EQUALS(initial_mhd, "none")) then
     init_mhd = init_mhd_none
  else if (CCTK_EQUALS(initial_mhd, "homo")) then
     init_mhd = init_mhd_homo
  else if (CCTK_EQUALS(initial_mhd, "turb")) then
     init_mhd = init_mhd_turb
  else if (CCTK_EQUALS(initial_mhd, "empty")) then
     init_mhd = init_mhd_empty
  else if (CCTK_EQUALS(initial_mhd, "alfven")) then
     init_mhd = init_mhd_alfven
  else if (CCTK_EQUALS(initial_mhd, "collision")) then
     init_mhd = init_mhd_collision
  else if (CCTK_EQUALS(initial_mhd, "compound")) then
     init_mhd = init_mhd_compound
  else if (CCTK_EQUALS(initial_mhd, "fast_rare")) then
     init_mhd = init_mhd_fast_rare
  else if (CCTK_EQUALS(initial_mhd, "fast_shock")) then
     init_mhd = init_mhd_fast_shock
  else if (CCTK_EQUALS(initial_mhd, "shock1")) then
     init_mhd = init_mhd_shock1
  else if (CCTK_EQUALS(initial_mhd, "shock2")) then
     init_mhd = init_mhd_shock2
  else if (CCTK_EQUALS(initial_mhd, "slow_rare")) then
     init_mhd = init_mhd_slow_rare
  else if (CCTK_EQUALS(initial_mhd, "slow_shock")) then
     init_mhd = init_mhd_slow_shock
  else if (CCTK_EQUALS(initial_mhd, "periodic")) then
     init_mhd = init_mhd_periodic
  else if (CCTK_EQUALS(initial_mhd, "spec")) then
     init_mhd = init_mhd_spec_init
  else
     call CCTK_WARN(0, "Illegal initial Plasma State selected")
  end if
  !
  !
  if (init_mhd /= init_mhd_none) then
     do k = kstart, cctk_lsh(3)
        do j = jstart, cctk_lsh(2)
           do i = istart, cctk_lsh(1)
               !
               ! get coordinates
               call getcoords(cctk_time, x,y,z, xx, i,j,k)
               call calc_position(shape,i,j,k,pos)
               !
               ! create initial matter field
               !
               select case (init_mhd)
               case (init_mhd_homo)
                  call makehomomhd(CCTK_PASS_FTOF,i,j,k,volume)
               case (init_mhd_turb)
                  call maketurbmhd(CCTK_PASS_FTOF,volume,radius_max,i,j,k)
               case (init_mhd_empty)
                  call makeempty(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_alfven)
                  call makealfven(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_collision)
                  call makecollision(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_compound)
                  call makecompound(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_fast_rare)
                  call makefastrarefaction(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_fast_shock)
                  call makefastshock(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_shock1)
                  call makeshocktube1(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_shock2)
                  call makeshocktube2(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_slow_rare)
                  call makeslowrarefaction(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_slow_shock)
                  call makeslowshock(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_periodic)
                  call makeperiodicshock(CCTK_PASS_FTOF,i,j,k)
               case (init_mhd_spec_init)
                  call make_spec_init(CCTK_PASS_FTOF, xx, itime, phase, freq_max, lx, ly, lz, i,j,k)
               case default
                  call CCTK_WARN(0, "internal error")
               end select
               !
               ! Initial guess for transition scale factors
    			scale_T1 = 6.97e11
    			scale_T2 = 1.38e15
               !
               call set_scalar(zero, psib0, pos)
               call set_scalar(zero, Qxx, pos)
               call set_scalar(zero, Qyy, pos)
               call set_scalar(zero, Qzz, pos)
               call set_scalar(zero, Qxy, pos)
               call set_scalar(zero, Qxz, pos)
               call set_scalar(zero, Qyz, pos)
               call set_scalar(zero, rho_star_dt, pos)
			   call set_scalar(zero, tau_dt, pos)
	           call set_scalar(zero, S_x_dt, pos)
	           call set_scalar(zero, S_y_dt, pos)
	           call set_scalar(zero, S_z_dt, pos)
	           call set_scalar(zero, bb_x_dt, pos)
	           call set_scalar(zero, bb_y_dt, pos)
	           call set_scalar(zero, bb_z_dt, pos)
	           call set_scalar(zero, psib_dt, pos)
	           call set_scalar(zero, rho_star_dark_dt, pos)
	           call set_scalar(zero, tau_dark_dt, pos)
	           call set_scalar(zero, SD_x_dt, pos)
	           call set_scalar(zero, SD_y_dt, pos)
	           call set_scalar(zero, SD_z_dt, pos)
	           call set_scalar(zero, alpha_dt, pos)
               call set_scalar(zero, betax_dt, pos)
               call set_scalar(zero, betay_dt, pos)
               call set_scalar(zero, betaz_dt, pos)
               call set_scalar(zero, bssn_bb_x_dt, pos)
               call set_scalar(zero, bssn_bb_y_dt, pos)
               call set_scalar(zero, bssn_bb_z_dt, pos)
	           call set_scalar(zero, phi_dt, pos)
	           call set_scalar(zero, bssn_gxx_dt, pos)
	           call set_scalar(zero, bssn_gyy_dt, pos)
	           call set_scalar(zero, bssn_gzz_dt, pos)
	           call set_scalar(zero, bssn_gxy_dt, pos)
	           call set_scalar(zero, bssn_gxz_dt, pos)
	           call set_scalar(zero, bssn_gyz_dt, pos)
	           call set_scalar(zero, bssn_Axx_dt, pos)
	           call set_scalar(zero, bssn_Ayy_dt, pos)
	           call set_scalar(zero, bssn_Azz_dt, pos)
	           call set_scalar(zero, bssn_Axy_dt, pos)
	           call set_scalar(zero, bssn_Axz_dt, pos)
	           call set_scalar(zero, bssn_Ayz_dt, pos)
	           call set_scalar(zero, bssn_K_dt, pos)
	           call set_scalar(zero, bssn_gamma_x_dt, pos)
	           call set_scalar(zero, bssn_gamma_y_dt, pos)
	           call set_scalar(zero, bssn_gamma_z_dt, pos)
               !
            end do
         end do
      end do
  end if
  !
  !
end subroutine grmhd_initial