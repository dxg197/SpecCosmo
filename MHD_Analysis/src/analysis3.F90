/*@@ Calculates the output data for Plasma Field @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine analysis3(CCTK_ARGUMENTS)
    use analytic
    use cactus_analysis
    use coords_analysis
    use constants_analysis
    use geodesic_analysis
    use geod_analysis
    use gowdylapse_analysis
    use gaugewave_analysis
    use noise_analysis
    use pointwise_analysis
    use scale_analysis
    use singular
    implicit none
	
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    
    integer, parameter :: calc_gauge_none        = 1
    integer, parameter :: calc_gauge_geodesic    = 2
    integer, parameter :: calc_gauge_gowdy       = 3
    integer, parameter :: calc_gauge_gwave       = 4
    integer, parameter :: calc_gauge_geod        = 5
    integer, parameter :: calc_gauge_scale       = 6
    
    CCTK_REAL xx(3)
    CCTK_REAL offset(3), delta_univ_age, HH_avg, HH_avg_p
    CCTK_REAL zero, zmax, zmin, lz, alph, beta(3)
    CCTK_REAL rho_od, b_x_od, b_y_od, b_z_od, temp_of, rho_otd, rho_c, rho_ratio, rho_of
    CCTK_REAL temp_od, pp_od, KijKij, vv_x_od, vv_y_od, vv_z_od, rho_otf
    CCTK_REAL gxx_d, gyy_d, gzz_d, gxy_d, gxz_d, gyz_d, h(3,3), rho_dark_od, rho_dark_of
    CCTK_REAL scale, results, igg3(3,3), kk3(3,3), HH_temp, Kijhij, scale_ratio, hijhij
    CCTK_REAL rho_avg, rhod_avg, pp_avg, b_x_avg, b_y_avg, b_z_avg, temp_avg, sum_lapse
    CCTK_REAL rho_avg2, rhod_avg2, rho_out_total_avg2, temp_avg2, sum_a, avg_density
    CCTK_REAL rho_avg3, rhod_avg3, rho_out_total_avg3, temp_avg3, sum_TrK, delta, delta_dot
    CCTK_REAL ienergy_avg, ienergy_avg2, ienergy_avg3, ienergy_d, ienergy_f 
    CCTK_REAL temp_of_rel, temp_od_rel, temp_of_nonrel, temp_od_nonrel
    CCTK_REAL temp_rel_avg, temp_nonrel_avg, temp_rel_avg2, temp_nonrel_avg2, temp_rel_avg3, temp_nonrel_avg3
    
    integer calc_gauge
    integer i,j,k,m,n,o,p
    integer shape(3), pos
    integer istart,jstart,kstart,iend,jend,kend
    integer ierr, handle, vindex
    
    !call cctk_coordrange(ierr, cctkgh, zmin, zmax, -1, "z", "cart3d")
    
    call CCTK_ReductionHandle(handle, "maximum")
    call CCTK_VarIndex(vindex, "grid::z")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmax, 1, vindex)
    
    call CCTK_ReductionHandle(handle, "minimum")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, zmin, 1, vindex)
    
    zmin = zmin + cctk_delta_space(3)*(cctk_nghostzones(3)-0.5)
    zmax = zmax - cctk_delta_space(3)*(cctk_nghostzones(3)-0.5)
    lz = zmax - zmin

    ! Set up shorthands
    ! -----------------

      zero = 0.0 

      istart = 1 
      jstart = 1 
      kstart = 1 

      iend = cctk_lsh(1) 
      jend = cctk_lsh(2) 
      kend = cctk_lsh(3) 
      
      shape(:) = cctk_lsh(:)
      
      n = 0
      m = 0
      sum_a = 0.0
      sum_TrK = 0.0
      sum_lapse = 0.0
      
      DO k = kstart, kend
       DO j = jstart, jend
        DO i = istart, iend
			
		  IF (TrK(i,j,k) < 0.0) THEN
 	      sum_a = sum_a + EXP(2.0*phi(i,j,k))
 	      sum_TrK = sum_TrK + TrK(i,j,k)
 	      n = n + 1
 	      END IF
 	      
 	    END DO
 	   END DO
 	  END DO 
 	       
      TrK_avg = sum_TrK/n
      scale = sum_a/n

      IF (isnan(TrK_avg)) THEN
      TrK_avg = TrK_avg_p
      aa_out_avg = aa_out_avg_p
      scale = aa_out_avg/aa0
      END IF
      
      aa_out_avg = aa0*scale
      
      IF (aa_out_avg_p.ne.0) THEN
      	scale_ratio = aa_out_avg/aa_out_avg_p
      else
      	scale_ratio = 1.0
      END IF
      
      aa_ratio_avg = scale_ratio
      
      scale_output = aa*scale
      
    call CCTK_ReductionHandle(handle, "average")
    results = 0.0

    call CCTK_VarIndex(vindex, "SpecGRMHD::rho")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_avg = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rho_avg failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::rho_dark")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rhod_avg = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rhod_avg failed!");
    ENDIF
   
    call CCTK_VarIndex(vindex, "SpecGRMHD::pp")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    pp_avg = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of pp_avg failed!");
    ENDIF

    call CCTK_VarIndex(vindex, "SpecGRMHD::b_x")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_x_avg = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of b_x_avg failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::b_y")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_y_avg = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of b_y_avg failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "SpecGRMHD::b_z")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    b_z_avg = results
    
    call CCTK_VarIndex(vindex, "MHD_Analysis::ienergy")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    ienergy_avg = results
    
    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of b_z_avg failed!");
    ENDIF

    ! -----------------
    ! Do the 
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend

          ! Get Position and Offsets using TGRTensor
 	      call calc_position(shape,i,j,k,pos)
          
          rho_od = rho(i,j,k) - rho_avg 
          rho_dark_od = rho_dark(i,j,k) - rhod_avg 
          rho_otd = rho_out_total(i,j,k) - rho_out_total_avg 
          temp_od = temp_out(i,j,k) - temp_out_avg 
          temp_od_rel = temp_out_rel(i,j,k) - temp_out_rel_avg 
          temp_od_nonrel = temp_out_nonrel(i,j,k) - temp_out_nonrel_avg 
          pp_od = pp(i,j,k) - pp_avg
          
          b_x_od = b_x(i,j,k) - b_x_avg  
          b_y_od = b_y(i,j,k) - b_y_avg  
          b_z_od = b_z(i,j,k) - b_z_avg  
          
          vv_x_od = Speed_Light*(vv_x(i,j,k) - vv_x_avg)  
          vv_y_od = Speed_Light*(vv_y(i,j,k) - vv_y_avg) 
          vv_z_od = Speed_Light*(vv_z(i,j,k) - vv_z_avg)  

          gxx_d = (gxx(i,j,k) - gxx_avg)/(aa_out_avg**2)
          gyy_d = (gyy(i,j,k) - gyy_avg)/(aa_out_avg**2)
          gzz_d = (gzz(i,j,k) - gzz_avg)/(aa_out_avg**2)
          gxy_d = (gxy(i,j,k) - gxy_avg)/(aa_out_avg**2)
          gxz_d = (gxz(i,j,k) - gxz_avg)/(aa_out_avg**2)
          gyz_d = (gyz(i,j,k) - gyz_avg)/(aa_out_avg**2)
          
          ienergy_d = ienergy(i,j,k) - ienergy_avg

          call set_scalar(rho_od/Density, rho_out_diff, pos)
          call set_scalar(rho_dark_od/Density, rho_dark_diff, pos)
          call set_scalar(rho_otd, rho_out_total_diff, pos)
          call set_scalar(temp_od, temp_out_diff, pos)
          call set_scalar(temp_od_rel, temp_out_rel_diff, pos)
          call set_scalar(temp_od_nonrel, temp_out_nonrel_diff, pos)
          call set_scalar(pp_od/Pressure, pp_out_diff, pos)
          call set_scalar(b_x_od/Gauss, b_x_out_diff, pos)
          call set_scalar(b_y_od/Gauss, b_y_out_diff, pos)
          call set_scalar(b_z_od/Gauss, b_z_out_diff, pos)
          call set_scalar(vv_x_od, vv_x_out_diff, pos)
          call set_scalar(vv_y_od, vv_y_out_diff, pos)
          call set_scalar(vv_z_od, vv_z_out_diff, pos)
              
          call set_scalar(gxx_d, gxx_diff, pos)
          call set_scalar(gyy_d, gyy_diff, pos)
          call set_scalar(gzz_d, gzz_diff, pos)
          call set_scalar(gxy_d, gxy_diff, pos)
          call set_scalar(gxz_d, gxz_diff, pos)
          call set_scalar(gyz_d, gyz_diff, pos)
          
          call set_scalar(ienergy_d, ienergy_diff, pos)
          
          END DO
       END DO
    END DO
    
    call CCTK_ReductionHandle(handle, "average")
    results = 0.0

    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_avg2 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rho_avg2 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_dark_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rhod_avg2 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rhod_avg2 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out_total_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_out_total_avg2 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rho_out_total_avg2 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_avg2 = results
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_rel_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_rel_avg2 = results
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_nonrel_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_nonrel_avg2 = results
    
    call CCTK_VarIndex(vindex, "MHD_Analysis::ienergy_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    ienergy_avg2 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of temp_avg2 failed!");
    ENDIF

    ! -----------------
    ! Do the 
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend

          ! Get Position and Offsets using TGRTensor
 	      call calc_position(shape,i,j,k,pos)

          rho_od = rho_out_diff(i,j,k) - rho_avg2 
          rho_dark_od = rho_dark_diff(i,j,k) - rhod_avg2 
          rho_otd = rho_out_total_diff(i,j,k) - rho_out_total_avg2
          temp_od = temp_out_diff(i,j,k) - temp_avg2  
          temp_od_rel = temp_out_rel_diff(i,j,k) - temp_rel_avg2 
          temp_od_nonrel = temp_out_nonrel_diff(i,j,k) - temp_nonrel_avg2  
          ienergy_d = ienergy_diff(i,j,k) - ienergy_avg2

          call set_scalar(rho_od, rho_out_diff, pos)
          call set_scalar(rho_dark_od, rho_dark_diff, pos)
          call set_scalar(rho_otd, rho_out_total_diff, pos)
          call set_scalar(temp_od, temp_out_diff, pos)
          call set_scalar(temp_od_rel, temp_out_rel_diff, pos)
          call set_scalar(temp_od_nonrel, temp_out_nonrel_diff, pos)
          call set_scalar(ienergy_d, ienergy_diff, pos)
          
          END DO
       END DO
    END DO
    
    call CCTK_ReductionHandle(handle, "average")
    results = 0.0

    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_avg3 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rho_avg3 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_dark_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rhod_avg3 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rhod_avg3 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::rho_out_total_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    rho_out_total_avg3 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of rho_out_total_avg3 failed!");
    ENDIF
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_avg3 = results
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_rel_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_rel_avg3 = results
    
    call CCTK_VarIndex(vindex, "mhd_analysis::temp_out_nonrel_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    temp_nonrel_avg3 = results
    
    call CCTK_VarIndex(vindex, "MHD_Analysis::ienergy_diff")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, results, 1, vindex)
    ienergy_avg3 = results

    IF (ierr.ne.0) THEN
       call CCTK_WARN(1,"Reduction of temp_avg3 failed!");
    ENDIF

    ! -----------------
    ! Do the 
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend

          ! Get Position and Offsets using TGRTensor
 	      call calc_position(shape,i,j,k,pos)

          rho_od = rho_out_diff(i,j,k) - rho_avg3 
          IF (rho_out_avg.ne.0) THEN
        	rho_of = rho_od/rho_out_avg
          END IF
          rho_dark_od = rho_dark_diff(i,j,k) - rhod_avg3 
          IF (rho_dark_avg.ne.0) THEN
        	rho_dark_of = rho_dark_od/rho_dark_avg
          END IF
          rho_otd = rho_out_total_diff(i,j,k) - rho_out_total_avg3
          IF (rho_out_total_avg.ne.0) THEN
        	rho_otf = rho_otd/rho_out_total_avg
          END IF
          temp_od = temp_out_diff(i,j,k) - temp_avg3  
          IF (temp_out_avg.ne.0) THEN
        	temp_of = temp_od/temp_out_avg
          END IF
          temp_od_rel = temp_out_rel_diff(i,j,k) - temp_rel_avg3  
          IF (temp_out_rel_avg.ne.0) THEN
        	temp_of_rel = temp_od_rel/temp_out_rel_avg
          END IF
          temp_od_nonrel = temp_out_nonrel_diff(i,j,k) - temp_nonrel_avg3  
          IF (temp_out_nonrel_avg.ne.0) THEN
        	temp_of_nonrel = temp_od_nonrel/temp_out_nonrel_avg
          END IF
          ienergy_d = ienergy_diff(i,j,k) - ienergy_avg3
          IF (ienergy_avg.ne.0) THEN
          	ienergy_f = ienergy_d/ienergy_avg
          END IF

          call set_scalar(rho_od, rho_out_diff, pos)
          call set_scalar(rho_dark_od, rho_dark_diff, pos)
          call set_scalar(rho_otd, rho_out_total_diff, pos)
          call set_scalar(rho_of, rho_out_frac, pos)
          call set_scalar(rho_dark_of, rho_dark_frac, pos)
          call set_scalar(rho_otf, rho_out_total_frac, pos)
          call set_scalar(temp_od, temp_out_diff, pos)
          call set_scalar(temp_of, temp_out_frac, pos)
          call set_scalar(temp_od_rel, temp_out_rel_diff, pos)
          call set_scalar(temp_of_rel, temp_out_rel_frac, pos)
          call set_scalar(temp_od_nonrel, temp_out_nonrel_diff, pos)
          call set_scalar(temp_of_nonrel, temp_out_nonrel_frac, pos)
          call set_scalar(rho_of-rho_dark_of, rho_diff_frac, pos)
          call set_scalar(ienergy_f, ienergy_frac, pos)
          
          END DO
       END DO
    END DO
    
    IF (calc_analytic.ne.0) THEN
    
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend
          call calc_position(shape,i,j,k,pos)
          call set_scalar(zero, rho_a_diff, pos)
          call set_scalar(zero, b_x_a_diff, pos)
          call set_scalar(zero, b_y_a_diff, pos)
          call set_scalar(zero, b_z_a_diff, pos)
          call set_scalar(zero, vv_x_a_diff, pos)
          call set_scalar(zero, vv_y_a_diff, pos)
          call set_scalar(zero, vv_z_a_diff, pos)
          call set_scalar(zero, rho_error, pos)
          call set_scalar(zero, b_x_error, pos)
          call set_scalar(zero, b_y_error, pos)
          call set_scalar(zero, b_z_error, pos)
          call set_scalar(zero, vv_x_error, pos)
          call set_scalar(zero, vv_y_error, pos)
          call set_scalar(zero, vv_z_error, pos)
          call set_scalar(zero, rho_error_percent, pos)
          call set_scalar(zero, b_x_error_percent, pos)
          call set_scalar(zero, b_y_error_percent, pos)
          call set_scalar(zero, b_z_error_percent, pos)
          call set_scalar(zero, vv_x_error_percent, pos)
          call set_scalar(zero, vv_y_error_percent, pos)
          call set_scalar(zero, vv_z_error_percent, pos)
          END DO
       END DO
    END DO
    
    call analytic_sol(CCTK_PASS_FTOF)
    
    END IF
   
    actual_time = (cctk_time - itime)/Second + itime
    
    time_step = cctk_delta_time
    
    HH_temp = -TrK_avg/3.0
    
    HH_out_avg = 3.085678e19*HH_temp*Second
    
    HH_avg = HH_temp*Second
    
    HH_avg_p = HH_out_avg_p/3.085678e19
    
    delta_univ_age = univ_age - univ_age_p
    
    ! Initial guess for transition scale factors
    ! scale_T1 = 6.97e11
    ! scale_T2 = 1.38e15
    
    ! Calculate Time
    
    IF ((rad_out_avg + rho_out_avg*epsilon(i,j,k)) >= (rho_out_avg + rho_dark_avg)) THEN
    	!univ_age = 1.0/(2.0*HH_temp*Second)
    	univ_age = univ_age_p + 2.0*itime*(scale_ratio-1.0)*(scale)**2.0
    	univ_age_T1 = univ_age
    	scale_T1 = scale
    	!PRINT *, 'Radiation Dominated Universe'
    else IF (((rho_out_avg + rho_dark_avg) > (rad_out_avg + rho_out_avg*epsilon(i,j,k))).and.((rho_out_avg + rho_dark_avg)*Density >= cosmo_constant/(8.0*pi))) THEN
    	!univ_age = 1.0/(1.5*HH_temp*Second)
    	univ_age = univ_age_p + 1.5*univ_age_T1*(scale_ratio-1.0)*(scale/scale_T1)**1.5
    	univ_age_T2 = univ_age
    	scale_T2 = scale
    	!PRINT *, 'Matter Dominated Universe'
    else IF (((rho_out_avg + rho_dark_avg) > (rad_out_avg + rho_out_avg*epsilon(i,j,k))).and.((rho_out_avg + rho_dark_avg)*Density < cosmo_constant/(8.0*pi))) THEN
    	!univ_age = 1.0/(HH_temp*Second)
    	univ_age = univ_age_p + (scale_ratio-1.0)*(scale/scale_T2)*(((HH_avg_p-HH_avg)/(delta_univ_age*HH_avg**2))*LOG(scale/scale_T2) + 1.0/(HH_avg*(scale/scale_T2)))
        !PRINT *, 'Dark Energy Dominated Universe'
    END IF
    
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend
          
          ! Get Position and Offsets using TGRTensor
 	      call calc_position(shape,i,j,k,pos)
   
          gw_energy_density(i,j,k) = 0.0
          KijKij = 0.0
          Kijhij = 0.0
          hijhij = 0.0
          
          igg3(1,1) = igxx(i,j,k)
    	  igg3(2,2) = igyy(i,j,k)
          igg3(3,3) = igzz(i,j,k)
          igg3(1,2) = igxy(i,j,k)
          igg3(1,3) = igxz(i,j,k)
          igg3(2,3) = igyz(i,j,k)
          igg3(2,1) = igxy(i,j,k)
          igg3(3,1) = igxz(i,j,k)
          igg3(3,2) = igyz(i,j,k)
          
          kk3(1,1) = kxx(i,j,k)
    	  kk3(2,2) = kyy(i,j,k)
    	  kk3(3,3) = kzz(i,j,k)
    	  kk3(1,2) = kxy(i,j,k)
    	  kk3(1,3) = kxz(i,j,k)
    	  kk3(2,3) = kyz(i,j,k)
    	  kk3(2,1) = kxy(i,j,k)
    	  kk3(3,1) = kxz(i,j,k)
    	  kk3(3,2) = kyz(i,j,k)
    	  
    	  h(1,1) = gxx_diff(i,j,k)
    	  h(2,2) = gyy_diff(i,j,k)
    	  h(3,3) = gzz_diff(i,j,k)
    	  h(1,2) = gxy_diff(i,j,k)
    	  h(1,3) = gxz_diff(i,j,k)
    	  h(2,3) = gyz_diff(i,j,k)
    	  h(2,1) = gxy_diff(i,j,k)
    	  h(3,1) = gxz_diff(i,j,k)
    	  h(3,2) = gyz_diff(i,j,k)
    	  
    	  rho_c = (TrK_avg**2.0)/(24.0*pi)
    	  
    	  total_rho(i,j,k) = alp_avg**2 * (Ttt_mhd(i,j,k) + Ttt_dark(i,j,k) + Ttt_vac(i,j,k))
    	  
          DO m = 1,3
          	DO n = 1,3
          		DO o = 1,3
          			DO p = 1,3
          				KijKij = KijKij + igg3(m,o)*igg3(n,p)*kk3(m,n)*kk3(o,p)
          				Kijhij = Kijhij + h(m,n)*igg3(m,o)*igg3(n,p)*kk3(o,p)
          				hijhij = hijhij + h(m,n)*igg3(m,o)*igg3(n,p)*h(o,p)
          			END DO
          		END DO
          	END DO
          END DO
          
          IF (HH_temp > 0.0) THEN
          	gw_energy_density(i,j,k) = (alp_avg/(3.0*HH_temp**2.0))*(alp_avg*KijKij - &
          	                            HH_temp*Kijhij*(aa_out_avg**2 + 1.0/aa_out_avg**2) - &
          	                            hijhij*HH_temp**2/alp_avg)
          END IF
          
          IF ((rho_c > 0.0).AND.(TrK(i,j,k) < 0.0)) THEN 
    		rho_ratio_gf(i,j,k) = total_rho(i,j,k)/rho_c
    	  ELSE 
    	  	rho_ratio_gf(i,j,k) = 1.0
    	  END IF
          
          call set_scalar(rho_c, rho_critical, pos)
          	
          END DO
       END DO
    END DO
    
    call CCTK_ReductionHandle(handle, "average")
    call CCTK_VarIndex(vindex, "MHD_Analysis::rho_ratio_gf")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, rho_ratio, 1, vindex)
    call CCTK_ReductionHandle(handle, "norm2")
    call CCTK_VarIndex(vindex, "MHD_Analysis::rho_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, rho_norm, 1, vindex)
    call CCTK_VarIndex(vindex, "MHD_Analysis::rho_dark_out")
    call CCTK_Reduce(ierr, cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, rho_dark_norm, 1, vindex)
    
    DO k = kstart, kend
       DO j = jstart, jend
          DO i = istart, iend
          
    IF ((keep_flat.ne.0).AND.(TrK(i,j,k) < 0.0).AND.(scale_output <= 1.0)) THEN 
          
    IF ((rho_ratio < min_rho_ratio).OR.(rho_ratio > max_rho_ratio)) THEN
    	  	TrK(i,j,k) = sqrt(abs(Flatness_Ratio*rho_ratio))*TrK(i,j,k)
    	  	TrK_p(i,j,k) = 0.5*(TrK(i,j,k) + TrK_p_p(i,j,k))
    END IF
    
    END IF
    
    	END DO
       END DO
    END DO
    
    delta = (rho_norm + rho_dark_norm)/2.80114e-27 - 1.0
    delta_dot = (rho_norm + rho_dark_norm - rho_norm_p - rho_dark_norm_p)/(2.80114e-27 * (univ_age-univ_age_p)) + 6.68231e-18
    
    HH_true = 3.085678e19*(HH_avg + (delta_dot/(3.0*(1.0+delta))))
    
    IF (CCTK_EQUALS(gauge_condition, "none")) THEN
     calc_gauge = calc_gauge_none
    else IF (CCTK_EQUALS(gauge_condition, "geodesic")) THEN
     calc_gauge = calc_gauge_geodesic
    else IF (CCTK_EQUALS(gauge_condition, "gowdy")) THEN
     calc_gauge = calc_gauge_gowdy
    else IF (CCTK_EQUALS(gauge_condition, "gwave")) THEN
     calc_gauge = calc_gauge_gwave
    else IF (CCTK_EQUALS(gauge_condition, "geod")) THEN
     calc_gauge = calc_gauge_geod
    else IF (CCTK_EQUALS(gauge_condition, "scale")) THEN
     calc_gauge = calc_gauge_scale
    else
     call CCTK_WARN(0, "Illegal initial gauge selected")
    END IF
    
    IF (calc_gauge /= calc_gauge_none) THEN
     DO k = 1, cctk_lsh(3)
        DO j = 1, cctk_lsh(2)
            DO i = 1, cctk_lsh(1)
               !
               ! get coordinates
               call getcoords(cctk_time, x,y,z, xx, i,j,k)
               call calc_position(shape,i,j,k,pos)
               !
               ! create gauge
               select case (calc_gauge)
               case (calc_gauge_geodesic)
                  call makegeodesic(alph, beta)
               case (calc_gauge_gowdy)
                  call glapse(cctk_time, xx(3), alph, beta, lz)
               case (calc_gauge_gwave)
                  call makegwavegauge(xx, cctk_time, lz, alph, beta)
               case (calc_gauge_geod)
                  call makegeod(phi(i,j,k),alph,beta)
               case (calc_gauge_scale)
                  call makescale(alph, beta, aa_out_avg)
               case default
                  call CCTK_WARN(0, "internal error")
               END select
               !  
               call getoffset(x,y,z, offset, i,j,k)
               beta = beta + offset
               !
               call makegaugenoise(alph, beta)
               !
               call set_scalar(alph, alpha, pos)
               call set_vector(beta, betax, betay, betaz, pos)
               !
            END DO
         END DO
      END DO
    END IF

END subroutine analysis3
      
