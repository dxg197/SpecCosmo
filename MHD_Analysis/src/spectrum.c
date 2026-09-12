/*@@ Calculates Spectrum of Variables @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

void Spectrum(CCTK_ARGUMENTS);
void Spec_Scalar( CCTK_ARGUMENTS, CCTK_REAL **scalar );
void Clean_Read_Scalar( CCTK_ARGUMENTS, CCTK_REAL *scalar );

/* Begin MHD_RHS */

void Spectrum(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
/*  Declare local variables */

  CCTK_INT i,j,k,m,n,p,four;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez, index, ierr, handle1, handle2;
  CCTK_REAL strainx, strainy, strainz, hdiffx, hdiffy, hdiffz;
  CCTK_REAL hsumx, hsumy, hsumz, gwspec_real, gwspec_imag;
  CCTK_REAL hplusx_spec_real, hplusy_spec_real, hplusz_spec_real;
  CCTK_REAL hcrossx_spec_real, hcrossy_spec_real, hcrossz_spec_real;
  CCTK_REAL hplusx_spec_imag, hplusy_spec_imag, hplusz_spec_imag;
  CCTK_REAL hcrossx_spec_imag, hcrossy_spec_imag, hcrossz_spec_imag;
  CCTK_REAL hplusx_spec_PSD, hplusy_spec_PSD, hplusz_spec_PSD;
  CCTK_REAL hcrossx_spec_PSD, hcrossy_spec_PSD, hcrossz_spec_PSD;
  CCTK_REAL rho_spec_real, temp_spec_real, temp1, rhod1, strain1, gwspec1;
  CCTK_REAL rho_spec_imag, temp_spec_imag, hubble1, gg[4][4], vv[4], bb[4];
  CCTK_REAL bb_mag_temp, vv_mag_avg, bb_mag_avg, scale_factor, utime, strain[4], scale;
  CCTK_REAL vv_mag_max, bb_mag_max, delta_rho_max, delta_dark_max, delta_temp_max, rho_avg;
  
  
/* Declare Arrays */

  CCTK_REAL **density,**temperature,**gwspec,**hplusx,**hplusy,**hplusz,
            **hcrossx,**hcrossy,**hcrossz,**vsquare;

/* Set up shorthands */

    sizex  = cctk_lsh[0]; 
    sizey  = cctk_lsh[1];
    sizez  = cctk_lsh[2];
	istart = cctk_nghostzones[0];
	jstart = cctk_nghostzones[1];
	kstart = cctk_nghostzones[2];	
	iend   = cctk_lsh[0] - cctk_nghostzones[0];
	jend   = cctk_lsh[1] - cctk_nghostzones[1];
	kend   = cctk_lsh[2] - cctk_nghostzones[2];
    
    four = 4;
      
    vsquare = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);  
    density = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    temperature = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    gwspec = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hplusx = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hplusy = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hplusz = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hcrossx = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hcrossy = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
    hcrossz = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
                  
    for(m=0; m < 4; m++) {
       vsquare[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       density[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       temperature[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       gwspec[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hplusx[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hplusy[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hplusz[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hcrossx[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hcrossy[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       hcrossz[m] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
    }
    
handle1 = CCTK_ReductionHandle("average");
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &rho_avg, 1, CCTK_VarIndex("MHD_Analysis::rho_out"));

/* Clean Scalars */

if (clean_scalars) {
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_out_diff); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_out_frac); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_dark_diff); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_dark_frac); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_out_total_diff); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, rho_out_total_frac);
	Clean_Read_Scalar(CCTK_PASS_CTOC, temp_out_diff); 
	Clean_Read_Scalar(CCTK_PASS_CTOC, temp_out_frac); 
}

for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
/* Define Arrays*/ 

		gg[1][1] = gxx[index];
        gg[2][2] = gyy[index];
        gg[3][3] = gzz[index];
        gg[1][2] = gxy[index];
        gg[1][3] = gxz[index];
        gg[2][3] = gyz[index];
        gg[2][1] = gxy[index];
        gg[3][1] = gxz[index];
        gg[3][2] = gyz[index];
        
        vv[1] = vv_x[index];
        vv[2] = vv_y[index];
        vv[3] = vv_z[index];
        
        bb[1] = b_x_out[index];
        bb[2] = b_y_out[index];
        bb[3] = b_z_out[index];
    	
    	vsquare[0][index] = 0.0;
    	
    	for(m=1; m < 4; m++) {
           	vsquare[0][index] += vv[m]*vv[m];
        }
    	
    	vv_mag[index] = pow(fabs(vsquare[0][index]),0.5);
        
        bb_mag_temp = 0.0;
        
        for(m=1; m < 4; m++) {
           	bb_mag_temp += bb[m]*bb[m];
        }
        
        bb_mag[index] = pow(fabs(bb_mag_temp),0.5);
        
        density[0][index] = rho_out_frac[index];
        
        temperature[0][index] = temp_out_frac[index];
        
        gwspec[0][index] = gw_energy_density[index];
        
        hplusx[0][index] = 0.5*(gyy_diff[index]-gzz_diff[index]);
        hplusy[0][index] = 0.5*(gzz_diff[index]-gxx_diff[index]);
        hplusz[0][index] = 0.5*(gxx_diff[index]-gyy_diff[index]);
        hcrossx[0][index] = gyz_diff[index];
        hcrossy[0][index] = gxz_diff[index];
        hcrossz[0][index] = gxy_diff[index];
        
        }
	}
}

/* Calc Averages */

ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &vv_mag_avg, 1, CCTK_VarIndex("MHD_Analysis::vv_mag"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &bb_mag_avg, 1, CCTK_VarIndex("MHD_Analysis::bb_mag"));

/* Calc Spectrum */ 

if (calc_spec) {
  
  Spec_Scalar(CCTK_PASS_CTOC, vsquare);
  Spec_Scalar(CCTK_PASS_CTOC, density);
  Spec_Scalar(CCTK_PASS_CTOC, temperature);
  Spec_Scalar(CCTK_PASS_CTOC, gwspec);
  Spec_Scalar(CCTK_PASS_CTOC, hplusx);
  Spec_Scalar(CCTK_PASS_CTOC, hplusy);
  Spec_Scalar(CCTK_PASS_CTOC, hplusz);
  Spec_Scalar(CCTK_PASS_CTOC, hcrossx);
  Spec_Scalar(CCTK_PASS_CTOC, hcrossy);
  Spec_Scalar(CCTK_PASS_CTOC, hcrossz);

}
  
for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
/* Define Arrays */
    	
    	gg[1][1] = gxx[index];
        gg[2][2] = gyy[index];
        gg[3][3] = gzz[index];
        gg[1][2] = gxy[index];
        gg[1][3] = gxz[index];
        gg[2][3] = gyz[index];
        gg[2][1] = gxy[index];
        gg[3][1] = gxz[index];
        gg[3][2] = gyz[index];
        
        rho_spec_real = density[1][index];
        temp_spec_real = temperature[1][index];
        gwspec_real = gwspec[1][index];
        hplusx_spec_real = hplusx[1][index];
        hplusy_spec_real = hplusy[1][index];
        hplusz_spec_real = hplusz[1][index];
        hcrossx_spec_real = hcrossx[1][index];
        hcrossy_spec_real = hcrossy[1][index];
        hcrossz_spec_real = hcrossz[1][index];
        
        rho_spec_imag = density[2][index];
        temp_spec_imag = temperature[2][index];
        gwspec_imag = gwspec[2][index];
        hplusx_spec_imag = hplusx[2][index];
        hplusy_spec_imag = hplusy[2][index];
		hplusz_spec_imag = hplusz[2][index];
        hcrossx_spec_imag = hcrossx[2][index];
        hcrossy_spec_imag = hcrossy[2][index];
        hcrossz_spec_imag = hcrossz[2][index];
        
        vv_squared_PSD[index] = vsquare[3][index];
        rho_spec_PSD[index] = density[3][index];
        temp_spec_PSD[index] = temperature[3][index];
        gwspec_PSD[index] = gwspec[3][index];
        hplusx_spec_PSD = hplusx[3][index];
        hplusy_spec_PSD = hplusy[3][index];
        hplusz_spec_PSD = hplusz[3][index];
        hcrossx_spec_PSD = hcrossx[3][index];
        hcrossy_spec_PSD = hcrossy[3][index];
        hcrossz_spec_PSD = hcrossz[3][index];
        
        vv_mag_diff[index] = vv_mag[index] - vv_mag_avg;
        bb_mag_diff[index] = bb_mag[index] - bb_mag_avg;
            
        hdiffx = 4.0*(hplusx_spec_imag*hcrossx_spec_real-hcrossx_spec_imag*hplusx_spec_real);
        hdiffy = 4.0*(hplusy_spec_imag*hcrossy_spec_real-hcrossy_spec_imag*hplusy_spec_real);
        hdiffz = 4.0*(hplusz_spec_imag*hcrossz_spec_real-hcrossz_spec_imag*hplusz_spec_real);
        
        hsumx = 2.0*(hplusx_spec_PSD + hcrossx_spec_PSD);
        hsumy = 2.0*(hplusy_spec_PSD + hcrossy_spec_PSD);
        hsumz = 2.0*(hplusz_spec_PSD + hcrossz_spec_PSD);
        
        strain[1] = pow(hplusx_spec_PSD + hcrossx_spec_PSD,0.5);
        strain[2] = pow(hplusy_spec_PSD + hcrossy_spec_PSD,0.5);
        strain[3] = pow(hplusz_spec_PSD + hcrossz_spec_PSD,0.5);
        
        geostrain[index] = 0.0;
        
        for(m=1; m < 4; m++) {
           	geostrain[index] += strain[m]*strain[m];
        }
        
        if (hsumx != 0.0) {
        Polarx[index] = hdiffx/hsumx;
        }
        if (hsumy != 0.0) {
        Polary[index] = hdiffy/hsumy;
        }
        if (hsumz != 0.0) {
        Polarz[index] = hdiffz/hsumz;
        }
        
        }
	}
}

ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &scale_factor, 1, CCTK_VarIndex("MHD_Analysis::scale_output"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &scale, 1, CCTK_VarIndex("MHD_Analysis::aa_out_avg"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &utime, 1, CCTK_VarIndex("MHD_Analysis::univ_age"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &hubble1, 1, CCTK_VarIndex("MHD_Analysis::HH_out_avg"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &temp1, 1, CCTK_VarIndex("MHD_Analysis::temp_out_avg"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &rhod1, 1, CCTK_VarIndex("MHD_Analysis::rho_dark_avg"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &strain1, 1, CCTK_VarIndex("MHD_Analysis::strain_avg"));
ierr = CCTK_Reduce(cctkGH, -1, handle1, 1, CCTK_VARIABLE_REAL, &gwspec1, 1, CCTK_VarIndex("MHD_Analysis::gwspec_avg"));

if (((scale_factor > 0.95*data_out_scale) && (scale_factor < 1.05*data_out_scale)) || ((scale_factor > 0.95) && (scale_factor < 1.05))) {

handle2 = CCTK_ReductionHandle("norm2");
ierr = CCTK_Reduce(cctkGH, -1, handle2, 1, CCTK_VARIABLE_REAL, &vv_mag_max, 1, CCTK_VarIndex("MHD_Analysis::vv_mag"));
ierr = CCTK_Reduce(cctkGH, -1, handle2, 1, CCTK_VARIABLE_REAL, &bb_mag_max, 1, CCTK_VarIndex("MHD_Analysis::bb_mag"));
ierr = CCTK_Reduce(cctkGH, -1, handle2, 1, CCTK_VARIABLE_REAL, &delta_rho_max, 1, CCTK_VarIndex("MHD_Analysis::rho_out_frac"));
ierr = CCTK_Reduce(cctkGH, -1, handle2, 1, CCTK_VARIABLE_REAL, &delta_dark_max, 1, CCTK_VarIndex("MHD_Analysis::rho_dark_frac"));
ierr = CCTK_Reduce(cctkGH, -1, handle2, 1, CCTK_VARIABLE_REAL, &delta_temp_max, 1, CCTK_VarIndex("MHD_Analysis::temp_out_frac"));

	printf("\n Scale Factor = %e \n",scale_factor);
	printf(" Universe Age = %e \n",utime);
	printf(" Hubble Parmeter = %e \n",hubble1);
	printf(" L2Norm B-field Mag= %e \n",bb_mag_max);
	printf(" L2Norm Velocity Mag= %e \n",vv_mag_max);
	printf(" L2Norm Density Perturbation = %e \n",delta_rho_max);
	printf(" L2Norm Dark Matter Density Perturbation = %e \n",delta_dark_max);
	printf(" L2Norm Temperature Perturbation = %e \n",delta_temp_max);
	printf(" Average Temperature = %e \n",temp1);
	printf(" Average Density = %e \n",rho_avg);
	printf(" Average Dark Matter Density = %e \n",rhod1);
	printf(" Average GW Strain = %e \n",strain1);
	printf(" Average GW Normalized Energy Density = %e \n \n",gwspec1);
	
}
    
    for(m=0; m < 4; m++) {
      free(vsquare[m]);
      free(density[m]);
      free(temperature[m]);
      free(gwspec[m]);
      free(hplusx[m]);
      free(hplusy[m]);
      free(hplusz[m]);
      free(hcrossx[m]);
      free(hcrossy[m]);
      free(hcrossz[m]);
    }
   
      free(vsquare);
      free(density);
      free(temperature);
      free(gwspec);
      free(hplusx);
      free(hplusy);
      free(hplusz);
      free(hcrossx);
      free(hcrossy);
      free(hcrossz);
      
}
    

