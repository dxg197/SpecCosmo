/*@@ Calculates Phase Transitions @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <time.h>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068
#define Speed_Light 2.99792458e8 
#define Meter 3.33564095198e-9
#define Second 1.0

void Transition(CCTK_ARGUMENTS);
void Add_Velocity(CCTK_ARGUMENTS, CCTK_REAL vv_add, CCTK_REAL aa_norm);


void Transition(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
/*  Declare local variables */

    CCTK_INT handle,ierr;
    CCTK_REAL aa_norm, scale_factor, Hubble, time, EW_scale_final, QCD_scale_final;
  	
  	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &scale_factor, 1, CCTK_VarIndex("MHD_Analysis::scale_output"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &Hubble, 1, CCTK_VarIndex("MHD_Analysis::HH_out_avg"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &time, 1, CCTK_VarIndex("MHD_Analysis::univ_age"));
	
	Hubble = Hubble/3.085678e19; 
	
	EW_scale_final = EW_scale*pow((time + trans_time_EW/Hubble)/time,0.5);
	
	QCD_scale_final = QCD_scale*pow((time + trans_time_QCD/Hubble)/time,0.5);
	
	aa_norm = scale_factor/aa;
	
	if (aa_norm == 0.0) {
		aa_norm = 1.0;
	}

if ((scale_factor >= EW_scale) && ((scale_factor <= EW_scale_final) || (CCTK_Equals(Biermann_mode,"none"))) && (time < 1.0)) {

	ierr = CCTK_ParameterSet("Biermann_mode", "SpecGRMHD", "EW");
	
	Add_Velocity(CCTK_PASS_CTOC, maxvel_EW, aa_norm);
	
	printf("EW Transition \n");

}
	
if ((scale_factor >= QCD_scale) && ((scale_factor <= QCD_scale_final) || (CCTK_Equals(Biermann_mode,"EW"))) && (time < 1.0)) {

	ierr = CCTK_ParameterSet("Biermann_mode", "SpecGRMHD", "QCD");
	
	ierr = CCTK_ParameterSet("eps_cof", "mhd_init", "30.875");
	
	Add_Velocity(CCTK_PASS_CTOC, maxvel_QCD, aa_norm);
	
	printf("QCD Transition \n");

}

if ((scale_factor > Nuc_scale) && (CCTK_Equals(Biermann_mode,"QCD")) && (time < 1000.0)) {
	
	ierr = CCTK_ParameterSet("Biermann_mode", "SpecGRMHD", "Nuc");
	
	//ierr = CCTK_ParameterSet("gam", "SpecGRMHD", "1.66666666666666666666667");
	
	ierr = CCTK_ParameterSet("eps_cof", "mhd_init", "1.68");
	
	printf("Primoridal Nucleosynthesis \n");

}
      
}

void Add_Velocity(CCTK_ARGUMENTS, CCTK_REAL vv_add, CCTK_REAL aa_norm)
{
	DECLARE_CCTK_ARGUMENTS;
  	DECLARE_CCTK_PARAMETERS;
  
/*  Declare local variables */

  	CCTK_INT i,j,k,m,n,p,four,handle,ierr,kx;
  	CCTK_INT istart,jstart,kstart,iend,jend,kend;
  	CCTK_INT sizex, sizey, sizez,index;
  	CCTK_INT npoints,nghost,freq_max, seed, seed1, seed2;
  	CCTK_REAL xmin, xmax, ymin, ymax, zmin, zmax;
  	CCTK_REAL lx,ly,lz,wcx,wcy,wcz,vel,phi0[10];
  	CCTK_REAL maxk3,uu4[4],uu3[4],igg[4][4],gg[4][4];
  	CCTK_REAL lapse, shift[4],conf,uu2;
  	CCTK_REAL beta[4],omega,shift2,w0,phase[10];
  	CCTK_REAL ghostx,ghosty,ghostz,vv3[4],rand1,rand2;
  	CCTK_REAL Tmunu[4][4], SS[4], HH_geo, hubble1, golden;
  	 	
  	istart = 0;
	jstart = 0;
	kstart = 0;	
	iend   = cctk_lsh[0];
	jend   = cctk_lsh[1];
	kend   = cctk_lsh[2];
	
	/*handle = CCTK_ReductionHandle("maximum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmax, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymax, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmax, 1, CCTK_VarIndex("grid::z"));
	
	handle = CCTK_ReductionHandle("minimum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmin, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymin, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmin, 1, CCTK_VarIndex("grid::z"));*/
	
	handle = CCTK_ReductionHandle("average");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &hubble1, 1, CCTK_VarIndex("MHD_Analysis::HH_out_avg"));

    //xmin = xmin0 + 0.5*CCTK_DELTA_SPACE(0)*(2*cctk_nghostzones[0]-1);
    //xmax = xmax0 - 0.5*CCTK_DELTA_SPACE(0)*(2*cctk_nghostzones[0]-1);
    lx = xmax0 - xmin0; 
    
	//ymin = ymin0 + 0.5*CCTK_DELTA_SPACE(1)*(2*cctk_nghostzones[1]-1);
    //ymax = ymax0 - 0.5*CCTK_DELTA_SPACE(1)*(2*cctk_nghostzones[1]-1);
    ly = ymax0 - ymin0; 
    
	//zmin = zmin0 + 0.5*CCTK_DELTA_SPACE(2)*(2*cctk_nghostzones[2]-1);
    //zmax = zmax0 - 0.5*CCTK_DELTA_SPACE(2)*(2*cctk_nghostzones[2]-1);
    lz = zmax0 - zmin0;  

    wcx = 2.0*pi/(lx*Meter);
    wcy = 2.0*pi/(ly*Meter);
    wcz = 2.0*pi/(lz*Meter);
    
    sizex = cctk_gsh[0];
    sizey = cctk_gsh[1];
    sizez = cctk_gsh[2];
    
    ghostx = cctk_nghostzones[0];
    ghosty = cctk_nghostzones[1];
    ghostz = cctk_nghostzones[2];
    
    maxk3 = 0.0;
    
    npoints = fmax(sizex,fmax(sizey,sizez));
    
    nghost = fmin(ghostx,fmin(ghosty,ghostz));
    
    HH_geo = hubble1/(Second*3.086e19);
  
    freq_max = (npoints - 2.0*nghost)/4;
    
    golden = 2.39996;
    
    for (kx = 1; kx <= freq_max; kx++) {
    	maxk3 += sqrt(2.0)*pow((kx*(wcx+wcy+wcz)/(3.0*HH_geo)),velocity_spec);
    }
    
    vel = vv_add/maxk3;
    
   for(m = 0; m < 10; m++) {
   	phi0[m] = 0.0;
   }
   
   //srand(time(NULL));
   if (CCTK_Equals(Biermann_mode,"none")) {
   		srand(1); //changed from 1
   }
   if (CCTK_Equals(Biermann_mode,"EW")) {
   		srand(1111); //changed from 11
   }
   if (CCTK_Equals(Biermann_mode,"QCD")) {
   		srand(292929); //changed from 29
   }
   if (CCTK_Equals(Biermann_mode,"Nuc")) {
   		srand(37373737); //changed from 37
   }
   
   for(i=0; i<10; i++) {
    	phase[i] = (rand() % 1000000000)/1000000000.0 - 0.5;
    	//printf("phase = %e \n",phase[i]);
   }
   
   if (vel != 0.0 ){
   
   for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
			
			for(m=0; m < 4; m++) {
  	 			uu4[m] = 0.0;
  	 			uu3[m] = 0.0;
  	 			SS[m] = 0.0;
			}
			
			lapse = alpha[index];
			
	if (random) {
		srand(time(NULL));
		uu4[1] = 2.0*vel*((rand() % 1000000000)/1000000000.0 - 0.5);
    	uu4[2] = 2.0*vel*((rand() % 1000000000)/1000000000.0 - 0.5);
    	uu4[3] = 2.0*vel*((rand() % 1000000000)/1000000000.0 - 0.5);
	}
    
   
    	seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/3.0))%1000; 
    	seed2 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/3.0))%1000; 
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((0.99e19)*pow(x[index]/(lx*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((0.99e19)*pow(y[index]/(ly*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed1) > 499) {seed1 = (int)pow(-1.0,seed1)*seed1;}
    	if (abs(seed2) > 499) {seed2 = (int)pow(-1.0,seed2)*seed2;}
    	srand(seed1);
    	rand1 = (rand() % 10000000000)/10000000000.0 - 0.5;
    	srand(seed2);
    	rand2 = (rand() % 10000000000)/10000000000.0 - 0.5;
    
   for(n = 1; n <= freq_max; n++) { 
    
    	w0 = n*wcz;
    	phi0[1] = (w0*(z[index]+lz*(phase[1]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[2] = (w0*(z[index]+lz*(phase[2]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[3] = (w0*(z[index]+lz*(phase[3]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	
    	if (!random) { 
    		uu4[1] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[1]);
    		uu4[2] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[2]); 
    		uu4[3] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[3]);
    	}
    
    }
    
    //printf("aa_norm = %e, x = %e, y = %e, z = %e, seed1 = %d, seed2 = %d, rand1 = %e, rand2 = %e, phase[1] = %e, phi0[1] = %e \n",aa_norm, x[index], y[index], z[index], seed1,seed2,rand1,rand2,phase[1],phi0[1]);
    
    	seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/3.0))%1000; 
    	seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/3.0))%1000; 
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(x[index]/(lx*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((0.99e19)*pow(x[index]/(lx*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((0.99e19)*pow(z[index]/(lz*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed1) > 499) {seed1 = (int)pow(-1.0,seed1)*seed1;}
    	if (abs(seed2) > 499) {seed2 = (int)pow(-1.0,seed2)*seed2;}
    	srand(seed1);
    	rand1 = (rand() % 10000000000)/10000000000.0 - 0.5;
    	srand(seed2);
    	rand2 = (rand() % 10000000000)/10000000000.0 - 0.5;
    	
    for(n = 1; n <= freq_max; n++) {  
    
    	w0 = n*wcy;
    	phi0[4] = (w0*(y[index]+ly*(phase[4]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[5] = (w0*(y[index]+ly*(phase[5]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[6] = (w0*(y[index]+ly*(phase[6]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	
    	if (!random) {
    		uu4[1] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[4]);
    		uu4[2] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[5]); 
    		uu4[3] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[6]);
    	}
    
    }
    
    	seed1 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/3.0))%1000; 
    	seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/3.0))%1000; 
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/5.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((1.0e19)*pow(y[index]/(ly*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((1.0e19)*pow(z[index]/(lz*aa_norm),1.0/7.0))%1000; }
    	if (abs(seed1) < 1) { seed1 = (long long int)((0.99e19)*pow(y[index]/(ly*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed2) < 1) { seed2 = (long long int)((0.99e19)*pow(z[index]/(lz*aa_norm),1.0/9.0))%1000; }
    	if (abs(seed1) > 499) {seed1 = (int)pow(-1.0,seed1)*seed1;}
    	if (abs(seed2) > 499) {seed2 = (int)pow(-1.0,seed2)*seed2;}
    	srand(seed1);
    	rand1 = (rand() % 10000000000)/10000000000.0 - 0.5;
    	srand(seed2);
    	rand2 = (rand() % 10000000000)/10000000000.0 - 0.5;
    	
    for(n = 1; n <= freq_max; n++) { 
    
    	w0 = n*wcx;
    	phi0[7] = (w0*(x[index]+lx*(phase[7]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[8] = (w0*(x[index]+lx*(phase[8]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	phi0[9] = (w0*(x[index]+lx*(phase[9]-0.5+rand1+rand2)))*Meter*aa_norm + golden*(n-1);
    	
    	if (!random) {
    		uu4[1] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[7]);
    		uu4[2] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[8]); 
    		uu4[3] += vel*(pow(w0/HH_geo,velocity_spec))*cos(phi0[9]);
    	}
    
    }
    	
        shift[1] = betax[index];
        shift[2] = betay[index];
        shift[3] = betaz[index];
			
		conf = exp(4.0*phi[index]);
  	        
        gg[0][0] = gtt[index];
        gg[0][1] = gtx[index];
        gg[0][2] = gty[index];
        gg[0][3] = gtz[index];
        gg[1][0] = gg[0][1];
        gg[2][0] = gg[0][2];
        gg[3][0] = gg[0][3];
        gg[1][1] = gxx[index];
        gg[2][2] = gyy[index];
        gg[3][3] = gzz[index];
        gg[1][2] = gxy[index];
        gg[1][3] = gxz[index];
        gg[2][3] = gyz[index];
        gg[2][1] = gxy[index];
        gg[3][1] = gxz[index];
        gg[3][2] = gyz[index];
        
        igg[0][0] = igtt[index];
        igg[0][1] = igtx[index];
        igg[0][2] = igty[index];
        igg[0][3] = igtz[index];
        igg[1][0] = igg[0][1];
        igg[2][0] = igg[0][2];
        igg[3][0] = igg[0][3];
        igg[1][1] = igxx[index];
        igg[2][2] = igyy[index];
        igg[3][3] = igzz[index];
        igg[1][2] = igxy[index];
        igg[1][3] = igxz[index];
        igg[2][3] = igyz[index];
        igg[2][1] = igg[1][2];
        igg[3][1] = igg[1][3];
        igg[3][2] = igg[2][3];
        
        uu3[1] = uu_x[index];
        uu3[2] = uu_y[index];
        uu3[3] = uu_z[index];
        
        
        for(m=1;m<4;m++) {
        	uu4[m] += uu3[m]; 
        }
        
    	uu2 = 0.0;
        
    	for(m=1;m<4;m++) {
      		for(n=1;n<4;n++) {
        		uu2  += gg[m][n]*uu4[m]*uu4[n];
      		}
    	}
    
    uu4[0] = sqrt(1.0 + uu2)/alpha[index];
    	
    	uu_t[index] = uu4[0];
    	uu_x[index] = uu4[1];
    	uu_y[index] = uu4[2];
    	uu_z[index] = uu4[3];
	
		beta[0] = beta_t[index];
		beta[1] = beta_x[index];
		beta[2] = beta_y[index];
		beta[3] = beta_z[index];
		
	for(m=1; m < 4; m++) {
  	 	vv3[m] = uu4[m]/(uu4[0]*alpha[index]) - shift[m];
	}
	
	vv_x[index] = vv3[1];
	vv_y[index] = vv3[2];
	vv_z[index] = vv3[3];
	
	
	rho_star[index] = rho[index]*uu4[0]*sqrtdetg[index];
	rho_star_p[index] = 0.5*(rho_star[index]+rho_star_p_p[index]);
	
    omega = rho[index] + pp[index] + rho[index]*epsilon[index] + trQ[index];
		
/* Calculate the Tmunu for MHD */
  		
  		Ttt_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_t[index] + (pp[index]+0.5*beta2[index])*igg[0][0] - 
  						 beta[0]*beta[0];
  		
  		Ttx_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_x[index] + (pp[index]+0.5*beta2[index])*igg[0][1] - 
  						 beta[0]*beta[1];
  		Tty_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_y[index] + (pp[index]+0.5*beta2[index])*igg[0][2] - 
  						 beta[0]*beta[2];
  		Ttz_mhd[index] = (omega+beta2[index])*uu_t[index]*uu_z[index] + (pp[index]+0.5*beta2[index])*igg[0][3] - 
  						 beta[0]*beta[3];
  		
        Txx_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_x[index] + (pp[index]+0.5*beta2[index])*igg[1][1] - 
                         beta[1]*beta[1] + Qxx[index];
  		Tyy_mhd[index] = (omega+beta2[index])*uu_y[index]*uu_y[index] + (pp[index]+0.5*beta2[index])*igg[2][2] - 
  		                 beta[2]*beta[2] + Qyy[index];
  		Tzz_mhd[index] = (omega+beta2[index])*uu_z[index]*uu_z[index] + (pp[index]+0.5*beta2[index])*igg[3][3] - 
  		                 beta[3]*beta[3] + Qzz[index];
  		Txy_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_y[index] + (pp[index]+0.5*beta2[index])*igg[1][2] -  
  		                 beta[1]*beta[2] + Qxy[index];
  		Txz_mhd[index] = (omega+beta2[index])*uu_x[index]*uu_z[index] + (pp[index]+0.5*beta2[index])*igg[1][3] - 
  		                 beta[1]*beta[3] + Qxz[index];
  		Tyz_mhd[index] = (omega+beta2[index])*uu_y[index]*uu_z[index] + (pp[index]+0.5*beta2[index])*igg[2][3] - 
  		                 beta[2]*beta[3] + Qyz[index];
  		
  	if (rho[index] <= 0.0) {
  	
    	Ttt_mhd[index] = beta2[index]*uu_t[index]*uu_t[index] - beta[0]*beta[0] + 0.5*beta2[index]*igg[0][0];
    	Ttx_mhd[index] = beta2[index]*uu_t[index]*uu_x[index] - beta[0]*beta[1] + 0.5*beta2[index]*igg[0][1];
  		Tty_mhd[index] = beta2[index]*uu_t[index]*uu_y[index] - beta[0]*beta[2] + 0.5*beta2[index]*igg[0][2];
  		Ttz_mhd[index] = beta2[index]*uu_t[index]*uu_z[index] - beta[0]*beta[3] + 0.5*beta2[index]*igg[0][3];
    	Txx_mhd[index] = beta2[index]*uu_x[index]*uu_x[index] - beta[1]*beta[1] + 0.5*beta2[index]*igg[1][1];
  		Tyy_mhd[index] = beta2[index]*uu_y[index]*uu_y[index] - beta[2]*beta[2] + 0.5*beta2[index]*igg[2][2];
  		Tzz_mhd[index] = beta2[index]*uu_z[index]*uu_z[index] - beta[3]*beta[3] + 0.5*beta2[index]*igg[3][3];
  		Txy_mhd[index] = beta2[index]*uu_x[index]*uu_y[index] - beta[1]*beta[2] + 0.5*beta2[index]*igg[1][2];
  		Txz_mhd[index] = beta2[index]*uu_x[index]*uu_z[index] - beta[1]*beta[3] + 0.5*beta2[index]*igg[1][3];
  		Tyz_mhd[index] = beta2[index]*uu_y[index]*uu_z[index] - beta[2]*beta[3] + 0.5*beta2[index]*igg[2][3];
  		
   }
   
   		Tmunu[0][0] = Ttt_mhd[index];
   		Tmunu[0][1] = Ttx_mhd[index];
   		Tmunu[0][2] = Tty_mhd[index];
   		Tmunu[0][3] = Ttz_mhd[index];
   		Tmunu[1][1] = Txx_mhd[index];
   		Tmunu[2][2] = Tyy_mhd[index];
   		Tmunu[3][3] = Tzz_mhd[index];
   		Tmunu[1][2] = Txy_mhd[index];
   		Tmunu[1][3] = Txz_mhd[index];
   		Tmunu[2][3] = Tyz_mhd[index];
   		
   		Tmunu[1][0] = Tmunu[0][1];
   		Tmunu[2][0] = Tmunu[0][2];
   		Tmunu[3][0] = Tmunu[0][3];
   		Tmunu[2][1] = Tmunu[1][2];
   		Tmunu[3][1] = Tmunu[1][3];
   		Tmunu[3][2] = Tmunu[2][3];
   
    for(m=1;m<4;m++) {
      for(n=1;n<4;n++) {
      	SS[m] += sqrtdetg[index]*gg[m][n]*Tmunu[0][n];
      }
    }
  		
        S_x[index] = SS[1];
  		S_y[index] = SS[2];
  		S_z[index] = SS[3];
  		
    	S_x_p[index] = 0.5*(S_x[index]+S_x_p_p[index]);
  		S_y_p[index] = 0.5*(S_y[index]+S_y_p_p[index]);
  		S_z_p[index] = 0.5*(S_z[index]+S_z_p_p[index]);
  		
        tau[index] = alpha[index]*sqrtdetg[index]*Tmunu[0][0]-rho_star[index];
        tau_p[index] = 0.5*(tau[index]+tau_p_p[index]);
   
   
   			}
   
   		}
   
   } 
   
   }

}
    
#undef pi 
#undef Speed_Light 
#undef Meter
#undef Second
