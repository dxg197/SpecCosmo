#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>

#include "mpi.h"
#include "Slab.h"
#include <fftw3-mpi.h>

#ifdef CCTK_MPI

/* Find out which driver to use */

#  include "cctk_DefineThorn.h"

#  if defined CARPET_CARPET
#    include "Carpet/Carpet/src/carpet_public.h"
#  endif
#  if defined CACTUSPUGH_PUGH
#    include "CactusPUGH/PUGH/src/include/pugh.h"
#  endif



/* Get MPI communicator */

static MPI_Comm get_mpi_comm (cGH const * restrict const cctkGH)
{
  /* CCTK_IsThornActive is an expensive function.  Cache the
     result. */
  static int initialised = 0;
  static int GetMPICommWorld_aliased;
  static int Carpet_active;
  static int PUGH_active;
  if (! initialised) {
    initialised = 1;
    GetMPICommWorld_aliased = CCTK_IsFunctionAliased ("GetMPICommWorld"); 
    Carpet_active =  CCTK_IsThornActive ("Carpet");
    PUGH_active = CCTK_IsThornActive ("PUGH");
  }
  if (GetMPICommWorld_aliased) {
  /*  return * (MPI_Comm const *) GetMPICommWorld (cctkGH); */
  }
#  if defined CARPET_CARPET
  if (Carpet_active) {
    return CarpetMPIComm ();
  }
#  endif
#  if defined CACTUSPUGH_PUGH
  if (PUGH_active) {
    return PUGH_pGH(cctkGH)->PUGH_COMM_WORLD;
  }
#  endif
  return MPI_COMM_WORLD;
}

#endif  /* CCTK_MPI */

void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_Write_Vector_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void SpecDeriv_SpecGalMethod_Vector_1d( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Vector_Derivative_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, iend, jend, kend, index;
	
	iend   = cctk_lsh[0];
	jend   = cctk_lsh[1];
	kend   = cctk_lsh[2];
    
    for(k=0; k < kend; k++)
	{
		for(j=0; j < jend; j++)
		{
			for(i=0; i < iend; i++)
			{
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                temp_phi0x[ index ] = vector[0][1][index];
                temp_phi0y[ index ] = vector[0][2][index];
                temp_phi0z[ index ] = vector[0][3][index];
			}
		}
	}
	SpecDeriv_SpecGalMethod_Vector_1d( CCTK_PASS_CTOC );
	SpecDeriv_Write_Vector_1d( CCTK_PASS_CTOC, vector );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_SpecGalMethod_Vector_1d( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t f_alloc_local, f_local_n0, f_local_0_start;
    ptrdiff_t b_alloc_local, b_local_n0, b_local_0_start;
    ptrdiff_t f_local_n1, f_local_1_end;
    ptrdiff_t b_local_n1, b_local_1_end;
    fftw_plan plan_forwardx;
    fftw_plan plan_forwardy;
    fftw_plan plan_forwardz;
    fftw_plan plan_backward3x;
    fftw_plan plan_backward3y;
    fftw_plan plan_backward3z;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementk;
    fftw_complex *Derivative_SGM3x, *Derivative_SGM3y, *Derivative_SGM3z;
    fftw_complex *phi_inx, *phi_iny, *phi_inz;
    fftw_complex IMAG_NUMBER = _Complex_I;
    fftw_complex ZERO_COMPLEXNO = (0.0,0.0);
	MPI_Comm comm;
  
    sizex  = cctk_gsh[0]-2*cctk_nghostzones[0]; 
    sizey  = cctk_gsh[1]-2*cctk_nghostzones[1];
    sizez  = cctk_gsh[2]-2*cctk_nghostzones[2]; 

    //comm = get_mpi_comm (cctkGH);
    fftw_mpi_init();
    
    f_alloc_local = fftw_mpi_local_size_1d(sizez, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_WISDOM_ONLY | FFTW_PATIENT, 
                    &f_local_n0, &f_local_0_start, &f_local_n1, &f_local_1_end);
    b_alloc_local = fftw_mpi_local_size_1d(sizez, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_WISDOM_ONLY | FFTW_PATIENT, 
                    &b_local_n0, &b_local_0_start, &b_local_n1, &b_local_1_end);
                    
    struct xferinfo info[3];
    
    for (d=0; d<3; ++d) {
      /* Source array: Cactus layout */
      info[d].src.gsh         = cctk_gsh[d];
      info[d].src.lbnd        = cctk_lbnd[d];
      info[d].src.lsh         = cctk_lsh[d];
      info[d].src.ash         = cctk_ash[d];
      info[d].src.lbbox       = cctk_bbox[2*d];
      info[d].src.ubbox       = cctk_bbox[2*d+1];
      info[d].src.nghostzones = cctk_nghostzones[d];
      
      /* Destination array: FFTW layout */
      info[d].dst.gsh         = cctk_gsh[d];
      info[0].dst.lbnd        = cctk_nghostzones[0];
      info[1].dst.lbnd        = cctk_nghostzones[1];
      info[2].dst.lbnd        = f_local_0_start+cctk_nghostzones[2];
      info[0].dst.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].dst.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].dst.lsh         = f_local_n0;
      info[0].dst.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].dst.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].dst.ash         = f_local_n0;
      info[d].dst.lbbox       = cctk_bbox[2*d];
      info[d].dst.ubbox       = cctk_bbox[2*d+1];
      info[d].dst.nghostzones = 0;
      
      /* Source slab: whole array */
      info[d].src.off = 0;
      info[d].src.len = cctk_gsh[d];
      info[d].src.str = 1;
      
      /* Destination slab: whole array */
      info[d].dst.off = 0;
      info[d].dst.len = cctk_gsh[d];
      info[d].dst.str = 1;
      
      /* No transformation */
      info[d].xpose = d;
      info[d].flip = 0;
    }
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0x,
       CCTK_VARIABLE_REAL, temp_phi0x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0y,
       CCTK_VARIABLE_REAL, temp_phi0y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0z,
       CCTK_VARIABLE_REAL, temp_phi0z);
    assert (! ierr);
    
    //  Begin Spectral-Galerkin Method
    //  Create phi, compute Forward FT
    
    phi_inx = fftw_alloc_complex(f_alloc_local);
    phi_iny = fftw_alloc_complex(f_alloc_local);
    phi_inz = fftw_alloc_complex(f_alloc_local);
    
    for(k = 0; k < f_local_n0; k++) 
    {
        phi_inx[k] = (fftw_complex)temp_phi0x[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_iny[k] = (fftw_complex)temp_phi0y[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inz[k] = (fftw_complex)temp_phi0z[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
    }
    
    plan_forwardx = fftw_mpi_plan_dft_1d(sizez,phi_inx,phi_inx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardy = fftw_mpi_plan_dft_1d(sizez,phi_iny,phi_iny,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardz = fftw_mpi_plan_dft_1d(sizez,phi_inz,phi_inz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    
    fftw_execute(plan_forwardx);
    fftw_execute(plan_forwardy);
    fftw_execute(plan_forwardz);
    
    //  Multiply i*k factor to FFTWave, create IK_FFTWavePadded
    Derivative_SGM3x = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3y = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3z = fftw_alloc_complex(b_alloc_local);
    
    for(k = 0; k < b_local_n0; k++) 
    {
                waveElementk = (CCTK_REAL)(k+b_local_0_start);  
                
                Fourier_3D_kwaveN3 = 0.0;
                
                if ( waveElementk < (sizez/2.0) ) 
                {
                Fourier_3D_kwaveN3 = waveElementk * Kgrid_res/(sizez*CCTK_DELTA_SPACE(2));
                }
                
                if ( waveElementk > (sizez/2.0) )
                {
                Fourier_3D_kwaveN3 = (waveElementk - sizez) * Kgrid_res/(sizez*CCTK_DELTA_SPACE(2));
                }
                
                Derivative_SGM3x[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inx[k];
                Derivative_SGM3y[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_iny[k];
                Derivative_SGM3z[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inz[k];
    }
    
    //  Create Derivative_SGM, compute Inverse FT 
    
    
    plan_backward3x = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3x,Derivative_SGM3x,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3y = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3y,Derivative_SGM3y,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3z = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3z,Derivative_SGM3z,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);

    fftw_execute(plan_backward3x);
    fftw_execute(plan_backward3y);
    fftw_execute(plan_backward3z);
    
    //   Store derivative into PHI 
    for(k = 0; k < b_local_n1; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1x[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2x[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3x[k]) / (sizez);
				temp_phi1y[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2y[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3y[k]) / (sizez);
				temp_phi1z[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2z[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3z[k]) / (sizez);
			}
		}
	}
	
	for (d=0; d<3; ++d) {
      /* Source array: FFTW layout */
      info[d].src.gsh         = cctk_gsh[d];
      info[0].src.lbnd        = cctk_nghostzones[0];
      info[1].src.lbnd        = cctk_nghostzones[1];
      info[2].src.lbnd        = b_local_1_end+cctk_nghostzones[2];
      info[0].src.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].src.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].src.lsh         = b_local_n1;
      info[0].src.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].src.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].src.ash         = b_local_n1;
      info[d].src.lbbox       = cctk_bbox[2*d];
      info[d].src.ubbox       = cctk_bbox[2*d+1];
      info[d].src.nghostzones = 0;
      
      /* Destination array: Cactus layout */
      info[d].dst.gsh         = cctk_gsh[d];
      info[d].dst.lbnd        = cctk_lbnd[d];
      info[d].dst.lsh         = cctk_lsh[d];
      info[d].dst.ash         = cctk_ash[d];
      info[d].dst.lbbox       = cctk_bbox[2*d];
      info[d].dst.ubbox       = cctk_bbox[2*d+1];
      info[d].dst.nghostzones = cctk_nghostzones[d];
      
      /* Source slab: whole array */
      info[d].src.off = 0;
      info[d].src.len = cctk_gsh[d];
      info[d].src.str = 1;
      
      /* Destination slab: whole array */
      info[d].dst.off = 0;
      info[d].dst.len = cctk_gsh[d];
      info[d].dst.str = 1;
      
      /* No transformation */
      info[d].xpose = d;
      info[d].flip = 0;
    }
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1x,
       CCTK_VARIABLE_REAL, temp_phi1x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2x,
       CCTK_VARIABLE_REAL, temp_phi2x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3x,
       CCTK_VARIABLE_REAL, temp_phi3x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1y,
       CCTK_VARIABLE_REAL, temp_phi1y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2y,
       CCTK_VARIABLE_REAL, temp_phi2y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3y,
       CCTK_VARIABLE_REAL, temp_phi3y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1z,
       CCTK_VARIABLE_REAL, temp_phi1z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2z,
       CCTK_VARIABLE_REAL, temp_phi2z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3z,
       CCTK_VARIABLE_REAL, temp_phi3z);
    assert (! ierr);
    
    //  Destroy/Free array pointers 

    fftw_destroy_plan(plan_forwardx);
    fftw_destroy_plan(plan_forwardy);
    fftw_destroy_plan(plan_forwardz);
    fftw_destroy_plan(plan_backward3x);
    fftw_destroy_plan(plan_backward3y);
    fftw_destroy_plan(plan_backward3z);
    fftw_free(phi_inx);
    fftw_free(phi_iny);
    fftw_free(phi_inz);
    fftw_free(Derivative_SGM3x);
    fftw_free(Derivative_SGM3y);
    fftw_free(Derivative_SGM3z);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Write_Vector_1d( CCTK_ARGUMENTS, CCTK_REAL ***vector )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, istart, jstart, kstart, iend, jend, kend, index;
	
	istart = cctk_nghostzones[0];
	jstart = cctk_nghostzones[1];
	kstart = cctk_nghostzones[2];
	iend   = cctk_lsh[ 0 ] - cctk_nghostzones[0];
	jend   = cctk_lsh[ 1 ] - cctk_nghostzones[1];
	kend   = cctk_lsh[ 2 ] - cctk_nghostzones[2];
    
    for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
            {
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                vector[1][1][index] = temp_phi1x[ index ];
                vector[2][1][index] = temp_phi2x[ index ];
                vector[3][1][index] = temp_phi3x[ index ];
                vector[1][2][index] = temp_phi1y[ index ];
                vector[2][2][index] = temp_phi2y[ index ];
                vector[3][2][index] = temp_phi3y[ index ];
                vector[1][3][index] = temp_phi1z[ index ];
                vector[2][3][index] = temp_phi2z[ index ];
                vector[3][3][index] = temp_phi3z[ index ];
			}
		}
	}
}
