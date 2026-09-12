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

void Calculate_A( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void Calculate_A_Write_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void Calculate_A_SpecGalMethod_Vector( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void Calculate_A( CCTK_ARGUMENTS, CCTK_REAL ***vector )
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
                temp_phi1x[ index ] = 0.0;
                temp_phi1y[ index ] = 0.0;
                temp_phi1z[ index ] = 0.0;
                temp_phi2x[ index ] = 0.0;
                temp_phi2y[ index ] = 0.0;
                temp_phi2z[ index ] = 0.0;
                temp_phi3x[ index ] = 0.0;
                temp_phi3y[ index ] = 0.0;
                temp_phi3z[ index ] = 0.0;
			}
		}
	}
	Calculate_A_SpecGalMethod_Vector( CCTK_PASS_CTOC );
	Calculate_A_Write_Vector( CCTK_PASS_CTOC, vector );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void Calculate_A_SpecGalMethod_Vector( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr, handle, rank;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    fftw_plan plan_forwardx;
    fftw_plan plan_forwardy;
    fftw_plan plan_forwardz;
    fftw_plan plan_backward1x;
    fftw_plan plan_backward2x;
    fftw_plan plan_backward3x;
    fftw_plan plan_backward1y;
    fftw_plan plan_backward2y;
    fftw_plan plan_backward3y;
    fftw_plan plan_backward1z;
    fftw_plan plan_backward2z;
    fftw_plan plan_backward3z;
    CCTK_REAL Fourier_3D_kwaveN1;
    CCTK_REAL Fourier_3D_kwaveN2;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementi, waveElementj, waveElementk, ksquare;
    CCTK_REAL lx, ly, lz, xmin, xmax, ymin, ymax, zmin, zmax;
    fftw_complex *phi_out1x, *phi_out1y, *phi_out1z;
    fftw_complex *phi_out2x, *phi_out2y, *phi_out2z;
    fftw_complex *phi_out3x, *phi_out3y, *phi_out3z;
    fftw_complex *Derivative_SGM1x, *Derivative_SGM1y, *Derivative_SGM1z;
    fftw_complex *Derivative_SGM2x, *Derivative_SGM2y, *Derivative_SGM2z;
    fftw_complex *Derivative_SGM3x, *Derivative_SGM3y, *Derivative_SGM3z;
    fftw_complex *phi_inx, *phi_iny, *phi_inz;
    fftw_complex *phi_outx, *phi_outy, *phi_outz;
    fftw_complex IMAG_NUMBER = _Complex_I;
    fftw_complex ZERO_COMPLEXNO = (0.0,0.0);
	MPI_Comm comm;
	
	handle = CCTK_ReductionHandle("maximum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmax, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymax, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmax, 1, CCTK_VarIndex("grid::z"));
	
	handle = CCTK_ReductionHandle("minimum");
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &xmin, 1, CCTK_VarIndex("grid::x"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &ymin, 1, CCTK_VarIndex("grid::y"));
	ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &zmin, 1, CCTK_VarIndex("grid::z"));

    xmin = xmin + 0.5*CCTK_DELTA_SPACE(0)*(2*cctk_nghostzones[0]-1);
    xmax = xmax - 0.5*CCTK_DELTA_SPACE(0)*(2*cctk_nghostzones[0]-1);
    lx = xmax - xmin; 
    
	ymin = ymin + 0.5*CCTK_DELTA_SPACE(1)*(2*cctk_nghostzones[1]-1);
    ymax = ymax - 0.5*CCTK_DELTA_SPACE(1)*(2*cctk_nghostzones[1]-1);
    ly = ymax - ymin; 
    
	zmin = zmin + 0.5*CCTK_DELTA_SPACE(2)*(2*cctk_nghostzones[2]-1);
    zmax = zmax - 0.5*CCTK_DELTA_SPACE(2)*(2*cctk_nghostzones[2]-1);
    lz = zmax - zmin;  
  
    sizex  = cctk_gsh[0]-2*cctk_nghostzones[0]; 
    sizey  = cctk_gsh[1]-2*cctk_nghostzones[1];
    sizez  = cctk_gsh[2]-2*cctk_nghostzones[2]; 

    //comm = get_mpi_comm (cctkGH);
    fftw_mpi_init();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) fftw_import_wisdom_from_filename("fftw.wis");
    fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD);
    
    alloc_local = fftw_mpi_local_size_3d(sizex, sizey, sizez, MPI_COMM_WORLD, &local_n0, &local_0_start);
    
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
      info[0].dst.lbnd        = 0.0;
      info[1].dst.lbnd        = 0.0;
      info[2].dst.lbnd        = local_0_start;
      info[0].dst.lsh         = cctk_gsh[0];
      info[1].dst.lsh         = cctk_gsh[1];
      info[2].dst.lsh         = local_n0;
      info[0].dst.ash         = cctk_gsh[0];
      info[1].dst.ash         = cctk_gsh[1];
      info[2].dst.ash         = local_n0;
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
    
    phi_inx = fftw_alloc_complex(alloc_local);
    phi_iny = fftw_alloc_complex(alloc_local);
    phi_inz = fftw_alloc_complex(alloc_local);
    phi_outx = fftw_alloc_complex(alloc_local);
    phi_outy = fftw_alloc_complex(alloc_local);
    phi_outz = fftw_alloc_complex(alloc_local);
    
    //Calculate the Vector Potential
    
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
              phi_inx[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0x[i + j*sizex + k*sizex*sizey];
              phi_iny[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0y[i + j*sizex + k*sizex*sizey];
              phi_inz[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0z[i + j*sizex + k*sizex*sizey];
              
              phi_outx[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outy[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outz[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
            }

        }

    }
    
    plan_forwardx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inx,phi_outx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_iny,phi_outy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inz,phi_outz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    
    fftw_execute(plan_forwardx);
    fftw_execute(plan_forwardy);
    fftw_execute(plan_forwardz);
    
    //  Multiply i*k factor to FFTWave, create IK_FFTWavePadded
    phi_out1x = fftw_alloc_complex(alloc_local);
    phi_out2x = fftw_alloc_complex(alloc_local);
    phi_out3x = fftw_alloc_complex(alloc_local);
    phi_out1y = fftw_alloc_complex(alloc_local);
    phi_out2y = fftw_alloc_complex(alloc_local);
    phi_out3y = fftw_alloc_complex(alloc_local);
    phi_out1z = fftw_alloc_complex(alloc_local);
    phi_out2z = fftw_alloc_complex(alloc_local);
    phi_out3z = fftw_alloc_complex(alloc_local);
    
    Derivative_SGM1x = fftw_alloc_complex(alloc_local);
    Derivative_SGM2x = fftw_alloc_complex(alloc_local);
    Derivative_SGM3x = fftw_alloc_complex(alloc_local);
    Derivative_SGM1y = fftw_alloc_complex(alloc_local);
    Derivative_SGM2y = fftw_alloc_complex(alloc_local);
    Derivative_SGM3y = fftw_alloc_complex(alloc_local);
    Derivative_SGM1z = fftw_alloc_complex(alloc_local);
    Derivative_SGM2z = fftw_alloc_complex(alloc_local);
    Derivative_SGM3z = fftw_alloc_complex(alloc_local);
    
    
    
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
                waveElementi = (CCTK_REAL)i;   
                waveElementj = (CCTK_REAL)j;   
                waveElementk = (CCTK_REAL)(k+local_0_start);  
                
                Fourier_3D_kwaveN1 = 0.0;
                Fourier_3D_kwaveN2 = 0.0;
                Fourier_3D_kwaveN3 = 0.0;
                
                if ( waveElementi < (sizex/2.0) )
                {
                Fourier_3D_kwaveN1 = waveElementi * Kgrid_res/lx;
                }
                
                if ( waveElementj < (sizey/2.0) )
                {
                Fourier_3D_kwaveN2 = waveElementj * Kgrid_res/ly;
                }
                
                if ( waveElementk < (sizez/2.0) ) 
                {
                Fourier_3D_kwaveN3 = waveElementk * Kgrid_res/lz;
                }
                
                if ( waveElementi > (sizex/2.0) )
                {
                Fourier_3D_kwaveN1 = (waveElementi - sizex) * Kgrid_res/lx;
                }
                
                if ( waveElementj > (sizey/2.0) )
                { 
                Fourier_3D_kwaveN2 = (waveElementj - sizey) * Kgrid_res/ly;
                }
                
                if ( waveElementk > (sizez/2.0) )
                {
                Fourier_3D_kwaveN3 = (waveElementk - sizez) * Kgrid_res/lz;
                }
                
                ksquare = Fourier_3D_kwaveN1*Fourier_3D_kwaveN1+Fourier_3D_kwaveN2*Fourier_3D_kwaveN2+Fourier_3D_kwaveN3*Fourier_3D_kwaveN3;
                
                if (ksquare != 0.0) {
                Derivative_SGM1x[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * (Fourier_3D_kwaveN2 
                                            * phi_outz[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                            
                Derivative_SGM2y[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * (Fourier_3D_kwaveN3 
                                            * phi_outx[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN1
                                            * phi_outz[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                            
                Derivative_SGM3z[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * (Fourier_3D_kwaveN1 
                                            * phi_outy[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                
                Derivative_SGM2x[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN2
                                            * (Fourier_3D_kwaveN2 
                                            * phi_outz[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                            
                Derivative_SGM3x[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN3
                                            * (Fourier_3D_kwaveN2 
                                            * phi_outz[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                            
                Derivative_SGM1y[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN1
                							* (Fourier_3D_kwaveN3 
                                            * phi_outx[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN1
                                            * phi_outz[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                            
                Derivative_SGM3y[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN3
                							* (Fourier_3D_kwaveN3 
                                            * phi_outx[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN1
                                            * phi_outz[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                                        
                Derivative_SGM1z[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN1
                							* (Fourier_3D_kwaveN1 
                                            * phi_outy[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                
                Derivative_SGM2z[i + j*sizex + k*sizex*sizey] = -Fourier_3D_kwaveN2
                							* (Fourier_3D_kwaveN1 
                                            * phi_outy[i + j*sizex + k*sizex*sizey] 
                                            - Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey])
                                            / ksquare;
                
            	} else {
            	Derivative_SGM1x[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM2x[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM3x[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM1y[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM2y[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM3y[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM1z[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM2z[i + j*sizex + k*sizex*sizey] = 0.0;
            	Derivative_SGM3z[i + j*sizex + k*sizex*sizey] = 0.0;
            	}
            	
            	// Dealias using truncation (2/3 rule)
                 
                 if (( waveElementi > (2.0*sizex/3.0) ) || ( waveElementi < (sizex/3.0) )) {
                 	Derivative_SGM1x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementj > (2.0*sizey/3.0) ) || ( waveElementj < (sizey/3.0) )) {
                 	Derivative_SGM2x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementk > (2.0*sizez/3.0) ) || ( waveElementk < (sizez/3.0) )) {
                 	Derivative_SGM3x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                                            
            }

        }

    }
    
    //  Create Derivative_SGM, compute Inverse FT 
    
    plan_backward1x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1x,phi_out1x,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2x,phi_out2x,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3x,phi_out3x,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1y,phi_out1y,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2y,phi_out2y,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3y,phi_out3y,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1z,phi_out1z,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2z,phi_out2z,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3z,phi_out3z,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);

    fftw_execute(plan_backward1x);
    fftw_execute(plan_backward2x);
    fftw_execute(plan_backward3x);
    fftw_execute(plan_backward1y);
    fftw_execute(plan_backward2y);
    fftw_execute(plan_backward3y);
    fftw_execute(plan_backward1z);
    fftw_execute(plan_backward2z);
    fftw_execute(plan_backward3z);
    
    //   Store derivative into PHI 
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1x[ i + j*sizex + k*sizex*sizey ] = creal(phi_out1x[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2x[ i + j*sizex + k*sizex*sizey ] = creal(phi_out2x[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3x[ i + j*sizex + k*sizex*sizey ] = creal(phi_out3x[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1y[ i + j*sizex + k*sizex*sizey ] = creal(phi_out1y[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2y[ i + j*sizex + k*sizex*sizey ] = creal(phi_out2y[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3y[ i + j*sizex + k*sizex*sizey ] = creal(phi_out3y[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1z[ i + j*sizex + k*sizex*sizey ] = creal(phi_out1z[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2z[ i + j*sizex + k*sizex*sizey ] = creal(phi_out2z[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3z[ i + j*sizex + k*sizex*sizey ] = creal(phi_out3z[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
			}
		}
	}
	
	for (d=0; d<3; ++d) {
      /* Source array: FFTW layout */
      info[d].src.gsh         = cctk_gsh[d];
      info[0].src.lbnd        = 0.0;
      info[1].src.lbnd        = 0.0;
      info[2].src.lbnd        = local_0_start;
      info[0].src.lsh         = cctk_gsh[0];
      info[1].src.lsh         = cctk_gsh[1];
      info[2].src.lsh         = local_n0;
      info[0].src.ash         = cctk_gsh[0];
      info[1].src.ash         = cctk_gsh[1];
      info[2].src.ash         = local_n0;
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
    fftw_destroy_plan(plan_backward1x);
    fftw_destroy_plan(plan_backward2x);
    fftw_destroy_plan(plan_backward3x);
    fftw_destroy_plan(plan_backward1y);
    fftw_destroy_plan(plan_backward2y);
    fftw_destroy_plan(plan_backward3y);
    fftw_destroy_plan(plan_backward1z);
    fftw_destroy_plan(plan_backward2z);
    fftw_destroy_plan(plan_backward3z);
    fftw_free(phi_inx);
    fftw_free(phi_iny);
    fftw_free(phi_inz);
    fftw_free(phi_outx);
    fftw_free(phi_outy);
    fftw_free(phi_outz);
    fftw_free(phi_out1x);
    fftw_free(phi_out1y);
    fftw_free(phi_out1z);
    fftw_free(phi_out2x);
    fftw_free(phi_out2y);
    fftw_free(phi_out2z);
    fftw_free(phi_out3x);
    fftw_free(phi_out3y);
    fftw_free(phi_out3z);
    fftw_free(Derivative_SGM1x);
    fftw_free(Derivative_SGM2x);
    fftw_free(Derivative_SGM3x);
    fftw_free(Derivative_SGM1y);
    fftw_free(Derivative_SGM2y);
    fftw_free(Derivative_SGM3y);
    fftw_free(Derivative_SGM1z);
    fftw_free(Derivative_SGM2z);
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
void Calculate_A_Write_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, istart, jstart, kstart, iend, jend, kend, index;
	
	istart = 0;
	jstart = 0;
	kstart = 0;
	iend   = cctk_lsh[ 0 ];
	jend   = cctk_lsh[ 1 ];
	kend   = cctk_lsh[ 2 ];
    
    for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
            {
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                vector[0][1][index] = temp_phi3y[ index ] - temp_phi2z[ index ];
                vector[0][2][index] = temp_phi1z[ index ] - temp_phi3x[ index ];
                vector[0][3][index] = temp_phi2x[ index ] - temp_phi1y[ index ];
                vector[1][1][index] = temp_phi1x[ index ];
                vector[2][2][index] = temp_phi2y[ index ];
                vector[3][3][index] = temp_phi3z[ index ];
                
			}
		}
	}
}
