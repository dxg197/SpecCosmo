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

void SpecDeriv_Vector_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_Write_Vector2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 );
void SpecDeriv_SpecGalMethod_Vector2( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Vector_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 )
{
	DECLARE_CCTK_ARGUMENTS
    // Declare and initialize variables
    CCTK_INT i, j, k, m, n, p, iend, jend, kend, index;
	
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
                temp_phi0x[ index ] = vector2[0][0][1][index];
                temp_phi0y[ index ] = vector2[0][0][2][index];
                temp_phi0z[ index ] = vector2[0][0][3][index];
                
                for(m=1; m < 4; m++) {
                	for(n=0; n < 4; n++) {
                		for(p=1; p < 4; p++) {
                			vector2[m][n][p][index] = 0.0;
                		}
                	}
                }
                
                for(m=0; m < 4; m++) {
                	for(n=1; n < 4; n++) {
                		for(p=1; p < 4; p++) {
                			vector2[m][n][p][index] = 0.0;
                		}
                	}
                }
			}
		}
	}
	SpecDeriv_SpecGalMethod_Vector2( CCTK_PASS_CTOC );
	SpecDeriv_Write_Vector2( CCTK_PASS_CTOC, vector2 );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_SpecGalMethod_Vector2( CCTK_ARGUMENTS )
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
    fftw_plan plan_backward11x;
    fftw_plan plan_backward22x;
    fftw_plan plan_backward33x;
    fftw_plan plan_backward11y;
    fftw_plan plan_backward22y;
    fftw_plan plan_backward33y;
    fftw_plan plan_backward11z;
    fftw_plan plan_backward22z;
    fftw_plan plan_backward33z;
    fftw_plan plan_backward12x;
    fftw_plan plan_backward13x;
    fftw_plan plan_backward23x;
    fftw_plan plan_backward12y;
    fftw_plan plan_backward13y;
    fftw_plan plan_backward23y;
    fftw_plan plan_backward12z;
    fftw_plan plan_backward13z;
    fftw_plan plan_backward23z;
    CCTK_REAL Fourier_3D_kwaveN1;
    CCTK_REAL Fourier_3D_kwaveN2;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementi, waveElementj, waveElementk;
    CCTK_REAL lx, ly, lz, xmin, xmax, ymin, ymax, zmin, zmax;
    fftw_complex *Derivative_SGM1x, *Derivative_SGM1y, *Derivative_SGM1z;
    fftw_complex *Derivative_SGM2x, *Derivative_SGM2y, *Derivative_SGM2z;
    fftw_complex *Derivative_SGM3x, *Derivative_SGM3y, *Derivative_SGM3z;
    fftw_complex *Derivative_SGM11x, *Derivative_SGM11y, *Derivative_SGM11z;
    fftw_complex *Derivative_SGM22x, *Derivative_SGM22y, *Derivative_SGM22z;
    fftw_complex *Derivative_SGM33x, *Derivative_SGM33y, *Derivative_SGM33z;
    fftw_complex *Derivative_SGM12x, *Derivative_SGM12y, *Derivative_SGM12z;
    fftw_complex *Derivative_SGM13x, *Derivative_SGM13y, *Derivative_SGM13z;
    fftw_complex *Derivative_SGM23x, *Derivative_SGM23y, *Derivative_SGM23z;
    fftw_complex *phi_inx, *phi_iny, *phi_inz;
    fftw_complex *Derivative_SGM1x_out, *Derivative_SGM1y_out, *Derivative_SGM1z_out;
    fftw_complex *Derivative_SGM2x_out, *Derivative_SGM2y_out, *Derivative_SGM2z_out;
    fftw_complex *Derivative_SGM3x_out, *Derivative_SGM3y_out, *Derivative_SGM3z_out;
    fftw_complex *Derivative_SGM11x_out, *Derivative_SGM11y_out, *Derivative_SGM11z_out;
    fftw_complex *Derivative_SGM22x_out, *Derivative_SGM22y_out, *Derivative_SGM22z_out;
    fftw_complex *Derivative_SGM33x_out, *Derivative_SGM33y_out, *Derivative_SGM33z_out;
    fftw_complex *Derivative_SGM12x_out, *Derivative_SGM12y_out, *Derivative_SGM12z_out;
    fftw_complex *Derivative_SGM13x_out, *Derivative_SGM13y_out, *Derivative_SGM13z_out;
    fftw_complex *Derivative_SGM23x_out, *Derivative_SGM23y_out, *Derivative_SGM23z_out;
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
    
    struct xferinfo info1[3];
    struct xferinfo info2[3];
    
    for (d=0; d<3; ++d) {
      /* Source array: Cactus layout */
      info1[d].src.gsh         = cctk_gsh[d];
      info1[d].src.lbnd        = cctk_lbnd[d];
      info1[d].src.lsh         = cctk_lsh[d];
      info1[d].src.ash         = cctk_ash[d];
      info1[d].src.lbbox       = cctk_bbox[2*d];
      info1[d].src.ubbox       = cctk_bbox[2*d+1];
      info1[d].src.nghostzones = cctk_nghostzones[d];
      
      /* Destination array: FFTW layout */
      info1[d].dst.gsh         = cctk_gsh[d];
      info1[0].dst.lbnd        = cctk_nghostzones[0];
      info1[1].dst.lbnd        = cctk_nghostzones[1];
      info1[2].dst.lbnd        = local_0_start+cctk_nghostzones[2];
      info1[0].dst.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info1[1].dst.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info1[2].dst.lsh         = local_n0;
      info1[0].dst.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info1[1].dst.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info1[2].dst.ash         = local_n0;
      info1[d].dst.lbbox       = cctk_bbox[2*d];
      info1[d].dst.ubbox       = cctk_bbox[2*d+1];
      info1[d].dst.nghostzones = 0; 
      
      /* Source slab: whole array */
      info1[d].src.off = 0;
      info1[d].src.len = cctk_gsh[d];
      info1[d].src.str = 1;
      
      /* Destination slab: whole array */
      info1[d].dst.off = 0;
      info1[d].dst.len = cctk_gsh[d];
      info1[d].dst.str = 1;
      
      /* No transformation */
      info1[d].xpose = d;
      info1[d].flip = 0;
    }
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0x,
       CCTK_VARIABLE_REAL, temp_phi0x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0y,
       CCTK_VARIABLE_REAL, temp_phi0y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
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
    Derivative_SGM1x = fftw_alloc_complex(alloc_local);
    Derivative_SGM2x = fftw_alloc_complex(alloc_local);
    Derivative_SGM3x = fftw_alloc_complex(alloc_local);
    Derivative_SGM1y = fftw_alloc_complex(alloc_local);
    Derivative_SGM2y = fftw_alloc_complex(alloc_local);
    Derivative_SGM3y = fftw_alloc_complex(alloc_local);
    Derivative_SGM1z = fftw_alloc_complex(alloc_local);
    Derivative_SGM2z = fftw_alloc_complex(alloc_local);
    Derivative_SGM3z = fftw_alloc_complex(alloc_local);
    Derivative_SGM11x = fftw_alloc_complex(alloc_local);
    Derivative_SGM22x = fftw_alloc_complex(alloc_local);
    Derivative_SGM33x = fftw_alloc_complex(alloc_local);
    Derivative_SGM11y = fftw_alloc_complex(alloc_local);
    Derivative_SGM22y = fftw_alloc_complex(alloc_local);
    Derivative_SGM33y = fftw_alloc_complex(alloc_local);
    Derivative_SGM11z = fftw_alloc_complex(alloc_local);
    Derivative_SGM22z = fftw_alloc_complex(alloc_local);
    Derivative_SGM33z = fftw_alloc_complex(alloc_local);
    Derivative_SGM12x = fftw_alloc_complex(alloc_local);
    Derivative_SGM13x = fftw_alloc_complex(alloc_local);
    Derivative_SGM23x = fftw_alloc_complex(alloc_local);
    Derivative_SGM12y = fftw_alloc_complex(alloc_local);
    Derivative_SGM13y = fftw_alloc_complex(alloc_local);
    Derivative_SGM23y = fftw_alloc_complex(alloc_local);
    Derivative_SGM12z = fftw_alloc_complex(alloc_local);
    Derivative_SGM13z = fftw_alloc_complex(alloc_local);
    Derivative_SGM23z = fftw_alloc_complex(alloc_local);
    
    Derivative_SGM1x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM11x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM22x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM33x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM11y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM22y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM33y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM11z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM22z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM33z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM12x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM13x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM23x_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM12y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM13y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM23y_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM12z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM13z_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM23z_out = fftw_alloc_complex(alloc_local);
    
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
                
                Derivative_SGM1x[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1
                                            * phi_outx[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2x[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3x[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1y[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1
                                            * phi_outy[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2y[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3y[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1z[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1
                                            * phi_outz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2z[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3z[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM11x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN1
                                            * phi_outx[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM22x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM33x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN3
                                            * Fourier_3D_kwaveN3
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM11y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1 
                                            * Fourier_3D_kwaveN1
                                            * phi_outy[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM22y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN2
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM33y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN3
                                            * Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM11z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN1
                                            * phi_outz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM22z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN2
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM33z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN3
                                            * Fourier_3D_kwaveN3
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM12x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN2
                                            * phi_outx[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM13x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN3
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM23x[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN3
                                            * phi_outx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM12y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN2
                                            * phi_outy[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM13y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM23y[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN3
                                            * phi_outy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM12z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN2
                                            * phi_outz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM13z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN1
                                            * Fourier_3D_kwaveN3
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM23z[i + j*sizex + k*sizex*sizey] = -1.0 
                                            * Fourier_3D_kwaveN2
                                            * Fourier_3D_kwaveN3
                                            * phi_outz[i + j*sizex + k*sizex*sizey];
                                            
                                            
                Derivative_SGM1x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                Derivative_SGM11x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM12x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM13x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM11y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM12y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM13y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM11z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM12z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM13z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                Derivative_SGM22x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM23x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM22y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM23y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM22z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM23z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                Derivative_SGM33x_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM33y_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM33z_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                // Dealias using truncation (2/3 rule)
                 
                 if (( waveElementi > (2.0*sizex/3.0) ) || ( waveElementi < (sizex/3.0) )) {
                 	Derivative_SGM1x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM11x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM11y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM11z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementj > (2.0*sizey/3.0) ) || ( waveElementj < (sizey/3.0) )) {
                 	Derivative_SGM2x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM22x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM22y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM22z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementk > (2.0*sizez/3.0) ) || ( waveElementk < (sizez/3.0) )) {
                 	Derivative_SGM3x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM33x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM33y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM33z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13z[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23x[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23y[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23z[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
            }

        }

    }
    
    //  Create Derivative_SGM, compute Inverse FT 
    
    
    plan_backward1x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1x,Derivative_SGM1x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2x,Derivative_SGM2x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3x,Derivative_SGM3x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1y,Derivative_SGM1y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2y,Derivative_SGM2y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3y,Derivative_SGM3y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1z,Derivative_SGM1z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2z,Derivative_SGM2z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3z,Derivative_SGM3z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward11x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM11x,Derivative_SGM11x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward22x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM22x,Derivative_SGM22x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward33x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM33x,Derivative_SGM33x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward11y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM11y,Derivative_SGM11y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward22y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM22y,Derivative_SGM22y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward33y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM33y,Derivative_SGM33y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward11z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM11z,Derivative_SGM11z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward22z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM22z,Derivative_SGM22z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward33z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM33z,Derivative_SGM33z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward12x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM12x,Derivative_SGM12x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward13x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM13x,Derivative_SGM13x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward23x = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM23x,Derivative_SGM23x_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward12y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM12y,Derivative_SGM12y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward13y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM13y,Derivative_SGM13y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward23y = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM23y,Derivative_SGM23y_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward12z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM12z,Derivative_SGM12z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward13z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM13z,Derivative_SGM13z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward23z = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM23z,Derivative_SGM23z_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);

    fftw_execute(plan_backward1x);
    fftw_execute(plan_backward2x);
    fftw_execute(plan_backward3x);
    fftw_execute(plan_backward1y);
    fftw_execute(plan_backward2y);
    fftw_execute(plan_backward3y);
    fftw_execute(plan_backward1z);
    fftw_execute(plan_backward2z);
    fftw_execute(plan_backward3z);
    fftw_execute(plan_backward11x);
    fftw_execute(plan_backward22x);
    fftw_execute(plan_backward33x);
    fftw_execute(plan_backward11y);
    fftw_execute(plan_backward22y);
    fftw_execute(plan_backward33y);
    fftw_execute(plan_backward11z);
    fftw_execute(plan_backward22z);
    fftw_execute(plan_backward33z);
    fftw_execute(plan_backward12x);
    fftw_execute(plan_backward13x);
    fftw_execute(plan_backward23x);
    fftw_execute(plan_backward12y);
    fftw_execute(plan_backward13y);
    fftw_execute(plan_backward23y);
    fftw_execute(plan_backward12z);
    fftw_execute(plan_backward13z);
    fftw_execute(plan_backward23z);
    
    //   Store derivative into PHI 
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi11x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM11x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi22x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM22x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi33x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi11y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM11y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi22y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM22y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi33y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi11z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM11z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi22z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM22z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi33z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi12x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM12x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi13x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM13x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi23x[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM23x_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi12y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM12y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi13y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM13y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi23y[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM23y_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi12z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM12z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi13z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM13z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi23z[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM23z_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
			}
		}
	}
	
	for (d=0; d<3; ++d) {
      /* Source array: FFTW layout */
      info2[d].src.gsh         = cctk_gsh[d];
      info2[0].src.lbnd        = cctk_nghostzones[0];
      info2[1].src.lbnd        = cctk_nghostzones[1];
      info2[2].src.lbnd        = local_0_start+cctk_nghostzones[2];
      info2[0].src.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info2[1].src.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info2[2].src.lsh         = local_n0;
      info2[0].src.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info2[1].src.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info2[2].src.ash         = local_n0;
      info2[d].src.lbbox       = cctk_bbox[2*d];
      info2[d].src.ubbox       = cctk_bbox[2*d+1];
      info2[d].src.nghostzones = 0;
      
      /* Destination array: Cactus layout */
      info2[d].dst.gsh         = cctk_gsh[d];
      info2[d].dst.lbnd        = cctk_lbnd[d];
      info2[d].dst.lsh         = cctk_lsh[d];
      info2[d].dst.ash         = cctk_ash[d];
      info2[d].dst.lbbox       = cctk_bbox[2*d];
      info2[d].dst.ubbox       = cctk_bbox[2*d+1];
      info2[d].dst.nghostzones = cctk_nghostzones[d];
      
      /* Source slab: whole array */
      info2[d].src.off = 0;
      info2[d].src.len = cctk_gsh[d];
      info2[d].src.str = 1;
      
      /* Destination slab: whole array */
      info2[d].dst.off = 0;
      info2[d].dst.len = cctk_gsh[d];
      info2[d].dst.str = 1;
      
      /* No transformation */
      info2[d].xpose = d;
      info2[d].flip = 0;
    }
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1x,
       CCTK_VARIABLE_REAL, temp_phi1x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2x,
       CCTK_VARIABLE_REAL, temp_phi2x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3x,
       CCTK_VARIABLE_REAL, temp_phi3x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1y,
       CCTK_VARIABLE_REAL, temp_phi1y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2y,
       CCTK_VARIABLE_REAL, temp_phi2y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3y,
       CCTK_VARIABLE_REAL, temp_phi3y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1z,
       CCTK_VARIABLE_REAL, temp_phi1z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2z,
       CCTK_VARIABLE_REAL, temp_phi2z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3z,
       CCTK_VARIABLE_REAL, temp_phi3z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi11x,
       CCTK_VARIABLE_REAL, temp_phi11x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi22x,
       CCTK_VARIABLE_REAL, temp_phi22x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi33x,
       CCTK_VARIABLE_REAL, temp_phi33x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi11y,
       CCTK_VARIABLE_REAL, temp_phi11y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi22y,
       CCTK_VARIABLE_REAL, temp_phi22y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi33y,
       CCTK_VARIABLE_REAL, temp_phi33y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi11z,
       CCTK_VARIABLE_REAL, temp_phi11z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi22z,
       CCTK_VARIABLE_REAL, temp_phi22z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi33z,
       CCTK_VARIABLE_REAL, temp_phi33z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi12x,
       CCTK_VARIABLE_REAL, temp_phi12x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi13x,
       CCTK_VARIABLE_REAL, temp_phi13x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi23x,
       CCTK_VARIABLE_REAL, temp_phi23x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi12y,
       CCTK_VARIABLE_REAL, temp_phi12y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi13y,
       CCTK_VARIABLE_REAL, temp_phi13y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi23y,
       CCTK_VARIABLE_REAL, temp_phi23y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi12z,
       CCTK_VARIABLE_REAL, temp_phi12z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi13z,
       CCTK_VARIABLE_REAL, temp_phi13z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi23z,
       CCTK_VARIABLE_REAL, temp_phi23z);
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
    fftw_destroy_plan(plan_backward11x);
    fftw_destroy_plan(plan_backward22x);
    fftw_destroy_plan(plan_backward33x);
    fftw_destroy_plan(plan_backward11y);
    fftw_destroy_plan(plan_backward22y);
    fftw_destroy_plan(plan_backward33y);
    fftw_destroy_plan(plan_backward11z);
    fftw_destroy_plan(plan_backward22z);
    fftw_destroy_plan(plan_backward33z);
    fftw_destroy_plan(plan_backward12x);
    fftw_destroy_plan(plan_backward13x);
    fftw_destroy_plan(plan_backward23x);
    fftw_destroy_plan(plan_backward12y);
    fftw_destroy_plan(plan_backward13y);
    fftw_destroy_plan(plan_backward23y);
    fftw_destroy_plan(plan_backward12z);
    fftw_destroy_plan(plan_backward13z);
    fftw_destroy_plan(plan_backward23z);  
    fftw_free(phi_inx);
    fftw_free(phi_iny);
    fftw_free(phi_inz);
    fftw_free(Derivative_SGM1x);
    fftw_free(Derivative_SGM2x);
    fftw_free(Derivative_SGM3x);
    fftw_free(Derivative_SGM1y);
    fftw_free(Derivative_SGM2y);
    fftw_free(Derivative_SGM3y);
    fftw_free(Derivative_SGM1z);
    fftw_free(Derivative_SGM2z);
    fftw_free(Derivative_SGM3z);
    fftw_free(Derivative_SGM11x);
    fftw_free(Derivative_SGM22x);
    fftw_free(Derivative_SGM33x);
    fftw_free(Derivative_SGM11y);
    fftw_free(Derivative_SGM22y);
    fftw_free(Derivative_SGM33y);
    fftw_free(Derivative_SGM11z);
    fftw_free(Derivative_SGM22z);
    fftw_free(Derivative_SGM33z);
    fftw_free(Derivative_SGM12x);
    fftw_free(Derivative_SGM13x);
    fftw_free(Derivative_SGM23x);
    fftw_free(Derivative_SGM12y);
    fftw_free(Derivative_SGM13y);
    fftw_free(Derivative_SGM23y);
    fftw_free(Derivative_SGM12z);
    fftw_free(Derivative_SGM13z);
    fftw_free(Derivative_SGM23z);
    fftw_free(phi_outx);
    fftw_free(phi_outy);
    fftw_free(phi_outz);
    fftw_free(Derivative_SGM1x_out);
    fftw_free(Derivative_SGM2x_out);
    fftw_free(Derivative_SGM3x_out);
    fftw_free(Derivative_SGM1y_out);
    fftw_free(Derivative_SGM2y_out);
    fftw_free(Derivative_SGM3y_out);
    fftw_free(Derivative_SGM1z_out);
    fftw_free(Derivative_SGM2z_out);
    fftw_free(Derivative_SGM3z_out);
    fftw_free(Derivative_SGM11x_out);
    fftw_free(Derivative_SGM22x_out);
    fftw_free(Derivative_SGM33x_out);
    fftw_free(Derivative_SGM11y_out);
    fftw_free(Derivative_SGM22y_out);
    fftw_free(Derivative_SGM33y_out);
    fftw_free(Derivative_SGM11z_out);
    fftw_free(Derivative_SGM22z_out);
    fftw_free(Derivative_SGM33z_out);
    fftw_free(Derivative_SGM12x_out);
    fftw_free(Derivative_SGM13x_out);
    fftw_free(Derivative_SGM23x_out);
    fftw_free(Derivative_SGM12y_out);
    fftw_free(Derivative_SGM13y_out);
    fftw_free(Derivative_SGM23y_out);
    fftw_free(Derivative_SGM12z_out);
    fftw_free(Derivative_SGM13z_out);
    fftw_free(Derivative_SGM23z_out);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Write_Vector2( CCTK_ARGUMENTS, CCTK_REAL ****vector2 )
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
    
	for(k=0; k < kstart; k++)
	{
		for(j=0; j < jstart; j++)
		{
			for(i=0; i < istart; i++)
            {
            	index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                vector2[1][0][1][index] = 0.0;
                vector2[2][0][1][index] = 0.0;
                vector2[3][0][1][index] = 0.0;
                vector2[1][0][2][index] = 0.0;
                vector2[2][0][2][index] = 0.0;
                vector2[3][0][2][index] = 0.0;
                vector2[1][0][3][index] = 0.0;
                vector2[2][0][3][index] = 0.0;
                vector2[3][0][3][index] = 0.0;
                
                vector2[0][1][1][index] = 0.0;
                vector2[0][2][1][index] = 0.0;
                vector2[0][3][1][index] = 0.0;
                vector2[0][1][2][index] = 0.0;
                vector2[0][2][2][index] = 0.0;
                vector2[0][3][2][index] = 0.0;
                vector2[0][1][3][index] = 0.0;
                vector2[0][2][3][index] = 0.0;
                vector2[0][3][3][index] = 0.0;
                
                vector2[1][1][1][index] = 0.0;
                vector2[2][2][1][index] = 0.0;
                vector2[3][3][1][index] = 0.0;
                vector2[1][1][2][index] = 0.0;
                vector2[2][2][2][index] = 0.0;
                vector2[3][3][2][index] = 0.0;
                vector2[1][1][3][index] = 0.0;
                vector2[2][2][3][index] = 0.0;
                vector2[3][3][3][index] = 0.0;
                
                vector2[1][2][1][index] = 0.0;
                vector2[1][3][1][index] = 0.0;
                vector2[2][3][1][index] = 0.0;
                vector2[1][2][2][index] = 0.0;
                vector2[1][3][2][index] = 0.0;
                vector2[2][3][2][index] = 0.0;
                vector2[1][2][3][index] = 0.0;
                vector2[1][3][3][index] = 0.0;
                vector2[2][3][3][index] = 0.0;
                
                vector2[2][1][1][index] = 0.0;
                vector2[3][1][1][index] = 0.0;
                vector2[3][2][1][index] = 0.0;
                vector2[2][1][2][index] = 0.0;
                vector2[3][1][2][index] = 0.0;
                vector2[3][2][2][index] = 0.0;
                vector2[2][1][3][index] = 0.0;
                vector2[3][1][3][index] = 0.0;
                vector2[3][2][3][index] = 0.0;
            }
        }
    }
    
    for(k=kend; k < cctk_lsh[2]; k++)
	{
		for(j=jend; j < cctk_lsh[1]; j++)
		{
			for(i=iend; i < cctk_lsh[0]; i++)
            {
            	index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                vector2[1][0][1][index] = 0.0;
                vector2[2][0][1][index] = 0.0;
                vector2[3][0][1][index] = 0.0;
                vector2[1][0][2][index] = 0.0;
                vector2[2][0][2][index] = 0.0;
                vector2[3][0][2][index] = 0.0;
                vector2[1][0][3][index] = 0.0;
                vector2[2][0][3][index] = 0.0;
                vector2[3][0][3][index] = 0.0;
                
                vector2[0][1][1][index] = 0.0;
                vector2[0][2][1][index] = 0.0;
                vector2[0][3][1][index] = 0.0;
                vector2[0][1][2][index] = 0.0;
                vector2[0][2][2][index] = 0.0;
                vector2[0][3][2][index] = 0.0;
                vector2[0][1][3][index] = 0.0;
                vector2[0][2][3][index] = 0.0;
                vector2[0][3][3][index] = 0.0;
                
                vector2[1][1][1][index] = 0.0;
                vector2[2][2][1][index] = 0.0;
                vector2[3][3][1][index] = 0.0;
                vector2[1][1][2][index] = 0.0;
                vector2[2][2][2][index] = 0.0;
                vector2[3][3][2][index] = 0.0;
                vector2[1][1][3][index] = 0.0;
                vector2[2][2][3][index] = 0.0;
                vector2[3][3][3][index] = 0.0;
                
                vector2[1][2][1][index] = 0.0;
                vector2[1][3][1][index] = 0.0;
                vector2[2][3][1][index] = 0.0;
                vector2[1][2][2][index] = 0.0;
                vector2[1][3][2][index] = 0.0;
                vector2[2][3][2][index] = 0.0;
                vector2[1][2][3][index] = 0.0;
                vector2[1][3][3][index] = 0.0;
                vector2[2][3][3][index] = 0.0;
                
                vector2[2][1][1][index] = 0.0;
                vector2[3][1][1][index] = 0.0;
                vector2[3][2][1][index] = 0.0;
                vector2[2][1][2][index] = 0.0;
                vector2[3][1][2][index] = 0.0;
                vector2[3][2][2][index] = 0.0;
                vector2[2][1][3][index] = 0.0;
                vector2[3][1][3][index] = 0.0;
                vector2[3][2][3][index] = 0.0;
            }
        }
    }
    
    for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
            {
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
                vector2[1][0][1][index] = temp_phi1x[ index ];
                vector2[2][0][1][index] = temp_phi2x[ index ];
                vector2[3][0][1][index] = temp_phi3x[ index ];
                vector2[1][0][2][index] = temp_phi1y[ index ];
                vector2[2][0][2][index] = temp_phi2y[ index ];
                vector2[3][0][2][index] = temp_phi3y[ index ];
                vector2[1][0][3][index] = temp_phi1z[ index ];
                vector2[2][0][3][index] = temp_phi2z[ index ];
                vector2[3][0][3][index] = temp_phi3z[ index ];
                
                vector2[0][1][1][index] = temp_phi1x[ index ];
                vector2[0][2][1][index] = temp_phi2x[ index ];
                vector2[0][3][1][index] = temp_phi3x[ index ];
                vector2[0][1][2][index] = temp_phi1y[ index ];
                vector2[0][2][2][index] = temp_phi2y[ index ];
                vector2[0][3][2][index] = temp_phi3y[ index ];
                vector2[0][1][3][index] = temp_phi1z[ index ];
                vector2[0][2][3][index] = temp_phi2z[ index ];
                vector2[0][3][3][index] = temp_phi3z[ index ];
                
                vector2[1][1][1][index] = temp_phi11x[ index ];
                vector2[2][2][1][index] = temp_phi22x[ index ];
                vector2[3][3][1][index] = temp_phi33x[ index ];
                vector2[1][1][2][index] = temp_phi11y[ index ];
                vector2[2][2][2][index] = temp_phi22y[ index ];
                vector2[3][3][2][index] = temp_phi33y[ index ];
                vector2[1][1][3][index] = temp_phi11z[ index ];
                vector2[2][2][3][index] = temp_phi22z[ index ];
                vector2[3][3][3][index] = temp_phi33z[ index ];
                
                vector2[1][2][1][index] = temp_phi12x[ index ];
                vector2[1][3][1][index] = temp_phi13x[ index ];
                vector2[2][3][1][index] = temp_phi23x[ index ];
                vector2[1][2][2][index] = temp_phi12y[ index ];
                vector2[1][3][2][index] = temp_phi13y[ index ];
                vector2[2][3][2][index] = temp_phi23y[ index ];
                vector2[1][2][3][index] = temp_phi12z[ index ];
                vector2[1][3][3][index] = temp_phi13z[ index ];
                vector2[2][3][3][index] = temp_phi23z[ index ];
                
                vector2[2][1][1][index] = temp_phi12x[ index ];
                vector2[3][1][1][index] = temp_phi13x[ index ];
                vector2[3][2][1][index] = temp_phi23x[ index ];
                vector2[2][1][2][index] = temp_phi12y[ index ];
                vector2[3][1][2][index] = temp_phi13y[ index ];
                vector2[3][2][2][index] = temp_phi23y[ index ];
                vector2[2][1][3][index] = temp_phi12z[ index ];
                vector2[3][1][3][index] = temp_phi13z[ index ];
                vector2[3][2][3][index] = temp_phi23z[ index ];
                
			}
		}
	}
}
