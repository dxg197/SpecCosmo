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

void SpecDeriv_Scalar_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_Write_Scalar2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 );
void SpecDeriv_SpecGalMethod_Scalar2( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Scalar_Derivative2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 )
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
                temp_phi0[index] = scalar2[0][0][index];
                scalar2[1][0][index] = 0.0;
                scalar2[2][0][index] = 0.0;
                scalar2[3][0][index] = 0.0;
                
                scalar2[0][1][index] = 0.0;
                scalar2[0][2][index] = 0.0;
                scalar2[0][3][index] = 0.0;
                
                scalar2[1][1][index] = 0.0;
                scalar2[2][2][index] = 0.0;
                scalar2[3][3][index] = 0.0;
                
                scalar2[1][2][index] = 0.0;
                scalar2[1][3][index] = 0.0;
                scalar2[2][3][index] = 0.0;
                
                scalar2[2][1][index] = 0.0;
                scalar2[3][1][index] = 0.0;
                scalar2[3][2][index] = 0.0;
			}
		}
	}
	
	SpecDeriv_SpecGalMethod_Scalar2( CCTK_PASS_CTOC );
	SpecDeriv_Write_Scalar2( CCTK_PASS_CTOC, scalar2 );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_SpecGalMethod_Scalar2( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr, handle, rank;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    fftw_plan plan_forward;
    fftw_plan plan_backward1;
    fftw_plan plan_backward2;
    fftw_plan plan_backward3;
    fftw_plan plan_backward11;
    fftw_plan plan_backward22;
    fftw_plan plan_backward33;
    fftw_plan plan_backward12;
    fftw_plan plan_backward13;
    fftw_plan plan_backward23;
    CCTK_REAL Fourier_3D_kwaveN1;
    CCTK_REAL Fourier_3D_kwaveN2;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementi, waveElementj, waveElementk;
    CCTK_REAL lx, ly, lz, xmin, xmax, ymin, ymax, zmin, zmax;
    fftw_complex *phi_in; 
    fftw_complex *Derivative_SGM1;
    fftw_complex *Derivative_SGM2;
    fftw_complex *Derivative_SGM3;  
    fftw_complex *Derivative_SGM11;
    fftw_complex *Derivative_SGM22;
    fftw_complex *Derivative_SGM33;
    fftw_complex *Derivative_SGM12;
    fftw_complex *Derivative_SGM13;
    fftw_complex *Derivative_SGM23;
    fftw_complex *phi_out; 
    fftw_complex *Derivative_SGM1_out;
    fftw_complex *Derivative_SGM2_out;
    fftw_complex *Derivative_SGM3_out;  
    fftw_complex *Derivative_SGM11_out;
    fftw_complex *Derivative_SGM22_out;
    fftw_complex *Derivative_SGM33_out;
    fftw_complex *Derivative_SGM12_out;
    fftw_complex *Derivative_SGM13_out;
    fftw_complex *Derivative_SGM23_out;
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
       CCTK_VARIABLE_REAL, temp_phi0,
       CCTK_VARIABLE_REAL, temp_phi0);
    assert (! ierr);
    
    //  Begin Spectral-Galerkin Method
    //  Create phi, compute Forward FT 
    
    phi_in = fftw_alloc_complex(alloc_local);
    phi_out = fftw_alloc_complex(alloc_local);
      
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
              phi_in[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0[i + j*sizex + k*sizex*sizey];
              phi_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
            }
        }
    }
    
    plan_forward = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_in,phi_out,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    fftw_execute(plan_forward);
    
    //  Multiply i*k factor to FFTWave, create IK_FFTWavePadded
    Derivative_SGM1 = fftw_alloc_complex(alloc_local);
    Derivative_SGM2 = fftw_alloc_complex(alloc_local);
    Derivative_SGM3 = fftw_alloc_complex(alloc_local);
    Derivative_SGM11 = fftw_alloc_complex(alloc_local);
    Derivative_SGM22 = fftw_alloc_complex(alloc_local);
    Derivative_SGM33 = fftw_alloc_complex(alloc_local);
    Derivative_SGM12 = fftw_alloc_complex(alloc_local);
    Derivative_SGM13 = fftw_alloc_complex(alloc_local);
    Derivative_SGM23 = fftw_alloc_complex(alloc_local);
    Derivative_SGM1_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM11_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM22_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM33_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM12_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM13_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM23_out = fftw_alloc_complex(alloc_local);
    
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
                
                Derivative_SGM1[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER * 
                                   Fourier_3D_kwaveN1 * 
                                   phi_out[i + j*sizex + k*sizex*sizey]; 
                  
                Derivative_SGM2[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER * 
                                   Fourier_3D_kwaveN2 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];
                  
                Derivative_SGM3[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER * 
                                   Fourier_3D_kwaveN3 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];
                  
                Derivative_SGM11[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN1 *
                                   Fourier_3D_kwaveN1 *
                                   phi_out[i + j*sizex + k*sizex*sizey]; 
                
                Derivative_SGM22[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN2 * 
                                   Fourier_3D_kwaveN2 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];
                  
                Derivative_SGM33[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN3 * 
                                   Fourier_3D_kwaveN3 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];
                
                Derivative_SGM12[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN1 * 
                                   Fourier_3D_kwaveN2 * 
                                   phi_out[i + j*sizex + k*sizex*sizey]; 
                
                Derivative_SGM13[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN1 * 
                                   Fourier_3D_kwaveN3 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];
                  
                Derivative_SGM23[i + j*sizex + k*sizex*sizey] = -1.0  * 
                                   Fourier_3D_kwaveN2 * 
                                   Fourier_3D_kwaveN3 * 
                                   phi_out[i + j*sizex + k*sizex*sizey];    
                                   
                Derivative_SGM1_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                  
                Derivative_SGM2_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                  
                Derivative_SGM3_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                  
                Derivative_SGM11_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                
                Derivative_SGM22_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                  
                Derivative_SGM33_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                Derivative_SGM12_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                
                Derivative_SGM13_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                  
                Derivative_SGM23_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;    
                
                // Dealias using truncation (2/3 rule)
                 
                 if (( waveElementi > (2.0*sizex/3.0) ) || ( waveElementi < (sizex/3.0) )) {
                 	Derivative_SGM1[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM11[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementj > (2.0*sizey/3.0) ) || ( waveElementj < (sizey/3.0) )) {
                 	Derivative_SGM2[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM22[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM12[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementk > (2.0*sizez/3.0) ) || ( waveElementk < (sizez/3.0) )) {
                 	Derivative_SGM3[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM33[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM13[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM23[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
            }

        }

    }
   
    //  Create Derivative_SGM, compute Inverse FT 
    
    plan_backward1 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1,Derivative_SGM1_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2,Derivative_SGM2_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3,Derivative_SGM3_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward11 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM11,Derivative_SGM11_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward22 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM22,Derivative_SGM22_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward33 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM33,Derivative_SGM33_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward12 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM12,Derivative_SGM12_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward13 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM13,Derivative_SGM13_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward23 = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM23,Derivative_SGM23_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);

    fftw_execute(plan_backward1);
    fftw_execute(plan_backward2);
    fftw_execute(plan_backward3);
    fftw_execute(plan_backward11);
    fftw_execute(plan_backward22);
    fftw_execute(plan_backward33);
    fftw_execute(plan_backward12);
    fftw_execute(plan_backward13);
    fftw_execute(plan_backward23);
    
    //   Store derivative into PHI 
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM1_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM2_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM3_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				
				temp_phi11[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM11_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi22[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM22_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi33[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM33_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				
				temp_phi12[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM12_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi13[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM13_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi23[i + j*sizex + k*sizex*sizey] = creal(Derivative_SGM23_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
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
       CCTK_VARIABLE_REAL, temp_phi1,
       CCTK_VARIABLE_REAL, temp_phi1);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2,
       CCTK_VARIABLE_REAL, temp_phi2);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3,
       CCTK_VARIABLE_REAL, temp_phi3);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi11,
       CCTK_VARIABLE_REAL, temp_phi11);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi22,
       CCTK_VARIABLE_REAL, temp_phi22);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi33,
       CCTK_VARIABLE_REAL, temp_phi33);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi12,
       CCTK_VARIABLE_REAL, temp_phi12);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi13,
       CCTK_VARIABLE_REAL, temp_phi13);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi23,
       CCTK_VARIABLE_REAL, temp_phi23);
    assert (! ierr);
    
    //  Destroy/Free array pointers 

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward1);
    fftw_destroy_plan(plan_backward2);
    fftw_destroy_plan(plan_backward3);
    fftw_destroy_plan(plan_backward11);
    fftw_destroy_plan(plan_backward22);
    fftw_destroy_plan(plan_backward33);
    fftw_destroy_plan(plan_backward12);
    fftw_destroy_plan(plan_backward13);
    fftw_destroy_plan(plan_backward23);
    fftw_free(phi_in);
    fftw_free(Derivative_SGM1);
    fftw_free(Derivative_SGM2);
    fftw_free(Derivative_SGM3);
    fftw_free(Derivative_SGM11);
    fftw_free(Derivative_SGM22);
    fftw_free(Derivative_SGM33);
    fftw_free(Derivative_SGM12);
    fftw_free(Derivative_SGM13);
    fftw_free(Derivative_SGM23);
    fftw_free(phi_out);
    fftw_free(Derivative_SGM1_out);
    fftw_free(Derivative_SGM2_out);
    fftw_free(Derivative_SGM3_out);
    fftw_free(Derivative_SGM11_out);
    fftw_free(Derivative_SGM22_out);
    fftw_free(Derivative_SGM33_out);
    fftw_free(Derivative_SGM12_out);
    fftw_free(Derivative_SGM13_out);
    fftw_free(Derivative_SGM23_out);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Write_Scalar2( CCTK_ARGUMENTS, CCTK_REAL ***scalar2 )
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
                scalar2[1][0][index] = 0.0;
                scalar2[2][0][index] = 0.0;
                scalar2[3][0][index] = 0.0;
                
                scalar2[0][1][index] = 0.0;
                scalar2[0][2][index] = 0.0;
                scalar2[0][3][index] = 0.0;
                
                scalar2[1][1][index] = 0.0;
                scalar2[2][2][index] = 0.0;
                scalar2[3][3][index] = 0.0;
                
                scalar2[1][2][index] = 0.0;
                scalar2[1][3][index] = 0.0;
                scalar2[2][3][index] = 0.0;
                
                scalar2[2][1][index] = 0.0;
                scalar2[3][1][index] = 0.0;
                scalar2[3][2][index] = 0.0;
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
                scalar2[1][0][index] = 0.0;
                scalar2[2][0][index] = 0.0;
                scalar2[3][0][index] = 0.0;
                
                scalar2[0][1][index] = 0.0;
                scalar2[0][2][index] = 0.0;
                scalar2[0][3][index] = 0.0;
                
                scalar2[1][1][index] = 0.0;
                scalar2[2][2][index] = 0.0;
                scalar2[3][3][index] = 0.0;
                
                scalar2[1][2][index] = 0.0;
                scalar2[1][3][index] = 0.0;
                scalar2[2][3][index] = 0.0;
                
                scalar2[2][1][index] = 0.0;
                scalar2[3][1][index] = 0.0;
                scalar2[3][2][index] = 0.0;
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
                scalar2[1][0][index] = temp_phi1[ index ];
                scalar2[2][0][index] = temp_phi2[ index ];
                scalar2[3][0][index] = temp_phi3[ index ];
                
                scalar2[0][1][index] = temp_phi1[ index ];
                scalar2[0][2][index] = temp_phi2[ index ];
                scalar2[0][3][index] = temp_phi3[ index ];
                
                scalar2[1][1][index] = temp_phi11[ index ];
                scalar2[2][2][index] = temp_phi22[ index ];
                scalar2[3][3][index] = temp_phi33[ index ];
                
                scalar2[1][2][index] = temp_phi12[ index ];
                scalar2[1][3][index] = temp_phi13[ index ];
                scalar2[2][3][index] = temp_phi23[ index ];
                
                scalar2[2][1][index] = temp_phi12[ index ];
                scalar2[3][1][index] = temp_phi13[ index ];
                scalar2[3][2][index] = temp_phi23[ index ];
			}
		}
	}
}
