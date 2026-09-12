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

void Clean_Read_Scalar( CCTK_ARGUMENTS, CCTK_REAL *scalar );
void Clean_Write_Scalar( CCTK_ARGUMENTS, CCTK_REAL *scalar );
void Clean_Scalar( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void Clean_Read_Scalar( CCTK_ARGUMENTS, CCTK_REAL *scalar )
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
                temp_phi0[ index ] = scalar[index];
			}
		}
	}
	
	Clean_Scalar( CCTK_PASS_CTOC );
	Clean_Write_Scalar( CCTK_PASS_CTOC, scalar );
}

/*******************************************************/
/*******************************************************/
/*                                                     */             
/*                                                     */
/*******************************************************/
/*******************************************************/
void Clean_Scalar( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr, rank;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    fftw_plan plan_forward;
    fftw_plan plan_backward;
    fftw_complex *phi_in, *phi_out, *phi_out2; 
    fftw_complex ZERO_COMPLEXNO = (0.0,0.0);
    MPI_Comm comm;
  
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
      info[0].dst.lbnd        = cctk_nghostzones[0];
      info[1].dst.lbnd        = cctk_nghostzones[1];
      info[2].dst.lbnd        = local_0_start+cctk_nghostzones[2];
      info[0].dst.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].dst.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].dst.lsh         = local_n0;
      info[0].dst.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].dst.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
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
       CCTK_VARIABLE_REAL, temp_phi0,
       CCTK_VARIABLE_REAL, temp_phi0);
    assert (! ierr);
    
    //  Begin Spectral-Galerkin Method
    //  Create phi, compute Forward FT 
    
    phi_in = fftw_alloc_complex(alloc_local);
    phi_out = fftw_alloc_complex(alloc_local);
    phi_out2 = fftw_alloc_complex(alloc_local);
    
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
              phi_in[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0[i + j*sizex + k*sizex*sizey];
              phi_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_out2[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
            }
        }
    }
    
    plan_forward = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_in,phi_out,MPI_COMM_WORLD,FFTW_FORWARD, FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    fftw_execute(plan_forward);

	for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {   
                // Dealias using truncation (2/3 rule)
                 
                 if (( i > (2.0*sizex/3.0) ) || ( i < (sizex/3.0) )) {
                 	phi_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                 }
                 
                 if (( j > (2.0*sizey/3.0) ) || ( j < (sizey/3.0) )) {
                 	phi_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                 }
                 
                 if (( k > (2.0*sizez/3.0) ) || ( k < (sizez/3.0) )) {
                 	phi_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                 }
            }

        }

    }     
    //  Create Derivative_SGM, compute Inverse FT 
    
    plan_backward = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_out,phi_out2,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    
    fftw_execute(plan_backward);
     
    //   Store derivative of PHI 
    
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1[i + j*sizex + k*sizex*sizey] = creal(phi_out2[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
			}
		}
	}
    
    for (d=0; d<3; ++d) {
      /* Source array: FFTW layout */
      info[d].src.gsh         = cctk_gsh[d];
      info[0].src.lbnd        = cctk_nghostzones[0];
      info[1].src.lbnd        = cctk_nghostzones[1];
      info[2].src.lbnd        = local_0_start+cctk_nghostzones[2];
      info[0].src.lsh         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].src.lsh         = cctk_gsh[1]-2*cctk_nghostzones[1];
      info[2].src.lsh         = local_n0;
      info[0].src.ash         = cctk_gsh[0]-2*cctk_nghostzones[0];
      info[1].src.ash         = cctk_gsh[1]-2*cctk_nghostzones[1];
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
       CCTK_VARIABLE_REAL, temp_phi1,
       CCTK_VARIABLE_REAL, temp_phi1);
    assert (! ierr);
	
    //  Destroy/Free array pointers 

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(phi_in);  
    fftw_free(phi_out);  
    fftw_free(phi_out2);  
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void Clean_Write_Scalar( CCTK_ARGUMENTS, CCTK_REAL *scalar )
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
                scalar[index] = temp_phi1[ index ];
			}
		}
	}
}
