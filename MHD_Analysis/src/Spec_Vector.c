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

void Spec_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void Spec_Write_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector );
void Spec_Analyze_Vector( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void Spec_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector )
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
                
                vector[1][1][index] = 0.0;
                vector[1][2][index] = 0.0;
                vector[1][3][index] = 0.0;
                vector[2][1][index] = 0.0;
                vector[2][2][index] = 0.0;
                vector[2][3][index] = 0.0;
                vector[3][1][index] = 0.0;
                vector[3][2][index] = 0.0;
                vector[3][3][index] = 0.0;
			}
		}
	}
	Spec_Analyze_Vector( CCTK_PASS_CTOC );
	Spec_Write_Vector( CCTK_PASS_CTOC, vector );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void Spec_Analyze_Vector( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr, rank;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    fftw_plan plan_forwardx;
    fftw_plan plan_forwardy;
    fftw_plan plan_forwardz;
    fftw_complex *phi_inx, *phi_iny, *phi_inz;
    fftw_complex *phi_outx, *phi_outy, *phi_outz;
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
    
    assert (cctk_dim <= 3);
    for (d=0; d<cctk_dim; ++d) {
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
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0x,
       CCTK_VARIABLE_REAL, temp_phi0x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0y,
       CCTK_VARIABLE_REAL, temp_phi0y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0z,
       CCTK_VARIABLE_REAL, temp_phi0z);
    assert (! ierr);
    
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
    
    //  compute Forward FT
    
    plan_forwardx = fftw_mpi_plan_dft_3d(sizez,sizey,sizex,phi_inx,phi_outx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardy = fftw_mpi_plan_dft_3d(sizez,sizey,sizex,phi_iny,phi_outy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardz = fftw_mpi_plan_dft_3d(sizez,sizey,sizex,phi_inz,phi_outz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    
    fftw_execute(plan_forwardx);
    fftw_execute(plan_forwardy);
    fftw_execute(plan_forwardz);
	
	for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1x[i + j*sizex + k*sizex*sizey] = creal(phi_outx[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
				temp_phi1y[i + j*sizex + k*sizex*sizey] = creal(phi_outy[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
				temp_phi1z[i + j*sizex + k*sizex*sizey] = creal(phi_outz[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
				temp_phi2x[i + j*sizex + k*sizex*sizey] = cimag(phi_outx[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
				temp_phi2y[i + j*sizex + k*sizex*sizey] = cimag(phi_outy[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
				temp_phi2z[i + j*sizex + k*sizex*sizey] = cimag(phi_outz[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez));
			    temp_phi3x[i + j*sizex + k*sizex*sizey] = (pow(creal(phi_outx[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0) +
				                                           pow(cimag(phi_outx[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0));
				temp_phi3y[i + j*sizex + k*sizex*sizey] = (pow(creal(phi_outy[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0) +
				                                           pow(cimag(phi_outy[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0));
				temp_phi3z[i + j*sizex + k*sizex*sizey] = (pow(creal(phi_outz[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0) +
				                                           pow(cimag(phi_outz[i + j*sizex + k*sizex*sizey]/(sizex*sizey*sizez)),2.0));
			}
		}
	}
	
	assert (cctk_dim <= 3);  
	for (d=0; d<cctk_dim; ++d) {
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
       CCTK_VARIABLE_REAL, temp_phi1x,
       CCTK_VARIABLE_REAL, temp_phi1x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1y,
       CCTK_VARIABLE_REAL, temp_phi1y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1z,
       CCTK_VARIABLE_REAL, temp_phi1z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2x,
       CCTK_VARIABLE_REAL, temp_phi2x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2y,
       CCTK_VARIABLE_REAL, temp_phi2y);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2z,
       CCTK_VARIABLE_REAL, temp_phi2z);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3x,
       CCTK_VARIABLE_REAL, temp_phi3x);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3y,
       CCTK_VARIABLE_REAL, temp_phi3y);
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
    fftw_free(phi_inx);
    fftw_free(phi_iny);
    fftw_free(phi_inz);
    fftw_free(phi_outx);
    fftw_free(phi_outy);
    fftw_free(phi_outz);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void Spec_Write_Vector( CCTK_ARGUMENTS, CCTK_REAL ***vector )
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
    //iend   = cctk_lsh[0];
	//jend   = cctk_lsh[1];
	//kend   = cctk_lsh[2];
    
    for(k=kstart; k < kend; k++)
	{
		for(j=jstart; j < jend; j++)
		{
			for(i=istart; i < iend; i++)
            {
				index = CCTK_GFINDEX3D( cctkGH, i, j, k );
				vector[1][1][index] = 0.0;
                vector[1][2][index] = 0.0;
                vector[1][3][index] = 0.0;
                vector[2][1][index] = 0.0;
                vector[2][2][index] = 0.0;
                vector[2][3][index] = 0.0;
                vector[3][1][index] = 0.0;
                vector[3][2][index] = 0.0;
                vector[3][3][index] = 0.0;
                
                if (!isnan(temp_phi1x[ index ])) {
                vector[1][1][index] = temp_phi1x[ index ];
                }
                if (!isnan(temp_phi1y[ index ])) {
                vector[1][2][index] = temp_phi1y[ index ];
                }
                if (!isnan(temp_phi1z[ index ])) {
                vector[1][3][index] = temp_phi1z[ index ];
                }
                if (!isnan(temp_phi2x[ index ])) {
                vector[2][1][index] = temp_phi2x[ index ];
                }
                if (!isnan(temp_phi2y[ index ])) {
                vector[2][2][index] = temp_phi2y[ index ];
                }
                if (!isnan(temp_phi2z[ index ])) {
                vector[2][3][index] = temp_phi2z[ index ];
                }
                if (!isnan(temp_phi3x[ index ])) {
                vector[3][1][index] = temp_phi3x[ index ];
                }
                if (!isnan(temp_phi3y[ index ])) {
                vector[3][2][index] = temp_phi3y[ index ];
                }
                if (!isnan(temp_phi3z[ index ])) {
                vector[3][3][index] = temp_phi3z[ index ];
                }
			}
		}
	}
}
