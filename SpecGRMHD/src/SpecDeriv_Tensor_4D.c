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

void SpecDeriv_Tensor_Derivative_4D( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D );
void SpecDeriv_Write_Tensor_4D( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D );
void SpecDeriv_SpecGalMethod_Tensor_4D( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative_4D( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D )
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
				
				for(m=1; m < 4; m++) {
                	for(n=0; n < 4; n++) {
                		for(p=0; p < 4; p++) {
                			tensor_4D[m][n][p][index] = 0.0;
                		}
                	}
                }
                
				temp_phi0tt[ index ] = tensor_4D[0][0][0][index];
				temp_phi0tx[ index ] = tensor_4D[0][0][1][index];
                temp_phi0ty[ index ] = tensor_4D[0][0][2][index];
                temp_phi0tz[ index ] = tensor_4D[0][0][3][index];
                temp_phi0xx[ index ] = tensor_4D[0][1][1][index];
                temp_phi0yy[ index ] = tensor_4D[0][2][2][index];
                temp_phi0zz[ index ] = tensor_4D[0][3][3][index];
                temp_phi0xy[ index ] = tensor_4D[0][1][2][index];
                temp_phi0xz[ index ] = tensor_4D[0][1][3][index];
                temp_phi0yz[ index ] = tensor_4D[0][2][3][index];
			}
		}
	}
	SpecDeriv_SpecGalMethod_Tensor_4D( CCTK_PASS_CTOC );
	SpecDeriv_Write_Tensor_4D( CCTK_PASS_CTOC, tensor_4D );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_SpecGalMethod_Tensor_4D( CCTK_ARGUMENTS )
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    //  Declare and initialize variables
    CCTK_INT i, j, k, d, ierr, handle, rank;
    ptrdiff_t sizex, sizey, sizez;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    fftw_plan plan_forwardtt;
    fftw_plan plan_forwardtx;
    fftw_plan plan_forwardty;
    fftw_plan plan_forwardtz;
    fftw_plan plan_forwardxx;
    fftw_plan plan_forwardyy;
    fftw_plan plan_forwardzz;
    fftw_plan plan_forwardxy;
    fftw_plan plan_forwardxz;
    fftw_plan plan_forwardyz;
    fftw_plan plan_backward1tt;
    fftw_plan plan_backward2tt;
    fftw_plan plan_backward3tt;
    fftw_plan plan_backward1tx;
    fftw_plan plan_backward2tx;
    fftw_plan plan_backward3tx;
    fftw_plan plan_backward1ty;
    fftw_plan plan_backward2ty;
    fftw_plan plan_backward3ty;
    fftw_plan plan_backward1tz;
    fftw_plan plan_backward2tz;
    fftw_plan plan_backward3tz;
    fftw_plan plan_backward1xx;
    fftw_plan plan_backward2xx;
    fftw_plan plan_backward3xx;
    fftw_plan plan_backward1yy;
    fftw_plan plan_backward2yy;
    fftw_plan plan_backward3yy;
    fftw_plan plan_backward1zz;
    fftw_plan plan_backward2zz;
    fftw_plan plan_backward3zz;
    fftw_plan plan_backward1xy;
    fftw_plan plan_backward2xy;
    fftw_plan plan_backward3xy;
    fftw_plan plan_backward1xz;
    fftw_plan plan_backward2xz;
    fftw_plan plan_backward3xz;
    fftw_plan plan_backward1yz;
    fftw_plan plan_backward2yz;
    fftw_plan plan_backward3yz;
    CCTK_REAL Fourier_3D_kwaveN1;
    CCTK_REAL Fourier_3D_kwaveN2;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementi, waveElementj, waveElementk;
    CCTK_REAL lx, ly, lz, xmin, xmax, ymin, ymax, zmin, zmax;
    fftw_complex *Derivative_SGM1xx, *Derivative_SGM1yy, *Derivative_SGM1zz;
    fftw_complex *Derivative_SGM2xx, *Derivative_SGM2yy, *Derivative_SGM2zz;
    fftw_complex *Derivative_SGM3xx, *Derivative_SGM3yy, *Derivative_SGM3zz;
    fftw_complex *Derivative_SGM1xy, *Derivative_SGM1yz, *Derivative_SGM1xz;
    fftw_complex *Derivative_SGM2xy, *Derivative_SGM2yz, *Derivative_SGM2xz;
    fftw_complex *Derivative_SGM3xy, *Derivative_SGM3yz, *Derivative_SGM3xz;
    fftw_complex *Derivative_SGM1tt, *Derivative_SGM1ty;
    fftw_complex *Derivative_SGM2tt, *Derivative_SGM2ty;
    fftw_complex *Derivative_SGM3tt, *Derivative_SGM3ty;
    fftw_complex *Derivative_SGM1tx, *Derivative_SGM1tz;
    fftw_complex *Derivative_SGM2tx, *Derivative_SGM2tz;
    fftw_complex *Derivative_SGM3tx, *Derivative_SGM3tz;
    fftw_complex *phi_intt, *phi_intx, *phi_inty, *phi_intz;
    fftw_complex *phi_inxx, *phi_inyy, *phi_inzz;
    fftw_complex *phi_inxy, *phi_inxz, *phi_inyz;
    fftw_complex *Derivative_SGM1xx_out, *Derivative_SGM1yy_out, *Derivative_SGM1zz_out;
    fftw_complex *Derivative_SGM2xx_out, *Derivative_SGM2yy_out, *Derivative_SGM2zz_out;
    fftw_complex *Derivative_SGM3xx_out, *Derivative_SGM3yy_out, *Derivative_SGM3zz_out;
    fftw_complex *Derivative_SGM1xy_out, *Derivative_SGM1yz_out, *Derivative_SGM1xz_out;
    fftw_complex *Derivative_SGM2xy_out, *Derivative_SGM2yz_out, *Derivative_SGM2xz_out;
    fftw_complex *Derivative_SGM3xy_out, *Derivative_SGM3yz_out, *Derivative_SGM3xz_out;
    fftw_complex *Derivative_SGM1tt_out, *Derivative_SGM1ty_out;
    fftw_complex *Derivative_SGM2tt_out, *Derivative_SGM2ty_out;
    fftw_complex *Derivative_SGM3tt_out, *Derivative_SGM3ty_out;
    fftw_complex *Derivative_SGM1tx_out, *Derivative_SGM1tz_out;
    fftw_complex *Derivative_SGM2tx_out, *Derivative_SGM2tz_out;
    fftw_complex *Derivative_SGM3tx_out, *Derivative_SGM3tz_out;
    fftw_complex *phi_outtt, *phi_outtx, *phi_outty, *phi_outtz;
    fftw_complex *phi_outxx, *phi_outyy, *phi_outzz;
    fftw_complex *phi_outxy, *phi_outxz, *phi_outyz;
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
       CCTK_VARIABLE_REAL, temp_phi0tt,
       CCTK_VARIABLE_REAL, temp_phi0tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0tx,
       CCTK_VARIABLE_REAL, temp_phi0tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0ty,
       CCTK_VARIABLE_REAL, temp_phi0ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0tz,
       CCTK_VARIABLE_REAL, temp_phi0tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0xx,
       CCTK_VARIABLE_REAL, temp_phi0xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0yy,
       CCTK_VARIABLE_REAL, temp_phi0yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0zz,
       CCTK_VARIABLE_REAL, temp_phi0zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0xy,
       CCTK_VARIABLE_REAL, temp_phi0xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0xz,
       CCTK_VARIABLE_REAL, temp_phi0xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info1, -1,
       CCTK_VARIABLE_REAL, temp_phi0yz,
       CCTK_VARIABLE_REAL, temp_phi0yz);
    assert (! ierr);
    
    //  Begin Spectral-Galerkin Method
    //  Create phi, compute Forward FT 
    
    phi_intt = fftw_alloc_complex(alloc_local);
    phi_intx = fftw_alloc_complex(alloc_local);
    phi_inty = fftw_alloc_complex(alloc_local);
    phi_intz = fftw_alloc_complex(alloc_local);
    phi_inxx = fftw_alloc_complex(alloc_local);
    phi_inyy = fftw_alloc_complex(alloc_local);
    phi_inzz = fftw_alloc_complex(alloc_local);
    phi_inxy = fftw_alloc_complex(alloc_local);
    phi_inxz = fftw_alloc_complex(alloc_local);
    phi_inyz = fftw_alloc_complex(alloc_local);
    
    phi_outtt = fftw_alloc_complex(alloc_local);
    phi_outtx = fftw_alloc_complex(alloc_local);
    phi_outty = fftw_alloc_complex(alloc_local);
    phi_outtz = fftw_alloc_complex(alloc_local);
    phi_outxx = fftw_alloc_complex(alloc_local);
    phi_outyy = fftw_alloc_complex(alloc_local);
    phi_outzz = fftw_alloc_complex(alloc_local);
    phi_outxy = fftw_alloc_complex(alloc_local);
    phi_outxz = fftw_alloc_complex(alloc_local);
    phi_outyz = fftw_alloc_complex(alloc_local);
    
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
              phi_intt[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0tt[i + j*sizex + k*sizex*sizey];
              phi_intx[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0tx[i + j*sizex + k*sizex*sizey];
              phi_inty[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0ty[i + j*sizex + k*sizex*sizey];
              phi_intz[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0tz[i + j*sizex + k*sizex*sizey];
              phi_inxx[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0xx[i + j*sizex + k*sizex*sizey];
              phi_inyy[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0yy[i + j*sizex + k*sizex*sizey];
              phi_inzz[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0zz[i + j*sizex + k*sizex*sizey];
              phi_inxy[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0xy[i + j*sizex + k*sizex*sizey];
              phi_inxz[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0xz[i + j*sizex + k*sizex*sizey];
              phi_inyz[i + j*sizex + k*sizex*sizey] = (fftw_complex)temp_phi0yz[i + j*sizex + k*sizex*sizey];
              
              phi_outtt[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outtx[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outty[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outtz[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outxx[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outyy[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outzz[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outxy[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outxz[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
              phi_outyz[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
            }

        }

    }
    
    plan_forwardtt = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_intt,phi_outtt,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardtx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_intx,phi_outtx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardty = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inty,phi_outty,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardtz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_intz,phi_outtz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardxx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inxx,phi_outxx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardyy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inyy,phi_outyy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardzz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inzz,phi_outzz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardxy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inxy,phi_outxy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardxz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inxz,phi_outxz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_forwardyz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,phi_inyz,phi_outyz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    
    fftw_execute(plan_forwardtt);
    fftw_execute(plan_forwardtx);
    fftw_execute(plan_forwardty);
    fftw_execute(plan_forwardtz);
    fftw_execute(plan_forwardxx);
    fftw_execute(plan_forwardyy);
    fftw_execute(plan_forwardzz);
    fftw_execute(plan_forwardxy);
    fftw_execute(plan_forwardxz);
    fftw_execute(plan_forwardyz);
    
    //  Multiply i*k factor to FFTWave, create IK_FFTWavePadded
    Derivative_SGM1tt = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tt = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tt = fftw_alloc_complex(alloc_local);
    Derivative_SGM1tx = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tx = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tx = fftw_alloc_complex(alloc_local);
    Derivative_SGM1ty = fftw_alloc_complex(alloc_local);
    Derivative_SGM2ty = fftw_alloc_complex(alloc_local);
    Derivative_SGM3ty = fftw_alloc_complex(alloc_local);
    Derivative_SGM1tz = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tz = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tz = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xx = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xx = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xx = fftw_alloc_complex(alloc_local);
    Derivative_SGM1yy = fftw_alloc_complex(alloc_local);
    Derivative_SGM2yy = fftw_alloc_complex(alloc_local);
    Derivative_SGM3yy = fftw_alloc_complex(alloc_local);
    Derivative_SGM1zz = fftw_alloc_complex(alloc_local);
    Derivative_SGM2zz = fftw_alloc_complex(alloc_local);
    Derivative_SGM3zz = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xy = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xy = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xy = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xz = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xz = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xz = fftw_alloc_complex(alloc_local);
    Derivative_SGM1yz = fftw_alloc_complex(alloc_local);
    Derivative_SGM2yz = fftw_alloc_complex(alloc_local);
    Derivative_SGM3yz = fftw_alloc_complex(alloc_local);  
    
    Derivative_SGM1tt_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tt_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tt_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1tx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1ty_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2ty_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3ty_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1tz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2tz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3tz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xx_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1yy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2yy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3yy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1zz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2zz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3zz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xy_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1xz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2xz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3xz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM1yz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM2yz_out = fftw_alloc_complex(alloc_local);
    Derivative_SGM3yz_out = fftw_alloc_complex(alloc_local);  
    
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
                
                Derivative_SGM1tt[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outtt[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2tt[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outtt[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3tt[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outtt[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1tx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outtx[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2tx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outtx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3tx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outtx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1ty[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outty[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2ty[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outty[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3ty[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outty[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1tz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outtz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2tz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outtz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3tz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outtz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1xx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outxx[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2xx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outxx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3xx[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outxx[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1yy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outyy[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2yy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outyy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3yy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outyy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1zz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outzz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2zz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outzz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3zz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outzz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1xy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outxy[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2xy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outxy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3xy[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outxy[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1xz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outxz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2xz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outxz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3xz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outxz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM1yz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN1 
                                            * phi_outyz[i + j*sizex + k*sizex*sizey]; 
                Derivative_SGM2yz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN2 
                                            * phi_outyz[i + j*sizex + k*sizex*sizey];
                Derivative_SGM3yz[i + j*sizex + k*sizex*sizey] = IMAG_NUMBER 
                                            * Fourier_3D_kwaveN3 
                                            * phi_outyz[i + j*sizex + k*sizex*sizey];
                                            
                
                
                Derivative_SGM1tt_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2tt_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3tt_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1tx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2tx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3tx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1ty_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2ty_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3ty_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1tz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2tz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3tz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1xx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2xx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3xx_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1yy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2yy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3yy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1zz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2zz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3zz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1xy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2xy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3xy_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1xz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2xz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3xz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM1yz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO; 
                Derivative_SGM2yz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                Derivative_SGM3yz_out[i + j*sizex + k*sizex*sizey] = ZERO_COMPLEXNO;
                
                // Dealias using truncation (2/3 rule)
                 
                 if (( waveElementi > (2.0*sizex/3.0) ) || ( waveElementi < (sizex/3.0) )) {
                 	Derivative_SGM1tt[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1tx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1ty[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1tz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1xx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1yy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1zz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1xy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1xz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM1yz[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementj > (2.0*sizey/3.0) ) || ( waveElementj < (sizey/3.0) )) {
                 	Derivative_SGM2tt[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2tx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2ty[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2tz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2xx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2yy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2zz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2xy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2xz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM2yz[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
                 if (( waveElementk > (2.0*sizez/3.0) ) || ( waveElementk < (sizez/3.0) )) {
                 	Derivative_SGM3tt[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3tx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3ty[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3tz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3xx[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3yy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3zz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3xy[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3xz[i + j*sizex + k*sizex*sizey] = 0.0;
                 	Derivative_SGM3yz[i + j*sizex + k*sizex*sizey] = 0.0;
                 }
                 
            }

        }

    }
    

    //  Create Derivative_SGM, compute Inverse FT 
    
    plan_backward1tt = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1tt,Derivative_SGM1tt_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2tt = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2tt,Derivative_SGM2tt_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3tt = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3tt,Derivative_SGM3tt_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1tx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1tx,Derivative_SGM1tx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2tx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2tx,Derivative_SGM2tx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3tx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3tx,Derivative_SGM3tx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1ty = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1ty,Derivative_SGM1ty_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2ty = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2ty,Derivative_SGM2ty_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3ty = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3ty,Derivative_SGM3ty_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1tz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1tz,Derivative_SGM1tz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2tz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2tz,Derivative_SGM2tz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3tz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3tz,Derivative_SGM3tz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1xx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1xx,Derivative_SGM1xx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2xx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2xx,Derivative_SGM2xx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3xx = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3xx,Derivative_SGM3xx_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1yy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1yy,Derivative_SGM1yy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2yy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2yy,Derivative_SGM2yy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3yy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3yy,Derivative_SGM3yy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1zz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1zz,Derivative_SGM1zz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2zz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2zz,Derivative_SGM2zz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3zz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3zz,Derivative_SGM3zz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1xy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1xy,Derivative_SGM1xy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2xy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2xy,Derivative_SGM2xy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3xy = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3xy,Derivative_SGM3xy_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1xz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1xz,Derivative_SGM1xz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2xz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2xz,Derivative_SGM2xz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3xz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3xz,Derivative_SGM3xz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward1yz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM1yz,Derivative_SGM1yz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward2yz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM2yz,Derivative_SGM2yz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);
    plan_backward3yz = fftw_mpi_plan_dft_3d(sizex,sizey,sizez,Derivative_SGM3yz,Derivative_SGM3yz_out,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_ESTIMATE);

    fftw_execute(plan_backward1tt);
    fftw_execute(plan_backward2tt);
    fftw_execute(plan_backward3tt);
    fftw_execute(plan_backward1tx);
    fftw_execute(plan_backward2tx);
    fftw_execute(plan_backward3tx);
    fftw_execute(plan_backward1ty);
    fftw_execute(plan_backward2ty);
    fftw_execute(plan_backward3ty);
    fftw_execute(plan_backward1tz);
    fftw_execute(plan_backward2tz);
    fftw_execute(plan_backward3tz);
    fftw_execute(plan_backward1xx);
    fftw_execute(plan_backward2xx);
    fftw_execute(plan_backward3xx);
    fftw_execute(plan_backward1yy);
    fftw_execute(plan_backward2yy);
    fftw_execute(plan_backward3yy);
    fftw_execute(plan_backward1zz);
    fftw_execute(plan_backward2zz);
    fftw_execute(plan_backward3zz);
    fftw_execute(plan_backward1xy);
    fftw_execute(plan_backward2xy);
    fftw_execute(plan_backward3xy);
    fftw_execute(plan_backward1xz);
    fftw_execute(plan_backward2xz);
    fftw_execute(plan_backward3xz);
    fftw_execute(plan_backward1yz);
    fftw_execute(plan_backward2yz);
    fftw_execute(plan_backward3yz);
    
    //   Store derivative into PHI 
    for(k = 0; k < local_n0; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1tt[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1tt_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2tt[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2tt_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3tt[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tt_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1tx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1tx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2tx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2tx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3tx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1ty[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1ty_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2ty[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2ty_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3ty[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3ty_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1tz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1tz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2tz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2tz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3tz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1xx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2xx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xx_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1yy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2yy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1zz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2zz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3zz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1xy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2xy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xy_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1xz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2xz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi1yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM1yz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi2yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM2yz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
				temp_phi3yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yz_out[i + j*sizex + k*sizex*sizey]) / (sizex*sizey*sizez);
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
       CCTK_VARIABLE_REAL, temp_phi1tt,
       CCTK_VARIABLE_REAL, temp_phi1tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2tt,
       CCTK_VARIABLE_REAL, temp_phi2tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3tt,
       CCTK_VARIABLE_REAL, temp_phi3tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1tx,
       CCTK_VARIABLE_REAL, temp_phi1tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2tx,
       CCTK_VARIABLE_REAL, temp_phi2tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3tx,
       CCTK_VARIABLE_REAL, temp_phi3tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1ty,
       CCTK_VARIABLE_REAL, temp_phi1ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2ty,
       CCTK_VARIABLE_REAL, temp_phi2ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3ty,
       CCTK_VARIABLE_REAL, temp_phi3ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1tz,
       CCTK_VARIABLE_REAL, temp_phi1tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2tz,
       CCTK_VARIABLE_REAL, temp_phi2tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3tz,
       CCTK_VARIABLE_REAL, temp_phi3tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1xx,
       CCTK_VARIABLE_REAL, temp_phi1xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2xx,
       CCTK_VARIABLE_REAL, temp_phi2xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3xx,
       CCTK_VARIABLE_REAL, temp_phi3xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1yy,
       CCTK_VARIABLE_REAL, temp_phi1yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2yy,
       CCTK_VARIABLE_REAL, temp_phi2yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3yy,
       CCTK_VARIABLE_REAL, temp_phi3yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1zz,
       CCTK_VARIABLE_REAL, temp_phi1zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2zz,
       CCTK_VARIABLE_REAL, temp_phi2zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3zz,
       CCTK_VARIABLE_REAL, temp_phi3zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1xy,
       CCTK_VARIABLE_REAL, temp_phi1xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2xy,
       CCTK_VARIABLE_REAL, temp_phi2xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3xy,
       CCTK_VARIABLE_REAL, temp_phi3xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1xz,
       CCTK_VARIABLE_REAL, temp_phi1xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2xz,
       CCTK_VARIABLE_REAL, temp_phi2xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3xz,
       CCTK_VARIABLE_REAL, temp_phi3xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi1yz,
       CCTK_VARIABLE_REAL, temp_phi1yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi2yz,
       CCTK_VARIABLE_REAL, temp_phi2yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info2, -1,
       CCTK_VARIABLE_REAL, temp_phi3yz,
       CCTK_VARIABLE_REAL, temp_phi3yz);
    assert (! ierr);
    
    //  Destroy/Free array pointers 

    fftw_destroy_plan(plan_forwardtt);
    fftw_destroy_plan(plan_forwardtx);
    fftw_destroy_plan(plan_forwardty);
    fftw_destroy_plan(plan_forwardtz);
    fftw_destroy_plan(plan_forwardxx);
    fftw_destroy_plan(plan_forwardyy);
    fftw_destroy_plan(plan_forwardzz);
    fftw_destroy_plan(plan_forwardxy);
    fftw_destroy_plan(plan_forwardxz);
    fftw_destroy_plan(plan_forwardyz);
    fftw_destroy_plan(plan_backward1tt);
    fftw_destroy_plan(plan_backward2tt);
    fftw_destroy_plan(plan_backward3tt);
    fftw_destroy_plan(plan_backward1tx);
    fftw_destroy_plan(plan_backward2tx);
    fftw_destroy_plan(plan_backward3tx);
    fftw_destroy_plan(plan_backward1ty);
    fftw_destroy_plan(plan_backward2ty);
    fftw_destroy_plan(plan_backward3ty);
    fftw_destroy_plan(plan_backward1tz);
    fftw_destroy_plan(plan_backward2tz);
    fftw_destroy_plan(plan_backward3tz);
    fftw_destroy_plan(plan_backward1xx);
    fftw_destroy_plan(plan_backward2xx);
    fftw_destroy_plan(plan_backward3xx);
    fftw_destroy_plan(plan_backward1yy);
    fftw_destroy_plan(plan_backward2yy);
    fftw_destroy_plan(plan_backward3yy);
    fftw_destroy_plan(plan_backward1zz);
    fftw_destroy_plan(plan_backward2zz);
    fftw_destroy_plan(plan_backward3zz);
    fftw_destroy_plan(plan_backward1xy);
    fftw_destroy_plan(plan_backward2xy);
    fftw_destroy_plan(plan_backward3xy);
    fftw_destroy_plan(plan_backward1xz);
    fftw_destroy_plan(plan_backward2xz);
    fftw_destroy_plan(plan_backward3xz);
    fftw_destroy_plan(plan_backward1yz);
    fftw_destroy_plan(plan_backward2yz);
    fftw_destroy_plan(plan_backward3yz);
    fftw_free(phi_intt);
    fftw_free(phi_intx);
    fftw_free(phi_inty);
    fftw_free(phi_intz);
    fftw_free(phi_inxx);
    fftw_free(phi_inyy);
    fftw_free(phi_inzz);
    fftw_free(phi_inxy);
    fftw_free(phi_inxz);
    fftw_free(phi_inyz);
    fftw_free(Derivative_SGM1tt);
    fftw_free(Derivative_SGM2tt);
    fftw_free(Derivative_SGM3tt);
    fftw_free(Derivative_SGM1tx);
    fftw_free(Derivative_SGM2tx);
    fftw_free(Derivative_SGM3tx);
    fftw_free(Derivative_SGM1ty);
    fftw_free(Derivative_SGM2ty);
    fftw_free(Derivative_SGM3ty);
    fftw_free(Derivative_SGM1tz);
    fftw_free(Derivative_SGM2tz);
    fftw_free(Derivative_SGM3tz);
    fftw_free(Derivative_SGM1xx);
    fftw_free(Derivative_SGM2xx);
    fftw_free(Derivative_SGM3xx);
    fftw_free(Derivative_SGM1yy);
    fftw_free(Derivative_SGM2yy);
    fftw_free(Derivative_SGM3yy);
    fftw_free(Derivative_SGM1zz);
    fftw_free(Derivative_SGM2zz);
    fftw_free(Derivative_SGM3zz);
    fftw_free(Derivative_SGM1xy);
    fftw_free(Derivative_SGM2xy);
    fftw_free(Derivative_SGM3xy);
    fftw_free(Derivative_SGM1xz);
    fftw_free(Derivative_SGM2xz);
    fftw_free(Derivative_SGM3xz);
    fftw_free(Derivative_SGM1yz);
    fftw_free(Derivative_SGM2yz);
    fftw_free(Derivative_SGM3yz);
    fftw_free(phi_outtt);
    fftw_free(phi_outtx);
    fftw_free(phi_outty);
    fftw_free(phi_outtz);
    fftw_free(phi_outxx);
    fftw_free(phi_outyy);
    fftw_free(phi_outzz);
    fftw_free(phi_outxy);
    fftw_free(phi_outxz);
    fftw_free(phi_outyz);
    fftw_free(Derivative_SGM1tt_out);
    fftw_free(Derivative_SGM2tt_out);
    fftw_free(Derivative_SGM3tt_out);
    fftw_free(Derivative_SGM1tx_out);
    fftw_free(Derivative_SGM2tx_out);
    fftw_free(Derivative_SGM3tx_out);
    fftw_free(Derivative_SGM1ty_out);
    fftw_free(Derivative_SGM2ty_out);
    fftw_free(Derivative_SGM3ty_out);
    fftw_free(Derivative_SGM1tz_out);
    fftw_free(Derivative_SGM2tz_out);
    fftw_free(Derivative_SGM3tz_out);
    fftw_free(Derivative_SGM1xx_out);
    fftw_free(Derivative_SGM2xx_out);
    fftw_free(Derivative_SGM3xx_out);
    fftw_free(Derivative_SGM1yy_out);
    fftw_free(Derivative_SGM2yy_out);
    fftw_free(Derivative_SGM3yy_out);
    fftw_free(Derivative_SGM1zz_out);
    fftw_free(Derivative_SGM2zz_out);
    fftw_free(Derivative_SGM3zz_out);
    fftw_free(Derivative_SGM1xy_out);
    fftw_free(Derivative_SGM2xy_out);
    fftw_free(Derivative_SGM3xy_out);
    fftw_free(Derivative_SGM1xz_out);
    fftw_free(Derivative_SGM2xz_out);
    fftw_free(Derivative_SGM3xz_out);
    fftw_free(Derivative_SGM1yz_out);
    fftw_free(Derivative_SGM2yz_out);
    fftw_free(Derivative_SGM3yz_out);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Write_Tensor_4D( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D )
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
                tensor_4D[1][0][0][index] = 0.0;
                tensor_4D[2][0][0][index] = 0.0;
                tensor_4D[3][0][0][index] = 0.0;
                tensor_4D[1][0][1][index] = 0.0;
                tensor_4D[2][0][1][index] = 0.0;
                tensor_4D[3][0][1][index] = 0.0;
                tensor_4D[1][0][2][index] = 0.0;
                tensor_4D[2][0][2][index] = 0.0;
                tensor_4D[3][0][2][index] = 0.0;
                tensor_4D[1][0][3][index] = 0.0;
                tensor_4D[2][0][3][index] = 0.0;
                tensor_4D[3][0][3][index] = 0.0;
                
                tensor_4D[1][1][0][index] = 0.0;
                tensor_4D[2][1][0][index] = 0.0;
                tensor_4D[3][1][0][index] = 0.0;
                tensor_4D[1][2][0][index] = 0.0;
                tensor_4D[2][2][0][index] = 0.0;
                tensor_4D[3][2][0][index] = 0.0;
                tensor_4D[1][3][0][index] = 0.0;
                tensor_4D[2][3][0][index] = 0.0;
                tensor_4D[3][3][0][index] = 0.0;
                
                tensor_4D[1][1][1][index] = 0.0;
                tensor_4D[2][1][1][index] = 0.0;
                tensor_4D[3][1][1][index] = 0.0;
                tensor_4D[1][2][2][index] = 0.0;
                tensor_4D[2][2][2][index] = 0.0;
                tensor_4D[3][2][2][index] = 0.0;
                tensor_4D[1][3][3][index] = 0.0;
                tensor_4D[2][3][3][index] = 0.0;
                tensor_4D[3][3][3][index] = 0.0;
                
                tensor_4D[1][1][2][index] = 0.0;
                tensor_4D[2][1][2][index] = 0.0;
                tensor_4D[3][1][2][index] = 0.0;
                tensor_4D[1][1][3][index] = 0.0;
                tensor_4D[2][1][3][index] = 0.0;
                tensor_4D[3][1][3][index] = 0.0;
                tensor_4D[1][2][3][index] = 0.0;
                tensor_4D[2][2][3][index] = 0.0;
                tensor_4D[3][2][3][index] = 0.0;
                
                tensor_4D[1][2][1][index] = 0.0;
                tensor_4D[2][2][1][index] = 0.0;
                tensor_4D[3][2][1][index] = 0.0;
                tensor_4D[1][3][1][index] = 0.0;
                tensor_4D[2][3][1][index] = 0.0;
                tensor_4D[3][3][1][index] = 0.0;
                tensor_4D[1][3][2][index] = 0.0;
                tensor_4D[2][3][2][index] = 0.0;
                tensor_4D[3][3][2][index] = 0.0;
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
                tensor_4D[1][0][0][index] = 0.0;
                tensor_4D[2][0][0][index] = 0.0;
                tensor_4D[3][0][0][index] = 0.0;
                tensor_4D[1][0][1][index] = 0.0;
                tensor_4D[2][0][1][index] = 0.0;
                tensor_4D[3][0][1][index] = 0.0;
                tensor_4D[1][0][2][index] = 0.0;
                tensor_4D[2][0][2][index] = 0.0;
                tensor_4D[3][0][2][index] = 0.0;
                tensor_4D[1][0][3][index] = 0.0;
                tensor_4D[2][0][3][index] = 0.0;
                tensor_4D[3][0][3][index] = 0.0;
                
                tensor_4D[1][1][0][index] = 0.0;
                tensor_4D[2][1][0][index] = 0.0;
                tensor_4D[3][1][0][index] = 0.0;
                tensor_4D[1][2][0][index] = 0.0;
                tensor_4D[2][2][0][index] = 0.0;
                tensor_4D[3][2][0][index] = 0.0;
                tensor_4D[1][3][0][index] = 0.0;
                tensor_4D[2][3][0][index] = 0.0;
                tensor_4D[3][3][0][index] = 0.0;
                
                tensor_4D[1][1][1][index] = 0.0;
                tensor_4D[2][1][1][index] = 0.0;
                tensor_4D[3][1][1][index] = 0.0;
                tensor_4D[1][2][2][index] = 0.0;
                tensor_4D[2][2][2][index] = 0.0;
                tensor_4D[3][2][2][index] = 0.0;
                tensor_4D[1][3][3][index] = 0.0;
                tensor_4D[2][3][3][index] = 0.0;
                tensor_4D[3][3][3][index] = 0.0;
                
                tensor_4D[1][1][2][index] = 0.0;
                tensor_4D[2][1][2][index] = 0.0;
                tensor_4D[3][1][2][index] = 0.0;
                tensor_4D[1][1][3][index] = 0.0;
                tensor_4D[2][1][3][index] = 0.0;
                tensor_4D[3][1][3][index] = 0.0;
                tensor_4D[1][2][3][index] = 0.0;
                tensor_4D[2][2][3][index] = 0.0;
                tensor_4D[3][2][3][index] = 0.0;
                
                tensor_4D[1][2][1][index] = 0.0;
                tensor_4D[2][2][1][index] = 0.0;
                tensor_4D[3][2][1][index] = 0.0;
                tensor_4D[1][3][1][index] = 0.0;
                tensor_4D[2][3][1][index] = 0.0;
                tensor_4D[3][3][1][index] = 0.0;
                tensor_4D[1][3][2][index] = 0.0;
                tensor_4D[2][3][2][index] = 0.0;
                tensor_4D[3][3][2][index] = 0.0;
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
                tensor_4D[1][0][0][index] = temp_phi1tt[ index ];
                tensor_4D[2][0][0][index] = temp_phi2tt[ index ];
                tensor_4D[3][0][0][index] = temp_phi3tt[ index ];
                tensor_4D[1][0][1][index] = temp_phi1tx[ index ];
                tensor_4D[2][0][1][index] = temp_phi2tx[ index ];
                tensor_4D[3][0][1][index] = temp_phi3tx[ index ];
                tensor_4D[1][0][2][index] = temp_phi1ty[ index ];
                tensor_4D[2][0][2][index] = temp_phi2ty[ index ];
                tensor_4D[3][0][2][index] = temp_phi3ty[ index ];
                tensor_4D[1][0][3][index] = temp_phi1tz[ index ];
                tensor_4D[2][0][3][index] = temp_phi2tz[ index ];
                tensor_4D[3][0][3][index] = temp_phi3tz[ index ];
                
                tensor_4D[1][1][0][index] = temp_phi1tx[ index ];
                tensor_4D[2][1][0][index] = temp_phi2tx[ index ];
                tensor_4D[3][1][0][index] = temp_phi3tx[ index ];
                tensor_4D[1][2][0][index] = temp_phi1ty[ index ];
                tensor_4D[2][2][0][index] = temp_phi2ty[ index ];
                tensor_4D[3][2][0][index] = temp_phi3ty[ index ];
                tensor_4D[1][3][0][index] = temp_phi1tz[ index ];
                tensor_4D[2][3][0][index] = temp_phi2tz[ index ];
                tensor_4D[3][3][0][index] = temp_phi3tz[ index ];
                
                tensor_4D[1][1][1][index] = temp_phi1xx[ index ];
                tensor_4D[2][1][1][index] = temp_phi2xx[ index ];
                tensor_4D[3][1][1][index] = temp_phi3xx[ index ];
                tensor_4D[1][2][2][index] = temp_phi1yy[ index ];
                tensor_4D[2][2][2][index] = temp_phi2yy[ index ];
                tensor_4D[3][2][2][index] = temp_phi3yy[ index ];
                tensor_4D[1][3][3][index] = temp_phi1zz[ index ];
                tensor_4D[2][3][3][index] = temp_phi2zz[ index ];
                tensor_4D[3][3][3][index] = temp_phi3zz[ index ];
                
                tensor_4D[1][1][2][index] = temp_phi1xy[ index ];
                tensor_4D[2][1][2][index] = temp_phi2xy[ index ];
                tensor_4D[3][1][2][index] = temp_phi3xy[ index ];
                tensor_4D[1][1][3][index] = temp_phi1xz[ index ];
                tensor_4D[2][1][3][index] = temp_phi2xz[ index ];
                tensor_4D[3][1][3][index] = temp_phi3xz[ index ];
                tensor_4D[1][2][3][index] = temp_phi1yz[ index ];
                tensor_4D[2][2][3][index] = temp_phi2yz[ index ];
                tensor_4D[3][2][3][index] = temp_phi3yz[ index ];
                
                tensor_4D[1][2][1][index] = temp_phi1xy[ index ];
                tensor_4D[2][2][1][index] = temp_phi2xy[ index ];
                tensor_4D[3][2][1][index] = temp_phi3xy[ index ];
                tensor_4D[1][3][1][index] = temp_phi1xz[ index ];
                tensor_4D[2][3][1][index] = temp_phi2xz[ index ];
                tensor_4D[3][3][1][index] = temp_phi3xz[ index ];
                tensor_4D[1][3][2][index] = temp_phi1yz[ index ];
                tensor_4D[2][3][2][index] = temp_phi2yz[ index ];
                tensor_4D[3][3][2][index] = temp_phi3yz[ index ];
			}
		}
	}
}
