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

void SpecDeriv_Tensor_Derivative2_4D_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2_4D );
void SpecDeriv_Write_Tensor2_4D_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2_4D );
void SpecDeriv_SpecGalMethod_Tensor2_4D_1d( CCTK_ARGUMENTS );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative2_4D_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2_4D )
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
                temp_phi0tt[ index ] = tensor2_4D[0][0][0][0][index];
				temp_phi0tx[ index ] = tensor2_4D[0][0][0][1][index];
                temp_phi0ty[ index ] = tensor2_4D[0][0][0][2][index];
                temp_phi0tz[ index ] = tensor2_4D[0][0][0][3][index];
                temp_phi0xx[ index ] = tensor2_4D[0][0][1][1][index];
                temp_phi0yy[ index ] = tensor2_4D[0][0][2][2][index];
                temp_phi0zz[ index ] = tensor2_4D[0][0][3][3][index];
                temp_phi0xy[ index ] = tensor2_4D[0][0][1][2][index];
                temp_phi0xz[ index ] = tensor2_4D[0][0][1][3][index];
                temp_phi0yz[ index ] = tensor2_4D[0][0][2][3][index];
			}
		}
	}
	SpecDeriv_SpecGalMethod_Tensor2_4D_1d( CCTK_PASS_CTOC );
	SpecDeriv_Write_Tensor2_4D_1d( CCTK_PASS_CTOC, tensor2_4D );
}

/*******************************************************/
/*******************************************************/
/*                                                     */
/*     Compute derivative of the periodic function     */
/*             SPECTRAL Galerkin Method                */
/*                                                     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_SpecGalMethod_Tensor2_4D_1d( CCTK_ARGUMENTS )
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
    fftw_plan plan_backward3tt;
    fftw_plan plan_backward3tx;
    fftw_plan plan_backward3ty;
    fftw_plan plan_backward3tz;
    fftw_plan plan_backward3xx;
    fftw_plan plan_backward3yy;
    fftw_plan plan_backward3zz;
    fftw_plan plan_backward3xy;
    fftw_plan plan_backward3xz;
    fftw_plan plan_backward3yz;
    fftw_plan plan_backward33tt;
    fftw_plan plan_backward33tx;
    fftw_plan plan_backward33ty;
    fftw_plan plan_backward33tz;
    fftw_plan plan_backward33xx;
    fftw_plan plan_backward33yy;
    fftw_plan plan_backward33zz;
    fftw_plan plan_backward33xy;
    fftw_plan plan_backward33xz;
    fftw_plan plan_backward33yz;
    CCTK_REAL Fourier_3D_kwaveN3;
    CCTK_REAL Kgrid_res = 8.0*atan(1.0);
    CCTK_REAL waveElementk;
    fftw_complex *Derivative_SGM3xx, *Derivative_SGM3yy, *Derivative_SGM3zz;
    fftw_complex *Derivative_SGM3xy, *Derivative_SGM3yz, *Derivative_SGM3xz;
    fftw_complex *Derivative_SGM3tt, *Derivative_SGM3ty;
    fftw_complex *Derivative_SGM3tx, *Derivative_SGM3tz;
    fftw_complex *Derivative_SGM33xx, *Derivative_SGM33xy;
    fftw_complex *Derivative_SGM33yy, *Derivative_SGM33yz;
    fftw_complex *Derivative_SGM33zz, *Derivative_SGM33xz;
    fftw_complex *Derivative_SGM33tt, *Derivative_SGM33tx;
    fftw_complex *Derivative_SGM33ty, *Derivative_SGM33tz;
    fftw_complex *phi_intt, *phi_intx, *phi_inty, *phi_intz;
    fftw_complex *phi_inxx, *phi_inyy, *phi_inzz;
    fftw_complex *phi_inxy, *phi_inxz, *phi_inyz;
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
       CCTK_VARIABLE_REAL, temp_phi0tt,
       CCTK_VARIABLE_REAL, temp_phi0tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0tx,
       CCTK_VARIABLE_REAL, temp_phi0tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0ty,
       CCTK_VARIABLE_REAL, temp_phi0ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0tz,
       CCTK_VARIABLE_REAL, temp_phi0tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0xx,
       CCTK_VARIABLE_REAL, temp_phi0xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0yy,
       CCTK_VARIABLE_REAL, temp_phi0yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0zz,
       CCTK_VARIABLE_REAL, temp_phi0zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0xy,
       CCTK_VARIABLE_REAL, temp_phi0xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0xz,
       CCTK_VARIABLE_REAL, temp_phi0xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, 3, info, -1,
       CCTK_VARIABLE_REAL, temp_phi0yz,
       CCTK_VARIABLE_REAL, temp_phi0yz);
    assert (! ierr);
    
    //  Begin Spectral-Galerkin Method
    //  Create phi, compute Forward FT 
    
    phi_intt = fftw_alloc_complex(f_alloc_local);
    phi_intx = fftw_alloc_complex(f_alloc_local);
    phi_inty = fftw_alloc_complex(f_alloc_local);
    phi_intz = fftw_alloc_complex(f_alloc_local);
    phi_inxx = fftw_alloc_complex(f_alloc_local);
    phi_inyy = fftw_alloc_complex(f_alloc_local);
    phi_inzz = fftw_alloc_complex(f_alloc_local);
    phi_inxy = fftw_alloc_complex(f_alloc_local);
    phi_inxz = fftw_alloc_complex(f_alloc_local);
    phi_inyz = fftw_alloc_complex(f_alloc_local);
    
    for(k = 0; k < f_local_n0; k++) 
    {
        phi_intt[k] = (fftw_complex)temp_phi0tt[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_intx[k] = (fftw_complex)temp_phi0tx[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inty[k] = (fftw_complex)temp_phi0ty[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_intz[k] = (fftw_complex)temp_phi0tz[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inxx[k] = (fftw_complex)temp_phi0xx[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inyy[k] = (fftw_complex)temp_phi0yy[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inzz[k] = (fftw_complex)temp_phi0zz[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inxy[k] = (fftw_complex)temp_phi0xy[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inxz[k] = (fftw_complex)temp_phi0xz[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
        phi_inyz[k] = (fftw_complex)temp_phi0yz[sizex-1 + (sizey-1)*sizex + k*sizex*sizey];
    }
    
    plan_forwardtt = fftw_mpi_plan_dft_1d(sizez,phi_intt,phi_intt,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardtx = fftw_mpi_plan_dft_1d(sizez,phi_intx,phi_intx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardty = fftw_mpi_plan_dft_1d(sizez,phi_inty,phi_inty,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardtz = fftw_mpi_plan_dft_1d(sizez,phi_intz,phi_intz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardxx = fftw_mpi_plan_dft_1d(sizez,phi_inxx,phi_inxx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardyy = fftw_mpi_plan_dft_1d(sizez,phi_inyy,phi_inyy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardzz = fftw_mpi_plan_dft_1d(sizez,phi_inzz,phi_inzz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardxy = fftw_mpi_plan_dft_1d(sizez,phi_inxy,phi_inxy,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardxz = fftw_mpi_plan_dft_1d(sizez,phi_inxz,phi_inxz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_forwardyz = fftw_mpi_plan_dft_1d(sizez,phi_inyz,phi_inyz,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    
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
    Derivative_SGM3tt = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3tx = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3ty = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3tz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3xx = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3yy = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3zz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3xy = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3xz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM3yz = fftw_alloc_complex(b_alloc_local); 
    
    Derivative_SGM33tt = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33tx = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33ty = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33tz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33xx = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33yy = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33zz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33xy = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33xz = fftw_alloc_complex(b_alloc_local);
    Derivative_SGM33yz = fftw_alloc_complex(b_alloc_local);
    
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
                
                Derivative_SGM3tt[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_intt[k];
                Derivative_SGM3tx[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_intx[k];
                Derivative_SGM3ty[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inty[k];
                Derivative_SGM3tz[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_intz[k];
                Derivative_SGM3xx[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inxx[k];
                Derivative_SGM3yy[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inyy[k];
                Derivative_SGM3zz[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inzz[k];
                Derivative_SGM3xy[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inxy[k];
                Derivative_SGM3xz[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inxz[k];
                Derivative_SGM3yz[k] = IMAG_NUMBER * Fourier_3D_kwaveN3 * phi_inyz[k];
                
                Derivative_SGM33tt[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_intt[k];
                Derivative_SGM33tx[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_intx[k];
                Derivative_SGM33ty[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inty[k];
                Derivative_SGM33tz[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_intz[k];
                
                Derivative_SGM33xx[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inxx[k];
                Derivative_SGM33yy[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inyy[k];
                Derivative_SGM33zz[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inzz[k];
                Derivative_SGM33xy[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inxy[k];
                Derivative_SGM33xz[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inxz[k];
                Derivative_SGM33yz[k] = -1.0 * Fourier_3D_kwaveN3 * Fourier_3D_kwaveN3 * phi_inyz[k];
    }
    
    //  Create Derivative_SGM, compute Inverse FT 
    
    plan_backward3tt = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3tt,Derivative_SGM3tt,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3tx = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3tx,Derivative_SGM3tx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3ty = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3ty,Derivative_SGM3ty,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3tz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3tz,Derivative_SGM3tz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3xx = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3xx,Derivative_SGM3xx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3yy = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3yy,Derivative_SGM3yy,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3zz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3zz,Derivative_SGM3zz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3xy = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3xy,Derivative_SGM3xy,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3xz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3xz,Derivative_SGM3xz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward3yz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM3yz,Derivative_SGM3yz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);

    plan_backward33tt = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33tt,Derivative_SGM33tt,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33tx = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33tx,Derivative_SGM33tx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33ty = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33ty,Derivative_SGM33ty,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33tz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33tz,Derivative_SGM33tz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33xx = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33xx,Derivative_SGM33xx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33yy = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33yy,Derivative_SGM33yy,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33zz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33zz,Derivative_SGM33zz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33xy = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33xy,Derivative_SGM33xy,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33xz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33xz,Derivative_SGM33xz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);
    plan_backward33yz = fftw_mpi_plan_dft_1d(sizez,Derivative_SGM33yz,Derivative_SGM33yz,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_WISDOM_ONLY | FFTW_PATIENT);

    fftw_execute(plan_backward3tt);
    fftw_execute(plan_backward3tx);
    fftw_execute(plan_backward3ty);
    fftw_execute(plan_backward3tz);
    fftw_execute(plan_backward3xx);
    fftw_execute(plan_backward3yy);
    fftw_execute(plan_backward3zz);
    fftw_execute(plan_backward3xy);
    fftw_execute(plan_backward3xz);
    fftw_execute(plan_backward3yz);
    
    fftw_execute(plan_backward33tt);
    fftw_execute(plan_backward33tx);
    fftw_execute(plan_backward33ty);
    fftw_execute(plan_backward33tz);
    fftw_execute(plan_backward33xx);
    fftw_execute(plan_backward33yy);
    fftw_execute(plan_backward33zz);
    fftw_execute(plan_backward33xy);
    fftw_execute(plan_backward33xz);
    fftw_execute(plan_backward33yz);
    
    //   Store derivative into PHI 
    for(k = 0; k < b_local_n1; k++) 
    {
        for(j = 0; j < sizey; j++) 
        {
            for(i = 0; i < sizex; i++) 
            {
				temp_phi1tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3tt[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tt[k]) / (sizez);
				temp_phi1tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3tx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tx[k]) / (sizez);
				temp_phi1ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3ty[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3ty[k]) / (sizez);
				temp_phi1tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3tz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3tz[k]) / (sizez);
				temp_phi1xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xx[k]) / (sizez);
				temp_phi1yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yy[k]) / (sizez);
				temp_phi1zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3zz[k]) / (sizez);
				temp_phi1xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xy[k]) / (sizez);
				temp_phi1xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xz[k]) / (sizez);
				temp_phi1yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yz[k]) / (sizez);
				
				temp_phi11tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33tt[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33tt[k]) / (sizez);
				temp_phi11tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33tx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33tx[k]) / (sizez);
				temp_phi11ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33ty[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33ty[k]) / (sizez);
				temp_phi11tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33tz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33tz[k]) / (sizez);
				temp_phi12tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23tt[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23tx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23ty[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23tz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				
				temp_phi1xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xx[k]) / (sizez);
				temp_phi1yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yy[k]) / (sizez);
				temp_phi1zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3zz[k]) / (sizez);
				temp_phi1xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xy[k]) / (sizez);
				temp_phi1xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3xz[k]) / (sizez);
				temp_phi1yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi2yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi3yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM3yz[k]) / (sizez);
				temp_phi11xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33xx[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33xx[k]) / (sizez);
				temp_phi11yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33yy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33yy[k]) / (sizez);
				temp_phi11zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33zz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33zz[k]) / (sizez);
				temp_phi11xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33xy[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33xy[k]) / (sizez);
				temp_phi11xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33xz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33xz[k]) / (sizez);
				temp_phi11yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi22yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi33yz[ i + j*sizex + k*sizex*sizey ] = creal(Derivative_SGM33yz[k]) / (sizez);
				temp_phi12xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23xx[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23yy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23zz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23xy[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23xz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi12yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi13yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
				temp_phi23yz[ i + j*sizex + k*sizex*sizey ] = 0.0;
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
       CCTK_VARIABLE_REAL, temp_phi1tt,
       CCTK_VARIABLE_REAL, temp_phi1tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2tt,
       CCTK_VARIABLE_REAL, temp_phi2tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3tt,
       CCTK_VARIABLE_REAL, temp_phi3tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1tx,
       CCTK_VARIABLE_REAL, temp_phi1tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2tx,
       CCTK_VARIABLE_REAL, temp_phi2tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3tx,
       CCTK_VARIABLE_REAL, temp_phi3tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1ty,
       CCTK_VARIABLE_REAL, temp_phi1ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2ty,
       CCTK_VARIABLE_REAL, temp_phi2ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3ty,
       CCTK_VARIABLE_REAL, temp_phi3ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1tz,
       CCTK_VARIABLE_REAL, temp_phi1tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2tz,
       CCTK_VARIABLE_REAL, temp_phi2tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3tz,
       CCTK_VARIABLE_REAL, temp_phi3tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1xx,
       CCTK_VARIABLE_REAL, temp_phi1xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2xx,
       CCTK_VARIABLE_REAL, temp_phi2xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3xx,
       CCTK_VARIABLE_REAL, temp_phi3xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1yy,
       CCTK_VARIABLE_REAL, temp_phi1yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2yy,
       CCTK_VARIABLE_REAL, temp_phi2yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3yy,
       CCTK_VARIABLE_REAL, temp_phi3yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1zz,
       CCTK_VARIABLE_REAL, temp_phi1zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2zz,
       CCTK_VARIABLE_REAL, temp_phi2zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3zz,
       CCTK_VARIABLE_REAL, temp_phi3zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1xy,
       CCTK_VARIABLE_REAL, temp_phi1xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2xy,
       CCTK_VARIABLE_REAL, temp_phi2xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3xy,
       CCTK_VARIABLE_REAL, temp_phi3xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1xz,
       CCTK_VARIABLE_REAL, temp_phi1xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2xz,
       CCTK_VARIABLE_REAL, temp_phi2xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3xz,
       CCTK_VARIABLE_REAL, temp_phi3xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi1yz,
       CCTK_VARIABLE_REAL, temp_phi1yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi2yz,
       CCTK_VARIABLE_REAL, temp_phi2yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi3yz,
       CCTK_VARIABLE_REAL, temp_phi3yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11tt,
       CCTK_VARIABLE_REAL, temp_phi11tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22tt,
       CCTK_VARIABLE_REAL, temp_phi22tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33tt,
       CCTK_VARIABLE_REAL, temp_phi33tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11tx,
       CCTK_VARIABLE_REAL, temp_phi11tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22tx,
       CCTK_VARIABLE_REAL, temp_phi22tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33tx,
       CCTK_VARIABLE_REAL, temp_phi33tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11ty,
       CCTK_VARIABLE_REAL, temp_phi11ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22ty,
       CCTK_VARIABLE_REAL, temp_phi22ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33ty,
       CCTK_VARIABLE_REAL, temp_phi33ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11tz,
       CCTK_VARIABLE_REAL, temp_phi11tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22tz,
       CCTK_VARIABLE_REAL, temp_phi22tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33tz,
       CCTK_VARIABLE_REAL, temp_phi33tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12tt,
       CCTK_VARIABLE_REAL, temp_phi12tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13tt,
       CCTK_VARIABLE_REAL, temp_phi13tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23tt,
       CCTK_VARIABLE_REAL, temp_phi23tt);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12tx,
       CCTK_VARIABLE_REAL, temp_phi12tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13tx,
       CCTK_VARIABLE_REAL, temp_phi13tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23tx,
       CCTK_VARIABLE_REAL, temp_phi23tx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12ty,
       CCTK_VARIABLE_REAL, temp_phi12ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13ty,
       CCTK_VARIABLE_REAL, temp_phi13ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23ty,
       CCTK_VARIABLE_REAL, temp_phi23ty);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12tz,
       CCTK_VARIABLE_REAL, temp_phi12tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13tz,
       CCTK_VARIABLE_REAL, temp_phi13tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23tz,
       CCTK_VARIABLE_REAL, temp_phi23tz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11xx,
       CCTK_VARIABLE_REAL, temp_phi11xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22xx,
       CCTK_VARIABLE_REAL, temp_phi22xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33xx,
       CCTK_VARIABLE_REAL, temp_phi33xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11yy,
       CCTK_VARIABLE_REAL, temp_phi11yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22yy,
       CCTK_VARIABLE_REAL, temp_phi22yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33yy,
       CCTK_VARIABLE_REAL, temp_phi33yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11zz,
       CCTK_VARIABLE_REAL, temp_phi11zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22zz,
       CCTK_VARIABLE_REAL, temp_phi22zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33zz,
       CCTK_VARIABLE_REAL, temp_phi33zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11xy,
       CCTK_VARIABLE_REAL, temp_phi11xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22xy,
       CCTK_VARIABLE_REAL, temp_phi22xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33xy,
       CCTK_VARIABLE_REAL, temp_phi33xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11xz,
       CCTK_VARIABLE_REAL, temp_phi11xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22xz,
       CCTK_VARIABLE_REAL, temp_phi22xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33xz,
       CCTK_VARIABLE_REAL, temp_phi33xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi11yz,
       CCTK_VARIABLE_REAL, temp_phi11yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi22yz,
       CCTK_VARIABLE_REAL, temp_phi22yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi33yz,
       CCTK_VARIABLE_REAL, temp_phi33yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12xx,
       CCTK_VARIABLE_REAL, temp_phi12xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13xx,
       CCTK_VARIABLE_REAL, temp_phi13xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23xx,
       CCTK_VARIABLE_REAL, temp_phi23xx);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12yy,
       CCTK_VARIABLE_REAL, temp_phi12yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13yy,
       CCTK_VARIABLE_REAL, temp_phi13yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23yy,
       CCTK_VARIABLE_REAL, temp_phi23yy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12zz,
       CCTK_VARIABLE_REAL, temp_phi12zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13zz,
       CCTK_VARIABLE_REAL, temp_phi13zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23zz,
       CCTK_VARIABLE_REAL, temp_phi23zz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12xy,
       CCTK_VARIABLE_REAL, temp_phi12xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13xy,
       CCTK_VARIABLE_REAL, temp_phi13xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23xy,
       CCTK_VARIABLE_REAL, temp_phi23xy);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12xz,
       CCTK_VARIABLE_REAL, temp_phi12xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13xz,
       CCTK_VARIABLE_REAL, temp_phi13xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23xz,
       CCTK_VARIABLE_REAL, temp_phi23xz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi12yz,
       CCTK_VARIABLE_REAL, temp_phi12yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi13yz,
       CCTK_VARIABLE_REAL, temp_phi13yz);
    assert (! ierr);
    
    ierr = Slab_Transfer
      (cctkGH, cctk_dim, info, -1,
       CCTK_VARIABLE_REAL, temp_phi23yz,
       CCTK_VARIABLE_REAL, temp_phi23yz);
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
    fftw_destroy_plan(plan_backward3tt);
    fftw_destroy_plan(plan_backward3tx);
    fftw_destroy_plan(plan_backward3ty);
    fftw_destroy_plan(plan_backward3tz);
    fftw_destroy_plan(plan_backward3xx);
    fftw_destroy_plan(plan_backward3yy);
    fftw_destroy_plan(plan_backward3zz);
    fftw_destroy_plan(plan_backward3xy);
    fftw_destroy_plan(plan_backward3xz);
    fftw_destroy_plan(plan_backward3yz);
    fftw_destroy_plan(plan_backward33tt);
    fftw_destroy_plan(plan_backward33tx);
    fftw_destroy_plan(plan_backward33ty);
    fftw_destroy_plan(plan_backward33tz);
    fftw_destroy_plan(plan_backward33xx);
    fftw_destroy_plan(plan_backward33yy);
    fftw_destroy_plan(plan_backward33zz);
    fftw_destroy_plan(plan_backward33xy);
    fftw_destroy_plan(plan_backward33xz);
    fftw_destroy_plan(plan_backward33yz);
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
    fftw_free(Derivative_SGM3tt);
    fftw_free(Derivative_SGM3tx);
    fftw_free(Derivative_SGM3ty);
    fftw_free(Derivative_SGM3tz);
    fftw_free(Derivative_SGM3xx);
    fftw_free(Derivative_SGM3yy);
    fftw_free(Derivative_SGM3zz);
    fftw_free(Derivative_SGM3xy);
    fftw_free(Derivative_SGM3xz);
    fftw_free(Derivative_SGM3yz);
    fftw_free(Derivative_SGM33tt);
    fftw_free(Derivative_SGM33tx);
    fftw_free(Derivative_SGM33ty);
    fftw_free(Derivative_SGM33tz);
    fftw_free(Derivative_SGM33xx);
    fftw_free(Derivative_SGM33yy);
    fftw_free(Derivative_SGM33zz);
    fftw_free(Derivative_SGM33xy);
    fftw_free(Derivative_SGM33xz);
    fftw_free(Derivative_SGM33yz);
    fftw_mpi_cleanup();
    fftw_cleanup();
}

/*******************************************************/
/*******************************************************/
/*        Function to Read GF into Output Array        */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Write_Tensor2_4D_1d( CCTK_ARGUMENTS, CCTK_REAL *****tensor2_4D )
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
                tensor2_4D[0][1][0][0][index] = temp_phi1tt[ index ];
                tensor2_4D[0][2][0][0][index] = temp_phi2tt[ index ];
                tensor2_4D[0][3][0][0][index] = temp_phi3tt[ index ];
                tensor2_4D[0][1][0][1][index] = temp_phi1tx[ index ];
                tensor2_4D[0][2][0][1][index] = temp_phi2tx[ index ];
                tensor2_4D[0][3][0][1][index] = temp_phi3tx[ index ];
                tensor2_4D[0][1][0][2][index] = temp_phi1ty[ index ];
                tensor2_4D[0][2][0][2][index] = temp_phi2ty[ index ];
                tensor2_4D[0][3][0][2][index] = temp_phi3ty[ index ];
                tensor2_4D[0][1][0][3][index] = temp_phi1tz[ index ];
                tensor2_4D[0][2][0][3][index] = temp_phi2tz[ index ];
                tensor2_4D[0][3][0][3][index] = temp_phi3tz[ index ];
                
                tensor2_4D[0][1][1][0][index] = temp_phi1tx[ index ];
                tensor2_4D[0][2][1][0][index] = temp_phi2tx[ index ];
                tensor2_4D[0][3][1][0][index] = temp_phi3tx[ index ];
                tensor2_4D[0][1][2][0][index] = temp_phi1ty[ index ];
                tensor2_4D[0][2][2][0][index] = temp_phi2ty[ index ];
                tensor2_4D[0][3][2][0][index] = temp_phi3ty[ index ];
                tensor2_4D[0][1][3][0][index] = temp_phi1tz[ index ];
                tensor2_4D[0][2][3][0][index] = temp_phi2tz[ index ];
                tensor2_4D[0][3][3][0][index] = temp_phi3tz[ index ];
                
                tensor2_4D[1][0][0][0][index] = temp_phi1tt[ index ];
                tensor2_4D[2][0][0][0][index] = temp_phi2tt[ index ];
                tensor2_4D[3][0][0][0][index] = temp_phi3tt[ index ];
                tensor2_4D[1][0][0][1][index] = temp_phi1tx[ index ];
                tensor2_4D[2][0][0][1][index] = temp_phi2tx[ index ];
                tensor2_4D[3][0][0][1][index] = temp_phi3tx[ index ];
                tensor2_4D[1][0][0][2][index] = temp_phi1ty[ index ];
                tensor2_4D[2][0][0][2][index] = temp_phi2ty[ index ];
                tensor2_4D[3][0][0][2][index] = temp_phi3ty[ index ];
                tensor2_4D[1][0][0][3][index] = temp_phi1tz[ index ];
                tensor2_4D[2][0][0][3][index] = temp_phi2tz[ index ];
                tensor2_4D[3][0][0][3][index] = temp_phi3tz[ index ];
                
                tensor2_4D[1][0][1][0][index] = temp_phi1tx[ index ];
                tensor2_4D[2][0][1][0][index] = temp_phi2tx[ index ];
                tensor2_4D[3][0][1][0][index] = temp_phi3tx[ index ];
                tensor2_4D[1][0][2][0][index] = temp_phi1ty[ index ];
                tensor2_4D[2][0][2][0][index] = temp_phi2ty[ index ];
                tensor2_4D[3][0][2][0][index] = temp_phi3ty[ index ];
                tensor2_4D[1][0][3][0][index] = temp_phi1tz[ index ];
                tensor2_4D[2][0][3][0][index] = temp_phi2tz[ index ];
                tensor2_4D[3][0][3][0][index] = temp_phi3tz[ index ];
                
                tensor2_4D[0][1][1][1][index] = temp_phi1xx[ index ];
                tensor2_4D[0][2][1][1][index] = temp_phi2xx[ index ];
                tensor2_4D[0][3][1][1][index] = temp_phi3xx[ index ];
                tensor2_4D[0][1][2][2][index] = temp_phi1yy[ index ];
                tensor2_4D[0][2][2][2][index] = temp_phi2yy[ index ];
                tensor2_4D[0][3][2][2][index] = temp_phi3yy[ index ];
                tensor2_4D[0][1][3][3][index] = temp_phi1zz[ index ];
                tensor2_4D[0][2][3][3][index] = temp_phi2zz[ index ];
                tensor2_4D[0][3][3][3][index] = temp_phi3zz[ index ];
                
                tensor2_4D[1][0][1][1][index] = temp_phi1xx[ index ];
                tensor2_4D[2][0][1][1][index] = temp_phi2xx[ index ];
                tensor2_4D[3][0][1][1][index] = temp_phi3xx[ index ];
                tensor2_4D[1][0][2][2][index] = temp_phi1yy[ index ];
                tensor2_4D[2][0][2][2][index] = temp_phi2yy[ index ];
                tensor2_4D[3][0][2][2][index] = temp_phi3yy[ index ];
                tensor2_4D[1][0][3][3][index] = temp_phi1zz[ index ];
                tensor2_4D[2][0][3][3][index] = temp_phi2zz[ index ];
                tensor2_4D[3][0][3][3][index] = temp_phi3zz[ index ];
                
                tensor2_4D[0][1][1][2][index] = temp_phi1xy[ index ];
                tensor2_4D[0][2][1][2][index] = temp_phi2xy[ index ];
                tensor2_4D[0][3][1][2][index] = temp_phi3xy[ index ];
                tensor2_4D[0][1][1][3][index] = temp_phi1xz[ index ];
                tensor2_4D[0][2][1][3][index] = temp_phi2xz[ index ];
                tensor2_4D[0][3][1][3][index] = temp_phi3xz[ index ];
                tensor2_4D[0][1][2][3][index] = temp_phi1yz[ index ];
                tensor2_4D[0][2][2][3][index] = temp_phi2yz[ index ];
                tensor2_4D[0][3][2][3][index] = temp_phi3yz[ index ];
                
                tensor2_4D[0][1][2][1][index] = temp_phi1xy[ index ];
                tensor2_4D[0][2][2][1][index] = temp_phi2xy[ index ];
                tensor2_4D[0][3][2][1][index] = temp_phi3xy[ index ];
                tensor2_4D[0][1][3][1][index] = temp_phi1xz[ index ];
                tensor2_4D[0][2][3][1][index] = temp_phi2xz[ index ];
                tensor2_4D[0][3][3][1][index] = temp_phi3xz[ index ];
                tensor2_4D[0][1][3][2][index] = temp_phi1yz[ index ];
                tensor2_4D[0][2][3][2][index] = temp_phi2yz[ index ];
                tensor2_4D[0][3][3][2][index] = temp_phi3yz[ index ];
                
                tensor2_4D[1][0][1][2][index] = temp_phi1xy[ index ];
                tensor2_4D[2][0][1][2][index] = temp_phi2xy[ index ];
                tensor2_4D[3][0][1][2][index] = temp_phi3xy[ index ];
                tensor2_4D[1][0][1][3][index] = temp_phi1xz[ index ];
                tensor2_4D[2][0][1][3][index] = temp_phi2xz[ index ];
                tensor2_4D[3][0][1][3][index] = temp_phi3xz[ index ];
                tensor2_4D[1][0][2][3][index] = temp_phi1yz[ index ];
                tensor2_4D[2][0][2][3][index] = temp_phi2yz[ index ];
                tensor2_4D[3][0][2][3][index] = temp_phi3yz[ index ];
                
                tensor2_4D[1][0][2][1][index] = temp_phi1xy[ index ];
                tensor2_4D[2][0][2][1][index] = temp_phi2xy[ index ];
                tensor2_4D[3][0][2][1][index] = temp_phi3xy[ index ];
                tensor2_4D[1][0][3][1][index] = temp_phi1xz[ index ];
                tensor2_4D[2][0][3][1][index] = temp_phi2xz[ index ];
                tensor2_4D[3][0][3][1][index] = temp_phi3xz[ index ];
                tensor2_4D[1][0][3][2][index] = temp_phi1yz[ index ];
                tensor2_4D[2][0][3][2][index] = temp_phi2yz[ index ];
                tensor2_4D[3][0][3][2][index] = temp_phi3yz[ index ];
                
                tensor2_4D[1][1][0][0][index] = temp_phi11tt[ index ];
                tensor2_4D[2][2][0][0][index] = temp_phi22tt[ index ];
                tensor2_4D[3][3][0][0][index] = temp_phi33tt[ index ];
                tensor2_4D[1][1][0][1][index] = temp_phi11tx[ index ];
                tensor2_4D[2][2][0][1][index] = temp_phi22tx[ index ];
                tensor2_4D[3][3][0][1][index] = temp_phi33tx[ index ];
                tensor2_4D[1][1][0][2][index] = temp_phi11ty[ index ];
                tensor2_4D[2][2][0][2][index] = temp_phi22ty[ index ];
                tensor2_4D[3][3][0][2][index] = temp_phi33ty[ index ];
                tensor2_4D[1][1][0][3][index] = temp_phi11tz[ index ];
                tensor2_4D[2][2][0][3][index] = temp_phi22tz[ index ];
                tensor2_4D[3][3][0][3][index] = temp_phi33tz[ index ];
                
                tensor2_4D[1][1][1][0][index] = temp_phi11tx[ index ];
                tensor2_4D[2][2][1][0][index] = temp_phi22tx[ index ];
                tensor2_4D[3][3][1][0][index] = temp_phi33tx[ index ];
                tensor2_4D[1][1][2][0][index] = temp_phi11ty[ index ];
                tensor2_4D[2][2][2][0][index] = temp_phi22ty[ index ];
                tensor2_4D[3][3][2][0][index] = temp_phi33ty[ index ];
                tensor2_4D[1][1][3][0][index] = temp_phi11tz[ index ];
                tensor2_4D[2][2][3][0][index] = temp_phi22tz[ index ];
                tensor2_4D[3][3][3][0][index] = temp_phi33tz[ index ];
                
                tensor2_4D[1][2][0][0][index] = temp_phi12tt[ index ];
                tensor2_4D[1][3][0][0][index] = temp_phi13tt[ index ];
                tensor2_4D[2][3][0][0][index] = temp_phi23tt[ index ];
                tensor2_4D[1][2][0][1][index] = temp_phi12tx[ index ];
                tensor2_4D[1][3][0][1][index] = temp_phi13tx[ index ];
                tensor2_4D[2][3][0][1][index] = temp_phi23tx[ index ];
                tensor2_4D[1][2][0][2][index] = temp_phi12ty[ index ];
                tensor2_4D[1][3][0][2][index] = temp_phi13ty[ index ];
                tensor2_4D[2][3][0][2][index] = temp_phi23ty[ index ];
                tensor2_4D[1][2][0][3][index] = temp_phi12tz[ index ];
                tensor2_4D[1][3][0][3][index] = temp_phi13tz[ index ];
                tensor2_4D[2][3][0][3][index] = temp_phi23tz[ index ];
                
                tensor2_4D[2][1][0][0][index] = temp_phi12tt[ index ];
                tensor2_4D[3][1][0][0][index] = temp_phi13tt[ index ];
                tensor2_4D[3][2][0][0][index] = temp_phi23tt[ index ];
                tensor2_4D[2][1][0][1][index] = temp_phi12tx[ index ];
                tensor2_4D[3][1][0][1][index] = temp_phi13tx[ index ];
                tensor2_4D[3][2][0][1][index] = temp_phi23tx[ index ];
                tensor2_4D[2][1][0][2][index] = temp_phi12ty[ index ];
                tensor2_4D[3][1][0][2][index] = temp_phi13ty[ index ];
                tensor2_4D[3][2][0][2][index] = temp_phi23ty[ index ];
                tensor2_4D[2][1][0][3][index] = temp_phi12tz[ index ];
                tensor2_4D[3][1][0][3][index] = temp_phi13tz[ index ];
                tensor2_4D[3][2][0][3][index] = temp_phi23tz[ index ];
                
                tensor2_4D[1][2][1][0][index] = temp_phi12tx[ index ];
                tensor2_4D[1][3][1][0][index] = temp_phi13tx[ index ];
                tensor2_4D[2][3][1][0][index] = temp_phi23tx[ index ];
                tensor2_4D[1][2][2][0][index] = temp_phi12ty[ index ];
                tensor2_4D[1][3][2][0][index] = temp_phi13ty[ index ];
                tensor2_4D[2][3][2][0][index] = temp_phi23ty[ index ];
                tensor2_4D[1][2][3][0][index] = temp_phi12tz[ index ];
                tensor2_4D[1][3][3][0][index] = temp_phi13tz[ index ];
                tensor2_4D[2][3][3][0][index] = temp_phi23tz[ index ];
                
                tensor2_4D[2][1][1][0][index] = temp_phi12tx[ index ];
                tensor2_4D[3][1][1][0][index] = temp_phi13tx[ index ];
                tensor2_4D[3][2][1][0][index] = temp_phi23tx[ index ];
                tensor2_4D[2][1][2][0][index] = temp_phi12ty[ index ];
                tensor2_4D[3][1][2][0][index] = temp_phi13ty[ index ];
                tensor2_4D[3][2][2][0][index] = temp_phi23ty[ index ];
                tensor2_4D[2][1][3][0][index] = temp_phi12tz[ index ];
                tensor2_4D[3][1][3][0][index] = temp_phi13tz[ index ];
                tensor2_4D[3][2][3][0][index] = temp_phi23tz[ index ];
                
                tensor2_4D[1][1][1][1][index] = temp_phi11xx[ index ];
                tensor2_4D[2][2][1][1][index] = temp_phi22xx[ index ];
                tensor2_4D[3][3][1][1][index] = temp_phi33xx[ index ];
                tensor2_4D[1][1][2][2][index] = temp_phi11yy[ index ];
                tensor2_4D[2][2][2][2][index] = temp_phi22yy[ index ];
                tensor2_4D[3][3][2][2][index] = temp_phi33yy[ index ];
                tensor2_4D[1][1][3][3][index] = temp_phi11zz[ index ];
                tensor2_4D[2][2][3][3][index] = temp_phi22zz[ index ];
                tensor2_4D[3][3][3][3][index] = temp_phi33zz[ index ];
                
                tensor2_4D[1][2][1][1][index] = temp_phi12xx[ index ];
                tensor2_4D[1][3][1][1][index] = temp_phi13xx[ index ];
                tensor2_4D[2][3][1][1][index] = temp_phi23xx[ index ];
                tensor2_4D[1][2][2][2][index] = temp_phi12yy[ index ];
                tensor2_4D[1][3][2][2][index] = temp_phi13yy[ index ];
                tensor2_4D[2][3][2][2][index] = temp_phi23yy[ index ];
                tensor2_4D[1][2][3][3][index] = temp_phi12zz[ index ];
                tensor2_4D[1][3][3][3][index] = temp_phi13zz[ index ];
                tensor2_4D[2][3][3][3][index] = temp_phi23zz[ index ];
                
                tensor2_4D[2][1][1][1][index] = temp_phi12xx[ index ];
                tensor2_4D[3][1][1][1][index] = temp_phi13xx[ index ];
                tensor2_4D[3][2][1][1][index] = temp_phi23xx[ index ];
                tensor2_4D[2][1][2][2][index] = temp_phi12yy[ index ];
                tensor2_4D[3][1][2][2][index] = temp_phi13yy[ index ];
                tensor2_4D[3][2][2][2][index] = temp_phi23yy[ index ];
                tensor2_4D[2][1][3][3][index] = temp_phi12zz[ index ];
                tensor2_4D[3][1][3][3][index] = temp_phi13zz[ index ];
                tensor2_4D[3][2][3][3][index] = temp_phi23zz[ index ];
                
                tensor2_4D[1][1][1][2][index] = temp_phi11xy[ index ];
                tensor2_4D[2][2][1][2][index] = temp_phi22xy[ index ];
                tensor2_4D[3][3][1][2][index] = temp_phi33xy[ index ];
                tensor2_4D[1][1][1][3][index] = temp_phi11xz[ index ];
                tensor2_4D[2][2][1][3][index] = temp_phi22xz[ index ];
                tensor2_4D[3][3][1][3][index] = temp_phi33xz[ index ];
                tensor2_4D[1][1][2][3][index] = temp_phi11yz[ index ];
                tensor2_4D[2][2][2][3][index] = temp_phi22yz[ index ];
                tensor2_4D[3][3][2][3][index] = temp_phi33yz[ index ];
                
                tensor2_4D[1][1][2][1][index] = temp_phi11xy[ index ];
                tensor2_4D[2][2][2][1][index] = temp_phi22xy[ index ];
                tensor2_4D[3][3][2][1][index] = temp_phi33xy[ index ];
                tensor2_4D[1][1][3][1][index] = temp_phi11xz[ index ];
                tensor2_4D[2][2][3][1][index] = temp_phi22xz[ index ];
                tensor2_4D[3][3][3][1][index] = temp_phi33xz[ index ];
                tensor2_4D[1][1][3][2][index] = temp_phi11yz[ index ];
                tensor2_4D[2][2][3][2][index] = temp_phi22yz[ index ];
                tensor2_4D[3][3][3][2][index] = temp_phi33yz[ index ];
                
                tensor2_4D[1][2][1][2][index] = temp_phi12xy[ index ];
                tensor2_4D[1][3][1][2][index] = temp_phi13xy[ index ];
                tensor2_4D[2][3][1][2][index] = temp_phi23xy[ index ];
                tensor2_4D[1][2][1][3][index] = temp_phi12xz[ index ];
                tensor2_4D[1][3][1][3][index] = temp_phi13xz[ index ];
                tensor2_4D[2][3][1][3][index] = temp_phi23xz[ index ];
                tensor2_4D[1][2][2][3][index] = temp_phi12yz[ index ];
                tensor2_4D[1][3][2][3][index] = temp_phi13yz[ index ];
                tensor2_4D[2][3][2][3][index] = temp_phi23yz[ index ];
                
                tensor2_4D[1][2][2][1][index] = temp_phi12xy[ index ];
                tensor2_4D[1][3][2][1][index] = temp_phi13xy[ index ];
                tensor2_4D[2][3][2][1][index] = temp_phi23xy[ index ];
                tensor2_4D[1][2][3][1][index] = temp_phi12xz[ index ];
                tensor2_4D[1][3][3][1][index] = temp_phi13xz[ index ];
                tensor2_4D[2][3][3][1][index] = temp_phi23xz[ index ];
                tensor2_4D[1][2][3][2][index] = temp_phi12yz[ index ];
                tensor2_4D[1][3][3][2][index] = temp_phi13yz[ index ];
                tensor2_4D[2][3][3][2][index] = temp_phi23yz[ index ];
                
                tensor2_4D[2][1][1][2][index] = temp_phi12xy[ index ];
                tensor2_4D[3][1][1][2][index] = temp_phi13xy[ index ];
                tensor2_4D[3][2][1][2][index] = temp_phi23xy[ index ];
                tensor2_4D[2][1][1][3][index] = temp_phi12xz[ index ];
                tensor2_4D[3][1][1][3][index] = temp_phi13xz[ index ];
                tensor2_4D[3][2][1][3][index] = temp_phi23xz[ index ];
                tensor2_4D[2][1][2][3][index] = temp_phi12yz[ index ];
                tensor2_4D[3][1][2][3][index] = temp_phi13yz[ index ];
                tensor2_4D[3][2][2][3][index] = temp_phi23yz[ index ];
                
                tensor2_4D[2][1][2][1][index] = temp_phi12xy[ index ];
                tensor2_4D[3][1][2][1][index] = temp_phi13xy[ index ];
                tensor2_4D[3][2][2][1][index] = temp_phi23xy[ index ];
                tensor2_4D[2][1][3][1][index] = temp_phi12xz[ index ];
                tensor2_4D[3][1][3][1][index] = temp_phi13xz[ index ];
                tensor2_4D[3][2][3][1][index] = temp_phi23xz[ index ];
                tensor2_4D[2][1][3][2][index] = temp_phi12yz[ index ];
                tensor2_4D[3][1][3][2][index] = temp_phi13yz[ index ];
                tensor2_4D[3][2][3][2][index] = temp_phi23yz[ index ];
			}
		}
	}
}
