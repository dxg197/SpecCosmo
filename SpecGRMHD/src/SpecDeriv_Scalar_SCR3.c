#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>


/*******************************************************/
/*******************************************************/
/*******************************************************/
/* 	        r=3 Shock Capturing Scheme                  */
/* 	       Written by Paul Smith May 2014               */
/*  Updated on:      10/19/2014                         */
/*  Notes:  updated e1 locations to avoid seg faults    */
/*                                                      */
/*******************************************************/
/*******************************************************/
/*******************************************************/

     
void SpecDeriv_Scalar_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL **scalar );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Scalar_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL **scalar )
{
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

    // Declare and initialize variables
    CCTK_INT i, j, k, m, istart, jstart, kstart, iend, jend, kend;
    CCTK_INT uu_index[4], vv_index[4],index;
    CCTK_INT uu2_index[4], vv2_index[4];
    CCTK_INT uu3_index[4], vv3_index[4];
    CCTK_INT uuindex, vvindex;
    CCTK_INT uu2index, vv2index, uu3index, vv3index;
	CCTK_REAL c01, c02, c03, c04, c11, c12, c13, c14, c21, c22, c23, c24;
	CCTK_REAL c01ru, c02ru, c03ru, c04ru, c11ru, c12ru, c13ru, c14ru, c21ru, c22ru, c23ru, c24ru;
	CCTK_REAL c01rv, c02rv, c03rv, c04rv, c11rv, c12rv, c13rv, c14rv, c21rv, c22rv, c23rv, c24rv;
	CCTK_REAL IS0, IS1, IS2, PM0, PM1, PM2, a0, a1, a2, fd, fdtg, tr1; 
	CCTK_REAL IS0ru, IS1ru, IS2ru, PM0ru, PM1ru, PM2ru, a0ru, a1ru, a2ru, fdru; 
	CCTK_REAL IS0rv, IS1rv, IS2rv, PM0rv, PM1rv, PM2rv, a0rv, a1rv, a2rv, fdrv;   
	CCTK_REAL R, Rvv, Ruu;   
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
				uu_index[1] = CCTK_GFINDEX3D( cctkGH, i+1, j, k );
				vv_index[1] = CCTK_GFINDEX3D( cctkGH, i-1, j, k );
				uu_index[2] = CCTK_GFINDEX3D( cctkGH, i, j+1, k );
				vv_index[2] = CCTK_GFINDEX3D( cctkGH, i, j-1, k );
				uu_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k+1 );
				vv_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k-1 );
				uu2_index[1] = CCTK_GFINDEX3D( cctkGH, i+2, j, k );
				vv2_index[1] = CCTK_GFINDEX3D( cctkGH, i-2, j, k );
				uu2_index[2] = CCTK_GFINDEX3D( cctkGH, i, j+2, k );
				vv2_index[2] = CCTK_GFINDEX3D( cctkGH, i, j-2, k );
				uu2_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k+2 );
				vv2_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k-2 );

				uu3_index[1] = CCTK_GFINDEX3D( cctkGH, i+3, j, k );
				vv3_index[1] = CCTK_GFINDEX3D( cctkGH, i-3, j, k );
				uu3_index[2] = CCTK_GFINDEX3D( cctkGH, i, j+3, k );
				vv3_index[2] = CCTK_GFINDEX3D( cctkGH, i, j-3, k );
				uu3_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k+3 );
				vv3_index[3] = CCTK_GFINDEX3D( cctkGH, i, j, k-3 );
				
                
				for (m = 1; m < 4; m++)
				{ 
					uuindex = uu_index[m];
					vvindex = vv_index[m];
					uu2index = uu2_index[m];
					vv2index = vv2_index[m];
					uu3index = uu3_index[m];
					vv3index = vv3_index[m];

					//*****************
					//*****SC Trigger
					//*****************
					fdtg = 0.0;
					if ((abs(scalar[0][index])+abs(scalar[0][uuindex]))!=0.0){
	
						tr1 = (scalar[0][index]-scalar[0][uuindex])/(abs(scalar[0][index]) + abs(scalar[0][uuindex]));
						fdtg = abs(tr1);
					
						} 
						
						//printf("fdtg for Tensor_4D: %13f\n", fdtg);
					
					
					if (fdtg>scp){
						//************************************************************************
						//****                         r=3 method                     ************
						//************************************************************************
					
						//****
						//   R  where Rj(x)=Rj(uj)
						//****
					
						//    scalar[0][index]
						//    CCTK_DELTA_SPACE(m-1)
					
						IS0 = ((1.0/2.0)*(pow((scalar[0][vvindex]-scalar[0][vv2index]),2)+pow((scalar[0][index]-scalar[0][vvindex]),2)))+(pow((scalar[0][index]-(2.0*scalar[0][vvindex])+scalar[0][vv2index]),2));
						IS1 = ((1.0/2.0)*(pow((scalar[0][index]-scalar[0][vvindex]),2)+pow((scalar[0][uuindex]-scalar[0][index]),2)))+(pow((scalar[0][uuindex]-(2.0*scalar[0][index])+scalar[0][vvindex]),2));
						IS2 = ((1.0/2.0)*(pow((scalar[0][uuindex]-scalar[0][index]),2)+pow((scalar[0][uu2index]-scalar[0][uuindex]),2)))+(pow((scalar[0][uu2index]-(2.0*scalar[0][uuindex])+scalar[0][index]),2));
					
						c01 = (scalar[0][index]-2.0*scalar[0][vvindex]+scalar[0][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c02 = (scalar[0][index]-scalar[0][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
						c03 = scalar[0][vvindex];
						c04 = (scalar[0][index]-2.0*scalar[0][vvindex]+scalar[0][vv2index])/(24.0);

						c11 = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c12 = (scalar[0][uuindex]-scalar[0][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
						c13 = scalar[0][index];
						c14 = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(24.0);
					
						c21 = (scalar[0][uu2index]-2.0*scalar[0][uuindex]+scalar[0][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c22 = (scalar[0][uu2index]-scalar[0][index])/(2.0*CCTK_DELTA_SPACE(m-1));
						c23 = scalar[0][uuindex];
						c24 = (scalar[0][uu2index]-2.0*scalar[0][uuindex]+scalar[0][index])/(24.0);
					
						PM0 = (c01*pow(scalar[0][index],2))+(c02*scalar[0][index])+c03+c04;
						PM1 = (c11*pow(scalar[0][index],2))+(c12*scalar[0][index])+c13+c14;
						PM2 = (c21*pow(scalar[0][index],2))+(c22*scalar[0][index])+c23+c24;
					
						//************************************************************************
						// r=3 method requires different alphas depending on the sign of the first derivative of scalar[0][index]
						//************************************************************************
					
						fd = (scalar[0][uuindex] - scalar[0][vvindex]) / (2.0*CCTK_DELTA_SPACE(m-1));

						if (fd>0.0){
						
							a0 = (1.0)/(12.0*pow((e1+IS0),3.0));
							a1 = (1.0)/(2.0*pow((e1+IS1),3.0));
							a2 = (1.0)/(4.0*pow((e1+IS2),3.0));
						
							}
						
						if(fd<=0.0){
						
							a0 = (1.0)/(4.0*pow((e1+IS0),3.0));
							a1 = (1.0)/(2.0*pow((e1+IS1),3.0));
							a2 = (1.0)/(12.0*pow((e1+IS2),3.0));
						
							}
					
						R=((a0/(a0+a1+a2))*PM0)+((a1/(a0+a1+a2))*PM1)+((a2/(a0+a1+a2))*PM2); //*****Updated for R=3
						//************************************************************************
						//************************************************************************

						//****
						//   Rj(x+1) where Rj(x+1)=Rj(u+1)
						//****
					
						IS0ru = ((1.0/2.0)*(pow((scalar[0][index]-scalar[0][vvindex]),2)+pow((scalar[0][uuindex]-scalar[0][index]),2)))+(pow((scalar[0][uuindex]-(2.0*scalar[0][index])+scalar[0][vvindex]),2));
						IS1ru = ((1.0/2.0)*(pow((scalar[0][uuindex]-scalar[0][index]),2)+pow((scalar[0][uu2index]-scalar[0][uuindex]),2)))+(pow((scalar[0][uu2index]-(2.0*scalar[0][uuindex])+scalar[0][index]),2));
						IS2ru = ((1.0/2.0)*(pow((scalar[0][uu2index]-scalar[0][uuindex]),2)+pow((scalar[0][uu3index]-scalar[0][uu2index]),2)))+(pow((scalar[0][uu3index]-(2.0*scalar[0][uu2index])+scalar[0][uuindex]),2));
					
						c01ru = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c02ru = (scalar[0][uuindex]-scalar[0][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
						c03ru = scalar[0][index];
						c04ru = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(24.0);

						c11ru = (scalar[0][uu2index]-2.0*scalar[0][uuindex]+scalar[0][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c12ru = (scalar[0][uu2index]-scalar[0][index])/(2.0*CCTK_DELTA_SPACE(m-1));
						c13ru = scalar[0][uuindex];
						c14ru = (scalar[0][uu2index]-2.0*scalar[0][uuindex]+scalar[0][index])/(24.0);
					
						c21ru = (scalar[0][uu3index]-2.0*scalar[0][uu2index]+scalar[0][uuindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c22ru = (scalar[0][uu3index]-scalar[0][uuindex])/(2.0*CCTK_DELTA_SPACE(m-1));
						c23ru = scalar[0][uu2index];
						c24ru = (scalar[0][uu3index]-2.0*scalar[0][uu2index]+scalar[0][uuindex])/(24.0);
					
					
						PM0ru = (c01ru*pow(scalar[0][uuindex],2))+(c02ru*scalar[0][uuindex])+c03ru+c04ru;
						PM1ru = (c11ru*pow(scalar[0][uuindex],2))+(c12ru*scalar[0][uuindex])+c13ru+c14ru;
						PM2ru = (c21ru*pow(scalar[0][uuindex],2))+(c22ru*scalar[0][uuindex])+c23ru+c24ru;
					
						//************************************************************************
						// r=3 method requires different alphas depending on the sign of the first derivative of scalar[0][index]
						//************************************************************************
					
						fdru = (scalar[0][uu2index] - scalar[0][index]) / (2.0*CCTK_DELTA_SPACE(m-1));

						if (fdru>0.0){
						
							a0ru = (1.0)/(12.0*pow((e1+IS0ru),3.0));
							a1ru = (1.0)/(2.0*pow((e1+IS1ru),3.0));
							a2ru = (1.0)/(4.0*pow((e1+IS2ru),3.0));
						
							}
						
						if(fdru<=0.0){
						
							a0ru = (1.0)/(4.0*pow((e1+IS0ru),3.0));
							a1ru = (1.0)/(2.0*pow((e1+IS1ru),3.0));
							a2ru = (1.0)/(12.0*pow((e1+IS2ru),3.0));
						
							}
					
						Ruu=((a0ru/(a0ru+a1ru+a2ru))*PM0ru)+((a1ru/(a0ru+a1ru+a2ru))*PM1ru)+((a2ru/(a0ru+a1ru+a2ru))*PM2ru); 
						//************************************************************************
						//************************************************************************

						//****
						//   Rj(x-1) where Rj(x-1)=Rj(u-1)
						//****
	
						IS0rv = ((1.0/2.0)*(pow((scalar[0][vv2index]-scalar[0][vv3index]),2)+pow((scalar[0][vvindex]-scalar[0][vv2index]),2)))+(pow((scalar[0][vvindex]-(2.0*scalar[0][vv2index])+scalar[0][vv3index]),2));
						IS1rv = ((1.0/2.0)*(pow((scalar[0][vvindex]-scalar[0][vv2index]),2)+pow((scalar[0][index]-scalar[0][vvindex]),2)))+(pow((scalar[0][index]-(2.0*scalar[0][vvindex])+scalar[0][vv2index]),2));
						IS2rv = ((1.0/2.0)*(pow((scalar[0][index]-scalar[0][vvindex]),2)+pow((scalar[0][uuindex]-scalar[0][index]),2)))+(pow((scalar[0][uuindex]-(2.0*scalar[0][index])+scalar[0][vvindex]),2));
					
						c01rv = (scalar[0][vvindex]-2.0*scalar[0][vv2index]+scalar[0][vv3index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c02rv = (scalar[0][vvindex]-scalar[0][vv3index])/(2.0*CCTK_DELTA_SPACE(m-1));
						c03rv = scalar[0][vv2index];
						c04rv = (scalar[0][vvindex]-2.0*scalar[0][vv2index]+scalar[0][vv3index])/(24.0);

						c11rv = (scalar[0][index]-2.0*scalar[0][vvindex]+scalar[0][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c12rv = (scalar[0][index]-scalar[0][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
						c13rv = scalar[0][vvindex];
						c14rv = (scalar[0][index]-2.0*scalar[0][vvindex]+scalar[0][vv2index])/(24.0);
					
						c21rv = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
						c22rv = (scalar[0][uuindex]-scalar[0][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
						c23rv = scalar[0][index];
						c24rv = (scalar[0][uuindex]-2.0*scalar[0][index]+scalar[0][vvindex])/(24.0);
					
						PM0rv = (c01rv*pow(scalar[0][vvindex],2))+(c02rv*scalar[0][vvindex])+c03rv+c04rv;
						PM1rv = (c11rv*pow(scalar[0][vvindex],2))+(c12rv*scalar[0][vvindex])+c13rv+c14rv;
						PM2rv = (c21rv*pow(scalar[0][vvindex],2))+(c22rv*scalar[0][vvindex])+c23rv+c24rv;
					
					
						//************************************************************************
						// r=3 method requires different alphas depending on the sign of the first derivative of scalar[0][index]
						//************************************************************************
					
						fdrv = (scalar[0][index] - scalar[0][vv2index]) / (2.0*CCTK_DELTA_SPACE(m-1));

						if (fdrv>0.0){

							a0rv = (1.0)/(12.0*pow((e1+IS0rv),3.0));
							a1rv = (1.0)/(2.0*pow((e1+IS1rv),3.0));
							a2rv = (1.0)/(4.0*pow((e1+IS2rv),3.0));
						
							}
						
						if(fdrv<=0.0){
						
							a0rv = (1.0)/(4.0*pow((e1+IS0rv),3.0));
							a1rv = (1.0)/(2.0*pow((e1+IS1rv),3.0));
							a2rv = (1.0)/(12.0*pow((e1+IS2rv),3.0));
						
							}

						Rvv = ((a0rv/(a0rv+a1rv+a2rv))*PM0rv)+((a1rv/(a0rv+a1rv+a2rv))*PM1rv)+((a2rv/(a0rv+a1rv+a2rv))*PM2rv); //*****Updated for R=3
					

						//************************************************************************
						//************************************************************************
					
						scalar[m][index]=(Ruu-Rvv)/(2*CCTK_DELTA_SPACE(m-1));  //***This needs to be updated for 4th order differential 

						//************************************************************************
						//************************************************************************
						
						}
					
					//****************************
					//** 4th order FD Method
					//****************************
					if (fdtg<=scp){
						scalar[m][index] = (scalar[0][vv2index] - 8.0*scalar[0][vvindex] + 8.0*scalar[0][uuindex] - scalar[0][uu2index]) / (12.0*CCTK_DELTA_SPACE(m-1));
						}
					//****************************
					//** 4th order  FD Method End
					//****************************
					
				  
				}

			}
		}
	}
}		
	


