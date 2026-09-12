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


void SpecDeriv_Tensor_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor )
{
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS
	
    // Declare and initialize variables
    CCTK_INT i, j, k, m, n, p, istart, jstart, kstart, iend, jend, kend;
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


			for (m = 1; m < 4; m++) {
				for (n = 1; n < 4; n++) {
					for (p = 1; p < 4; p++) {
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
						if ((abs(tensor[0][n][p][index])+abs(tensor[0][n][p][uuindex]))!=0.0){
	
							tr1 = (tensor[0][n][p][index]-tensor[0][n][p][uuindex])/(abs(tensor[0][n][p][index]) + abs(tensor[0][n][p][uuindex]));
							fdtg = abs(tr1);
						 
							} 
							
						//printf("fdtg for Tensor_4D: %13f\n", fdtg);

					
						if (fdtg>scp){


							//************************************************************************
							//****                         r=3 method                     ************
							//************************************************************************
					
							//****
							//   R where Rj(x)=Rj(uj)
							//****
					
							//  tensor[0][n][p]
							//  CCTK_DELTA_SPACE(m-1)
					
							IS0 = ((1.0/2.0)*(pow((tensor[0][n][p][vvindex]-tensor[0][n][p][vv2index]),2)+pow((tensor[0][n][p][index]-tensor[0][n][p][vvindex]),2)))+(pow((tensor[0][n][p][index]-(2.0*tensor[0][n][p][vvindex])+tensor[0][n][p][vv2index]),2));
							IS1 = ((1.0/2.0)*(pow((tensor[0][n][p][index]-tensor[0][n][p][vvindex]),2)+pow((tensor[0][n][p][uuindex]-tensor[0][n][p][index]),2)))+(pow((tensor[0][n][p][uuindex]-(2.0*tensor[0][n][p][index])+tensor[0][n][p][vvindex]),2));
							IS2 = ((1.0/2.0)*(pow((tensor[0][n][p][uuindex]-tensor[0][n][p][index]),2)+pow((tensor[0][n][p][uu2index]-tensor[0][n][p][uuindex]),2)))+(pow((tensor[0][n][p][uu2index]-(2.0*tensor[0][n][p][uuindex])+tensor[0][n][p][index]),2));
					
							c01 = (tensor[0][n][p][index]-2.0*tensor[0][n][p][vvindex]+tensor[0][n][p][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c02 = (tensor[0][n][p][index]-tensor[0][n][p][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03 = tensor[0][n][p][vvindex];
							c04 = (tensor[0][n][p][index]-2.0*tensor[0][n][p][vvindex]+tensor[0][n][p][vv2index])/(24.0);

							c11 = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c12 = (tensor[0][n][p][uuindex]-tensor[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13 = tensor[0][n][p][index];
							c14 = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(24.0);
					
							c21 = (tensor[0][n][p][uu2index]-2.0*tensor[0][n][p][uuindex]+tensor[0][n][p][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c22 = (tensor[0][n][p][uu2index]-tensor[0][n][p][index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23 = tensor[0][n][p][uuindex];
							c24 = (tensor[0][n][p][uu2index]-2.0*tensor[0][n][p][uuindex]+tensor[0][n][p][index])/(24.0);
					
							PM0 = (c01*pow(tensor[0][n][p][index],2))+(c02*tensor[0][n][p][index])+c03+c04;
							PM1 = (c11*pow(tensor[0][n][p][index],2))+(c12*tensor[0][n][p][index])+c13+c14;
							PM2 = (c21*pow(tensor[0][n][p][index],2))+(c22*tensor[0][n][p][index])+c23+c24;
					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor[0][n][p][index]
							//************************************************************************
					
							fd = (tensor[0][n][p][uuindex] - tensor[0][n][p][vvindex]) / (2.0*CCTK_DELTA_SPACE(m-1));

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
					
							IS0ru = ((1.0/2.0)*(pow((tensor[0][n][p][index]-tensor[0][n][p][vvindex]),2)+pow((tensor[0][n][p][uuindex]-tensor[0][n][p][index]),2)))+(pow((tensor[0][n][p][uuindex]-(2.0*tensor[0][n][p][index])+tensor[0][n][p][vvindex]),2));
							IS1ru = ((1.0/2.0)*(pow((tensor[0][n][p][uuindex]-tensor[0][n][p][index]),2)+pow((tensor[0][n][p][uu2index]-tensor[0][n][p][uuindex]),2)))+(pow((tensor[0][n][p][uu2index]-(2.0*tensor[0][n][p][uuindex])+tensor[0][n][p][index]),2));
							IS2ru = ((1.0/2.0)*(pow((tensor[0][n][p][uu2index]-tensor[0][n][p][uuindex]),2)+pow((tensor[0][n][p][uu3index]-tensor[0][n][p][uu2index]),2)))+(pow((tensor[0][n][p][uu3index]-(2.0*tensor[0][n][p][uu2index])+tensor[0][n][p][uuindex]),2));
					
							c01ru = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c02ru = (tensor[0][n][p][uuindex]-tensor[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03ru = tensor[0][n][p][index];
							c04ru = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(24.0);

							c11ru = (tensor[0][n][p][uu2index]-2.0*tensor[0][n][p][uuindex]+tensor[0][n][p][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c12ru = (tensor[0][n][p][uu2index]-tensor[0][n][p][index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13ru = tensor[0][n][p][uuindex];
							c14ru = (tensor[0][n][p][uu2index]-2.0*tensor[0][n][p][uuindex]+tensor[0][n][p][index])/(24.0);
					
							c21ru = (tensor[0][n][p][uu3index]-2.0*tensor[0][n][p][uu2index]+tensor[0][n][p][uuindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c22ru = (tensor[0][n][p][uu3index]-tensor[0][n][p][uuindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23ru = tensor[0][n][p][uu2index];
							c24ru = (tensor[0][n][p][uu3index]-2.0*tensor[0][n][p][uu2index]+tensor[0][n][p][uuindex])/(24.0);
					
					
							PM0ru = (c01ru*pow(tensor[0][n][p][uuindex],2))+(c02ru*tensor[0][n][p][uuindex])+c03ru+c04ru;
							PM1ru = (c11ru*pow(tensor[0][n][p][uuindex],2))+(c12ru*tensor[0][n][p][uuindex])+c13ru+c14ru;
							PM2ru = (c21ru*pow(tensor[0][n][p][uuindex],2))+(c22ru*tensor[0][n][p][uuindex])+c23ru+c24ru;
					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor[0][n][p][index]
							//************************************************************************
					
							fdru = (tensor[0][n][p][uu2index] - tensor[0][n][p][index]) / (2.0*CCTK_DELTA_SPACE(m-1));

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
	
							IS0rv = ((1.0/2.0)*(pow((tensor[0][n][p][vv2index]-tensor[0][n][p][vv3index]),2)+pow((tensor[0][n][p][vvindex]-tensor[0][n][p][vv2index]),2)))+(pow((tensor[0][n][p][vvindex]-(2.0*tensor[0][n][p][vv2index])+tensor[0][n][p][vv3index]),2));
							IS1rv = ((1.0/2.0)*(pow((tensor[0][n][p][vvindex]-tensor[0][n][p][vv2index]),2)+pow((tensor[0][n][p][index]-tensor[0][n][p][vvindex]),2)))+(pow((tensor[0][n][p][index]-(2.0*tensor[0][n][p][vvindex])+tensor[0][n][p][vv2index]),2));
							IS2rv = ((1.0/2.0)*(pow((tensor[0][n][p][index]-tensor[0][n][p][vvindex]),2)+pow((tensor[0][n][p][uuindex]-tensor[0][n][p][index]),2)))+(pow((tensor[0][n][p][uuindex]-(2.0*tensor[0][n][p][index])+tensor[0][n][p][vvindex]),2));
					
							c01rv = (tensor[0][n][p][vvindex]-2.0*tensor[0][n][p][vv2index]+tensor[0][n][p][vv3index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c02rv = (tensor[0][n][p][vvindex]-tensor[0][n][p][vv3index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03rv = tensor[0][n][p][vv2index];
							c04rv = (tensor[0][n][p][vvindex]-2.0*tensor[0][n][p][vv2index]+tensor[0][n][p][vv3index])/(24.0);

							c11rv = (tensor[0][n][p][index]-2.0*tensor[0][n][p][vvindex]+tensor[0][n][p][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c12rv = (tensor[0][n][p][index]-tensor[0][n][p][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13rv = tensor[0][n][p][vvindex];
							c14rv = (tensor[0][n][p][index]-2.0*tensor[0][n][p][vvindex]+tensor[0][n][p][vv2index])/(24.0);
					
							c21rv = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2));
							c22rv = (tensor[0][n][p][uuindex]-tensor[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23rv = tensor[0][n][p][index];
							c24rv = (tensor[0][n][p][uuindex]-2.0*tensor[0][n][p][index]+tensor[0][n][p][vvindex])/(24.0);
					
							PM0rv = (c01rv*pow(tensor[0][n][p][vvindex],2))+(c02rv*tensor[0][n][p][vvindex])+c03rv+c04rv;
							PM1rv = (c11rv*pow(tensor[0][n][p][vvindex],2))+(c12rv*tensor[0][n][p][vvindex])+c13rv+c14rv;
							PM2rv = (c21rv*pow(tensor[0][n][p][vvindex],2))+(c22rv*tensor[0][n][p][vvindex])+c23rv+c24rv;
					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor[0][n][p][index]
							//************************************************************************
					
							fdrv = (tensor[0][n][p][index] - tensor[0][n][p][vv2index]) / (2.0*CCTK_DELTA_SPACE(m-1));

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
					
							tensor[m][n][p][index]=(Ruu-Rvv)/(2*CCTK_DELTA_SPACE(m-1));  //***This needs to be updated for 4th order differential

							//************************************************************************
							//************************************************************************
						}				

						//****************************
						//** 4th order FD Method
						//****************************
						if (fdtg<=scp){
							tensor[m][n][p][index] = (tensor[0][n][p][vv2index] - 8.0*tensor[0][n][p][vvindex] + 8.0*tensor[0][n][p][uuindex] - tensor[0][n][p][uu2index]) / (12.0*CCTK_DELTA_SPACE(m-1));
							}
						//****************************
						//** 4th order  FD Method End
						//****************************
					
						
                      } 
                   }
                }
			}
		}
	}
}

