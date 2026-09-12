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


void SpecDeriv_Tensor_Derivative_4D_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D );

/*******************************************************/
/*******************************************************/
/*           Function to GF from Input Array           */
/*     (assume same # of modes in x-y-z direction)     */
/*******************************************************/
/*******************************************************/
void SpecDeriv_Tensor_Derivative_4D_SCR3( CCTK_ARGUMENTS, CCTK_REAL ****tensor_4D )
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
	CCTK_REAL IS0, IS1, IS2, PM0, PM1, PM2, a0, a1, a2, fd, fdtg; 
	CCTK_REAL IS0ru, IS1ru, IS2ru, PM0ru, PM1ru, PM2ru, a0ru, a1ru, a2ru, fdru; 
	CCTK_REAL IS0rv, IS1rv, IS2rv, PM0rv, PM1rv, PM2rv, a0rv, a1rv, a2rv, fdrv;   
	CCTK_REAL R, Rvv, Ruu, tr, tr1;   
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
				for (n = 0; n < 4; n++) {
					for (p = 0; p < 4; p++) {
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
						if ((abs(tensor_4D[0][n][p][index])+abs(tensor_4D[0][n][p][uuindex]))!=0.0){
	
							tr1 = (tensor_4D[0][n][p][index]-tensor_4D[0][n][p][uuindex])/(abs(tensor_4D[0][n][p][index]) + abs(tensor_4D[0][n][p][uuindex]));
							fdtg = abs(tr1);
						 
							} 
							
						//printf("fdtg for Tensor_4D: %13f\n", fdtg);

					
						if (fdtg>scp){

							//printf("R=3 Method Triggered\n");
							//************************************************************************
							//****                         r=3 method                     ************
							//************************************************************************
					
							//****
							//   R  where Rj(x)=Rj(uj)
							//****
					
							//  tensor_4D[0][n][p]
							//  CCTK_DELTA_SPACE(m-1)
					
							IS0 = ((1.0/2.0)*(pow((tensor_4D[0][n][p][vvindex]-tensor_4D[0][n][p][vv2index]),2.0)+pow((tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vvindex]),2.0)))+(pow((tensor_4D[0][n][p][index]-(2.0*tensor_4D[0][n][p][vvindex])+tensor_4D[0][n][p][vv2index]),2.0));
							IS1 = ((1.0/2.0)*(pow((tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vvindex]),2.0)+pow((tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][index]),2.0)))+(pow((tensor_4D[0][n][p][uuindex]-(2.0*tensor_4D[0][n][p][index])+tensor_4D[0][n][p][vvindex]),2.0));
							IS2 = ((1.0/2.0)*(pow((tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][index]),2.0)+pow((tensor_4D[0][n][p][uu2index]-tensor_4D[0][n][p][uuindex]),2.0)))+(pow((tensor_4D[0][n][p][uu2index]-(2.0*tensor_4D[0][n][p][uuindex])+tensor_4D[0][n][p][index]),2.0));
					
							c01 = (tensor_4D[0][n][p][index]-2.0*tensor_4D[0][n][p][vvindex]+tensor_4D[0][n][p][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c02 = (tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03 = tensor_4D[0][n][p][vvindex];
							c04 = (tensor_4D[0][n][p][index]-2.0*tensor_4D[0][n][p][vvindex]+tensor_4D[0][n][p][vv2index])/(24.0);

							c11 = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c12 = (tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13 = tensor_4D[0][n][p][index];
							c14 = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(24.0);
					
							c21 = (tensor_4D[0][n][p][uu2index]-2.0*tensor_4D[0][n][p][uuindex]+tensor_4D[0][n][p][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c22 = (tensor_4D[0][n][p][uu2index]-tensor_4D[0][n][p][index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23 = tensor_4D[0][n][p][uuindex];
							c24 = (tensor_4D[0][n][p][uu2index]-2.0*tensor_4D[0][n][p][uuindex]+tensor_4D[0][n][p][index])/(24.0);
					
							PM0 = (c01*pow(tensor_4D[0][n][p][index],2.0))+(c02*tensor_4D[0][n][p][index])+c03+c04;
							PM1 = (c11*pow(tensor_4D[0][n][p][index],2.0))+(c12*tensor_4D[0][n][p][index])+c13+c14;
							PM2 = (c21*pow(tensor_4D[0][n][p][index],2.0))+(c22*tensor_4D[0][n][p][index])+c23+c24;
						
							//printf("IS0 for Tensor_4D: %13f\n", IS0);
							//printf("IS1 for Tensor_4D: %13f\n", IS1);
							//printf("IS2 for Tensor_4D: %13f\n", IS2);
						
							//printf("c01 for Tensor_4D: %13f\n", c01);
							//printf("c02 for Tensor_4D: %13f\n", c02);
							//printf("c03 for Tensor_4D: %13f\n", c03);
							//printf("c04 for Tensor_4D: %13f\n", c04);
						
							//printf("c11 for Tensor_4D: %13f\n", c11);
							//printf("c12 for Tensor_4D: %13f\n", c12);
							//printf("c13 for Tensor_4D: %13f\n", c13);
							//printf("c14 for Tensor_4D: %13f\n", c14);
						
							//printf("c21 for Tensor_4D: %13f\n", c21);
							//printf("c22 for Tensor_4D: %13f\n", c22);
							//printf("c23 for Tensor_4D: %13f\n", c23);
							//printf("c24 for Tensor_4D: %13f\n", c24);
						
							//printf("PM0 for Tensor_4D: %13f\n", PM0);
							//printf("PM1 for Tensor_4D: %13f\n", PM1);
							//printf("PM2 for Tensor_4D: %13f\n", PM2);
						

					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor_4D[0][n][p][index]
							//************************************************************************
					
							fd = (tensor_4D[0][n][p][uuindex] - tensor_4D[0][n][p][vvindex]) / (2.0*CCTK_DELTA_SPACE(m-1));
							//printf("fd for Tensor_4D: %13f\n", fd);

							if (fd>0.0){
						
								a0 = (1.0)/(12.0*pow((e1+IS0),3.0));
								a1 = (1.0)/(2.0*pow((e1+IS1),3.0));
								a2 = (1.0)/(4.0*pow((e1+IS2),3.0));
								//printf("if fd for Tensor_4D >0.0 a0: %13f\n", a0);
								//printf("if fd for Tensor_4D >0.0 a1: %13f\n", a1);
								//printf("if fd for Tensor_4D >0.0 a2: %13f\n", a2);
							
								}
						
							if(fd<=0.0){
						
								a0 = (1.0)/(4.0*pow((e1+IS0),3.0));
								a1 = (1.0)/(2.0*pow((e1+IS1),3.0));
								a2 = (1.0)/(12.0*pow((e1+IS2),3.0));
							
								//printf("if fd for Tensor_4D <0.0 a0: %13f\n", a0);
								//printf("if fd for Tensor_4D <0.0 a1: %13f\n", a1);
								//printf("if fd for Tensor_4D <0.0 a2: %13f\n", a2);

						
								}
					
							R=((a0/(a0+a1+a2))*PM0)+((a1/(a0+a1+a2))*PM1)+((a2/(a0+a1+a2))*PM2); //*****Updated for R=3
							//printf("R for Tensor_4D: %13f\n", R);
							//************************************************************************
							//************************************************************************

							//****
							//   Rj(x+1) where Rj(x+1)=Rj(u+1)
							//****
					
							IS0ru = ((1.0/2.0)*(pow((tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vvindex]),2.0)+pow((tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][index]),2.0)))+(pow((tensor_4D[0][n][p][uuindex]-(2.0*tensor_4D[0][n][p][index])+tensor_4D[0][n][p][vvindex]),2.0));
							IS1ru = ((1.0/2.0)*(pow((tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][index]),2.0)+pow((tensor_4D[0][n][p][uu2index]-tensor_4D[0][n][p][uuindex]),2.0)))+(pow((tensor_4D[0][n][p][uu2index]-(2.0*tensor_4D[0][n][p][uuindex])+tensor_4D[0][n][p][index]),2.0));
							IS2ru = ((1.0/2.0)*(pow((tensor_4D[0][n][p][uu2index]-tensor_4D[0][n][p][uuindex]),2.0)+pow((tensor_4D[0][n][p][uu3index]-tensor_4D[0][n][p][uu2index]),2.0)))+(pow((tensor_4D[0][n][p][uu3index]-(2.0*tensor_4D[0][n][p][uu2index])+tensor_4D[0][n][p][uuindex]),2.0));
					
							c01ru = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c02ru = (tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03ru = tensor_4D[0][n][p][index];
							c04ru = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(24.0);

							c11ru = (tensor_4D[0][n][p][uu2index]-2.0*tensor_4D[0][n][p][uuindex]+tensor_4D[0][n][p][index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c12ru = (tensor_4D[0][n][p][uu2index]-tensor_4D[0][n][p][index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13ru = tensor_4D[0][n][p][uuindex];
							c14ru = (tensor_4D[0][n][p][uu2index]-2.0*tensor_4D[0][n][p][uuindex]+tensor_4D[0][n][p][index])/(24.0);
					
							c21ru = (tensor_4D[0][n][p][uu3index]-2.0*tensor_4D[0][n][p][uu2index]+tensor_4D[0][n][p][uuindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c22ru = (tensor_4D[0][n][p][uu3index]-tensor_4D[0][n][p][uuindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23ru = tensor_4D[0][n][p][uu2index];
							c24ru = (tensor_4D[0][n][p][uu3index]-2.0*tensor_4D[0][n][p][uu2index]+tensor_4D[0][n][p][uuindex])/(24.0);
					
					
							PM0ru = (c01ru*pow(tensor_4D[0][n][p][uuindex],2.0))+(c02ru*tensor_4D[0][n][p][uuindex])+c03ru+c04ru;
							PM1ru = (c11ru*pow(tensor_4D[0][n][p][uuindex],2.0))+(c12ru*tensor_4D[0][n][p][uuindex])+c13ru+c14ru;
							PM2ru = (c21ru*pow(tensor_4D[0][n][p][uuindex],2.0))+(c22ru*tensor_4D[0][n][p][uuindex])+c23ru+c24ru;
						
							//printf("IS0ru for Tensor_4D: %13f\n", IS0ru);
							//printf("IS1ru for Tensor_4D: %13f\n", IS1ru);
							//printf("IS2ru for Tensor_4D: %13f\n", IS2ru);
						
							//printf("c01ru for Tensor_4D: %13f\n", c01ru);
							//printf("c02ru for Tensor_4D: %13f\n", c02ru);
							//printf("c03ru for Tensor_4D: %13f\n", c03ru);
							//printf("c04ru for Tensor_4D: %13f\n", c04ru);
						
							//printf("c11ru for Tensor_4D: %13f\n", c11ru);
							//printf("c12ru for Tensor_4D: %13f\n", c12ru);
							//printf("c13ru for Tensor_4D: %13f\n", c13ru);
							//printf("c14ru for Tensor_4D: %13f\n", c14ru);
						
							//printf("c21ru for Tensor_4D: %13f\n", c21ru);
							//printf("c22ru for Tensor_4D: %13f\n", c22ru);
							//printf("c23ru for Tensor_4D: %13f\n", c23ru);
							//printf("c24ru for Tensor_4D: %13f\n", c24ru);
						
							//printf("PM0ru for Tensor_4D: %13f\n", PM0ru);
							//printf("PM1ru for Tensor_4D: %13f\n", PM1ru);
							//printf("PM2ru for Tensor_4D: %13f\n", PM2ru);

					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor_4D[0][n][p][index]
							//************************************************************************
					
							fdru = (tensor_4D[0][n][p][uu2index] - tensor_4D[0][n][p][index]) / (2.0*CCTK_DELTA_SPACE(m-1));
							//printf("fdru for Tensor_4D: %13f\n", fdru);

							if (fdru>0.0){
						
								a0ru = (1.0)/(12.0*pow((e1+IS0ru),3.0));
								a1ru = (1.0)/(2.0*pow((e1+IS1ru),3.0));
								a2ru = (1.0)/(4.0*pow((e1+IS2ru),3.0));
								//printf("if fd for Tensor_4D >0.0 a0ru: %13f\n", a0ru);
								//printf("if fd for Tensor_4D >0.0 a1ru: %13f\n", a1ru);
								//printf("if fd for Tensor_4D >0.0 a2ru: %13f\n", a2ru);

						
								}
						
							if(fdru<=0.0){
						
								a0ru = (1.0)/(4.0*pow((e1+IS0ru),3.0));
								a1ru = (1.0)/(2.0*pow((e1+IS1ru),3.0));
								a2ru = (1.0)/(12.0*pow((e1+IS2ru),3.0));
								//printf("if fd for Tensor_4D <0.0 a0ru: %13f\n", a0ru);
								//printf("if fd for Tensor_4D <0.0 a1ru: %13f\n", a1ru);
								//printf("if fd for Tensor_4D <0.0 a2ru: %13f\n", a2ru);

						
								}
					
							Ruu=((a0ru/(a0ru+a1ru+a2ru))*PM0ru)+((a1ru/(a0ru+a1ru+a2ru))*PM1ru)+((a2ru/(a0ru+a1ru+a2ru))*PM2ru);
							//printf("Ruu for Tensor_4D: %13f\n", Ruu); 
							//************************************************************************
							//************************************************************************

							//****
							//   Rj(x-1) where Rj(x-1)=Rj(u-1)
							//****
	
							IS0rv = ((1.0/2.0)*(pow((tensor_4D[0][n][p][vv2index]-tensor_4D[0][n][p][vv3index]),2.0)+pow((tensor_4D[0][n][p][vvindex]-tensor_4D[0][n][p][vv2index]),2.0)))+(pow((tensor_4D[0][n][p][vvindex]-(2.0*tensor_4D[0][n][p][vv2index])+tensor_4D[0][n][p][vv3index]),2.0));
							IS1rv = ((1.0/2.0)*(pow((tensor_4D[0][n][p][vvindex]-tensor_4D[0][n][p][vv2index]),2.0)+pow((tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vvindex]),2.0)))+(pow((tensor_4D[0][n][p][index]-(2.0*tensor_4D[0][n][p][vvindex])+tensor_4D[0][n][p][vv2index]),2.0));
							IS2rv = ((1.0/2.0)*(pow((tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vvindex]),2.0)+pow((tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][index]),2.0)))+(pow((tensor_4D[0][n][p][uuindex]-(2.0*tensor_4D[0][n][p][index])+tensor_4D[0][n][p][vvindex]),2.0));
					
							c01rv = (tensor_4D[0][n][p][vvindex]-2.0*tensor_4D[0][n][p][vv2index]+tensor_4D[0][n][p][vv3index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c02rv = (tensor_4D[0][n][p][vvindex]-tensor_4D[0][n][p][vv3index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c03rv = tensor_4D[0][n][p][vv2index];
							c04rv = (tensor_4D[0][n][p][vvindex]-2.0*tensor_4D[0][n][p][vv2index]+tensor_4D[0][n][p][vv3index])/(24.0);

							c11rv = (tensor_4D[0][n][p][index]-2.0*tensor_4D[0][n][p][vvindex]+tensor_4D[0][n][p][vv2index])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c12rv = (tensor_4D[0][n][p][index]-tensor_4D[0][n][p][vv2index])/(2.0*CCTK_DELTA_SPACE(m-1));
							c13rv = tensor_4D[0][n][p][vvindex];
							c14rv = (tensor_4D[0][n][p][index]-2.0*tensor_4D[0][n][p][vvindex]+tensor_4D[0][n][p][vv2index])/(24.0);
					
							c21rv = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(2.0*pow((CCTK_DELTA_SPACE(m-1)),2.0));
							c22rv = (tensor_4D[0][n][p][uuindex]-tensor_4D[0][n][p][vvindex])/(2.0*CCTK_DELTA_SPACE(m-1));
							c23rv = tensor_4D[0][n][p][index];
							c24rv = (tensor_4D[0][n][p][uuindex]-2.0*tensor_4D[0][n][p][index]+tensor_4D[0][n][p][vvindex])/(24.0);
					
							PM0rv = (c01rv*pow(tensor_4D[0][n][p][vvindex],2.0))+(c02rv*tensor_4D[0][n][p][vvindex])+c03rv+c04rv;
							PM1rv = (c11rv*pow(tensor_4D[0][n][p][vvindex],2.0))+(c12rv*tensor_4D[0][n][p][vvindex])+c13rv+c14rv;
							PM2rv = (c21rv*pow(tensor_4D[0][n][p][vvindex],2.0))+(c22rv*tensor_4D[0][n][p][vvindex])+c23rv+c24rv;
						
							//printf("IS0rv for Tensor_4D: %13f\n", IS0rv);
							//printf("IS1rv for Tensor_4D: %13f\n", IS1rv);
							//printf("IS2rv for Tensor_4D: %13f\n", IS2rv);
					
							//printf("c01rv for Tensor_4D: %13f\n", c01rv);
							//printf("c02rv for Tensor_4D: %13f\n", c02rv);
							//printf("c03rv for Tensor_4D: %13f\n", c03rv);
							//printf("c04rv for Tensor_4D: %13f\n", c04rv);
					
							//printf("c11rv for Tensor_4D: %13f\n", c11rv);
							//printf("c12rv for Tensor_4D: %13f\n", c12rv);
							//printf("c13rv for Tensor_4D: %13f\n", c13rv);
							//printf("c14rv for Tensor_4D: %13f\n", c14rv);
						
							//printf("c21rv for Tensor_4D: %13f\n", c21rv);
							//printf("c22rv for Tensor_4D: %13f\n", c22rv);
							//printf("c23rv for Tensor_4D: %13f\n", c23rv);
							//printf("c24rv for Tensor_4D: %13f\n", c24rv);
						
							//printf("PM0rv for Tensor_4D: %13f\n", PM0rv);
							//printf("PM1rv for Tensor_4D: %13f\n", PM1rv);
							//printf("PM2rv for Tensor_4D: %13f\n", PM2rv);

					
					
					
							//************************************************************************
							// r=3 method requires different alphas depending on the sign of the first derivative of tensor_4D[0][n][p][index]
							//************************************************************************
					
							fdrv = (tensor_4D[0][n][p][index] - tensor_4D[0][n][p][vv2index]) / (2.0*CCTK_DELTA_SPACE(m-1));

							if (fdrv>0.0){

								a0rv = (1.0)/(12.0*pow((e1+IS0rv),3.0));
								a1rv = (1.0)/(2.0*pow((e1+IS1rv),3.0));
								a2rv = (1.0)/(4.0*pow((e1+IS2rv),3.0));
								//printf("if fd for Tensor_4D >0.0 a0rv: %13f\n", a0rv);
								//printf("if fd for Tensor_4D >0.0 a1rv: %13f\n", a1rv);
								//printf("if fd for Tensor_4D >0.0 a2rv: %13f\n", a2rv);

						
								}
						
							if(fdrv<=0.0){
						
								a0rv = (1.0)/(4.0*pow((e1+IS0rv),3.0));
								a1rv = (1.0)/(2.0*pow((e1+IS1rv),3.0));
								a2rv = (1.0)/(12.0*pow((e1+IS2rv),3.0));
								//printf("if fd for Tensor_4D <0.0 a0rv: %13f\n", a0rv);
								//printf("if fd for Tensor_4D <0.0 a1rv: %13f\n", a1rv);
								//printf("if fd for Tensor_4D <0.0 a2rv: %13f\n", a2rv);

						
								}

							Rvv = ((a0rv/(a0rv+a1rv+a2rv))*PM0rv)+((a1rv/(a0rv+a1rv+a2rv))*PM1rv)+((a2rv/(a0rv+a1rv+a2rv))*PM2rv); //*****Updated for R=3
							//printf("Rvv for Tensor_4D: %13f\n", Rvv); 
							//printf("CCTK_DELTA_SPACE(m-1): %13f\n", CCTK_DELTA_SPACE(m-1));
							//************************************************************************
							//************************************************************************
					
							tensor_4D[m][n][p][index]=(Ruu-Rvv)/(2.0*CCTK_DELTA_SPACE(m-1));
							//printf("tensor_4D[m][n][p][index] for Tensor_4D: %13f\n", tensor_4D[m][n][p][index]);

							//************************************************************************
							//************************************************************************
						}

						//****************************
						//** 4th order FD Method
						//****************************
						if (fdtg<=scp){
							////printf("4th Order FD Method Triggered\n");
							tensor_4D[m][n][p][index] = (tensor_4D[0][n][p][vv2index] - 8.0*tensor_4D[0][n][p][vvindex] + 8.0*tensor_4D[0][n][p][uuindex] - tensor_4D[0][n][p][uu2index]) / (12.0*CCTK_DELTA_SPACE(m-1));
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
