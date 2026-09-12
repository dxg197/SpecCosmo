#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h> 

void Calculate_A(CCTK_ARGUMENTS, CCTK_REAL ***vector );

void B_Clean(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT i,j,k,m,n,p,four;
  CCTK_INT istart,jstart,kstart,iend,jend,kend;
  CCTK_INT sizex, sizey, sizez, index;
  
  /* Declare Arrays */

  CCTK_REAL ***bb3;

/* Set up shorthands */

    sizex  = cctk_lsh[0]; 
    sizey  = cctk_lsh[1];
    sizez  = cctk_lsh[2];
	istart = cctk_nghostzones[0];
	jstart = cctk_nghostzones[1];
	kstart = cctk_nghostzones[2];	
	iend   = cctk_lsh[0] - 2.0*cctk_nghostzones[0];
	jend   = cctk_lsh[1] - 2.0*cctk_nghostzones[1];
	kend   = cctk_lsh[2] - 2.0*cctk_nghostzones[2];
    
    four = 4;
      
    bb3 = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **)*four);
                  
    for(m=0; m < 4; m++) {
       bb3[m] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *)*four);
       
       for(n=0; n < 4; n++) {
            bb3[m][n] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*sizex*sizey*sizez);
       }
    }
    
	for(k=kstart; k < sizez; k++)
	{
		for(j=jstart; j < sizey; j++)
		{
			for(i=istart; i < sizex; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );
        
        	bb3[0][1][index] = bb_x[index];
    		bb3[0][2][index] = bb_y[index];
    		bb3[0][3][index] = bb_z[index];
        
        	}
		}
	}
    
    Calculate_A(CCTK_PASS_CTOC, bb3);
    
	for(k=0; k < sizez; k++)
	{
		for(j=0; j < sizey; j++)
		{
			for(i=0; i < sizex; i++)
			{
  
             index = CCTK_GFINDEX3D( cctkGH, i, j, k );
             
    		bb_x[index] = bb3[0][1][index];
 			bb_y[index] = bb3[0][2][index];
 			bb_z[index] = bb3[0][3][index];
 			
 			aa_x[index] = bb3[1][1][index];
 			aa_y[index] = bb3[2][2][index];
 			aa_z[index] = bb3[3][3][index];
 			
 			}
		}
	}

	for(m=0; m < 4; m++) {
      for(n=0; n < 4; n++) {
        free(bb3[m][n]);
      }
    }
    
    for(m=0; m < 4; m++) {
      free(bb3[m]);
    }
   
      free(bb3);
      	
}