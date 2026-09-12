#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ExternalVariables.h"

void MoL_Adams_Adapt(CCTK_ARGUMENTS);
void MoL_AMAdd(CCTK_ARGUMENTS);
void Adams_setstep(CCTK_ARGUMENTS);
void MoL_RK_Adapt(CCTK_ARGUMENTS);
void RK_setstep(CCTK_ARGUMENTS);
void RK_ReduceAdaptiveError(CCTK_ARGUMENTS);

/* Adams-Moulton */

static void order1 (CCTK_REAL* restrict const OldVar,
                    CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
#pragma omp parallel for simd
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += dt * RHSVar0[index];
  }
}

static void order2 (CCTK_REAL* restrict const OldVar,
                    CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const f0 = + (1.0/2.0) * dt;
  CCTK_REAL const f1 = + (1.0/2.0) * dt;
#pragma omp parallel for simd
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += f0 * RHSVar0[index] + f1 * RHSVar1[index];
  }
}

static void order3 (CCTK_REAL* restrict const OldVar,
                    CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const f0 = + (5.0/12.0) * dt;
  CCTK_REAL const f1 = + (8.0/12.0) * dt;
  CCTK_REAL const f2 = - (1.0/12.0) * dt;
#pragma omp parallel for simd
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index];
  }
}

static void order4 (CCTK_REAL* restrict const OldVar,
                    CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const* restrict const RHSVar3 = RHSVars[3];
  CCTK_REAL const f0 = + ( 9.0/24.0) * dt;
  CCTK_REAL const f1 = + (19.0/24.0) * dt;
  CCTK_REAL const f2 = - ( 5.0/24.0) * dt;
  CCTK_REAL const f3 = + ( 1.0/24.0) * dt;
#pragma omp parallel for simd
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index] +
      f3 * RHSVar3[index];
  }
}

static void order5 (CCTK_REAL* restrict const OldVar,
                    CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const* restrict const RHSVar3 = RHSVars[3];
  CCTK_REAL const* restrict const RHSVar4 = RHSVars[4];
  CCTK_REAL const f0 = + (251.0/720.0) * dt;
  CCTK_REAL const f1 = + (646.0/720.0) * dt;
  CCTK_REAL const f2 = - (264.0/720.0) * dt;
  CCTK_REAL const f3 = + (106.0/720.0) * dt;
  CCTK_REAL const f4 = - ( 19.0/720.0) * dt;
#pragma omp parallel for simd
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index] +
      f3 * RHSVar3[index] + f4 * RHSVar4[index];
  }
}

/* Array of function pointers */
static
void (* const orders[]) (CCTK_REAL* restrict const OldVar,
                         CCTK_REAL* restrict const UpdateVar,
                         CCTK_REAL const* restrict *restrict const RHSVars,
                         CCTK_REAL const dt,
                         int const totalsize)
= { order1, order2, order3, order4, order5 };
static int const max_order = sizeof orders / sizeof *orders;

void MoL_AMAdd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL const dt = CCTK_DELTA_TIME;
  
  /* Determine the order of accuracy */
  int order;
  if (CCTK_EQUALS(AB_Type,"1")) {
    order = 1 ;
  } else if (CCTK_EQUALS(AB_Type,"2")) {
    order = 2;
  } else if (CCTK_EQUALS(AB_Type,"3")) {
    order = 3;
  } else if (CCTK_EQUALS(AB_Type,"4")) {
    order = 4;
  } else if (CCTK_EQUALS(AB_Type,"5")) {
    order = 5;
  } else {
    abort();
  }
  if (AB_initially_reduce_order) {
    /* Reduce the order for the first time steps */
    int const iteration = cctk_iteration; 
    /*int iteration = 1 + lrint((cctk_time - cctk_initial_time) / dt);*/
    /*if ((order > iteration) && (interation <= max_order) && (iteration >= 1)) {*/
    if (order > iteration) {
      order = iteration;
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Reducing Adams order to %d", order);
    }
  }
  //printf ("MoL Adams: iter=%d, order=%d\n", cctk_iteration, order);
  assert (order >= 1 && order <= max_order);
  
  int totalsize = 1;
  for (int arraydim = 0; arraydim < cctk_dim; arraydim++) {
    totalsize *= cctk_ash[arraydim];
  }
  
  for (int var = 0; var < MoLNumEvolvedVariables; var++) {
    CCTK_REAL* restrict const OldVar =
      CCTK_VarDataPtrI(cctkGH, 1, EvolvedVariableIndex[var]);
    CCTK_REAL* restrict const UpdateVar =
      CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
    CCTK_REAL const* restrict RHSVars[order];
    for (int tl = 0; tl < order; tl++) {
      RHSVars[tl] = CCTK_VarDataPtrI(cctkGH, tl, RHSVariableIndex[var]);
    }
    
    /* Add RHS */
    (orders[order-1]) (OldVar, UpdateVar, RHSVars, dt, totalsize);
  }/* var */
  
  for (int var = 0; var < MoLNumEvolvedArrayVariables; var++) {
    CCTK_REAL* restrict const OldVar =
      CCTK_VarDataPtrI(cctkGH, 1, EvolvedArrayVariableIndex[var]);
    CCTK_REAL* restrict const UpdateVar =
      CCTK_VarDataPtrI(cctkGH, 0, EvolvedArrayVariableIndex[var]);
    CCTK_REAL const* restrict RHSVars[order];
    for (int tl = 0; tl < order; tl++) {
      RHSVars[tl] = CCTK_VarDataPtrI(cctkGH, tl, RHSArrayVariableIndex[var]);
    }
    
    int const groupindex =
      CCTK_GroupIndexFromVarI(EvolvedArrayVariableIndex[var]);
    cGroupDynamicData arraydata;
    int const ierr = CCTK_GroupDynamicData(cctkGH, groupindex, &arraydata);
    if (ierr) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, 
                  "The driver does not return group information for group '%s'.",
                  CCTK_GroupName(groupindex));
    }
    int arraytotalsize = 1;
    for (int arraydim = 0; arraydim < arraydata.dim; arraydim++) {
      arraytotalsize *= arraydata.ash[arraydim];
    }
    
    /* Add RHS */
    (orders[order-1]) (OldVar, UpdateVar, RHSVars, dt, totalsize);
  }/* var */
}

void MoL_Adams_Adapt(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr, handle, index, i, j, k;  
  CCTK_REAL max_error = 0.0;
  CCTK_REAL const dt = CCTK_DELTA_TIME;
  CCTK_REAL chg_fac;
  CCTK_REAL diff, power; 
  
  /* potentially change dt here */
  
  if (CCTK_EQUALS(AB_Type,"1")) {
    power = 1.0/2.0;
  } else if (CCTK_EQUALS(AB_Type,"2")) {
    power = 1.0/3.0;
  } else if (CCTK_EQUALS(AB_Type,"3")) {
    power = 1.0/4.0;
  } else if (CCTK_EQUALS(AB_Type,"4")) {
    power = 1.0/5.0;
  } else if (CCTK_EQUALS(AB_Type,"5")) {
    power = 1.0/6.0;
  } else {
    abort();
  }
  
  if (CCTK_Equals(AB_Step_Type,"Average")) {
  	handle = CCTK_ReductionHandle("average");
  }
  if (CCTK_Equals(AB_Step_Type,"Maximum")) {
  	handle = CCTK_ReductionHandle("maximum");
  }
  if (CCTK_Equals(AB_Step_Type,"Minimum")) {
  	handle = CCTK_ReductionHandle("minimum");
  }
  if (CCTK_Equals(AB_Step_Type,"Norm")) {
  	handle = CCTK_ReductionHandle("norm2");
  }

    
  for (int var = 0; var < MoLNumEvolvedVariables; var++) {

    CCTK_REAL* restrict const OldVar = CCTK_VarDataPtrI(cctkGH, 1, EvolvedVariableIndex[var]);
    CCTK_REAL* restrict const UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);

    for(k=cctk_nghostzones[2]; k < cctk_lsh[2]-cctk_nghostzones[2]; k++)
	{
		for(j=cctk_nghostzones[1]; j < cctk_lsh[1]-cctk_nghostzones[1]; j++)
		{
			for(i=cctk_nghostzones[0]; i < cctk_lsh[0]-cctk_nghostzones[0]; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );

			if ((OldVar[index] != 0.0) || (UpdateVar[index] != 0.0)) {
				temp_diff[index] = (UpdateVar[index] - OldVar[index])/(fabs(OldVar[index]) + fabs(UpdateVar[index]));
			} else {
				temp_diff[index] = 0.0;
			}

			}
		}
	}
 

    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &diff, 1, CCTK_VarIndex("MoL::temp_diff"));
    
      if (diff > max_error) {
    	  max_error = diff;
      }

  }/* var */
  
  for (int var = 0; var < MoLNumEvolvedArrayVariables; var++) {
    
    CCTK_REAL* restrict const OldVar = CCTK_VarDataPtrI(cctkGH, 1, EvolvedArrayVariableIndex[var]);
    CCTK_REAL* restrict const UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedArrayVariableIndex[var]);

    for(k=cctk_nghostzones[2]; k < cctk_lsh[2]-cctk_nghostzones[2]; k++)
	{
		for(j=cctk_nghostzones[1]; j < cctk_lsh[1]-cctk_nghostzones[1]; j++)
		{
			for(i=cctk_nghostzones[0]; i < cctk_lsh[0]-cctk_nghostzones[0]; i++)
			{
			index = CCTK_GFINDEX3D( cctkGH, i, j, k );

			if ((OldVar[index] != 0.0) || (UpdateVar[index] != 0.0)) {
				temp_diff[index] = (UpdateVar[index] - OldVar[index])/(fabs(OldVar[index]) + fabs(UpdateVar[index]));
			} else {
				temp_diff[index] = 0.0;
			}

			}
		}
	}

    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &diff, 1, CCTK_VarIndex("MoL::temp_diff"));
    
      if (diff > max_error) {
    	  max_error = diff;
      }
      
  }/* var */
  
 if (max_error != 0.0) {
    chg_fac = pow(maximum_relative_error/max_error,power);
    
if (limit_AB_steps) {
    if (chg_fac > maximum_increase) {
        chg_fac = maximum_increase;
    }
    
    if (chg_fac < maximum_decrease) {
        chg_fac = maximum_decrease;
    }
}
    
if (big_AB_steps) {
    if ((chg_fac < maximum_increase) && (chg_fac > maximum_decrease)) {
    	chg_fac = 1.0;
    }
}
    
    *NewDT = chg_fac*CCTK_DELTA_TIME;
 }
 
 if (print_adams_step) {
     printf ("max_error = %e, time = %e, time step = %e\n", max_error, cctk_time, *NewDT);
 }

}

void Adams_setstep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if ((*NewDT < maximum_stepsize) && (*NewDT > minimum_stepsize)){
  cctkGH->cctk_delta_time = *NewDT;
  }
  //CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", (double)cctkGH->cctk_delta_time);
  
}

void MoL_RK_Adapt(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT i, j, k, index;  
  CCTK_REAL total;
  
  assert (cctk_dim <= 3);
  int imin[3], imax[3];
  for (int d = 0; d < cctk_dim; d++)
  {
    imin[d] = cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - cctk_nghostzones[d];
  }
  
  /* potentially change dt here */
  

    
  for (int var = 0; var < MoLNumEvolvedVariables; var++) {

    CCTK_REAL const * restrict const OldVar = CCTK_VarDataPtrI(cctkGH, 1, EvolvedVariableIndex[var]);
    CCTK_REAL const * restrict const UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
    CCTK_REAL const * restrict const ErrorVar = CCTK_VarDataPtrI(cctkGH, var, CCTK_FirstVarIndex("MOL::ERRORESTIMATE"));

    for(k=imin[2]; k < imax[2]; k++)
	{
		for(j=imin[1]; j < imax[1]; j++)
		{
			for(i=imin[0]; i < imax[0]; i++)
			{
			index = CCTK_GFINDEX3D(cctkGH, i, j, k);
			
			total = fabs(OldVar[index]) + fabs(UpdateVar[index]);

			if ((total != 0.0) && (fabs(ErrorVar[index]) < total)) {
				temp_diff[index] = ErrorVar[index]/total;
			} else {
				temp_diff[index] = 0.0;
			}
			
			//if (fabs(temp_diff[index]) > 1.0e-20) {
			//printf ("temp_diff = %4.3e, ErrorVar = %4.3e, OldVar = %4.3e, UpdateVar = %4.3e\n", temp_diff[index], ErrorVar[index], OldVar[index], UpdateVar[index]);
			//}

		  }
	   }
	}

  }/* var */
  
  for (int var = 0; var < MoLNumEvolvedArrayVariables; var++) {
    
    CCTK_REAL const * restrict const OldVar = CCTK_VarDataPtrI(cctkGH, 1, EvolvedArrayVariableIndex[var]);
    CCTK_REAL const * restrict const UpdateVar = CCTK_VarDataPtrI(cctkGH, 0, EvolvedArrayVariableIndex[var]);
    CCTK_REAL const * restrict const ErrorVar = CCTK_VarDataPtrI(cctkGH, var, CCTK_FirstVarIndex("MOL::ERRORESTIMATE"));

    for(k=imin[2]; k < imax[2]; k++)
	{
		for(j=imin[1]; j < imax[1]; j++)
		{
			for(i=imin[0]; i < imax[0]; i++)
			{
			index = CCTK_GFINDEX3D(cctkGH, i, j, k);
			
			total = fabs(OldVar[index]) + fabs(UpdateVar[index]);

			if ((total != 0.0) && (fabs(ErrorVar[index]) < total)) {
				temp_diff[index] = (ErrorVar[index])/total;
			} else {
				temp_diff[index] = 0.0;
			}

		   }
	    }
	}
  }
 
 }
 
 void RK_ReduceAdaptiveError(CCTK_ARGUMENTS)
 {
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr, handle;  
  CCTK_REAL max_error = 0.0;
  CCTK_REAL chg_fac, total;
  CCTK_REAL diff, power; 
  
  if (CCTK_EQUALS(ODE_Method,"RK45")) {
    power = 1.0/5.0;
  } else if (CCTK_EQUALS(ODE_Method,"RK45CK")) {
    power = 1.0/5.0;
  } else if (CCTK_EQUALS(ODE_Method,"RK65")) {
    power = 1.0/6.0;
  } else if (CCTK_EQUALS(ODE_Method,"RK87")) {
    power = 1.0/8.0;
  } else {
    abort();
  }
    
  if (CCTK_Equals(AB_Step_Type,"Average")) {
  	handle = CCTK_ReductionHandle("average");
  }
  if (CCTK_Equals(AB_Step_Type,"Maximum")) {
  	handle = CCTK_ReductionHandle("maximum");
  }
  if (CCTK_Equals(AB_Step_Type,"Minimum")) {
  	handle = CCTK_ReductionHandle("minimum");
  }
  if (CCTK_Equals(AB_Step_Type,"Norm")) {
  	handle = CCTK_ReductionHandle("norm2");
  }
 

    ierr = CCTK_Reduce(cctkGH, -1, handle, 1, CCTK_VARIABLE_REAL, &diff, 1, CCTK_VarIndex("MoL::temp_diff"));
    
      if (diff > max_error) {
    	  max_error = diff;
      }
      
  /* var */
 
 chg_fac = 1.0; 
  
 if (max_error != 0.0) {
    chg_fac = pow(maximum_relative_error/max_error,power);
    
if (limit_AB_steps) {
    if (chg_fac > maximum_increase) {
        chg_fac = maximum_increase;
    }
    
    if (chg_fac < maximum_decrease) {
        chg_fac = maximum_decrease;
    }
}
    
if (big_AB_steps) {
    if ((chg_fac < maximum_increase) && (chg_fac > maximum_decrease)) {
    	chg_fac = 1.0;
    }
}	
    *NewDT = fabs(chg_fac*CCTK_DELTA_TIME);
 }
 
 if (print_adams_step) {
     printf ("max_error = %4.3e, old time step = %4.3e, new time step = %4.3e, chg_fac = %4.3e\n", max_error, CCTK_DELTA_TIME, *NewDT, chg_fac);
 }

}

void RK_setstep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if ((*NewDT < maximum_stepsize) && (*NewDT > minimum_stepsize)){
  cctkGH->cctk_delta_time = *NewDT;
  }
  //CCTK_VInfo (CCTK_THORNSTRING, "Setting time step to %g", (double)cctkGH->cctk_delta_time);
  
}


