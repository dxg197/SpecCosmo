 /*@@
   @file      MoLRegister.c
   @date      Wed Jan  9 23:01:36 2002
   @author    
   @desc 
   Routine to register the variables with the MoL thorn.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/* Prototypes */
void GRMHD_MoLRegister(CCTK_ARGUMENTS);

void GRMHD_MoLRegister(CCTK_ARGUMENTS) 
{

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr = 0,  metricg, inv_metricg, curvg, var, var_Misc, cons, temp1, temp2;

  if (! fix_lapse) {
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::alpha"), CCTK_VarIndex("SpecGRMHD::alpha_dt"));
  }
  
  if (! fix_shift) {
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::betax"), CCTK_VarIndex("SpecGRMHD::betax_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::betay"), CCTK_VarIndex("SpecGRMHD::betay_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::betaz"), CCTK_VarIndex("SpecGRMHD::betaz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_bb_x"), CCTK_VarIndex("SpecGRMHD::bssn_bb_x_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_bb_y"), CCTK_VarIndex("SpecGRMHD::bssn_bb_y_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_bb_z"), CCTK_VarIndex("SpecGRMHD::bssn_bb_z_dt"));
  }
  
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::phi"), CCTK_VarIndex("SpecGRMHD::phi_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gxx"), CCTK_VarIndex("SpecGRMHD::bssn_gxx_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gxy"), CCTK_VarIndex("SpecGRMHD::bssn_gxy_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gxz"), CCTK_VarIndex("SpecGRMHD::bssn_gxz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gyy"), CCTK_VarIndex("SpecGRMHD::bssn_gyy_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gyz"), CCTK_VarIndex("SpecGRMHD::bssn_gyz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gzz"), CCTK_VarIndex("SpecGRMHD::bssn_gzz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Axx"), CCTK_VarIndex("SpecGRMHD::bssn_Axx_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Axy"), CCTK_VarIndex("SpecGRMHD::bssn_Axy_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Axz"), CCTK_VarIndex("SpecGRMHD::bssn_Axz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Ayy"), CCTK_VarIndex("SpecGRMHD::bssn_Ayy_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Ayz"), CCTK_VarIndex("SpecGRMHD::bssn_Ayz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_Azz"), CCTK_VarIndex("SpecGRMHD::bssn_Azz_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::TrK"), CCTK_VarIndex("SpecGRMHD::bssn_K_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gamma_x"), CCTK_VarIndex("SpecGRMHD::bssn_gamma_x_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gamma_y"), CCTK_VarIndex("SpecGRMHD::bssn_gamma_y_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bssn_gamma_z"), CCTK_VarIndex("SpecGRMHD::bssn_gamma_z_dt"));
                          
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::rho_star"), CCTK_VarIndex("SpecGRMHD::rho_star_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::tau"), CCTK_VarIndex("SpecGRMHD::tau_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::S_x"), CCTK_VarIndex("SpecGRMHD::S_x_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::S_y"), CCTK_VarIndex("SpecGRMHD::S_y_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::S_z"), CCTK_VarIndex("SpecGRMHD::S_z_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bb_x"), CCTK_VarIndex("SpecGRMHD::bb_x_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bb_y"), CCTK_VarIndex("SpecGRMHD::bb_y_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::bb_z"), CCTK_VarIndex("SpecGRMHD::bb_z_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::psib0"), CCTK_VarIndex("SpecGRMHD::psib_dt"));
  
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::rho_star_dark"), CCTK_VarIndex("SpecGRMHD::rho_star_dark_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::tau_dark"), CCTK_VarIndex("SpecGRMHD::tau_dark_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::SD_x"), CCTK_VarIndex("SpecGRMHD::SD_x_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::SD_y"), CCTK_VarIndex("SpecGRMHD::SD_y_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::SD_z"), CCTK_VarIndex("SpecGRMHD::SD_z_dt"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("SpecGRMHD::EErad"), CCTK_VarIndex("SpecGRMHD::EErad_dt"));
  
  metricg = CCTK_GroupIndex("SpecGRMHD::metric");
  inv_metricg = CCTK_GroupIndex("SpecGRMHD::inv_metric");
  curvg = CCTK_GroupIndex("SpecGRMHD::curv");
  var = CCTK_GroupIndex("SpecGRMHD::calcvars");
  var_Misc = CCTK_GroupIndex("SpecGRMHD::bssn_variables_Misc");
  cons = CCTK_GroupIndex("SpecGRMHD::adm_bssn_cons");
  temp1 = CCTK_GroupIndex("SpecGRMHD::temp_variables1");
  temp2 = CCTK_GroupIndex("SpecGRMHD::temp_variables2");
   
  if (CCTK_IsFunctionAliased("MoLRegisterConstrained"))
  {
    //ierr += MoLRegisterConstrainedGroup(var);
    //ierr += MoLRegisterConstrainedGroup(var_Misc);
    //ierr += MoLRegisterConstrainedGroup(cons);
  }
  else
  {
    CCTK_WARN(0, "MoL function MoLRegisterConstrained not aliased");
    ierr++;
  } 
  
  if (CCTK_IsFunctionAliased("MoLRegisterSaveAndRestoreGroup"))
  {
    //ierr += MoLRegisterSaveAndRestoreGroup(temp1);
    //ierr += MoLRegisterSaveAndRestoreGroup(temp2);
  }
  else
  {
    CCTK_WARN(0, "MoL function MoLRegisterSaveAndRestoreGroup not aliased");
    ierr++;
  }
  
  if (ierr) CCTK_WARN(0,"Problems registering variables with MoL");
}
