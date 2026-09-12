/*@@ Sets the symmetries for GRMHD @@*/

#include "cctk.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

static const char *rcsid = "$Header: $";

CCTK_FILEVERSION(SpecCosmo_SpecGRMHD_InitSymBound_c)

void SpecGRMHD_InitSymBound(CCTK_ARGUMENTS);

void SpecGRMHD_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
      
  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"SpecGRMHD::rho_star");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::rho_star_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::alpha");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::tau");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::tau_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::psib0");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::phi");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::EErad");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::TrK");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Axx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Ayy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Azz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betax1");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betay2");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betaz3");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Speed_Sound");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::sqrtdetg");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::epsilon");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::rho");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::pp");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uu_t");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::div_vv");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::beta_t");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::divb");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::rho_total");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttt_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uud_t");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttt_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::rho_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::trQ"); 
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txx_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tyy_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tzz_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txx_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tyy_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tzz_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_xx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_yy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_zz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_xx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_yy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_zz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::ham");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::AATF");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::detg1");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Axy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betax2");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betay1");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txy_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txy_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_xy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_xy");
  
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = -1;

  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Axz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betax3");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betaz1");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txz_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Txz_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_xz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_xz");
  
  sym[0] = 1;
  sym[1] = -1;
  sym[2] = -1;

  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_Ayz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betay3");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betaz2");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tyz_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tyz_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Qyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ricci_yz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Sij_yz");
      
  sym[0] = -1;
  sym[1] = 1;
  sym[2] = 1;
  
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::S_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::SD_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bb_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::alpha_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betax");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_bb_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dk_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gamma_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uu_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vv_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::b_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::beta_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttx_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uud_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vvd_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttx_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yyx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zzx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Si_x");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::momx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::dGx");

  sym[0] = 1;
  sym[1] = -1;
  sym[2] = 1;
  
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::S_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::SD_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bb_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::alpha_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betay");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_bb_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dk_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gamma_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uu_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vv_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::b_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::beta_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tty_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uud_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vvd_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Tty_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zzy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xyx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Si_y");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::momy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::dGy");

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::S_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::SD_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bb_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::alpha_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::betaz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_bb_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dk_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_gamma_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uu_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vv_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::b_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::beta_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttz_mhd");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::uud_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::vvd_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Ttz_dark");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zxx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zyy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zzz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xzx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yzy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Si_z");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::momz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::dGz");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_zxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_yxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::bssn_dkij_xyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zxy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yxz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xyz");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_xzy");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_yzx");
  SetCartSymVN(cctkGH, sym,"SpecGRMHD::Chris_zyx");

  return;
}
