! MHD_Anlaysis

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



module singular
  use coords_analysis
  use constants_analysis
  implicit none
  private
  public makesingular
  !
contains
  !
  subroutine makesingular(alph, alph_x, alph_y, alph_z, beta, oldlapse)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in) :: oldlapse
    CCTK_REAL, intent(out) :: alph, alph_x, alph_y, alph_z, beta(3)
    !
    alph = oldlapse
    alph_x = 0.0
    alph_y = 0.0
    alph_z = 0.0
    !
    ! betai = 0
    beta(:) = 0
    !
  end subroutine makesingular
  !
end module singular
