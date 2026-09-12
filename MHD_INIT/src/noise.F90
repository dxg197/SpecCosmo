! MHD_Init
! noise: Add noise to initial data
! Add some noise to the metric and extrinsic curvature.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module noise
  implicit none
  private
  public makenoise
  public makegaugenoise
  !
contains
  !
  subroutine makenoise(g, kk)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(inout) :: g(3,3), kk(3,3)
    !
    CCTK_REAL rand
    integer i,j
    !
    if (noise_amplitude_metric > 0) then
       do i=1,3
          do j=i,3! only for j>=i
             call RANDOM_NUMBER(rand)
             g(i,j) = g(i,j) + 2 * noise_amplitude_metric * (rand - 0.5)
          end do
       end do
       !
       ! symmetries
       g(2,1) = g(1,2)
       g(3,1) = g(1,3)
       g(3,2) = g(2,3)
    end if
    !
    if (noise_amplitude_curv > 0) then
       do i=1,3
          do j=i,3! only for j>=i
             call RANDOM_NUMBER(rand)
             kk(i,j) = kk(i,j) + 2 * noise_amplitude_curv * (rand - 0.5)
          end do
       end do
       !
       ! symmetries
       kk(2,1) = kk(1,2)
       kk(3,1) = kk(1,3)
       kk(3,2) = kk(2,3)
    end if
    !
  end subroutine makenoise
  !
  subroutine makegaugenoise(alph, beta)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(inout) :: alph, beta(3)
    !
    CCTK_REAL rand
    integer i
    !
    if (noise_amplitude_lapse > 0) then
       call RANDOM_NUMBER(rand)
       alph = alph + 2 * noise_amplitude_lapse * (rand - 0.5)
    end if
    !
    if (noise_amplitude_shift > 0) then
       do i=1,3
          call RANDOM_NUMBER(rand)
          beta(i) = beta(i) + 2 * noise_amplitude_shift * (rand - 0.5)
       end do
    end if
    !
  end subroutine makegaugenoise
  !
end module noise
