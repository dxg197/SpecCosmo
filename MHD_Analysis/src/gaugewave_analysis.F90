! Simple_Init
! gaugewave: Analytic solution for Goiel's gravitational waves:
! Provide initial data and gauge conditions for nonlinear gauge gravitational waves.
 
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module gaugewave_analysis
  implicit none
  private
  public makegaugewave
  public makegwavegauge
  public calcgwavegaugerhs
  !
contains
  !
  subroutine makegaugewave (xx, time, lz, g, kk)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)  :: xx(3), time, lz
    CCTK_REAL, intent(out) :: g(3,3), kk(3,3)
    CCTK_REAL wc, pi
    integer i
    !
    pi = 4.0d0*atan(1.0d0)
    wc = nw*2.0d0*pi/lz
    !
    g(:,:) = 0.0d0 
    forall (i=1:3) g(i,i) = 1.0d0 
    !
    kk(:,:) = 0.0d0
    !
    g(1,1) = exp(Amp*sin(wc*(time+xx(3))))
    kk(1,1) = -0.5d0*wc*exp(0.5d0*Amp*sin(wc*(time+xx(3))))*Amp*cos(wc*(time+xx(3)))
    !
  end subroutine makegaugewave
  !
  subroutine makegwavegauge (xx, time, lz, alph, beta)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)  :: xx(3), time, lz
    CCTK_REAL, intent(out) :: alph, beta(3)
    CCTK_REAL wc, pi
    !
    pi = 4.0d0*atan(1.0d0)
    wc = nw*2.0d0*pi/lz
    !
    alph = exp(0.5d0*Amp*sin(wc*(time+xx(3))))
    !alph_z = 0.5d0*Amp*wc*cos(wc*(time+xx(3)))*exp(0.5d0*Amp*sin(wc*(time+xx(3))))
    !
    ! betai = 0
    beta(:) = 0.d0
    !
  end subroutine makegwavegauge
  !
  subroutine calcgwavegaugerhs (xx, time, lz, alph_dot, beta_dot)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)  :: xx(3), time, lz
    CCTK_REAL, intent(out) :: alph_dot, beta_dot(3)
    CCTK_REAL wc, pi
    !
    pi = 4.0d0*atan(1.0d0)
    wc = nw*2.0d0*pi/lz
    !
    alph_dot = 0.5d0*wc*exp(0.5d0*Amp*sin(pi/2.0d0+wc*(time+xx(3))))*Amp*cos(pi/2.0d0+wc*(time+xx(3)))
    !
    ! betai,0 = betai,t = 0
    beta_dot(:) = 0.0d0
  end subroutine calcgwavegaugerhs
    !
end module gaugewave_analysis
