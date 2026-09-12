#include "cctk.h"

module constants_init
  implicit none
  private
  
  public delta2, delta3, delta4
  public eta4
  public epsilon2, epsilon3, epsilon4
  
  public pi
  
  public Grav_Const, Speed_Light, Rad_Const, Boltz_Const, Eps_Naught, Solar_Mass, Proton_Mass
  public Kilogram, Second, Meter, Kelvin, Density, SIE, Pressure, Ampere, Tesla, Gauss, E_Field
  public G_Kilogram, G_Second, G_Meter, G_Kelvin, G_Density, G_SIE, G_Pressure, G_Amp, G_Tesla, G_Gauss, G_E_Field
  public cc, GC, kB, cfC
  
  CCTK_REAL, parameter :: zero = 0
  
  integer, parameter :: rk = kind(zero)
  
  
  CCTK_REAL, parameter :: delta2(2,2) &
       = reshape((/ 1,0,   0,1 /), (/2,2/))
  
  CCTK_REAL, parameter :: delta3(3,3) &
       = reshape((/ 1,0,0,   0,1,0,   0,0,1 /), (/3,3/))
  
  CCTK_REAL, parameter :: delta4(0:3,0:3) &
       = reshape((/ 1,0,0,0,   0,1,0,0,   0,0,1,0,   0,0,0,1 /), (/4,4/))
  
  
  
  CCTK_REAL, parameter :: eta4(0:3,0:3) &
       = reshape((/ -1,0,0,0,   0,1,0,0,   0,0,1,0,   0,0,0,1 /), (/4,4/))
  
  
  
  CCTK_REAL, parameter :: epsilon2(2,2) &
       = reshape ((/ 0,+1,   -1,0 /), (/2,2/))
  
  CCTK_REAL, parameter :: epsilon3(3,3,3) &
       = reshape ((/ 0,0,0,    0,0,+1,   0,-1,0, &
       &             0,0,-1,   0,0,0,    +1,0,0, &
       &             0,+1,0,   -1,0,0,   0,0,0   /), (/3,3,3/))
  
  CCTK_REAL, parameter :: epsilon4(0:3,0:3,0:3,0:3) &
       = reshape ((/ 0,0,0,0,    0,0,0,0,    0,0,0,0,   0,0,0,0,  &
       &             0,0,0,0,    0,0,0,0,    0,0,0,+1,  0,0,-1,0, &
       &             0,0,0,0,    0,0,0,-1,   0,0,0,0,   0,+1,0,0, &
       &             0,0,0,0,    0,0,+1,0,   0,-1,0,0,  0,0,0,0,  &
       
       &             0,0,0,0,    0,0,0,0,    0,0,0,-1,  0,0,+1,0, &
       &             0,0,0,0,    0,0,0,0,    0,0,0,0,   0,0,0,0,  &
       &             0,0,0,+1,   0,0,0,0,    0,0,0,0,   -1,0,0,0, &
       &             0,0,-1,0,   0,0,0,0,    +1,0,0,0,  0,0,0,0,  &
       
       &             0,0,0,0,    0,0,0,+1,   0,0,0,0,   0,-1,0,0, &
       &             0,0,0,-1,   0,0,0,0,    0,0,0,0,   +1,0,0,0, &
       &             0,0,0,0,    0,0,0,0,    0,0,0,0,   0,0,0,0,  &
       &             0,+1,0,0,   -1,0,0,0,   0,0,0,0,   0,0,0,0,  &
       
       &             0,0,0,0,    0,0,-1,0,   0,+1,0,0,  0,0,0,0,  &
       &             0,0,+1,0,   0,0,0,0,    -1,0,0,0,  0,0,0,0,  &
       &             0,-1,0,0,   +1,0,0,0,   0,0,0,0,   0,0,0,0,  &
       &             0,0,0,0,    0,0,0,0,    0,0,0,0,   0,0,0,0   /), &
       &          (/4,4,4,4/))
  
  
  
  CCTK_REAL, parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068_rk
  
! List constants for conversions

  CCTK_REAL, parameter :: Grav_Const = 6.6742867D-11       !  units are meter^3 kilogram^-1 second^-2

  CCTK_REAL, parameter :: Speed_Light = 2.99792458D8       !  units are meter second^-1

  CCTK_REAL, parameter :: Rad_Const   = 7.565767D-16       !  units are Joules meter^-3 Kelvin^-4

  CCTK_REAL, parameter :: Boltz_Const = 1.380650424D-23    !  units are kilogram meter^2 second^-2 Kelvin^-1

  CCTK_REAL, parameter :: Eps_Naught = 8.854187817D-12     !  units are Amps^2 second^4 kilogram^-1 meter^-3

  CCTK_REAL, parameter :: Solar_Mass = 1.9989225D30        !  units are kilogram
  
  CCTK_REAL, parameter :: Proton_Mass = 1.67262163783D-27  !  units are kilogram
  
!/************************************************************************************************************
!* SI Variables in terms of Geometrized Units; by J. Hamilton and C. Ballard                                 *
!*	for more information on actual variable manipulation                                                     * 
!*	see C. Ballard Fall 2005 Research Paper                                                                  *
!*	 (Calculated by Mathematica)                                                                             *
!************************************************************************************************************/

  CCTK_REAL, parameter :: Kilogram = 2.4770954945D-36      !  units are seconds

  CCTK_REAL, parameter :: Second = 1.0                     !  units are seconds

  CCTK_REAL, parameter :: Meter = 3.33564095198D-9         !  units are seconds

  CCTK_REAL, parameter :: Kelvin = 3.8052664687D-76        !  units are seconds

  CCTK_REAL, parameter :: Density = 6.6742867D-11          !  units are seconds^-2
             
                                                           !  Specific Internal Energy
  CCTK_REAL, parameter :: SIE = 1.11265006D-17             !  unitless in Geometrized Units

  CCTK_REAL, parameter :: Pressure = 7.4261454669D-28      !  units are seconds^-2

  CCTK_REAL, parameter :: Ampere = 2.87449146D-26          !  unitless in Geometrized Units

  CCTK_REAL, parameter :: Tesla = 8.6175086029D-11         !  units are seconds^-1
  
  CCTK_REAL, parameter :: Gauss = 8.6175086029D-7          !  units are seconds^-1

  CCTK_REAL, parameter :: E_Field = 2.874491463D-19        !  units are seconds^-1
  
  
!/************************************************************************************************************
!* Geometrized Variables in terms of SI units; by C. Ballard and J.Hamilton                                  *
!*	for more information on actual variable manipulation                                                     * 
!*	see C. Ballard Fall 2005 Research Paper                                                                  *
!*	 (Calculated by Mathematica)                                                                             *
!************************************************************************************************************/

  CCTK_REAL, parameter :: G_Kilogram = 4.036986068D35      !  units are kilogram

  CCTK_REAL, parameter :: G_Second = 1.0                   !  units are seconds

  CCTK_REAL, parameter :: G_Meter = 2.997924580D8          !  units are meters

  CCTK_REAL, parameter :: G_Kelvin = 2.627936856D75        !  units are Kelvin

  CCTK_REAL, parameter :: G_Density = 1.498287450D10       !  units are kilogram meter^-3
                                      
                                                           !  Specific Internal Energy
  CCTK_REAL, parameter :: G_SIE = 8.987551787D16           !  units are meter^2 second ^-2 Kelvin^-1

  CCTK_REAL, parameter :: G_Pressure = 1.346593606D27      !  units are kilogram meter^-1 seconds^-2

  CCTK_REAL, parameter :: G_Amp = 3.478876225D25           !  units are sqrt(kilogram meter^3 second^-4)

  CCTK_REAL, parameter :: G_Tesla = 1.160428200D10         !  units are kilogram Amp^-1 second^-2
  
  CCTK_REAL, parameter :: G_Gauss = 1.160428200D6          !  units are kilogram Amp^-1 second^-2

  CCTK_REAL, parameter :: G_E_Field = 3.478876222D18       !  units are kilogram meter Amp^-1 second^-3
  
      
! CCTK_REAL, parameter :: CONSTANTS                        !  currently not reserved for anything

  CCTK_REAL, parameter :: cc = 1                           !  speed of light is 299792458 m/s

  CCTK_REAL, parameter :: GC = 1                           !  Newton's Gravitational Constant 6.67e-11 kg m^3/s^2

  CCTK_REAL, parameter :: kB = 1                           !  Boltzmann Constant 1.38e-23 (kg m^2)/(s^2 K)
  
!/*********************************************************************************************************************
!*Coulomb's Force Constant 1/4*pi*epsilon_naught                                                                      *
!*	{epsilon_naught 8.85D-12 (s^4*A^2)/(kg*m^3)}                                                                      *
!*********************************************************************************************************************/

  CCTK_REAL, parameter :: cfC = 1/(4*pi)                   !  Coloumb Force Constant
  
end module constants_init

