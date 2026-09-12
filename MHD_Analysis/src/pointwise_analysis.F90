#include "cctk.h"

module pointwise_analysis
  implicit none
  private
  public calc_position
  public calc_offsets
  public get_scalar
  public get_vector
  public get_tensor
  public get_connections
  public set_scalar
  public set_vector
  public set_tensor
  public set_connections
contains
  subroutine calc_position(shape, i,j,k, pos)
    integer, intent(in)  :: shape(3)
    integer, intent(in)  :: i,j,k
    integer, intent(out) :: pos
    pos = 1 + i-1 + shape(1) * (j-1 + shape(2) * (k-1))
  end subroutine calc_position
  subroutine calc_offsets (shape, off)
    integer, intent(in)  :: shape(3)
    integer, intent(out) :: off(3)
    off(1) = 1
    off(2) = shape(1)
    off(3) = shape(1) * shape(2)
  end subroutine calc_offsets
  subroutine get_scalar(a, f, pos)
    CCTK_REAL, intent(in)  :: a(*)
    CCTK_REAL, intent(out) :: f
    integer,   intent(in)  :: pos
    f = a(pos)
  end subroutine get_scalar
  subroutine get_vector(ax,ay,az, f, pos)
    CCTK_REAL, intent(in)  :: ax(*),ay(*),az(*)
    CCTK_REAL, intent(out) :: f(3)
    integer,   intent(in)  :: pos
    f(1) = ax(pos)
    f(2) = ay(pos)
    f(3) = az(pos)
  end subroutine get_vector
  subroutine get_tensor(axx,axy,axz,ayy,ayz,azz, f, pos)
    CCTK_REAL, intent(in)  :: axx(*),axy(*),axz(*),ayy(*),ayz(*),azz(*)
    CCTK_REAL, intent(out) :: f(3,3)
    integer,   intent(in)  :: pos
    f(1,1) = axx(pos)
    f(1,2) = axy(pos)
    f(1,3) = axz(pos)
    f(2,1) = axy(pos)
    f(2,2) = ayy(pos)
    f(2,3) = ayz(pos)
    f(3,1) = axz(pos)
    f(3,2) = ayz(pos)
    f(3,3) = azz(pos)
  end subroutine get_tensor
  subroutine get_connections(&
       gammaxxx,gammaxxy,gammaxxz,gammaxyy,gammaxyz,gammaxzz, &
       gammayxx,gammayxy,gammayxz,gammayyy,gammayyz,gammayzz, &
       gammazxx,gammazxy,gammazxz,gammazyy,gammazyz,gammazzz, gamma, pos)
    CCTK_REAL, intent(in)  :: gammaxxx(*),gammaxxy(*),gammaxxz(*)
    CCTK_REAL, intent(in)  :: gammaxyy(*),gammaxyz(*),gammaxzz(*)
    CCTK_REAL, intent(in)  :: gammayxx(*),gammayxy(*),gammayxz(*)
    CCTK_REAL, intent(in)  :: gammayyy(*),gammayyz(*),gammayzz(*)
    CCTK_REAL, intent(in)  :: gammazxx(*),gammazxy(*),gammazxz(*)
    CCTK_REAL, intent(in)  :: gammazyy(*),gammazyz(*),gammazzz(*)
    CCTK_REAL, intent(out) :: gamma(3,3,3)
    integer,   intent(in)  :: pos
    gamma(1,1,1) = gammaxxx(pos)
    gamma(1,1,2) = gammaxxy(pos)
    gamma(1,1,3) = gammaxxz(pos)
    gamma(1,2,1) = gammaxxy(pos)
    gamma(1,2,2) = gammaxyy(pos)
    gamma(1,2,3) = gammaxyz(pos)
    gamma(1,3,1) = gammaxxz(pos)
    gamma(1,3,2) = gammaxyz(pos)
    gamma(1,3,3) = gammaxzz(pos)
    gamma(2,1,1) = gammayxx(pos)
    gamma(2,1,2) = gammayxy(pos)
    gamma(2,1,3) = gammayxz(pos)
    gamma(2,2,1) = gammayxy(pos)
    gamma(2,2,2) = gammayyy(pos)
    gamma(2,2,3) = gammayyz(pos)
    gamma(2,3,1) = gammayxz(pos)
    gamma(2,3,2) = gammayyz(pos)
    gamma(2,3,3) = gammayzz(pos)
    gamma(3,1,1) = gammazxx(pos)
    gamma(3,1,2) = gammazxy(pos)
    gamma(3,1,3) = gammazxz(pos)
    gamma(3,2,1) = gammazxy(pos)
    gamma(3,2,2) = gammazyy(pos)
    gamma(3,2,3) = gammazyz(pos)
    gamma(3,3,1) = gammazxz(pos)
    gamma(3,3,2) = gammazyz(pos)
    gamma(3,3,3) = gammazzz(pos)
  end subroutine get_connections
  subroutine set_scalar(f, a, pos)
    CCTK_REAL, intent(in) :: f
    CCTK_REAL             :: a(*)
    integer,   intent(in) :: pos
    a(pos) = f
  end subroutine set_scalar
  subroutine set_vector(f, ax,ay,az, pos)
    CCTK_REAL, intent(in) :: f(3)
    CCTK_REAL             :: ax(*),ay(*),az(*)
    integer,   intent(in) :: pos
    ax(pos) = f(1)
    ay(pos) = f(2)
    az(pos) = f(3)
  end subroutine set_vector
  subroutine set_tensor(f, axx,axy,axz,ayy,ayz,azz, pos)
    CCTK_REAL, intent(in) :: f(3,3)
    CCTK_REAL             :: axx(*),axy(*),axz(*),ayy(*),ayz(*),azz(*)
    integer,   intent(in) :: pos
    axx(pos) = f(1,1)
    axy(pos) = f(1,2)
    axz(pos) = f(1,3)
    ayy(pos) = f(2,2)
    ayz(pos) = f(2,3)
    azz(pos) = f(3,3)
  end subroutine set_tensor
  subroutine set_connections(gamma, &
       gammaxxx,gammaxxy,gammaxxz,gammaxyy,gammaxyz,gammaxzz, &
       gammayxx,gammayxy,gammayxz,gammayyy,gammayyz,gammayzz, &
       gammazxx,gammazxy,gammazxz,gammazyy,gammazyz,gammazzz, pos)
    CCTK_REAL, intent(in) :: gamma(3,3,3)
    CCTK_REAL             :: gammaxxx(*),gammaxxy(*),gammaxxz(*)
    CCTK_REAL             :: gammaxyy(*),gammaxyz(*),gammaxzz(*)
    CCTK_REAL             :: gammayxx(*),gammayxy(*),gammayxz(*)
    CCTK_REAL             :: gammayyy(*),gammayyz(*),gammayzz(*)
    CCTK_REAL             :: gammazxx(*),gammazxy(*),gammazxz(*)
    CCTK_REAL             :: gammazyy(*),gammazyz(*),gammazzz(*)
    integer,   intent(in)  :: pos
    gammaxxx(pos) = gamma(1,1,1)
    gammaxxy(pos) = gamma(1,1,2)
    gammaxxz(pos) = gamma(1,1,3)
    gammaxyy(pos) = gamma(1,2,2)
    gammaxyz(pos) = gamma(1,2,3)
    gammaxzz(pos) = gamma(1,3,3)
    gammayxx(pos) = gamma(2,1,1)
    gammayxy(pos) = gamma(2,1,2)
    gammayxz(pos) = gamma(2,1,3)
    gammayyy(pos) = gamma(2,2,2)
    gammayyz(pos) = gamma(2,2,3)
    gammayzz(pos) = gamma(2,3,3)
    gammazxx(pos) = gamma(3,1,1)
    gammazxy(pos) = gamma(3,1,2)
    gammazxz(pos) = gamma(3,1,3)
    gammazyy(pos) = gamma(3,2,2)
    gammazyz(pos) = gamma(3,2,3)
    gammazzz(pos) = gamma(3,3,3)
  end subroutine set_connections
end module pointwise_analysis
