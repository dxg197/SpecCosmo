! MHD_Analysis
! coords: Coordinate calculations:
! Calculate the current coordinates, depending on the time, grid position,
! and shift offsets.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module coords_analysis
  implicit none
  private
  public getoffset, getcoords, getscalar, getvector, gettensor, setscalar, setvector, settensor
contains
  !
  subroutine getoffset(x,y,z, offset, i,j,k)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)  :: x(:,:,:), y(:,:,:), z(:,:,:)
    CCTK_REAL, intent(out) :: offset(3)
    integer,   intent(in)  :: i, j, k
    !
    CCTK_REAL xx(3), rr
    !
    call getvector(x, y, z, xx, i,j,k)
    if(any(xx/=0)) then
       rr = sqrt(sum(xx**2))
       offset =(/ shift_vx, shift_vy, shift_vz /) + shift_vr * xx / rr
    else
       offset = 0
    end if
  end subroutine getoffset
  !
  subroutine getcoords(time, x,y,z, xx, i,j,k)
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(in)  :: time
    CCTK_REAL, intent(in)  :: x(:,:,:), y(:,:,:), z(:,:,:)
    CCTK_REAL, intent(out) :: xx(3)
    integer,   intent(in)  :: i, j, k
    !
    CCTK_REAL offset(3)
    CCTK_REAL rr
    !
    call getvector(x, y, z, xx, i,j,k)
    offset =(/ shift_vx, shift_vy, shift_vz /)
    if(any(xx/=0)) then
       rr = sqrt(sum(xx**2))
       offset = offset + shift_vr * xx / rr
    end if
    xx = xx + time * offset
  end subroutine getcoords
  !
  subroutine getscalar(arr, s, i,j,k)
    CCTK_REAL, intent(in)  :: arr(:,:,:)
    CCTK_REAL, intent(out) :: s
    integer,   intent(in)  :: i,j,k
    !
    ! get scalar from array
    s = arr(i,j,k)
  end subroutine getscalar
  !
  subroutine getvector(arrx, arry, arrz, v, i,j,k)
    CCTK_REAL, intent(in)  :: arrx(:,:,:), arry(:,:,:), arrz(:,:,:)
    CCTK_REAL, intent(out) :: v(3)
    integer,   intent(in)  :: i,j,k
    !
    ! get scalars from arrays
    v(1) = arrx(i,j,k)
    v(2) = arry(i,j,k)
    v(3) = arrz(i,j,k)
  end subroutine getvector
  !
  subroutine gettensor(arrxx, arrxy, arrxz, arryy, arryz, arrzz, t, i,j,k)
    CCTK_REAL, intent(in)  :: arrxx(:,:,:), arrxy(:,:,:), arrxz(:,:,:), &
         &                    arryy(:,:,:), arryz(:,:,:), arrzz(:,:,:)
    CCTK_REAL, intent(out) :: t(3,3)
    integer,   intent(in)  :: i,j,k
    !
    ! get scalars from arrays
    t(1,1) = arrxx(i,j,k)
    t(1,2) = arrxy(i,j,k)
    t(1,3) = arrxz(i,j,k)
    t(2,2) = arryy(i,j,k)
    t(2,3) = arryz(i,j,k)
    t(3,3) = arrzz(i,j,k)
    !
    ! symmetries
    t(2,1) = t(1,2)
    t(3,1) = t(1,3)
    t(3,2) = t(2,3)
  end subroutine gettensor
  !
  subroutine setscalar(s, arr, i,j,k)
    CCTK_REAL, intent(in)    :: s
    CCTK_REAL, intent(inout) :: arr(:,:,:)
    integer,   intent(in)    :: i,j,k
    !
    ! set array from scalar
    arr(i,j,k) = s
  end subroutine setscalar
  !
  subroutine setvector(v, arrx, arry, arrz, i,j,k)
    CCTK_REAL, intent(in)    :: v(3)
    CCTK_REAL, intent(inout) :: arrx(:,:,:), arry(:,:,:), arrz(:,:,:)
    integer,   intent(in)    :: i,j,k
    !
    ! set arrays from scalars
    arrx(i,j,k) = v(1)
    arry(i,j,k) = v(2)
    arrz(i,j,k) = v(3)
  end subroutine setvector
  !
  subroutine settensor(t, arrxx, arrxy, arrxz, arryy, arryz, arrzz, i,j,k)
    CCTK_REAL, intent(in)    :: t(3,3)
    CCTK_REAL, intent(inout) :: arrxx(:,:,:), arrxy(:,:,:), arrxz(:,:,:), &
         &                      arryy(:,:,:), arryz(:,:,:), arrzz(:,:,:)
    integer,   intent(in)    :: i,j,k
    !
    ! set arrays from scalars
    arrxx(i,j,k) = t(1,1)
    arrxy(i,j,k) = t(1,2)
    arrxz(i,j,k) = t(1,3)
    arryy(i,j,k) = t(2,2)
    arryz(i,j,k) = t(2,3)
    arrzz(i,j,k) = t(3,3)
  end subroutine settensor
  !
end module coords_analysis
