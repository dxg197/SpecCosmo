! MHD_Init
! cactus: Cactus interfaces:
! Fortran 90 interface definition for routines in Cactus.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

module cactus
  implicit none
  public
  !
  interface
     !
     ! from Cactus:
     !
     subroutine cctk_barrier(cctkgh)
       implicit none
       CCTK_POINTER cctkgh
     end subroutine cctk_barrier
     !
     subroutine cctk_coordrange(ierr, cctkgh, lower, upper, dir, name, &
          &                      systemname)
       implicit none
       integer ierr
       CCTK_POINTER  cctkgh
       CCTK_REAL     lower
       CCTK_REAL     upper
       integer       dir
       character*(*) name
       character*(*) systemname
     end subroutine cctk_coordrange
     !
     function cctk_equals(arg1, arg2)
       implicit none
       integer       cctk_equals
       CCTK_POINTER  arg1
       character*(*) arg2
     end function cctk_equals
     !
     subroutine cctk_groupbboxgn(ierr, cctkgh, nbbox, bbox, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       nbbox
       integer       bbox(nbbox)
       character*(*) groupname
     end subroutine cctk_groupbboxgn
     !
     subroutine cctk_groupgshgn(ierr, cctkgh, ngsh, gsh, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       ngsh
       integer       gsh(ngsh)
       character*(*) groupname
     end subroutine cctk_groupgshgn
     !
     subroutine cctk_grouplbndgn(ierr, cctkgh, nlbnd, lbnd, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       nlbnd
       integer       lbnd(nlbnd)
       character*(*) groupname
     end subroutine cctk_grouplbndgn
     !
     subroutine cctk_grouplshgn(ierr, cctkgh, nlsh, lsh, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       nlsh
       integer       lsh(nlsh)
       character*(*) groupname
     end subroutine cctk_grouplshgn
     !
     subroutine cctk_groupnghostzonesgn &
          (ierr, cctkgh, nnghostzones, nghostzones, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       nnghostzones
       integer       nghostzones(nnghostzones)
       character*(*) groupname
     end subroutine cctk_groupnghostzonesgn
     !
     subroutine cctk_groupubndgn(ierr, cctkgh, nubnd, ubnd, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       nubnd
       integer       ubnd(nubnd)
       character*(*) groupname
     end subroutine cctk_groupubndgn
     !
     subroutine cctk_info(thorn, message)
       implicit none
       character*(*) thorn
       character*(*) message
     end subroutine cctk_info
     !
     function cctk_isthornactive(name)
       implicit none
       integer       cctk_isthornactive
       character*(*) name
     end function cctk_isthornactive
     !
     function cctk_myproc(cctkgh)
       implicit none
       integer      cctk_myproc
       CCTK_POINTER cctkgh
     end function cctk_myproc
     !
     function cctk_nprocs(cctkgh)
       implicit none
       integer      cctk_nprocs
       CCTK_POINTER cctkgh
     end function cctk_nprocs
     !
     subroutine cctk_paramwarn(thorn, message)
       implicit none
       character*(*) thorn
       character*(*) message
     end subroutine cctk_paramwarn
     !
     subroutine cctk_reductionarrayhandle(handle, reduction)
       implicit none
       integer       handle
       character*(*) reduction
     end subroutine cctk_reductionarrayhandle
     !
     subroutine cctk_registerbanner(ierr, banner)
       implicit none
       integer       ierr
       character*(*) banner
     end subroutine cctk_registerbanner
     !
     subroutine cctk_syncgroup(cctkgh, groupname)
       implicit none
       CCTK_POINTER  cctkgh
       character*(*) groupname
     end subroutine cctk_syncgroup
     !
     subroutine cctk_warn(level, line, file, thorn, message)
       implicit none
       integer       level
       integer       line
       character*(*) file
       character*(*) thorn
       character*(*) message
     end subroutine cctk_warn
     !
     !
     !
     ! from thorn CactusBase/Boundary:
     subroutine bndcopygn(ierr, cctkgh, sw, togroupname, fromgroupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       character*(*) togroupname
       character*(*) fromgroupname
     end subroutine bndcopygn
     !
     subroutine bndflatgn(ierr, cctkgh, sw, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       character*(*) groupname
     end subroutine bndflatgn
     !
     ! from thorn CactusBase/Boundary:
     subroutine bndradiativevn(ierr, cctkgh, sw, var0, v0, &
          &                     tovarname, fromvarname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       CCTK_REAL     var0
       CCTK_REAL     v0
       character*(*) tovarname
       character*(*) fromvarname
     end subroutine bndradiativevn
     !
     ! from thorn CactusBase/Boundary:
     subroutine bndrobinvn(ierr, cctkgh, sw, finf, npow, varname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       CCTK_REAL     finf
       integer       npow
       character*(*) varname
     end subroutine bndrobinvn
     !
     subroutine bndscalargn(ierr, cctkgh, sw, var0, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       CCTK_REAL     var0
       character*(*) groupname
     end subroutine bndscalargn
     !
     !
     !
     ! from thorn CactusBase/CartGrid3D:
     !
     subroutine cartsymgn(ierr, cctkgh, groupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       character*(*) groupname
     end subroutine cartsymgn
     !
     subroutine setcartsymvn(ierr, cctkgh, sw, varname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       character*(*) varname
     end subroutine setcartsymvn
     !
  end interface
  !
end module cactus
