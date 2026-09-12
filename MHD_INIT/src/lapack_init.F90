#include "cctk.h"

module lapack_init
  implicit none
  public
  
  ! ev: eigenvalues
  ! ge (general), gg (generalised eigenvalues), sb (symmetric banded),
  ! sp (symmetric packed), st (symmetric tridiagonal), sy (symmetric)
  
  interface geev
     subroutine sgeev (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, &
          &            ldvr, work, lwork, info)
       character jobvl, jobvr
       integer   info, lda, ldvl, ldvr, lwork, n
       real      a(lda,n), vl(ldvl,n), vr(ldvr,n), &
            &    wi(n), work(lwork), wr(n)
     end subroutine sgeev
     subroutine dgeev (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, &
          &            ldvr, work, lwork, info)
       character        jobvl, jobvr
       integer          info, lda, ldvl, ldvr, lwork, n
       double precision a(lda,n), vl(ldvl,n), vr(ldvr,n), &
            &           wi(n), work(lwork), wr(n)
     end subroutine dgeev
  end interface
  
  interface syev
     subroutine ssyev (jobz, uplo, n, a, lda, w, work, lwork, info)
       character jobz, uplo
       integer   info, lda, lwork, n
       real      a(lda,n), w(n), work(lwork)
     end subroutine ssyev
     subroutine dsyev (jobz, uplo, n, a, lda, w, work, lwork, info)
       character        jobz, uplo
       integer          info, lda, lwork, n
       double precision a(lda,n), w(n), work(lwork)
     end subroutine dsyev
  end interface
  
  
  
  ! sv: solve
  ! gb (general banded), ge (general), gt (general tridiagonal),
  ! pb (symmetric positive banded), po (symmetric positive),
  ! pp (symmetric positive packed), pt (symmetric positive tridiagonal),
  ! sp (symmetric packed), sy (symmetric)
  
  interface gesv
     subroutine sgesv (n, nrhs, a, lda, ipiv, b, ldb, info)
       implicit none
       integer n
       integer nrhs
       integer lda
       real    a(lda,n)
       integer ipiv(n)
       integer ldb
       real    b(ldb,nrhs)
       integer info
     end subroutine sgesv
     subroutine dgesv (n, nrhs, a, lda, ipiv, b, ldb, info)
       implicit none
       integer          n
       integer          nrhs
       integer          lda
       double precision a(lda,n)
       integer          ipiv(n)
       integer          ldb
       double precision b(ldb,nrhs)
       integer          info
     end subroutine dgesv
  end interface

  interface posv
     subroutine sposv (uplo, n, nrhs, a, lda, b, ldb, info)
       implicit none
       character uplo
       integer   n
       integer   nrhs
       integer   lda
       real      a(lda,n)
       integer   ldb
       real      b(ldb,nrhs)
       integer   info
     end subroutine sposv
     subroutine dposv (uplo, n, nrhs, a, lda, b, ldb, info)
       implicit none
       character        uplo
       integer          n
       integer          nrhs
       integer          lda
       double precision a(lda,n)
       integer          ldb
       double precision b(ldb,nrhs)
       integer          info
     end subroutine dposv
  end interface

  interface sysv
     subroutine ssysv (uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
       implicit none
       character uplo
       integer   n
       integer   nrhs
       integer   lda
       real      a(lda,n)
       integer   ipiv(n)
       integer   ldb
       real      b(ldb,nrhs)
       integer   lwork
       real      work(lwork)
       integer   info
     end subroutine ssysv
     subroutine dsysv (uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
       implicit none
       character        uplo
       integer          n
       integer          nrhs
       integer          lda
       double precision a(lda,n)
       integer          ipiv(n)
       integer          ldb
       double precision b(ldb,nrhs)
       integer          lwork
       double precision work(lwork)
       integer          info
     end subroutine dsysv
  end interface
  
  
  
  ! trf: LU factorisation
  ! gb (general banded), ge (general), gt (general tridiagonal),
  ! pb (symmetric positive banded), po (symmetric positive),
  ! pp (symmetric positive packed), pt (symmetric positive tridiagonal),
  ! sp (symmetric packed), sy (symmetric)
  
  interface sytrf
     subroutine ssytrf (m, n, a, lda, ipiv, info)
       implicit none
       integer m
       integer n
       integer lda
       real    a(lda,n)
       integer ipiv(n)
       integer info
     end subroutine ssytrf
     subroutine dsytrf (m, n, a, lda, ipiv, info)
       implicit none
       integer          m
       integer          n
       integer          lda
       double precision a(lda,n)
       integer          ipiv(n)
       integer          info
     end subroutine dsytrf
  end interface
  
end module lapack_init
