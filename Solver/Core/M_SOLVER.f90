!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License. 
!
!   Contributors list:
!   - G. Delhommeau
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh 
!   - A. Babarit 
!   2022 02 10
!   R. Kurnia, added LU decomp and GMRES solver
!
!--------------------------------------------------------------------------------------
MODULE M_SOLVER

  IMPLICIT NONE

  PUBLIC :: GAUSSZ
 
  CONTAINS

  !---------------------------------------------------------------------------!

  SUBROUTINE GAUSSZ(A,N, Ainv)

    ! Input
    INTEGER,                  INTENT(IN)   :: N
    COMPLEX, DIMENSION(N, N), INTENT(IN)   :: A
    ! Output
    COMPLEX, DIMENSION(N, N), INTENT(OUT)  :: Ainv

    ! Local variables
    INTEGER :: I, J, K, L, IL
    COMPLEX :: C, P
    COMPLEX, DIMENSION(N,2*N) :: WS ! Working matrix
    REAL, PARAMETER :: EPS=1E-20

    ! Initialize working matrix
    WS(:, 1:N)     = A
    WS(:, N+1:2*N) = CMPLX(0., 0.)
    DO I = 1, N
      WS(I, N+I) = CMPLX(1., 0.)
    END DO


    ! Gauss pivot inversion
    DO J = 1, N-1
      K = J
      DO I = J+1, N
        IF ((ABS(WS(K, J))-ABS(WS(I, J))) <= 0.0) K = I
      END DO
      IF (K <= J) THEN 
        DO L = J, 2*N
          C = WS(J, L)
          WS(J, L) = WS(K, L)
          WS(K, L) = C
        END DO
      ELSE 
        IF (ABS(WS(J, J)) <= EPS) THEN
          WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
          STOP
        END IF
      END IF
      DO K = J+1, 2*N
        P = WS(J, K)/WS(J, J)
        DO I = J+1, N
          WS(I, K) = WS(I, K) - WS(I, J)*P
        END DO
      END DO
    END DO
    IF (ABS(WS(N, N)) < EPS) THEN 
      WRITE(*, '(A,E16.6)') 'PIVOT INFERIEUR A ', EPS
      STOP
    END IF
    DO IL = N+1, 2*N
      DO J = N, 1, -1
        WS(J, IL) = WS(J, IL)/WS(J, J)
        DO I = 1, J-1
          WS(I, IL) = WS(I, IL) - WS(I, J)*WS(J, IL)
        END DO
      END DO
    END DO
   
    ! Extract output
     Ainv(:, :) = WS(:, N+1:2*N)
    
     RETURN

  END SUBROUTINE

!   SUBROUTINE LU_INVERS_MATRIX(Ainv,M,N,ID_DP)
!    INTEGER:: M,N,I,J,INFO,K, ID_DP
!    COMPLEX, DIMENSION(M,N)  :: Ainv
!    INTEGER, DIMENSION(M)    :: IPIV
!    COMPLEX, DIMENSION(M)    :: WORK  ! work array for LAPACK
!    
!    ! ZGETRF is for complex double variable, use CGETRF for complex variable
!    ! Note that all of variables in this codes are defined as complex/ real (single precision)
!    ! but in the compile process,it is forced to be double precision with '-r8', see the makefile
!    ! if the -r8 is removed, then CGETRF and CGETRI must be used.
!    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
!    ! using partial pivoting with row interchanges.
!    
!    IF (ID_DP.EQ.1) THEN
!      CALL ZGETRF(M,N,Ainv,M,IPIV,INFO) !LAPACK function
!    ELSE
!      CALL CGETRF(M,N,Ainv,M,IPIV,INFO) !LAPACK function
!    END IF
!
!    IF (INFO /= 0) THEN
!     stop 'Matrix is numerically singular!'
!    END IF
!
!    ! ZGETRI is for complex double variable, use CGETRI for complex variable
!    ! ZGETRI computes the inverse of a matrix using the LU factorization
!    ! computed by ZGETRF.
!     
!      IF (ID_DP.EQ.1) THEN
!      call ZGETRI(M, Ainv, M, IPIV, WORK, M, info)
!      ELSE
!      call CGETRI(M, Ainv, M, IPIV, WORK, M, info)
!      END IF
!    IF (INFO /= 0) THEN
!         STOP 'Matrix inversion failed!'
!    END IF
!
!    RETURN
!       
!    END SUBROUTINE
!
!    SUBROUTINE LU_SOLVER(A,B,M,N,ID_DP)
!    INTEGER:: M,N,I,J,INFO,K, ID_DP
!    COMPLEX, DIMENSION(M,N)  :: A
!    INTEGER, DIMENSION(M)    :: IPIV
!    COMPLEX, DIMENSION(M)    :: B
!    
!    ! ZGETRF is for complex double variable, use CGETRF for complex variable
!    ! Note that all of variables in this codes are defined as complex/ real (single precision)
!    ! but in the compile process,it is forced to be double precision with '-r8', see the makefile
!    ! if the -r8 is removed, then CGETRF and CGETRI must be used.
!    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
!    ! using partial pivoting with row interchanges.
!    
!    IF (ID_DP.EQ.1) THEN
!      CALL ZGETRF(M,N,A,M,IPIV,INFO) !LAPACK function
!    ELSE
!      CALL CGETRF(M,N,A,M,IPIV,INFO) !LAPACK function
!    END IF
!
!    IF (INFO /= 0) THEN
!     stop 'Matrix is numerically singular!'
!    END IF
!
!    ! ZGETRS is for complex double variable, use CGETRI for complex variable
!    ! ZGETRS solves the lineary system using the LU factorization
!    ! computed by ZGETRF.
!     
!      IF (ID_DP.EQ.1) THEN
!      call ZGETRS('N',M, 1, A, M, IPIV, B, M, info)
!      ELSE
!      call CGETRS('N',M, 1, A, M, IPIV, B, M, info)
!      END IF
!    IF (INFO /= 0) THEN
!         STOP 'Solution procedure failed!'
!    END IF
!
!    RETURN
!       
!    END SUBROUTINE
!   
!
!    SUBROUTINE GMRES_SOLVER(A,B,nD,mD,ZOL_GMRES,ID_DP,TOLGMRES,NITERGMRES)
!!   For GMRES variables and parameters   
!    integer mD,nD  ! dimension for othonormal basis mD, and the matrix size nD x nD
!    integer i,j,lda, lwork, ldstrt, ID_DP,NITERGMRES
!    integer revcom, colx, coly, colz, nbscal
!    integer irc(5), icntl(8), info(3)
!    integer matvec, precondLeft, precondRight, dotProd
!    parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
!    integer nout
!    complex B(nD),A(nD,nD)
!    complex, dimension(:), allocatable :: work
!    real  cntl(5), rinfo(2),TOLGMRES
!    complex ZERO, ONE
!    complex ZOL_GMRES(nD)
!    parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))
!    
!    lda=nD
!    ldstrt = mD
!    !lwork = ldstrt**2 + ldstrt*(lda+5) + 5*lda + 1
!    lwork = ldstrt**2 + ldstrt*(lda+5) + 6*lda + 2      !if Icntl(5)=0 or 1 and Icntl(8)=1
!    allocate(work(lwork))
!    !write(*,*) lwork
!
!
!    !if (nD.gt.lda) then
!    !    STOP
!    !endif
!
!!***************************************
!! setting up the RHS
!    DO i=1,nD
!     work(nD+i)=B(i)
!    END DO
!    
!!*****************************************
!!** Reverse communication implementation
!!*
!
!!*******************************************************
!!** Initialize the control parameters to default value
!!*******************************************************
!!*
!IF (ID_DP.EQ.1) THEN
!      call init_zgmres(icntl,cntl)
!ELSE
!      call init_cgmres(icntl,cntl)
!ENDIF
!!*
!!*************************
!!*c Tune some parameters
!!*************************
!!*The tolerance
!      icntl(1)=TOLGMRES
!!* Save the convergence history on standard output
!      icntl(3) = 40
!!* Maximum number of iterations
!      icntl(7) = NITERGMRES
!!*
!!* preconditioner location
!      icntl(4) = 1
!!* orthogonalization scheme
!      icntl(5)=0
!!* initial guess
!      icntl(6) = 0
!!* residual calculation strategy at restart
!      icntl(8) = 1
!!****************************************
!!*
!
!IF (ID_DP.EQ.1) THEN
!10     call drive_zgmres(nD,nD,mD,lwork,work,irc,icntl,cntl,info,rinfo)
!               revcom = irc(1)
!               colx   = irc(2)
!               coly   = irc(3)
!               colz   = irc(4)
!               nbscal = irc(5)
!        !*
!        IF (revcom.eq.matvec) then
!        !* perform the matrix vector product
!        !*        work(colz) <-- A * work(colx)
!                call zgemv('N',nD,nD,ONE,A,lda,work(colx),1,ZERO,work(colz),1)
!                goto 10
!        !*
!        ELSE IF (revcom.eq.precondLeft) then
!        !* perform the left preconditioning
!        !*         work(colz) <-- M^{-1} * work(colx)
!        
!                 call zcopy(nD,work(colx),1,work(colz),1)
!                 call ztrsm('L','L','N','N',nD,1,ONE,A,lda,work(colz),nD)
!                 goto 10
!        !*
!        ELSE IF (revcom.eq.precondRight) then
!        !* perform the right preconditioning
!                 
!                 call zcopy(nD,work(colx),1,work(colz),1)
!                 call ztrsm('L','U','N','N',nD,1,ONE,A,lda,work(colz),nD)
!                 goto 10
!        !*
!        ELSE IF (revcom.eq.dotProd) then
!        !*      perform the scalar product
!        !*      work(colz) <-- work(colx) work(coly)
!        !*
!        
!                 call zgemv('C',nD,nbscal,ONE,work(colx),nD,work(coly),1,ZERO,work(colz),1)
!                 goto 10
!        END IF
!ELSE
!12     call drive_cgmres(nD,nD,mD,lwork,work,irc,icntl,cntl,info,rinfo)
!               revcom = irc(1)
!               colx   = irc(2)
!               coly   = irc(3)
!               colz   = irc(4)
!               nbscal = irc(5)
!        !*
!        IF (revcom.eq.matvec) then
!        !* perform the matrix vector product
!        !*        work(colz) <-- A * work(colx)
!                call cgemv('N',nD,nD,ONE,A,lda,work(colx),1,ZERO,work(colz),1)
!                goto 12
!        !*
!        ELSE IF (revcom.eq.precondLeft) then
!        !* perform the left preconditioning
!        !*         work(colz) <-- M^{-1} * work(colx)
!        
!                 call ccopy(nD,work(colx),1,work(colz),1)
!                 call ctrsm('L','L','N','N',nD,1,ONE,A,lda,work(colz),nD)
!                 goto 12
!        !*
!        ELSE IF (revcom.eq.precondRight) then
!        !* perform the right preconditioning
!                 
!                 call ccopy(nD,work(colx),1,work(colz),1)
!                 call ctrsm('L','U','N','N',nD,1,ONE,A,lda,work(colz),nD)
!                 goto 12
!        !*
!        ELSE IF (revcom.eq.dotProd) then
!        !*      perform the scalar product
!        !*      work(colz) <-- work(colx) work(coly)
!        !*
!        
!                 call cgemv('C',nD,nbscal,ONE,work(colx),nD,work(coly),1,ZERO,work(colz),1)
!                 goto 12
!        END IF
!END IF
!!*
!        if (info(1).eq.0) then
!        write(*,*) ' Normal exit'
!        write(*,*) ' Convergence after ', info(2),' iterations'
!       ! write(*,*) ' Backward error - preconditioned system', rinfo(1)
!       ! write(*,*) ' Backward error - unpreconditioned system', rinfo(2)
!       ! write(*,*) ' Solution : '
!        do j=1,nD
!       ! write(*,*) work(j)
!        ZOL_GMRES(j)=work(j)
!        end do
!       ! write(*,*) ' Optimal size for workspace ', info(3)
!        else if (info(1).eq.-1) then
!        write(*,*) ' Bad value of n'
!        else if (info(1).eq.-2) then
!        write(*,*) ' Bad value of m'
!        else if (info(1).eq.-3) then
!        write(*,*) ' Too small workspace. '
!        write(*,*) ' Minimal value should be ', info(2)
!        else if (info(1).eq.-4) then
!        write(*,*) ' No convergence after ', icntl(7), ' iterations'
!        write(*,*) ' switched to LU decomposition solver'
!        WRITE(100,*) 'GMRES do not converge, switched to LU decomposition solver!'
!        CALL LU_SOLVER(A,B,nD,nD,ID_DP)
!        do j=1,nD
!        ZOL_GMRES(j)=B(j)
!        end do
!        else if (info(1).eq.-5) then
!        write(*,*) ' Type of preconditioner not specified'
!        endif
!!*******************************
!       deallocate(work)
!    END SUBROUTINE
END MODULE 
