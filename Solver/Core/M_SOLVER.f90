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
!
!--------------------------------------------------------------------------------------
MODULE M_SOLVER

  IMPLICIT NONE
 
  CONTAINS
!---------------------------------------------------------------------------!
    SUBROUTINE GAUSSZ(A,NMAX,N,M)
    INTEGER:: N,M,NMAX,I,J,K,L,IL
    !COMPLEX A(NMAX,1),C,P   
    COMPLEX A(NMAX,2*NMAX),C,P    !RK: size matrix A is changed to be consistent
    REAL:: EPS
 
    EPS=1.E-20
	DO J=1,N-1
	    K=J
	    DO I=J+1,N
	        IF  ((ABS(A(K,J))-ABS(A(I,J))).LE.0.) K=I
        END DO
	    IF ((K-J).LE.0.) THEN 
            DO L=J,M
	            C=A(J,L)
	            A(J,L)=A(K,L)
                A(K,L)=C
            END DO
        ELSE 
            IF ((ABS(A(J,J))-EPS).LE.0.) THEN
                WRITE(*,'(A,E16.6)') 'PIVOT INFERIEUR A ',EPS
                STOP
            END IF
        END IF
        DO K=J+1,M
	        P=A(J,K)/A(J,J)
	        DO I=J+1,N
                A(I,K)=A(I,K)-A(I,J)*P
            END DO
        END DO
    END DO
	IF ((ABS(A(N,N))-EPS).LE.0.) THEN 
	    WRITE(*,'(A,E16.6)') 'PIVOT INFERIEUR A ',EPS
        STOP
    END IF
    DO IL=N+1,M
	    DO J=N,1,-1
	        A(J,IL)=A(J,IL)/A(J,J)
	        DO I=1,J-1
                A(I,IL)=A(I,IL)-A(I,J)*A(J,IL)
            END DO
        END DO
    END DO
	RETURN
        
!	EPS=1.E-20
!	DO 10 J=1,N-1
!	K=J
!	DO 20 I=J+1,N
!	IF(ABS(A(K,J))-ABS(A(I,J)))30,20,20
!    30 K=I
!    20 CONTINUE
!	IF(K-J)50,40,50
!    50 DO 60 L=J,M
!	C=A(J,L)
!	A(J,L)=A(K,L)
!    60 A(K,L)=C
!    40 IF(ABS(A(J,J))-EPS)120,120,70
!    70 DO 80 K=J+1,M
!	P=A(J,K)/A(J,J)
!	DO 80 I=J+1,N
!    80 A(I,K)=A(I,K)-A(I,J)*P
!    10 CONTINUE
!	IF(ABS(A(N,N))-EPS)120,120,90
!    90 DO 100 IL=N+1,M
!	DO 100 J=N,1,-1
!	A(J,IL)=A(J,IL)/A(J,J)
!	DO 100 I=1,J-1
!    100 A(I,IL)=A(I,IL)-A(I,J)*A(J,IL)
!	RETURN
!    120 WRITE(*,500)EPS
!    500 FORMAT(5X,'PIVOT INFERIEUR A ',1P,E16.6)
!	STOP
	
    END SUBROUTINE

   SUBROUTINE LU_INVERS_MATRIX(Ainv,M,N)
    INTEGER:: M,N,I,J,INFO,K
    COMPLEX, DIMENSION(M,N)  :: Ainv
    INTEGER, DIMENSION(M)    :: IPIV
    COMPLEX, DIMENSION(M)    :: WORK  ! work array for LAPACK
    
    ! ZGETRF is for complex double variable, use CGETRF for complex variable
    ! Note that all of variables in this codes are defined as complex/ real (single precision)
    ! but in the compile process,it is forced to be double precision with '-r8', see the makefile
    ! if the -r8 is removed, then CGETRF and CGETRI must be used.
    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL ZGETRF(M,N,Ainv,M,IPIV,INFO) !LAPACK function
    IF (INFO /= 0) THEN
     stop 'Matrix is numerically singular!'
    END IF

    ! ZGETRI is for complex double variable, use CGETRI for complex variable
    ! ZGETRI computes the inverse of a matrix using the LU factorization
    ! computed by ZGETRF.
     call ZGETRI(M, Ainv, M, IPIV, WORK, M, info)

    IF (INFO /= 0) THEN
         STOP 'Matrix inversion failed!'
    END IF

    RETURN
       
    END SUBROUTINE

END MODULE 
