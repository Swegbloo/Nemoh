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
!
!--------------------------------------------------------------------------------------
MODULE Elementary_functions

  IMPLICIT NONE

  PUBLIC :: GG, CIH, SIH, CROSS_PRODUCT,PL2,PL5

CONTAINS

  COMPLEX FUNCTION GG(Z, CEX)
    ! Estimation of ∫_z^∞ exp(-t)/t dt
    ! See p.367 of G. Delhommeau thesis (referenced as [Del]).

    COMPLEX, INTENT(IN) :: Z, CEX
    COMPLEX             :: Y

    IF (REAL(Z) < -16.0) THEN                                      ! Case 1 p. 368 in [Del]
      Y = 1./Z
      GG = Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))
    ELSE IF (ABS(AIMAG(Z)) > 10.0) THEN                            ! Case 3 p. 368 in [Del]
      GG = 0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE IF (REAL(Z) > -0.5) THEN                                  ! Case 2 p. 368 in [Del]
      GG = -(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+    &
        Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+                    &
        Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+                  &
        Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+                     &
        Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01 &
        +Z*(-0.8700861E-03+Z*0.2989204E-03))))))
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE                                                           ! Case 4 p. 369 in [Del]
      IF (AIMAG(Z) < 0) THEN
        GG = ((((((( (1.000000, 1.3935496E-06)*Z+ (15.82958, -20.14222))  &
          *Z+ (-70.52863, -227.9511))*Z+ (-985.4221, -226.6272))*Z        &
          + (-1202.318, 1580.907))*Z+ (953.2441, 1342.447))*Z             &
          + (417.3716, -196.6665))*Z+ (-9.881266, -24.24952))/            &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, -20.14481))*Z  &
          + (-55.66969, -248.1167))*Z+ (-1068.640, -434.4707))*Z          &
          + (-2082.250, 1522.471))*Z+ (383.3455, 2730.378))*Z             &
          + (1216.791, 351.7189))*Z+ (115.3926, -161.2647))*Z             &
          + (-3.777369, -4.510900))-(0., 3.14159265)*CEX
      ELSE
        GG = ((((((( (1.000000, -1.3935496E-06)*Z+ (15.82958, 20.14222))  &
          *Z+ (-70.52863, 227.9511))*Z+ (-985.4221, 226.6272))*Z          &
          + (-1202.318, -1580.907))*Z+ (953.2441, -1342.447))*Z           &
          + (417.3716, 196.6665))*Z+ (-9.881266, 24.24952))/              &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, 20.14481))*Z   &
          + (-55.66969, 248.1167))*Z+ (-1068.640, 434.4707))*Z            &
          + (-2082.250, -1522.471))*Z+ (383.3455, -2730.378))*Z           &
          + (1216.791, -351.7189))*Z+ (115.3926, 161.2647))*Z             &
          + (-3.777369, 4.510900))+(0., 3.14159265)*CEX
      END IF
    END IF
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION CIH(AK,Z,H)
    REAL, INTENT(IN) :: AK,Z,H

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      CIH=COSH(AK*(Z+H))/COSH(AK*H)
    ELSE
      CIH=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION SIH(AK,Z,H)
    REAL, INTENT(IN) :: AK,Z,H

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      SIH=SINH(AK*(Z+H))/COSH(AK*H)
    ELSE
      SIH=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  FUNCTION CIH_Vect(AK,Z,H,NZ) result(CIHV)
    INTEGER           , INTENT(IN) :: NZ
    REAL              , INTENT(IN) :: AK,H
    REAL,DIMENSION(NZ), INTENT(IN) :: Z  
    REAL,DIMENSION(NZ)             :: CIHV

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      CIHV=COSH(AK*(Z+H))/COSH(AK*H)
    ELSE
      CIHV=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

   FUNCTION SIH_Vect(AK,Z,H,NZ) result(SIHV)
    INTEGER           , INTENT(IN) :: NZ
    REAL              , INTENT(IN) :: AK,H
    REAL,DIMENSION(NZ), INTENT(IN) :: Z  
    REAL,DIMENSION(NZ)             :: SIHV

    IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
      SIHV=SINH(AK*(Z+H))/COSH(AK*H)
    ELSE
      SIHV=EXP(AK*Z)
    ENDIF
    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  FUNCTION CROSS_PRODUCT(A, B)
    REAL, DIMENSION(3) :: CROSS_PRODUCT
    REAL, DIMENSION(3), INTENT(IN) :: A, B

    CROSS_PRODUCT(1) = A(2)*B(3) - A(3)*B(2)
    CROSS_PRODUCT(2) = A(3)*B(1) - A(1)*B(3)
    CROSS_PRODUCT(3) = A(1)*B(2) - A(2)*B(1)
  END FUNCTION CROSS_PRODUCT

  FUNCTION CROSS_PRODUCT_COMPLEX(A, B) RESULT(CROSS_PRODUCT)
    COMPLEX, DIMENSION(3) :: CROSS_PRODUCT
    COMPLEX, DIMENSION(3), INTENT(IN) :: A, B

    CROSS_PRODUCT(1) = A(2)*B(3) - A(3)*B(2)
    CROSS_PRODUCT(2) = A(3)*B(1) - A(1)*B(3)
    CROSS_PRODUCT(3) = A(1)*B(2) - A(2)*B(1)
  END FUNCTION CROSS_PRODUCT_COMPLEX

  !-------------------------------------------------------------------------------!

  REAL FUNCTION X0(Y)
    ! Solve x * tanh(x) = y by dichotomy

    REAL, INTENT(IN) :: Y
    REAL             :: X_LOW, X_MEAN, X_UP
    REAL             :: EPS, STEP

    X_LOW = 0.0

    ! Find upper bound for x
    X_UP = 0.0
    STEP = MAX(Y, SQRT(Y))
    DO WHILE (X_UP*TANH(X_UP) < Y)
      X_UP = X_UP + STEP
    END DO

    ! Dichotomy
    EPS = 5.E-6
    DO WHILE (X_UP - X_LOW > EPS*X_UP)
      X_MEAN = (X_LOW + X_UP)/2
      IF (X_MEAN*TANH(X_MEAN) < Y) THEN
        X_LOW = X_MEAN
      ELSE
        X_UP = X_MEAN
      END IF
    END DO
  
    X0 = X_MEAN

    RETURN
  END FUNCTION

  !-------------------------------------------------------------------------------!

  REAL FUNCTION PL2(U1,U2,U3,XU)
    REAL::U1,U2,U3,XU
    PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
    RETURN
  END FUNCTION
  
  REAL FUNCTION PL5(U1,U2,U3,U4,U5,XU)
    REAL::U1,U2,U3,U4,U5,XU
    PL5=((XU-U1)*(XU-U2)*(XU-U3)*(XU-U4))/&
      ((U5-U1)*(U5-U2)*(U5-U3)*(U5-U4))
    RETURN
  END FUNCTION

  FUNCTION Fun_closest(N,w,wobs) result(Iwobs)
    INTEGER, INTENT(IN):: N
    REAL, INTENT(IN)   :: w(N),wobs
    INTEGER            :: Iwobs,Iw,Nw
    REAL               :: mindist,mindistN
    
    mindist=w(N)-w(1)
    DO Iw=1,N
        IF (abs(w(Iw)-wobs)<=mindist)THEN
                Iwobs=Iw
                mindist=abs(w(Iw)-wobs)
        ENDIF
    ENDDO
    RETURN     
  END FUNCTION

  FUNCTION Fun_MIN(N,vect) result(minvalue)
    INTEGER,            INTENT(IN):: N
    REAL,DIMENSION(N),  INTENT(IN):: vect
    REAL                          :: minvalue
    INTEGER                       ::I
    minvalue=vect(1)
    DO I=1,N
    minvalue=MIN(vect(I),minvalue)
    ENDDO    
  END FUNCTION

  FUNCTION Fun_MAX(N,vect) result(maxvalue)
    INTEGER,            INTENT(IN):: N
    REAL,DIMENSION(N),  INTENT(IN):: vect
    REAL                          :: maxvalue
    INTEGER                       :: I
    maxvalue=vect(1)
    DO I=1,N
    maxvalue=MAX(vect(I),maxvalue)
    ENDDO    
  END FUNCTION

  FUNCTION Fun_KronDelta(m,n) result(delta)
        INTEGER, INTENT(IN) :: m,n
        INTEGER             :: delta
        IF (m.EQ.n) delta=1
        IF (m.NE.n) delta=0
  END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Bessel functions added by RK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION fun_BESSJ (N,X) RESULT(BESSJ)

!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      IMPLICIT NONE
      INTEGER, PARAMETER :: IACC = 40
      REAL, PARAMETER :: BIGNO = 1.D10, BIGNI = 1.D-10
      INTEGER M, N, J, JSUM
      REAL  X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N.EQ.0) THEN
      BESSJ = fun_BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = fun_BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = fun_BESSJ0(X)
      BJ  = fun_BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END FUNCTION

      FUNCTION fun_BESSJ0 (X) RESULT(BESSJ0)
      IMPLICIT NONE
      REAL X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL  Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END FUNCTION
! ---------------------------------------------------------------------------
      FUNCTION fun_BESSJ1 (X) RESULT(BESSJ1)
      IMPLICIT NONE
      REAL X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL  Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
      
      END FUNCTION

      

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Elementary_functions
