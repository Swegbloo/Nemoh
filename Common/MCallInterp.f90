MODULE MCallInterp

USE linear_interpolation_module

IMPLICIT NONE 

CONTAINS

 FUNCTION FUN_INTERP1_REAL(X,VAR,NX,XOUT,NXOUT) RESULT(VAROUT)
   INTEGER,            INTENT(IN) :: NX,NXOUT
   REAL,DIMENSION(NX), INTENT(IN) :: X,VAR
   REAL,DIMENSION(NXOUT)          :: XOUT,VAROUT
   Type(linear_interp_1d)         :: interp1
   INTEGER                        :: iflag,I
   
   CALL interp1%initialize(X,VAR,iflag)
   DO I=1,NXOUT
     CALL interp1%evaluate(XOUT(I), VAROUT(I))
   ENDDO
   CALL interp1%destroy()
 END FUNCTION

 FUNCTION FUN_INTERP1_COMPLEX(X,VAR,NX,XOUT,NXOUT) RESULT(VAROUT)
   INTEGER,               INTENT(IN) :: NX,NXOUT
   REAL,DIMENSION(NX),    INTENT(IN) :: X
   COMPLEX,DIMENSION(NX), INTENT(IN) :: VAR

   REAL,DIMENSION(NXOUT)             :: XOUT
   COMPLEX,DIMENSION(NXOUT)          :: VAROUT
   Type(linear_interp_1d)            :: interp1
   INTEGER                           :: iflag
   REAL,DIMENSION(NXOUT)             :: VAROUT_R,VAROUT_I

   VAROUT_R=FUN_INTERP1_REAL(X,REAL(VAR),NX,XOUT,NXOUT)
   VAROUT_I=FUN_INTERP1_REAL(X,AIMAG(VAR),NX,XOUT,NXOUT)
   VAROUT=CMPLX(VAROUT_R,VAROUT_I)
 END FUNCTION

END MODULE
