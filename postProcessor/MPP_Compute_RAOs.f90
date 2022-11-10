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
!   - R. Kurnia (2022)
!
!--------------------------------------------------------------------------------------
MODULE MPP_Compute_RAOs

    USE CONSTANTS,          ONLY:II,PI
    USE MResults
    USE MPP_ReadInputFiles, ONLY:TMech 
    USE M_SOLVER,           ONLY:LU_INVERS_MATRIX
    USE MNemohCal,          ONLY:TNemCal
!
    IMPLICIT NONE


CONTAINS
    SUBROUTINE Compute_RAOs(RAOS,Results,MechCoef)
!
!
!   Inputs/outputs
    TYPE(TResults),             INTENT(IN) :: Results
    TYPE(TMech),                INTENT(IN) :: MechCoef
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,Results%Nbeta) :: RAOs
!   Locals
    INTEGER :: Iw,Ibeta
    REAL    :: w
    COMPLEX,DIMENSION(Results%Nradiation,Results%Nintegration):: MAT_A,invMAT_A
    COMPLEX,DIMENSION(Results%Nintegration)                   :: ExcitForce
!
   DO Iw=1,Results%Nw
      w= Results%w(Iw)
      MAT_A=-(MechCoef%MassMat+Results%AddedMass(Iw,:,:))*w*w                   &
            -II*w*(Results%RadiationDamping(Iw,:,:)+MechCoef%DampCoefMat_EXT)   &
            +MechCoef%StiffMat+MechCoef%StiffMat_EXT
      CALL LU_INVERS_MATRIX(MAT_A,Results%Nradiation, invMAT_A)
      DO Ibeta=1,Results%Nbeta
       ExcitForce=Results%DiffractionForce(Iw,Ibeta,:)+Results%FroudeKrylovForce(Iw,Ibeta,:)
       RAOs(:,Iw,Ibeta)=MATMUL(invMAT_A,ExcitForce) 
      ENDDO
   ENDDO

    END SUBROUTINE Compute_RAOs

    SUBROUTINE SAVE_RAO(RAOs,w,beta,Nintegration,Nw,Nbeta,IndxForce,dirname,filename,InpNEMOHCAL)
    CHARACTER(LEN=*),                   INTENT(IN) :: filename,dirname
    TYPE(TNemCal),                      INTENT(IN) :: InpNEMOHCAL
    INTEGER,                            INTENT(IN) :: Nintegration,Nw, Nbeta
    REAL, DIMENSION(Nw),                INTENT(IN) :: w
    REAL, DIMENSION(Nbeta),             INTENT(IN) :: beta
    INTEGER,DIMENSION(Nintegration),    INTENT(IN) :: IndxForce
    COMPLEX,DIMENSION(Nintegration,Nw,Nbeta),&
                                        INTENT(IN) :: RAOs
    INTEGER :: Iw, Ibeta,Iinteg,u,Ninteg2
    REAL,DIMENSION(Nintegration*2) ::wRAOS 
    REAL                           ::RAOSphase 
    CHARACTER(LEN=23) :: FreqVar_text,FreqVar_text2
    REAL,DIMENSION(Nw)::freqVar
    
    IF (InpNEMOHCAL%OptOUTPUT%FreqType==1.OR.InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr) THEN
        FreqVar_text='VARIABLES="w (rad/s)"'
        freqVar=w
        FreqVar_text2='Number of pulsation= '
    ELSEIF (InpNEMOHCAL%OptOUTPUT%FreqType==2) THEN
        FreqVar_text='VARIABLES="f (Hz)"'
        freqVar=w/2/PI
        FreqVar_text2='Number of frequency= '
    ELSEIF (InpNEMOHCAL%OptOUTPUT%FreqType==3) THEN
        FreqVar_text='VARIABLES="T (s)"'
        freqVar=2*PI/w
        FreqVar_text2='Number of periode= '
    ENDIF


    CALL make_directory(dirname)
    Ninteg2=2*Nintegration
    OPEN(NEWUNIT=u, FILE=TRIM(dirname)//TRIM(filename),ACTION='WRITE')
    WRITE(u,'(13(A,X))') FreqVar_text, '|X| (m/m)','|Y| (m/m)',' |Z| (m/m)',  &
                '|phi| (deg)',' |theta| (deg)',' |psi| (deg)', &
                 'ang(x) (deg)',' ang(y) (deg)',' ang(z) (deg)', &
                 'ang(phi) (deg)',' ang(theta) (deg)',' ang(psi) (deg)'
    WRITE(u,*) 'Number of column (Nvariables*Nbody)=',Nintegration 
    DO Ibeta=1,Nbeta
    WRITE(u,'((A,X),F10.3,2(A,X),I4)') 'beta=',beta(Ibeta)*180/PI,'(deg),',FreqVar_text2,Nw
        DO Iw=1,Nw
             DO Iinteg=1,Nintegration
                wRAOS(Iinteg)=ABS(RAOs(Iinteg,Iw,Ibeta))
                IF (IndxForce(Iinteg).GT.3) wRAOS(Iinteg)=wRAOS(Iinteg)*180/PI
                RAOSphase=ATAN2(AIMAG(RAOs(Iinteg,Iw,Ibeta)),REAL(RAOs(Iinteg,Iw,Ibeta)))
                wRAOS(Iinteg+Nintegration)=180/PI*RAOSphase
             ENDDO
             WRITE(u,'((F10.3,X),<Ninteg2>(E14.7,X))') freqVar(Iw),(wRAOS(Iinteg),Iinteg=1,Ninteg2)
        ENDDO
    ENDDO
    CLOSE(u)

    END SUBROUTINE

   SUBROUTINE  make_directory(dirname)
          CHARACTER(LEN=*),       INTENT(IN) ::dirname
          LOGICAL                            ::existdir
          INQUIRE (DIRECTORY=dirname, EXIST=existdir)       
          IF (.NOT.existdir) CALL SYSTEM('mkdir '//dirname)
   END SUBROUTINE

END MODULE
