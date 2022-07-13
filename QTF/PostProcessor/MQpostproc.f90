!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - Postprocessing module
!   MQpostproc
!   Description:  read and write data
!   Contributors list:
!    - Ruddy Kurnia (ECN)
!--------------------------------------------------------------------------------------
!
MODULE MQpostproc

USE MNemohCal,          ONLY: TNemCal
USE CONSTANTS

IMPLICIT NONE
INTEGER, PARAMETER      :: IDDUOK=1,IDHASBO=1,IDHASFS=1,IDQTFP=1
!INTEGER, PARAMETER      :: STORE_HALF_DIAGONAL=1
CONTAINS

  SUBROUTINE READWRITE_QTFDATA(inpNEMCAL,wd)
          !Input variable
          TYPE(TNemCal),                 INTENT(IN) :: inpNEMCAL
          CHARACTER(LEN=*),              INTENT(IN) :: wd              
          !local variables
          REAL                                      :: betai,betaj
          REAL, DIMENSION(7)                        :: lineDUOK,lineHASBO,lineHASFS, &
                                                       lineASYMP
          REAL,ALLOCATABLE,DIMENSION(:,:)           :: QTFtotR,QTFtotI
          REAL,ALLOCATABLE,DIMENSION(:)             :: w
          INTEGER                                   :: NwQ,Nbeta,I,J,IDDOF,ID_MP
          INTEGER                                   :: IwQ,Iw1,Iw2
          INTEGER                                   :: uo_m             !unit write file
          INTEGER                                   :: u1,u2,u3,u4      !unit read file
          INTEGER                                   :: IDCONTRIB(3)
          INTEGER,DIMENSION(inpNEMCAL%Nintegtot)    :: DOF
          INTEGER                                   :: Ninteg,IintegS,Iinteg, &
                                                       Ibeta1,Ibeta2
          CHARACTER*1                               :: strDOF
          CHARACTER*1,DIMENSION(2)                  :: str_MP
                  
          !CALL CHECK_QTF_DATA_EXIST(inpNEMCAL,wd)
          IDCONTRIB(1:3)=0
          DOF(:)=1      !DOF to be write
          !DOF(1)=1
          !DOF(3)=1
          !DOF(5)=1       
          Nbeta= inpNEMCAL%waveinput%NBeta
          NwQ  = inpNEMCAL%qtfinput%omega(1)
          Ninteg=inpNEMCAL%Nintegtot
          str_MP(1)='M'
          str_MP(2)='P'
          DO ID_MP=1,InpNEMCAL%qtfinput%switch_QTFP+1
            !open output file 
            OPEN(NEWUNIT=uo_m, FILE=TRIM(wd)//'/results/QTF/OUT_QTF'//str_MP(ID_MP)//'_N.dat',&
                     ACTION='WRITE')
            WRITE(uo_m,'(9(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [deg]','beta2[deg]','DOF',&
                    'MOD(QTF)/rho/g','PHASE(QTF)[deg]','Re(QTF)/rho/g','Im(QTF)/rho/g'

              DO IintegS=1,Ninteg  
              !open QTF_DUOK file
                IF (InpNEMCAL%qtfinput%switch_qtfduok==IDDUOK) THEN
                   OPEN(NEWUNIT=u1, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_DUOK.dat',&
                           STATUS='UNKNOWN', ACTION='READ')
                   READ(u1,*) 
                   IDCONTRIB(1)=1
                ENDIF

                !open QTF_HASBO file
                IF (InpNEMCAL%qtfinput%switch_qtfhasbo==IDHASBO) THEN
                  OPEN(NEWUNIT=u2, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASBO.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u2,*)
                  IDCONTRIB(2)=1
                ENDIF
                !open QTF_HASFS file

                IF (InpNEMCAL%qtfinput%switch_qtfhasfs==IDHASFS) THEN
                  OPEN(NEWUNIT=u3, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASFS.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u3,*)
                  OPEN(NEWUNIT=u4, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_ASYMP.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u4,*)
                  IDCONTRIB(3)=1
                ENDIF

               
                DO Ibeta1=1,Nbeta
                    DO Ibeta2=1,Nbeta
                    ALLOCATE(w(NwQ))
                    ALLOCATE(QTFtotR(NwQ,NwQ))
                    ALLOCATE(QTFtotI(NwQ,NwQ))
                      DO IwQ=0,NwQ-1
                          DO Iw1=IwQ+1,NwQ
                           Iw2=Iw1-IwQ
                           DO Iinteg=1,Ninteg 
                            IF (Iinteg==IintegS) THEN
                             QTFtotR(Iw1,Iw2)=0
                             QTFtotI(Iw1,Iw2)=0
                             IF (IDCONTRIB(1)==1) THEN
                              READ(u1,*) lineDUOK(:)
                              w(Iw1)=lineDUOK(1)
                              betai=lineDUOK(3)*180/PI
                              betaj=lineDUOK(4)*180/PI
                              QTFtotR(Iw1,Iw2)=QtFtotR(Iw1,Iw2)+lineDUOK(6)
                              QTFtotI(Iw1,Iw2)=QtFtotI(Iw1,Iw2)+lineDUOK(7)
                             ENDIF
                             IF (IDCONTRIB(2)==1) THEN
                               READ(u2,*) lineHASBO(:)
                               w(Iw1)=lineHASBO(1)
                               betai=lineHASBO(3)*180/PI
                               betaj=lineHASBO(4)*180/PI
                               QTFtotR(Iw1,Iw2)=QtFtotR(Iw1,Iw2)+lineHASBO(6)
                               QTFtotI(Iw1,Iw2)=QtFtotI(Iw1,Iw2)+lineHASBO(7)
                             ENDIF
                             IF (IDCONTRIB(3)==1) THEN
                               READ(u3,*) lineHASFS(:)
                               READ(u4,*) lineASYMP(:)
                               w(Iw1)=lineHASFS(1)
                               betai=lineHASFS(3)*180/PI
                               betaj=lineHASFS(4)*180/PI
                               QTFtotR(Iw1,Iw2)=QtFtotR(Iw1,Iw2)+lineHASFS(6)
                               QTFtotI(Iw1,Iw2)=QtFtotI(Iw1,Iw2)+lineHASFS(7)
                               QTFtotR(Iw1,Iw2)=QtFtotR(Iw1,Iw2)+lineASYMP(6)
                               QTFtotI(Iw1,Iw2)=QtFtotI(Iw1,Iw2)+lineASYMP(7)
                             ENDIF
                            ELSE
                               IF (IDCONTRIB(1)==1) READ(u1,*)
                               IF (IDCONTRIB(2)==1) READ(u2,*)
                               IF (IDCONTRIB(3)==1) READ(u3,*)
                               IF (IDCONTRIB(3)==1) READ(u4,*)
                            ENDIF
                           ENDDO
                          ENDDO
                      ENDDO

                      DO Iw1=1,NwQ
                        DO Iw2=1,Iw1
                          QTFtotR(Iw1,Iw2)=QTFtotR(Iw1,Iw2)/inpNEMCAL%Env%RHO/inpNEMCAL%Env%G
                          QTFtotI(Iw1,Iw2)=QTFtotI(Iw1,Iw2)/inpNEMCAL%Env%RHO/inpNEMCAL%Env%G

                          WRITE(uo_m,'(4(F10.3,X),I2,4(X,E14.7))')              &
                                w(Iw1),w(Iw2),betai,betaj,IintegS,              &
                                SQRT(QTFtotR(Iw1,Iw2)**2+QTFtotI(Iw1,Iw2)**2),  &
                                ATAN(QTFtotI(Iw1,Iw2)/QTFtotR(Iw1,Iw2))*180/PI, &
                                QTFtotR(Iw1,Iw2),QTFtotI(Iw1,Iw2)
                        ENDDO
                      ENDDO
                    DEALLOCATE(w)
                    DEALLOCATE(QTFtotR)
                    DEALLOCATE(QTFtotI)
                    ENDDO
                ENDDO
                IF (IDCONTRIB(1)==1) CLOSE(u1)
                IF (IDCONTRIB(2)==1) CLOSE(u2)
                IF (IDCONTRIB(3)==1) CLOSE(u3)
                IF (IDCONTRIB(3)==1) CLOSE(u4)

            END DO
            print*,'results/OUT_QTF',str_MP(ID_MP),'.dat saved!'
            CLOSE(uo_m)
          END DO
         
  END SUBROUTINE

  LOGICAL FUNCTION  exist_file(inp,filename)
      TYPE(TNemCal),          INTENT(IN) :: inp
      CHARACTER(LEN=*),       INTENT(IN) :: filename
          INQUIRE (FILE=filename, EXIST=exist_file)       
          IF (.NOT.exist_file) THEN
               PRINT*,filename,' data is missing!'
          ENDIF
  END FUNCTION
END MODULE

