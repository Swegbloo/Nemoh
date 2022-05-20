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
INTEGER, PARAMETER      :: STORE_HALF_DIAGONAL=1
CONTAINS

  SUBROUTINE READWRITE_QTFDATA(inpNEMCAL,wd)
          !Input variable
          TYPE(TNemCal),          INTENT(IN) :: inpNEMCAL
          CHARACTER(LEN=*),       INTENT(IN) :: wd              
          !local variables
          REAL, ALLOCATABLE,DIMENSION(:)  :: w
          REAL, ALLOCATABLE,DIMENSION(:)  :: lineDUOKR,lineHASBOR,lineHASFSR, &
                                             lineDUOKI,lineHASBOI,lineHASFSI, &
                                             lineASYMPR,lineASYMPI
          REAL, ALLOCATABLE,DIMENSION(:)  :: QTFtotR,QTFtotI
          INTEGER                         :: Nw,I,J,IDDOF,ID_MP
          INTEGER                         :: uo_m             !unit write file
          INTEGER                         :: u1R,u2R,u3R,u4R  !unit read file
          INTEGER                         :: u1I,u2I,u3I,u4I  !unit read file
          INTEGER                         :: IDCONTRIB(3)
          INTEGER                         :: DOF(6)
          CHARACTER*1                     :: strDOF
          CHARACTER*1,DIMENSION(2)        :: str_MP        
          !CALL CHECK_QTF_DATA_EXIST(inpNEMCAL,wd)
          IDCONTRIB(:)=0
          DOF(:)=0      !DOF to be write
          DOF(1)=1
          DOF(3)=1
          DOF(5)=1        
          Nw=inpNEMCAL%waveinput%NFreq
          str_MP(1)='M'
          str_MP(2)='P'
          DO ID_MP=1,InpNEMCAL%qtfinput%switch_QTFP+1
            !open output file 
            OPEN(NEWUNIT=uo_m, FILE=TRIM(wd)//'/results/QTF/OUT_QTF'//str_MP(ID_MP)//'.dat',&
                     ACTION='WRITE')
            WRITE(uo_m,'(9(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [deg]','beta2[deg]','DOF',&
                    'MOD(QTF)/rho/g','PHASE(QTF)[deg]','Re(QTF)/rho/g','Im(QTF)/rho/g'
            DO IDDOF=1,6
              ALLOCATE(QTFtotR(Nw),QTFtotI(Nw),w(Nw))
              ALLOCATE(lineDUOKR(Nw+1),lineHASBOR(Nw+1),lineHASFSR(Nw+1),lineASYMPR(Nw+1))
              ALLOCATE(lineDUOKI(Nw+1),lineHASBOI(Nw+1),lineHASFSI(Nw+1),lineASYMPI(Nw+1))
              
              !open QTF_DUOK file
              IF (DOF(IDDOF).EQ.1) THEN
              WRITE(strDOF,'(I0)') IDDOF
              print*,'QTF',str_MP(ID_MP),'DOF',strDOF
                IF (InpNEMCAL%qtfinput%switch_qtfduok==IDDUOK) THEN
                   OPEN(NEWUNIT=u1R, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_DUOK_DOF_'//strDOF//'_PR.dat',&
                           STATUS='UNKNOWN', ACTION='READ')
                   READ(u1R,*) lineDUOKR(:)
                   w=lineDUOKR(2:Nw+1)
                   OPEN(NEWUNIT=u1I, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_DUOK_DOF_'//strDOF//'_PI.dat',&
                           STATUS='UNKNOWN', ACTION='READ')
                   READ(u1I,*) lineDUOKI(:)
                   IDCONTRIB(1)=1
                ENDIF

                !open QTF_HASBO file
                IF (InpNEMCAL%qtfinput%switch_qtfhasbo==IDHASBO) THEN
                  OPEN(NEWUNIT=u2R, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASBO_DOF_'//strDOF//'_PR.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u2R,*) lineHASBOR(:)
                  w=lineHASBOR(2:Nw+1)
                  OPEN(NEWUNIT=u2I, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASBO_DOF_'//strDOF//'_PI.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u2I,*) lineHASBOI(:)
                   IDCONTRIB(2)=1
                ENDIF
                !open QTF_HASFS file
                IF (InpNEMCAL%qtfinput%switch_qtfhasfs==IDHASFS) THEN
                  OPEN(NEWUNIT=u3R, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASFS_DOF_'//strDOF//'_PR.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u3R,*) lineHASFSR(:)
                  w=lineHASFSR(2:Nw+1)
                  OPEN(NEWUNIT=u3I, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_HASFS_DOF_'//strDOF//'_PI.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u3I,*) lineHASFSI(:)
                   OPEN(NEWUNIT=u4R, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_ASYMP_DOF_'//strDOF//'_PR.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u4R,*) lineASYMPR(:)
                  OPEN(NEWUNIT=u4I, FILE=TRIM(wd)//'/results/QTF/QTF'//str_MP(ID_MP)//'_ASYMP_DOF_'//strDOF//'_PI.dat',&
                          STATUS='UNKNOWN', ACTION='READ')
                  READ(u4I,*) lineASYMPI(:)
                  IDCONTRIB(3)=1
                ENDIF
         
                DO I=1,Nw
                     QTFtotR(:)=0
                     QTFtotI(:)=0

                    IF (IDCONTRIB(1)==1) THEN
                        READ(u1R,*) lineDUOKR(:)
                        READ(u1I,*) lineDUOKI(:)
                        QTFtotR(:)=QTFtotR(:)+lineDUOKR(2:Nw+1)
                        QTFtotI(:)=QTFtotI(:)+lineDUOKI(2:Nw+1)
                    ENDIF
                    
                    IF (IDCONTRIB(2)==1) THEN
                        READ(u2R,*) lineHASBOR(:)
                        READ(u2I,*) lineHASBOI(:)
                        QTFtotR=QTFtotR+lineHASBOR(2:Nw+1)
                        QTFtotI=QTFtotI+lineHASBOI(2:Nw+1)
                    ENDIF
                    IF (IDCONTRIB(3)==1) THEN
                        READ(u3R,*) lineHASFSR(:)
                        READ(u3I,*) lineHASFSI(:)
                        READ(u4R,*) lineASYMPR(:)
                        READ(u4I,*) lineASYMPI(:)
                        QTFtotR=QTFtotR+lineHASFSR(2:Nw+1)+lineASYMPR(2:Nw+1)
                        QTFtotI=QTFtotI+lineHASFSI(2:Nw+1)+lineASYMPI(2:Nw+1)
                    ENDIF
                    QTFtotR=QTFtotR/inpNEMCAL%Env%RHO/inpNEMCAL%Env%G
                    QTFtotI=QTFtotI/inpNEMCAL%Env%RHO/inpNEMCAL%Env%G
                    DO J=I**STORE_HALF_DIAGONAL,Nw
                            WRITE(uo_m,'(4(F10.3,X),I2,4(X,E14.7))') w(I),w(J), 0.,0.,IDDOF,         &
                               SQRT(QTFtotR(J)**2+QTFtotI(J)**2),ATAN(QTFtotI(J)/QTFtotR(J))*180/PI, &
                                      QTFtotR(J),QTFtotI(J)
                    END DO
                END DO
                CLOSE(u1R)
                CLOSE(u1I)
                CLOSE(u2R)
                CLOSE(u3I)
              ENDIF
              DEALLOCATE(QTFtotR,QTFtotI,w)
              DEALLOCATE(lineDUOKR,lineHASBOR,lineHASFSR,lineASYMPR)
              DEALLOCATE(lineDUOKI,lineHASBOI,lineHASFSI,lineASYMPI)
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

