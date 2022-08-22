MODULE MQSolverOutputFiles

USE  MFileDirectoryList!, ONLY:OutQTFDir,OutFileDM,OutFileDP,OutFileHBM,&
                       !       OutFileHBP,make_directory 
IMPLICIT NONE

CONTAINS

        SUBROUTINE INITIALIZE_OUTPUT_FILES(workdir)
             CHARACTER(len=*)   :: workdir
             INTEGER            :: u,Iterm
             CHARACTER(len=1)   :: strT
             
             CALL make_directory(workdir//OutQTFDir)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileDM, ACTION='WRITE')
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileDP, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileHBM, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileHBP, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             DO Iterm=1,6
             WRITE(strT,'(I0.1)') Iterm
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileDM_term//strT//'.dat', &
                     ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileDP_term//strT//'.dat', &
                     ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileHBM_term//strT//'.dat', &
                     ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutQTFDir//OutFileHBP_term//strT//'.dat', &
                     ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             ENDDO
       END

       SUBROUTINE WRITE_QTF_DATA(wd,FileM,FileP,Ninteg,w1,w2,beta1,beta2,QTFdat)
          CHARACTER(LEN=*),    INTENT(IN)::wd,FileM,FileP
          INTEGER,             INTENT(IN)::Ninteg
          REAL,                INTENT(IN)::w1,w2,beta1,beta2
          COMPLEX,DIMENSION(Ninteg,2),INTENT(IN):: QTFdat
          Integer :: Iinteg,u1,u2

          OPEN(NEWUNIT=u1, FILE=wd//OutQTFDir//FileM, ACTION='WRITE',POSITION='APPEND')
          OPEN(NEWUNIT=u2, FILE=wd//OutQTFDir//FileP, ACTION='WRITE',POSITION='APPEND')

          DO Iinteg=1,Ninteg
             WRITE(u1,'(4(F10.3,X),I3,2(X,E14.7))')               & 
                  w1,w2,beta1,beta2,Iinteg,                       &
                  REAL(QTFdat(Iinteg,1)),AIMAG(QTFdat(Iinteg,1))
             WRITE(u2,'(4(F10.3,X),I3,2(X,E14.7))')               & 
                  w1,w2,beta1,beta2,Iinteg,                       &
                  REAL(QTFdat(Iinteg,2)),AIMAG(QTFdat(Iinteg,2))
          ENDDO
          CLOSE(u1)
          CLOSE(u2)
       END SUBROUTINE

END MODULE
