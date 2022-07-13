MODULE MQSolverOutputFiles

IMPLICIT NONE
CHARACTER(LEN=*),PARAMETER ::  OutFileDM='/results/QTF/QTFM_DUOK.dat'
CHARACTER(LEN=*),PARAMETER ::  OutFileDP='/results/QTF/QTFP_DUOK.dat'
CHARACTER(LEN=*),PARAMETER ::  OutFileHBM='/results/QTF/QTFM_HASBO.dat'
CHARACTER(LEN=*),PARAMETER ::  OutFileHBP='/results/QTF/QTFP_HASBO.dat'
  !Output variables
  !DUOK      : Quadratic QTF
  !HASBO     : Potential QTF-Haskind on Body term
  !HASFS     : Potential QTF-Haskind on Free-Surface in Finite domain
  !HASFS_ASYM: Potential QTF-Haskind on Free-Surface in Asymptotic domain

CONTAINS

        SUBROUTINE INITIALIZE_OUTPUT_FILES(workdir)
             CHARACTER(len=*)   ::workdir
             INTEGER            :: u

             OPEN(NEWUNIT=u,FILE=workdir//OutFileDM, ACTION='WRITE')
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutFileDP, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutFileHBM, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
             OPEN(NEWUNIT=u,FILE=workdir//OutFileHBP, ACTION='WRITE') 
             WRITE(u,'(7(A,X))') 'w1[rad/s]','w2[rad/s]','beta1 [rad]','beta2[rad]',&
                                  'DOF','Re(QTF)','Im(QTF)'
             CLOSE(u)
       END

       SUBROUTINE WRITE_QTF_DATA(wd,FileM,FileP,Ninteg,w1,w2,beta1,beta2,QTFdat)
          CHARACTER(LEN=*),    INTENT(IN)::wd,FileM,FileP
          INTEGER,             INTENT(IN)::Ninteg
          REAL,                INTENT(IN)::w1,w2,beta1,beta2
          COMPLEX,DIMENSION(Ninteg,2),INTENT(IN):: QTFdat
          Integer :: Iinteg,u1,u2

          OPEN(NEWUNIT=u1, FILE=wd//FileM, ACTION='WRITE',POSITION='APPEND')
          OPEN(NEWUNIT=u2, FILE=wd//FileP, ACTION='WRITE',POSITION='APPEND')

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
