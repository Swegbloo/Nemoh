!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------

MODULE MQTFWrite

! write data of QTF solver in files

CONTAINS

    SUBROUTINE WriteFa(ID,Environment,Nw,w,Nbeta,beta,M,Kh,RAID,AMOR)
    
      USE MIdentification
      USE MEnvironment
      USE QTFCOM_VAR
    
      IMPLICIT NONE  
      
      TYPE(TID) :: ID                     ! Calculation identification data
      TYPE(TEnvironment) :: Environment   ! Environment data   
      INTEGER :: Nw
      REAL,DIMENSION(:),ALLOCATABLE :: w
      INTEGER :: Nbeta
      REAL,DIMENSION(:),ALLOCATABLE :: beta
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: M,Kh
      REAL,DIMENSION(:,:),ALLOCATABLE :: RAID,AMOR
      
      INTEGER :: i,j,k
      
      REAL :: T,AMZ,ILIN,VA,WR,AM0
      INTEGER :: IMIN,NM
      

      !!!!!!! unused parameters
      VA = 0
      T = 2.0*PI/w(1)
      IMIN = 1
      AMZ = 0
      ILIN = 0
      WR = w(1)
      AM0 = 0
      
      OPEN(UNIT=12,FILE=ID%ID(1:ID%lID)//'/QTF/FA.RES',ACCESS='DIRECT',STATUS='UNKNOWN',RECL=4*4*NFA)
      
      WRITE(12,REC=1)NC,NCO,NSYMY,NP,IMX,IXX,Environment%XEFF,Environment%YEFF,ZEFF,(IMXC(I),I=1,NCO),(IPOS(I),I=1,NC),ZER1,T,IMIN,Environment%Depth,AMZ,ILIN,Environment%RHO,VA,beta(1),WR,AM0,NIN,(AIND(I),I=1,NIN)
      WRITE(12,REC=2)(X(I),I=1,NP)
      WRITE(12,REC=3)(Y(I),I=1,NP)
      WRITE(12,REC=4)(Z(I),I=1,NP)
      WRITE(12,REC=5)(M1(I),I=1,IMX),(M2(I),I=1,IMX),(M3(I),I=1,IMX),(M4(I),I=1,IMX)
      WRITE(12,REC=6)(P(I),I=1,IMX),(Q(I),I=1,IMX),(R(I),I=1,IMX)
      WRITE(12,REC=7)(XM(I),I=1,IMX),(YM(I),I=1,IMX),(ZM(I),I=1,IMX)
      WRITE(12,REC=8)(AIRE(I),I=1,IMX),(TDIS(I),I=1,IMX),  (TDIS(I),I=1,IMX)!  (DIST(I),I=1,IMX)   AC: We write again TDIS instead of DIST that is was not precomputed
      WRITE(12,REC=9)(XM(I),I=IMX+1,IXX),(YM(I),I=IMX+1,IXX),(ALF(I),I=1,NFFL),(IND(I),I=1,NFFL)
      
      !AC: NOTE: This should not be in fa.RES because it is not geometrical properties
      !AC: NOTE: Only Inertia and Hydrostatic matrix are used for QTF. Additional stiffness and damping are needed by mechanical solver only (so not here).
      NM=6*NC
      WRITE(12,REC=10)(((M(I,J,K),I=1,6),J=1,6),((Kh(I,J,K),I=1,6),J=1,6),K=1,NC),(( RAID(I,J) ,I=1,NM),J=1,NM),(( AMOR(I,J) ,I=1,NM),J=1,NM)  ! ,(XCDG(K),YCDG(K),ZCDG(K),K=1,NC) AC:unused
      
      CLOSE(UNIT=12)
      
    END SUBROUTINE WriteFa
    
    
    SUBROUTINE QTFWriteSing(ID,iw,Nw,Nradiation,Nbeta,Npanels, ZIGB,ZIGS)
      
      USE MIdentification
      USE QTFCOM_VAR
    
      IMPLICIT NONE  
      
      TYPE(TID) :: ID                     ! Calculation identification data
      INTEGER :: Nw,Nradiation,Nbeta,Npanels,k,i
      COMPLEX,DIMENSION(:,:),ALLOCATABLE :: ZIGB,ZIGS
      COMPLEX :: chgfmt1,chgfmt2
      CHARACTER*5 :: str
      INTEGER :: iw 
      
      WRITE(str,'(I5)') iw
      
      OPEN(UNIT=10,FILE=ID%ID(1:ID%lID)//'/QTF/sing.wat',ACCESS='DIRECT', STATUS='REPLACE',RECL=4*4*NFA)
!       OPEN(UNIT=11,FILE=ID%ID(1:ID%lID)//'/QTF/sing'//str//'.dat')
      
      do k=1,(Nradiation+Nbeta)
	WRITE(10,REC=k)(ZIGB(i,k),i=1,Npanels),(ZIGS(i,k),i=1,Npanels)
! 	WRITE(11,*)(ZIGB(i,k),i=1,Npanels),(ZIGS(i,k),i=1,Npanels)
! 	WRITE(*,*) k, (ZIGB(i,k),i=1,Npanels),(ZIGS(i,k),i=1,Npanels)
      end do
      
!       WRITE(*,*) "-------"
!       
!       CLOSE(11)
      CLOSE(10)
      
    END SUBROUTINE 
    
END MODULE 

