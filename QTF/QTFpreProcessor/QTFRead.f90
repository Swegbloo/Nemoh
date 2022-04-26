!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------

MODULE MQTFRead

! read data needed by QTF solver in files

CONTAINS

  SUBROUTINE Readwbeta(ID,Nbodies,Nw,w,Nbeta,beta,Nradiation,Nintegration)

    USE MIdentification
    USE QTFCOM_VAR
    
    IMPLICIT NONE

    TYPE(TID) :: ID                     ! Calculation identification data
    INTEGER :: Nbodies,Nw,Nradiation,Nbeta,Nintegration
    REAL :: wmin,wmax,betamax,betamin,dw,dbeta, dwtemp
    REAL,DIMENSION(:),ALLOCATABLE :: w
    REAL,DIMENSION(:),ALLOCATABLE :: beta
    
    INTEGER :: M,c,i,j, FreqType

    Nradiation=0
    Nintegration=0

    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
    DO c=1,7
	READ(10,*)
    END DO

    DO c=1,Nbodies
	DO i=1,3
	    READ(10,*)
	END DO
	READ(10,*) M
	Nradiation=Nradiation+M
	DO i=1,M
	    READ(10,*)
	END DO
	READ(10,*) M
	Nintegration=Nintegration+M
	DO i=1,M
	    READ(10,*)
	END DO
	READ(10,*) M
	DO i=1,M
	    READ(10,*)
	END DO
    END DO
    READ(10,*)
    READ(10,*) FreqType,Nw,wmin,wmax
    READ(10,*) Nbeta,betamin,betamax

    CLOSE(10)
    IF (FreqType.Eq.2) THEN
        wmin=2*PI*wmin
        wmax=2*PI*wmax

    ELSEIF (FreqType.Eq.3) THEN
         WRITE(*,*) 'Input periode is not allow in QTF, please changed in NEMOH.CAL' 
         STOP
    END IF

         
  !   ----- Fill w and beta ----------
    
    ALLOCATE(w(Nw))
    IF (Nw.GT.1) THEN
	
	dw = (wmax-wmin)/(Nw-1)
        dwtemp = wmax/Nw
        IF (abs(dwtemp-dw).LT.0.001) THEN  ! added by RK
        dw =wmin  
        END IF

	DO j=1,Nw
	    w(j)=wmin+dw*(j-1)
	END DO
	
	IF(dw.NE.wmin) THEN
	  WRITE(*,*) "ERROR: minimum frequency must be equal to frequency step!"
	END IF
	
    ELSE
	WRITE(*,*) "ERROR: QTF can not be run with 1 frequency only!"
    END IF

    ALLOCATE(beta(Nbeta))
    IF (Nbeta.GT.1) THEN

    !         DO j=1,Nbeta
    !             beta(j)=(betamin+(betamax-betamin)*(j-1)/(Nbeta-1))*PI/180.
    !         END DO
	WRITE(*,*) "ERROR: QTF not available yet for several wave directions!"
	
    ELSE
	beta(1)=betamin*PI/180.
    END IF

  END SUBROUTINE

  SUBROUTINE ReadInertia(ID,Nbodies,M)

    USE MIdentification

    IMPLICIT NONE

    TYPE(TID) :: ID                     ! Calculation identification data
    INTEGER :: Nbodies
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: M
    REAL :: mass
    INTEGER :: i,j,c
    
    ALLOCATE(M(6,6,1))
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Inertia.dat')
    
    DO c=1,Nbodies
    
      DO i=1,6
! 	DO j=1,6
	
! 	  M(i,j,c) = 0.0
	  READ(10,*) (M(i,j,c), j=1,6)
	
! 	END DO
      END DO
      
!       READ(10,*) mass
      
!       M(1,1,c) = mass
!       M(2,2,c) = mass
!       M(3,3,c) = mass
      
!       DO i=1,3
! 	READ(10,*) (M(i+3,j+3,c), j=1,3)
!       END DO
      
    END DO
    
    CLOSE(10)

  END SUBROUTINE ReadInertia

  SUBROUTINE ReadKh(ID,Nbodies,Kh)

    USE MIdentification

    IMPLICIT NONE

    TYPE(TID) :: ID                     ! Calculation identification data
    INTEGER :: Nbodies
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: Kh
    INTEGER :: i,j,c
    
    ALLOCATE(Kh(6,6,Nbodies))
    
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Kh.dat')
    
    DO c=1,Nbodies
      
      DO i=1,6
	READ(10,*) (Kh(i,j,c), j=1,6)
      END DO
      
    END DO
    
    CLOSE(10)

  END SUBROUTINE ReadKh

  SUBROUTINE ReadRaidEXT(ID,Nbodies,RAID)

    USE MIdentification

    IMPLICIT NONE

    TYPE(TID) :: ID                     ! Calculation identification data
    INTEGER :: Nbodies
    REAL,DIMENSION(:,:),ALLOCATABLE :: RAID
    INTEGER :: i,j,c
    
    ALLOCATE(RAID(6*Nbodies,6*Nbodies))
    
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Km.dat')
    DO i=1,6*Nbodies
	READ(10,*) (RAID(i,j), j=1,6*Nbodies)
    END DO
    CLOSE(10)

  END SUBROUTINE ReadRaidEXT

  SUBROUTINE ReadAmorEXT(ID,Nbodies,AMOR)

    USE MIdentification

    IMPLICIT NONE

    TYPE(TID) :: ID                     ! Calculation identification data
    INTEGER :: Nbodies
    REAL,DIMENSION(:,:),ALLOCATABLE :: AMOR
    INTEGER :: i,j,c
    
    ALLOCATE(AMOR(6*Nbodies,6*Nbodies))
    
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mechanics/Km.dat') !???? Checked the file name
    DO i=1,6*Nbodies
	READ(10,*) (AMOR(i,j), j=1,6*Nbodies)
    END DO
    CLOSE(10)
    
    CLOSE(10)

  END SUBROUTINE ReadAmorEXT
  
  SUBROUTINE QTFReadFirstOrderLoad(ID, Nw,Nradiation,Nbeta,Nintegration, CA0, CM0, Fe)
  
    USE MIdentification
    USE QTFCOM_VAR

    IMPLICIT NONE
    
    TYPE(TID) :: ID 
    INTEGER :: Nw, Nbeta, Nradiation, Nintegration
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: CA0,CM0
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: Fe
    REAL :: Amp(6), Phase(6)
    
    REAL :: w0
    INTEGER :: j,k,l
    
    ALLOCATE(CA0(Nw,Nradiation,Nradiation))
    ALLOCATE(CM0(Nw,Nradiation,Nradiation))
    ALLOCATE(Fe(Nw,Nintegration,Nbeta))
    
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/CA.dat')
    READ(10,*) 
    DO l=1,Nw
      READ(10,*) ! TODO : check if frequencies are consistent between 1st and 2nd order 
      DO j=1,Nradiation
	  READ(10,'(6(X,E13.6))') ( CA0(l,j,k),k=1,Nradiation )
      END DO
    END DO
    CLOSE(10)

    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/CM.dat')
    READ(10,*) 
    DO l=1,Nw
      READ(10,*) ! TODO : check if frequencies are consistent between 1st and 2nd order 
      DO j=1,Nradiation
	  READ(10,'(6(X,E13.6))') ( CM0(l,j,k),k=1,Nradiation )
      END DO
    END DO
    CLOSE(10)

    !!!!!!!!!!   Correction temporaire pulsation irregulaires !!!!!!!!!!
    ! JND : 09/2015
    ! pulsation irreguliere : 39 et 32
!~     DO j =1,Nradiation 
!~ 	DO k=1,Nradiation
!~ 	    CA0(39,j,k)=(  CA0(37,j,k)+2*CA0(40,j,k))/3.
!~ 	    CM0(39,j,k)=(  CM0(37,j,k)+2*CM0(40,j,k))/3.
!~ 	    CA0(38,j,k)=(2*CA0(37,j,k)+  CA0(40,j,k))/3.
!~ 	    CM0(38,j,k)=(2*CM0(37,j,k)+  CM0(40,j,k))/3.
!~ 	    
!~ 	    CA0(32,j,k)=(3*CA0(30,j,k)+  CA0(34,j,k))/4.
!~ 	    CM0(32,j,k)=(3*CM0(30,j,k)+  CM0(34,j,k))/4.
!~ 	    CA0(31,j,k)=(2*CA0(30,j,k)+2*CA0(34,j,k))/4.
!~ 	    CM0(31,j,k)=(2*CM0(30,j,k)+2*CM0(34,j,k))/4.
!~ 	    CA0(33,j,k)=(  CA0(30,j,k)+3*CA0(34,j,k))/4.
!~ 	    CM0(33,j,k)=(  CM0(30,j,k)+3*CM0(34,j,k))/4.
!~ 	    
!~ 	ENDDO
!~     ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/Fe.dat')
    READ(10,*) 
    READ(10,*)
    DO j=1,Nbeta
	READ(10,*)
	DO l=1,Nw
!         TODO: CHECK PHASES CONVENTION !!!!!!!!!!!  
	  READ(10,'(F7.4,6(X,E13.6),6(X,F7.2))') w0,(Amp(k),k=1,Nintegration),(Phase(k),k=1,Nintegration)
	  
	  DO k=1,Nintegration
	    
	    Phase(k) = Phase(k)*PI/180.0  ! pass to radians
	    
	    !A.C. : WARNING: phase convention is modified here
	   ! Phase(k) = Phase(k)-PI/2.0   ! closed by RK
	    
	    Fe(l,k,j) = Amp(k)*CEXP(CMPLX(0.,1.)*Phase(k))
	  
	  END DO
	  
	END DO
    END DO
    CLOSE(10)
  
    !!!!!!!!!!   Correction temporaire pulsation irregulaires !!!!!!!!!!
    ! JND : 09/2015
    ! pulsation irreguliere : 39 et 32
!~     DO j =1,Nbeta
!~ 	DO k=1,Nintegration
!~ 	    Fe(39,k,j)=(  Fe(37,k,j)+2*Fe(40,k,j))/3.
!~ 	    Fe(38,k,j)=(2*Fe(37,k,j)+  Fe(40,k,j))/3.
!~ 	    
!~ 	    Fe(32,k,j)=(3*Fe(30,k,j)+  Fe(34,k,j))/4.
!~ 	    Fe(31,k,j)=(2*Fe(30,k,j)+2*Fe(34,k,j))/4.
!~ 	    Fe(33,k,j)=(  Fe(30,k,j)+3*Fe(34,k,j))/4.
!~ 	    
!~ 	ENDDO
!~     ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  
  END SUBROUTINE QTFReadFirstOrderLoad
  
  SUBROUTINE QTFReadGrin(ID)
    
    USE MIdentification
    USE QTFCOM_VAR

    IMPLICIT NONE
    
    TYPE(TID) :: ID 
    INTEGER :: IR,JZ,I,J
    
    OPEN(UNIT=44,FILE=ID%ID(1:ID%lID)//'/QTF/GRIN.QAT',FORM='UNFORMATTED',STATUS='OLD')
    READ(44)IR,JZ,(XR(I),I=1,IR),(XZ(J),J=1,JZ)
    Do I=1,IR
        Do J=1,JZ
        write(*,'(I5,I5,E15.4,E15.4)') IR,JZ,XR(I),XZ(J)
        Enddo
    Enddo
    DO 1515 J=1,JZ
    READ(44)(APD1X(I,J),I=1,IR),(APD1Z(I,J),I=1,IR), (APD2X(I,J),I=1,IR),(APD2Z(I,J),I=1,IR)
    DO I=1,IR
    write(*,'(E15.4,E15.4,E15.4,E15.4)') APD1X(I,J),APD1Z(I,J),APD2X(I,J),APD2Z(I,J)
    ENDDO
    1515 CONTINUE
    CLOSE(UNIT=44)
    read(*,*)
  END SUBROUTINE QTFReadGrin
  
  SUBROUTINE QTFReadMotion(ID,Nw,Nintegration,Nbeta,RAO)
    
    USE MIdentification
    USE QTFCOM_VAR

    IMPLICIT NONE
    
    TYPE(TID) :: ID
    INTEGER :: Nw, Nbeta, Nintegration
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: RAO
    REAL :: Amp(6), Phase(6), w0
    INTEGER :: j,k,l
    
    ALLOCATE(RAO(Nw,Nintegration,Nbeta))
    
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Motion/RAO.dat')
    READ(10,*) 
    READ(10,*)
    DO j=1,Nbeta
	READ(10,*)
	DO l=1,Nw
!         TODO: CHECK PHASES CONVENTION !!!!!!!!!!!  
	  READ(10,*) w0,(Amp(k),k=1,Nintegration),(Phase(k),k=1,Nintegration)
	  
	  DO k=1,Nintegration
	    
	    Phase(k) = Phase(k)*PI/180.0  ! pass to radians
	    
	    !A.C. : WARNING: phase convention is modified here
 	    ! Phase(k) = Phase(k)-PI/2.0 ! close by RK   
	    IF (k>3) THEN
		Amp(k)=Amp(k)*PI/180.0	  ! JND :degre de liberte en rotation en deg/m --> rad/m 
	    ENDIF
	    
	    RAO(l,k,j) =Amp(k)*CEXP(CMPLX(0.,1.)*Phase(k))
	    
	  END DO

	END DO
    !!!!!!!!!!   Correction temporaire pulsation irregulaires !!!!!!!!!!
    ! JND : 09/2015
    ! pulsation irreguliere : 39 et 32
!~ 	DO k=1,Nintegration
!~ 	    RAO(31,k,j) = (3*RAO(30,k,j)+  RAO(34,k,j))/4.
!~ 	    RAO(32,k,j) = (2*RAO(30,k,j)+2*RAO(34,k,j))/4.
!~ 	    RAO(33,k,j) = (  RAO(30,k,j)+3*RAO(34,k,j))/4.
!~ 	    
!~ 	    RAO(38,k,j) = (2*RAO(37,k,j)+  RAO(40,k,j))/3.
!~ 	    RAO(39,k,j) = (  RAO(37,k,j)+2*RAO(40,k,j))/3.
!~ 	ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO
    CLOSE(10)
    
  END SUBROUTINE QTFReadMotion
  
  SUBROUTINE QTFReadSources(ID,iw,Nw,Nradiation,Nbeta,Npanels, ZIGB,ZIGS)
  
    USE MIdentification
    USE QTFCOM_VAR
    
    IMPLICIT NONE
    
    TYPE(TID) :: ID
    INTEGER :: iw,Nw,Nradiation,Nbeta,Npanels, Pbnumber, idiffrad,i,k
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: ZIGB,ZIGS
    COMPLEX,DIMENSION(Npanels,Nradiation+Nbeta) :: ZIGBl,ZIGSl,ZIGBh,ZIGSh
    CHARACTER*5 :: str
    REAL :: RE,IM
    
    idiffrad = 0
    ! read radiation sources at frequency w
    DO k = 1,Nradiation
    
      idiffrad = idiffrad + 1
      Pbnumber = (Nbeta+Nradiation)*(iw-1)+Nbeta+k
    
      !WRITE(str,'(I5)') Pbnumber 
      WRITE(str, '(I0.5)') Pbnumber
      WRITE(*,'(I5,A,I5)') Pbnumber, 'iw=',iw
      OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources/sources.'//str//'.dat')
   
      do i=1,Npanels
	READ(11,*) RE,IM 
       ZIGB(I,idiffrad) = CMPLX(RE,IM)

      end do
      
      do i=1,Npanels
	READ(11,*) RE,IM 
       ZIGS(I,idiffrad) = CMPLX(RE,IM)
      end do
     
      CLOSE(11)
    
    END DO
    
    DO k=1,Nbeta
    
      idiffrad = idiffrad + 1
      Pbnumber = (Nbeta+Nradiation)*(iw-1)+k
      
      !WRITE(str,'(I5)') Pbnumber
      WRITE(str, '(I0.5)') Pbnumber
      OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources/sources.'//str//'.dat')
   
      do i=1,IMX
	READ(11,*) RE,IM 
 	ZIGB(I,idiffrad) = CMPLX(RE,IM)        ! Now convention same as the potential in first order NEMOH  
	! Change convention for diffraction problems
	! ZIGB(I,idiffrad) = CMPLX(IM,-RE)  ! comments RK: it is multiply by [-i] also compensated in PhiI (missing i) 
      end do
      
      do i=1,IMX
	READ(11,*) RE,IM 
 	ZIGS(I,idiffrad) = CMPLX(RE,IM)
	! Change convention for difffraction problems
	!ZIGS(I,idiffrad) = CMPLX(IM,-RE)
	
      end do
    
      CLOSE(11)
    
    END DO
    !!!!!!!!!!   Correction temporaire pulsation irregulaires !!!!!!!!!!
    ! JND : 09/2015
    ! pulsation irreguliere : 39 et 32
!~     IF (iw<32+2 .AND. iw>32-2) THEN
!~ 	idiffrad = 0
!~ 	! read radiation sources at frequency w
!~ 	DO k = 1,Nradiation
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(30-1)+Nbeta+k
!~ 	
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGBl(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	  
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGSl(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	 
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	
!~ 	DO k=1,Nbeta
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(30-1)+k
!~ 	  
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGB(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGBl(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	  
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGS(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGSl(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	idiffrad = 0
!~ 	! read radiation sources at frequency w
!~ 	DO k = 1,Nradiation
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(34-1)+Nbeta+k
!~ 	
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGBh(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	  
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGSh(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	 
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	
!~ 	DO k=1,Nbeta
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(34-1)+k
!~ 	  
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGB(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGBh(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	  
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGS(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGSh(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	
!~ 	idiffrad = 0
!~ 	DO k=1,Nradiation
!~ 	    idiffrad = idiffrad + 1
!~ 	    DO i=1,Npanels
!~ 		ZIGB(I,idiffrad)=((iw-30)*ZIGBh(I,idiffrad)+(34-iw)*ZIGBl(I,idiffrad))/4.
!~ 		ZIGS(I,idiffrad)=((iw-30)*ZIGSh(I,idiffrad)+(34-iw)*ZIGSl(I,idiffrad))/4.
!~ 	    ENDDO
!~ 	ENDDO
!~ 	DO k=1,Nbeta
!~ 	    idiffrad = idiffrad + 1
!~ 	    DO i=1,IMX
!~ 		ZIGB(I,idiffrad)=((iw-30)*ZIGBh(I,idiffrad)+(34-iw)*ZIGBl(I,idiffrad))/4.
!~ 		ZIGS(I,idiffrad)=((iw-30)*ZIGSh(I,idiffrad)+(34-iw)*ZIGSl(I,idiffrad))/4.
!~ 	    ENDDO
!~ 	ENDDO
!~     ELSEIF (iw<40 .AND. iw>37) THEN
!~ 	idiffrad = 0
!~ 	! read radiation sources at frequency w
!~ 	DO k = 1,Nradiation
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(37-1)+Nbeta+k
!~ 	
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGBl(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	  
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGSl(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	 
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	
!~ 	DO k=1,Nbeta
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(37-1)+k
!~ 	  
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGB(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGBl(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	  
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGS(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGSl(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	idiffrad = 0
!~ 	! read radiation sources at frequency w
!~ 	DO k = 1,Nradiation
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(40-1)+Nbeta+k
!~ 	
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGBh(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	  
!~ 	  do i=1,Npanels
!~ 	    READ(11,*) RE,IM 
!~ 	    ZIGSh(I,idiffrad) = CMPLX(RE,IM)
!~ 	  end do
!~ 	 
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO
!~ 	
!~ 	DO k=1,Nbeta
!~ 	
!~ 	  idiffrad = idiffrad + 1
!~ 	  Pbnumber = (Nbeta+Nradiation)*(40-1)+k
!~ 	  
!~ 	  WRITE(str,'(I5)') Pbnumber
!~ 	  OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!~        
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGB(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGBh(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	  
!~ 	  do i=1,IMX
!~ 	    READ(11,*) RE,IM 
!~     ! 	ZIGS(I,idiffrad) = CMPLX(RE,IM)
!~ 	    ! Change convention for difffraction problems
!~ 	    ZIGSh(I,idiffrad) = CMPLX(IM,-RE)
!~ 	    
!~ 	  end do
!~ 	
!~ 	  CLOSE(11)
!~ 	
!~ 	END DO    
!~ 	
!~ 	
!~ 	
!~ 	idiffrad = 0
!~ 	DO k=1,Nradiation
!~ 	    idiffrad = idiffrad + 1
!~ 	    DO i=1,Npanels
!~ 		ZIGB(I,idiffrad)=((iw-37)*ZIGBh(I,idiffrad)+(40-iw)*ZIGBl(I,idiffrad))/3.
!~ 		ZIGS(I,idiffrad)=((iw-37)*ZIGSh(I,idiffrad)+(40-iw)*ZIGSl(I,idiffrad))/3.
!~ 	    ENDDO
!~ 	ENDDO
!~ 	DO k=1,Nbeta
!~ 	    idiffrad = idiffrad + 1
!~ 	    DO i=1,IMX
!~ 		ZIGB(I,idiffrad)=((iw-37)*ZIGBh(I,idiffrad)+(40-iw)*ZIGBl(I,idiffrad))/3.
!~ 		ZIGS(I,idiffrad)=((iw-37)*ZIGSh(I,idiffrad)+(40-iw)*ZIGSl(I,idiffrad))/3.
!~ 	    ENDDO
!~ 	ENDDO
!~     ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE QTFReadSources
  
  
   
END MODULE 
