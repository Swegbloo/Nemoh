!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------

PROGRAM Main

  ! Init of variables and loop over frequencies are made here

  USE MIdentification
  USE MEnvironment
  USE MMesh
  USE QTFCOM_VAR   ! Check QTFCOM_VAR to see the global variables
  USE MQTFRead     ! Functions reading Nemoh 1st order output and other input files
  USE MQTFWrite    ! Functions writing files for 2nd order solver
  USE MQTFGeometry ! The functions doing geometrical preprocessing on the body mesh for 2nd order (contour)  
  USE MQTFInit     ! The actual preprocessing job for QTF 
  USE MQTFBaseFunctions

  IMPLICIT NONE

  TYPE(TID) :: ID                     ! Calculation identification data
  TYPE(TMesh) :: Mesh                 ! Mesh data
  TYPE(TEnvironment) :: Environment   ! Environment data   
  INTEGER :: Nw
  REAL,DIMENSION(:),ALLOCATABLE :: w
  INTEGER :: Nradiation, Nintegration
  INTEGER :: Nbeta
  REAL,DIMENSION(:),ALLOCATABLE :: beta
  
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: M,Kh
  REAL,DIMENSION(:,:),ALLOCATABLE :: AMOR,RAID
  INTEGER :: i,j,k,l,c, ipb
  REAL :: mass, dw
  
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: CA0,CM0
  COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: Fe0
  COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: RAO
  
  COMPLEX,DIMENSION(:,:),ALLOCATABLE :: ZIGB,ZIGS
  
   INTEGER :: IPrintPOT,NPRINTW1, IDOFPRINT  !! added by rk to print  potentials on body elements
   CHARACTER*50 FMTV,FMT1
    IPrintPOT=1
    NPRINTW1=11
    IDOFPRINT=1
   
  
  !   --- read input datas -------

  CALL ReadTID(ID)
  CALL ReadTMesh(Mesh, TRIM(ID%ID)//'/mesh/')                             ! read mesh data & construct a structural variable MESH 
  CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/Nemoh.cal')        ! Airy potential, velocity, pressure
  CALL Readwbeta(ID,Mesh%Nbodies,Nw,w,Nbeta,beta,Nradiation,Nintegration) ! define w for given w interval input in NEMOH.Cal
  CALL ReadInertia(ID,Mesh%Nbodies,M)                                     ! Read Mass inertia matrix (NdofxNdof) from NEMOH (preProc)
  CALL ReadKh(ID,Mesh%Nbodies,Kh)                                         ! Read stifness matrix from NEMOH (preProc)
  CALL ReadRaidEXT(ID,Mesh%Nbodies,RAID)                                  ! Read additional stif mat i.e mooring (user specified)
  CALL ReadAmorEXT(ID,Mesh%Nbodies,AMOR)                                  ! Read additional ???? check the variable and the data
  
  ! Clean old results !
  CALL SYSTEM('rm '//ID%ID(1:ID%lID)//'/QTF/*.RES 2> /dev/null')
  CALL SYSTEM('rm '//ID%ID(1:ID%lID)//'/QTF/*.wat 2> /dev/null') 
  
  
  dw =w(2)-w(1) ! (w(Nw)-w(1))/(Nw-1) changed by RK
  
  !   ----- Check some limitations due to the current state of the code -------

  IF(Mesh%Nbodies.NE.1) THEN
    WRITE(*,*) "ERROR: QTF working only for 1 body"
  ENDIF

  IF(Nradiation.NE.6) THEN
    WRITE(*,*) "ERROR: QTF working only for 1 body with convential 6 DOF"
  ENDIF
  
  IF(Nintegration.NE.6) THEN
    WRITE(*,*) "ERROR: QTF working only for 1 body with convential 6 DOF"
  ENDIF
  
  IF(dw.NE.w(1)) THEN
    WRITE(*,*) "ERROR: QTF working only if wmin = dw"
  ENDIF

  !     ---------- Print summary --------------

  WRITE(*,*) ' '
  WRITE(*,*) ' Summary of calculation'
  WRITE(*,*) ' '
  WRITE(*,*) ' Working directory: ', ID%ID(1:ID%lID)
  WRITE(*,*) ' '
  IF (Environment%Depth.GT.0.) THEN
      WRITE(*,'(A,F7.2,A)') '  -> Water depth = ',Environment%Depth,' m'
  ELSE
      WRITE(*,'(A)') '  -> Infinite water depth'
  END IF
  WRITE(*,'(A,I5,A,F7.4,A,F7.4,A,F7.4)') '  ->',Nw,' wave frequencies from ',w(1),' to ',w(Nw), " with freq step = ", dw
  WRITE(*,'(A,I5)') '  -> Number of faces on body =',Mesh%Npanels
  WRITE(*,*) 'Press enter to continue'
  READ(*,*)
  !!!!!!!!!  Geometrical preProcessing !!!!!!!!!!!!!!!
  CALL QTFGEOM(ID,Mesh)       !  it is updating the global variables

  !!!!!!!!! Initialization of quantities for all frequencies !!!!!!!!!
  
  !init some global variables
  G = Environment%G
  XEFF = Environment%XEFF      ! coord of point where incident wave is measured
  YEFF = Environment%YEFF      ! the coord is userdefined in NEMOH.cal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !ZEFF = -71.6 ! it was always 0, this value should be read in NEMOH.cal !Centre of Gravity
  ZEFF = 0 ! it was always 0, this value should be read in NEMOH.cal !Centre of Gravity
               ! this is treated as COG for the QTF Calculation 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NC = 1 ! 1 body
  NCO = NC  ! number of body
  IPOS(1) = 0 ! TODO: check this value.
  IMXC(1) = IMX ! TODO: check this value.  IMX is number of panels on body mesh
  AIND(1) = beta(1)  ! wave direction (only one direction is available)
  NIN = 1
  WRITE(*,*) "beta = ",beta(1)
  !Write FA.RES
  CALL WriteFa(ID,Environment,Nw,w,Nbeta,beta,M,Kh,RAID,AMOR)
  
  ! Read added mass, damping and excitation force
  CALL QTFReadFirstOrderLoad(ID, Nw,Nradiation,Nbeta,Nintegration, CA0, CM0, Fe0)

  ! Read (precise) green functions
  ! CALL QTFReadGrin(ID)                       !it is updating global variable of green functions
    CALL CREK(ID) !preparing second term green function
  ! Read motion RAO
  CALL QTFReadMotion(ID,Nw,Nintegration,Nbeta,RAO)
  
  ALLOCATE(ZIGB(Mesh%Npanels,Nradiation+Nbeta))
  ALLOCATE(ZIGS(Mesh%Npanels,Nradiation+Nbeta))


  IF (IPrintPOT==1) THEN
        WRITE(FMT1,'(I2.2)') NPRINTW1
        WRITE(FMTV,'(I1)') IDOFPRINT
        OPEN(188,FILE=ID%ID(1:ID%lID)//'/QTF/PLOTBODYPOT_W'//TRIM(FMT1)//'_DOF'//TRIM(FMTV)//'.dat')   
        WRITE(188,*) 'IDpanel','XM     ','YM     ','ZM  ','RE(PHI1)  ','IM(PHI1)     ','PHI1_X  .._Y  .._Z'
        WRITE(188,*) 'IDpanel','XM     ','YM     ','ZM  ','RE(PSIM)  ', 'IM(PSI)     ','PSIM_X  .._Y  ..Z'
  ENDIF

  
  !!!!!!!!!  LOOP ON FREQUENCIES !!!!!!!!!!!!!!!
  DO i=1,Nw
  
    !Reformat data for QTF preProcessor
    do k=1,Nradiation
      do l=1,Nradiation
	CM(k,l,1,1) = CM0(i,k,l)
	CA(k,l,1,1) = CA0(i,k,l)
      end do
    end do
    
    do j=1,Nbeta
      do k =1,Nintegration
	ZF1(k,1,j) = Fe0(i,k,j)
	!NOTE: Motion ZA should depend on wave direction
	ZA(k,1) = RAO(i,k,j) 
      end do
    end do

    CALL QTFReadSources(ID,i,Nw,Nradiation,Nbeta,Mesh%Npanels, ZIGB,ZIGS)
    
!     WRITE(*,*) "---sing.wat for w = ----", w(i)
      
    CALL QTFWriteSing(ID,i,Nw,Nradiation,Nbeta,Mesh%Npanels, ZIGB,ZIGS)
    
    ! Call the preProcessor for this frequency
    CALL QTFInit(ID,DPI/w(i),beta(1),i,Environment%Depth,IPrintPOT,NPRINTW1,IDOFPRINT)
    
  END DO
    CLOSE(188)
    DEALLOCATE(ZIGB)
    DEALLOCATE(ZIGS)
    DEALLOCATE(w)
    DEALLOCATE(beta)
    DEALLOCATE(M)
    DEALLOCATE(Kh)
    DEALLOCATE(CA0)
    DEALLOCATE(CM0)
    DEALLOCATE(Fe0)
    DEALLOCATE(RAO)
    

END PROGRAM Main
