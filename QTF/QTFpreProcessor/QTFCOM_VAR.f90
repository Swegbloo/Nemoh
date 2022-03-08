!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!   - Gerard Delhommeau
!--------------------------------------------------------------------------------------

  MODULE QTFCOM_VAR

    ! This data is initialized by the Main function 
    ! And passed along to the QTFinit function that is doind the actual job of the QTF preproc
  
    IMPLICIT NONE
    
    INTEGER,PARAMETER :: LE = 6! Output in terminal
    INTEGER,PARAMETER :: LN = 1        ! Maximum number of body (1 for now)
    INTEGER,PARAMETER :: NI = 1        ! Maximum number of wave direction (1 for now)
    INTEGER,PARAMETER :: NPT = 10000   ! Maximum number of points in body mesh (counting contour) TODO: make dynamic
    INTEGER,PARAMETER :: NFA = 5000   ! Maximum number of panels in body mesh (counting contour) TODO: make dynamic
    REAL, PARAMETER :: PI4 = ATAN(1.0), PI = 4.0*PI4, DPI = 2.0*PI,QPI = 4.0*PI
    
    !Number of bodies!
    INTEGER :: NC,NCO
    
    ! Environment data !
    REAL :: G,XEFF,YEFF,ZEFF
    INTEGER :: NIN          !Number of incidence = 1
    
    ! Mesh data !
    REAL:: X(NPT),Y(NPT),Z(NPT) ! Points of body mesh + contour
    REAL :: TDIS(NFA) ! Caracteristic distance of panels  ! ,DIST(NFA)     AC: only TDIS is used 
    INTEGER :: M1(NFA),M2(NFA),M3(NFA),M4(NFA) ! body mesh connectivities
    REAL :: P(NFA),Q(NFA),R(NFA),AIRE(NFA) ! P,Q,R: unit mesh normals   AIRE: surface of panels
    REAL:: XM(2*NFA),YM(2*NFA),ZM(2*NFA) !XM, YM, ZM: center of panels
    REAL :: ALF(NFA)
    INTEGER :: IND(NFA) ! ALF = length of the contour lines   IND = panel to which the contour line belongs
    
    INTEGER :: NP, IMX, NSYMY !Number of points on body mesh, number of panels on body mesh, symmetry Y
    INTEGER :: IXX, NFFL  ! Total number of panels (counting the body+contour) on the body and number of lines on the contour
    
    ! Note: IXX = IMX + NFFL
    
    REAL :: ZER1       ! Geometrical tolerance
    
    ! 1st order load
    REAL :: CM(6,6,LN,LN),CA(6,6,LN,LN)      ! Added mass and damping for 1 frequency and 1 body
    COMPLEX :: ZF1(6,LN,NI)   ! Excitation force (complex) for 1 body and 1 wave direction

!   Data containing green functions (GRIN.QAT)
    REAL :: XR(700),XZ(130),APD1X(700,130),APD1Z(700,130), APD2X(700,130),APD2Z(700,130)

!   Complex RAO of motion
    COMPLEX :: ZA(6,LN)
    
!   Sources ZIGB, ZIGS computed by Nemoh 1st order are not defined here. They are passed to QTFinit through the file "sing.wat"
      
    ! For QTFInit and QTFBaseFunctions
    INTEGER,PARAMETER :: NEXR=31,NPIN=101!in NEMOH1 251
    INTEGER :: IPOS(LN),IMXC(LN),KK(5),AIND(NI)
    REAL :: FSP(NFA),FSM(NFA)
    REAL :: VSXP(NFA),VSYP(NFA),VSZP(NFA)
    REAL :: VSXM(NFA),VSYM(NFA),VSZM(NFA)
    REAL :: FS1(NFA,2),FS2(NFA,2)                        
    REAL :: VSX1(NFA,2),VSY1(NFA,2),VSZ1(NFA,2)                                
    REAL :: VSX2(NFA,2),VSY2(NFA,2),VSZ2(NFA,2)                                
    REAL :: AMBDA(NEXR),AR(NEXR) 
    REAL :: XJAC(16,NFA),XGA(16,NFA),YGA(16,NFA),ZGA(16,NFA)
    INTEGER :: NG
    REAL :: CQ(NPIN),QQ(NPIN)
    REAL :: VSXP1(NFA),VSYP1(NFA),VSZP1(NFA),VSXM1(NFA),VSYM1(NFA),VSZM1(NFA),VSXP2(NFA),VSYP2(NFA),VSZP2(NFA),VSXM2(NFA),VSYM2(NFA),VSZM2(NFA),S1B(NFA),S1S(NFA),S2B(NFA),S2S(NFA),SP1(NFA),SM1(NFA),SP2(NFA),SM2(NFA)                                      
    
    
  END MODULE QTFCOM_VAR
