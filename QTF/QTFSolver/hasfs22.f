!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!    - Fabien Robaux (EDF/INNOSEA)
!    - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)    
!    - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!   version 2022
!   R. Kurnia
!   HASFS for calculating the potentials on the free-surface mesh for each input wave-frequency
!   User mesh input is now only the free-surface mesh, the discretized
!   contour is not needed anymore, the Green-fromulation is not used in
!   Hasfscalc.f90
!   User-mesh should be named as SF_L12_2.dat and placed in the mesh
!   folder
!--------------------------------------------------------------------------------------

 
      PROGRAM SFLBO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!PROGRAMME SFLBO 								  !!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! CE PROGRAMME CALCULE LES POTENTIELS ET LEURS DERIVEES EN LES POINTS DE LA !!!!!!!!!
    !!!!! SURFACE LIBRE. LE MAILLAGE UTILISEE EST ENREGISTRE DANS                   !!!!!!!!!
    !!!!! LE FICHIER SF_L12.dat. LES POTENTIELS SONT CALCULES SUR LES POINTS        !!!!!!!!!
    !!!!! CENTRAUX DES FACETTES DE CE MAILLAGE ET SONT ENREGISTRES DANS             !!!!!!!!!
    !!!!! POTENTIELS_SL.RES, DANS LEQUEL NOUS TRAVAILLONS SOUS FORME                !!!!!!!!!
    !!!!! D'ENREGISTREMENT.                                                         !!!!!!!!!
    !!!!! IL VIENT APRES DUOK ET HASBO MEME SI IL EN EST TOUT A FAIT INDEPENDANT    !!!!!!!!!
    !!!!! IL VIENT AVANT HASFSCALC QUI, LUI, EFFECTUE LE CALCUL A PARTIR DES        !!!!!!!!!
    !!!!! QUANTITES CALCULEES ICI                                                   !!!!!!!!!
    !!!!!                                                                           !!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    
      USE PARA
      COMMON X(NPT),Y(NPT),Z(NPT),DIST(NFA),M1(NFA),M2(NFA)
      COMMON M3(NFA),M4(NFA),P(NFA),Q(NFA),R(NFA),AIRE(NFA),TDIS(NFA)
      COMMON XM(2*NFA),YM(2*NFA),ZM(2*NFA),IND(NFA)
      COMMON IPOS(LN),IMXC(LN),KK(5),AIND(NI)
      COMMON FSP(NFA),FSM(NFA)
      COMMON SP1(NFA),SM1(NFA),SP2(NFA),SM2(NFA)
      COMMON VSXP(NFA),VSYP(NFA),VSZP(NFA)
      COMMON VSXM(NFA),VSYM(NFA),VSZM(NFA)
      COMMON VSXP1(NFA),VSYP1(NFA),VSZP1(NFA)
      COMMON VSXM1(NFA),VSYM1(NFA),VSZM1(NFA)
      COMMON VSXP2(NFA),VSYP2(NFA),VSZP2(NFA)
      COMMON VSXM2(NFA),VSYM2(NFA),VSZM2(NFA)
      COMMON FS1(NFA,2),FS2(NFA,2)
      COMMON VSX1(NFA,2),VSY1(NFA,2),VSZ1(NFA,2)
      COMMON VSX2(NFA,2),VSY2(NFA,2),VSZ2(NFA,2)
      COMMON AMBDA(NEXR),AR(NEXR),XL(5),YL(5),ZL(5)
      COMMON RR(5),DRX(5),DRY(5),DRZ(5)
      COMMON/CST/NCS,NCO,NSYMY,NP,IMX,IXX,NIN,XEFF,YEFF,ZEFF,
     &PI4,DPI,QPI,CQ(NPIN),QQ(NPIN)
      COMMON/SLM/XMSL(NFASL),YMSL(NFASL),NTOT
      COMMON/JAC/XJAC(16,NFA),XGA(16,NFA),YGA(16,NFA),ZGA(16,NFA),NG
      COMMON/FIC/XR(700),XZ(130),APD1X(700,130),APD1Z(700,130),
     &APD2X(700,130),APD2Z(700,130)
      CHARACTER*10 ID,IDEN
      CHARACTER*5 FMT9
      INTEGER :: lID
      COMMON/CAR/IDEN
      DIMENSION XSL(NFASL),YSL(NFASL)
      DIMENSION M1SL(NFASL),M2SL(NFASL),M3SL(NFASL),M4SL(NFASL)
      COMPLEX ZSB(NFA),ZSS(NFA)
      COMPLEX ZPB,ZPS,ZVXB,ZVXS,ZVYB,ZVYS,ZVZB,ZVZS
      REAL,DIMENSION(6,6,LN):: AIN,ASH
      REAL,DIMENSION(6*LN,6*LN):: AMOR,RAID
      REAL,DIMENSION(2,NFA):: DNSRADR,DNSRADM  !DENSITE DE SOURCE RADIE
      REAL,DIMENSION(2,NFA):: DNSPERR,DNSPERM  !DENSITE DE SOURCE PERTURBE
      DIMENSION LK(5),CN(2*NFA,6),CNSL(3,2*NFASL),NT(NPIN),
     & AIRESL(NFASL)
      LOGICAL :: FICHEY
      REAL :: PTPROCHE(2)
      CHARACTER(20) :: NOMF4,file8
      ChARACTER(21) :: file9
      CHARACTER(16) :: NOMB1,NOMB2,NOMB3,NOMB4,NOMB5
      CHARACTER(16) :: NOMH0,NOMH1,NOMH2,NOMH3,NOMH4,NOMH5
      COMPLEX S1(6),S2(6),S3(6),S4(6),S5(6),S6(6),S7(6),
     &E1(6),E2(6),E3(6),E4(6)
      COMPLEX,DIMENSION(:,:),ALLOCATABLE ::POTPERSLB,POTINCSLB,POTTOTSLB
      COMPLEX,DIMENSION(:,:),ALLOCATABLE ::POTPERSLS,POTINCSLS,POTTOTSLS
      COMPLEX,DIMENSION(:,:,:), ALLOCATABLE ::POTRADSLB,POTRADSLS
      COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::VPERSLB,VPERSLS
      COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::VINCSLB,VINCSLS
      COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::VTOTSLB,VTOTSLS
      COMPLEX,DIMENSION(:,:,:,:),ALLOCATABLE::VRADSLB,VRADSLS
      REAL::S1BPER(NFA),S1SPER(NFA),S2BPER(NFA),S2SPER(NFA)
      REAL::S1BRAD(NFA), FMT,FMT1
      REAL::RSL,BEXP,RT,NPASR0,RCEREXT0,RIND
      REAL::X1,X2,X3,X4,Y1,Y2,Y3,Y4,H1,H2,G1,G2,ZSL
      INTEGER::NTRAN,IOP,UNIT9,ITERTEMP
      COMPLEX::PSL
      LOGICAL::dir_e  !to check existing files or directories
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      OPEN(UNIT=7,FILE="ID.dat")
      READ(7,*) lID
      READ(7,*) ID
      CLOSE(7)
    
    



    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
    !!!!!     OUVERTURE ET ENREGISTREMENT DE GEOMETRIE COMPLETE        !!!!!!
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
      NOMF4=ID(1:lID)//'/QTF/FA.RES'
      OPEN(UNIT=12,FILE=NOMF4,ACCESS='DIRECT',
     #STATUS='UNKNOWN',RECL=4*4*NFA)
      READ(12,REC=1)NC,NCO,NSYMY,NP,IMX,IXX,XEFF,YEFF,ZEFF,
     1(IMXC(I),I=1,NCO),(IPOS(I),I=1,NC),ZER,T,IMIN,H,AMZ,ILIN,
     #RHO,VA,BETA1,W1,AM0,NIN,(AIND(I),I=1,NIN)
      JQ=6*NC
      I1=IMX+1
      READ(12,REC=2)(X(I),I=1,NP)
      READ(12,REC=3)(Y(I),I=1,NP)
      READ(12,REC=4)(Z(I),I=1,NP)
      READ(12,REC=5)(M1(I),I=1,IMX),(M2(I),I=1,IMX),
     #(M3(I),I=1,IMX),(M4(I),I=1,IMX)
      READ(12,REC=6)(CN(I,1),I=1,IMX),(CN(I,2),I=1,IMX),
     #(CN(I,3),I=1,IMX)
      READ(12,REC=7)(XM(I),I=1,IMX),(YM(I),I=1,IMX),(ZM(I),I=1,IMX)
      READ(12,REC=8)(AIRE(I),I=1,IMX),(TDIS(I),I=1,IMX),
     #(DIST(I),I=1,IMX)
      IXMX=IXX-IMX
      READ(12,REC=9)(XM(I),I=I1,IXX),(YM(I),I=I1,IXX),
     #(IND(I),I=1,IXMX)
      READ(12,REC=10)(((AIN(I,J,K),I=1,6),J=1,6),
     #((ASH(I,J,K),I=1,6),J=1,6),K=1,NC),((RAID(I,J),I=1,JQ),J=1,JQ),
     #((AMOR(I,J),I=1,JQ),J=1,JQ)


      NFAC=IMX
      NFFL=IXX-NFAC!NOMBRE DE LIGNE DE LA SURFACE DE FLOTAISON
      WRITE(*,*)'NFAC=',NFAC,'NFFL=',NFFL,'IXX=',IXX
      NL=NFAC+NFFL!NL=IXX
      XCDG=XEFF
      YCDG=YEFF
      ZCDG=ZEFF
      PI4=ATAN(1.)
      DPI=2.*PI
      QPI=4.*PI
      NJ=NSYMY+1
      NG=1
      DO I=1,IMX
        CALL SOMGO(X,Y,Z,M1,M2,M3,M4,I,XJAC,XGA,YGA,ZGA,NG)
        IF(NG.EQ.1) THEN
          XGA(1,I)=XM(I)
          YGA(1,I)=YM(I)
          ZGA(1,I)=ZM(I)
        ENDIF
        DO L=1,NG
          XJAC(L,I)=XJAC(L,I)/AIRE(I)
        END DO
      END DO
      PRINT *,'NG = ',NG
    !
    ! ---- FICHIERS DES QUANTITES DU PREMIER ORDRE  -----
      NOMB1=ID(1:lID)//'/QTF/B1.RES'
      NOMB2=ID(1:lID)//'/QTF/B2.RES'
      NOMB3=ID(1:lID)//'/QTF/B3.RES'
      NOMB4=ID(1:lID)//'/QTF/B4.RES'
    ! ---- FICHIERS PROBLEMES ADDITIONNELS DE RADIATION  -----
      NOMH1=ID(1:lID)//'/QTF/H1.RES'
      NOMH2=ID(1:lID)//'/QTF/H2.RES'
      NOMH3=ID(1:lID)//'/QTF/H3.RES'
      NOMH4=ID(1:lID)//'/QTF/H4.RES'
      NR1=3+24*NC+72*NC**2+NFFL*2*NJ
      NR3=NL*2*NJ

      OPEN(21,FILE=NOMB1,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR1)
      OPEN(22,FILE=NOMB2,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
      OPEN(23,FILE=NOMB3,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
      OPEN(24,FILE=NOMB4,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)

      NOMB5=ID(1:lID)//'/QTF/AH.RES'
      OPEN(UNIT=40,FILE=NOMB5,STATUS='UNKNOWN')


c~       PRINT*,'LONGUEUR D''ADIMENSIONNALISATION ?'
c~         READ*,XLAD
      XLAD=1
      AD=RHO*G*XLAD
    ! SOMMATION SUR 1 SEUL CORPS
      KNC=1
    !
    ! PROCEDURE VALABLE SI W(NT(I)-NT(I-1))=W(I)-W(I-1)
c~       PRINT *,'NOMBRE DE PULSATIONS ETUDIEES'
      READ(LL,*)N
c~       WRITE(LE,*)"ENTRER LES NUMEROS D''ENREGISTREMENTS DANS L'ORDRE"
c~       WRITE(LE,*)"DES PULSATIONS CROISSANTES (PERIODES DECROISSANTES)"
      READ(LL,*)(NT(I),I=1,N)
      READ(LL,*) LQTFP
      
!!!!!!!!!!!!!! version matrice carre pour les QTF+!!!!!!!!!!!!!!!!!!!!!!
!~        IF (LHASK>1 .AND. LQTFP==1) THEN                               
! SI LES QTF+ SONT CALCULES, IL FAUT REDUIRE LE NOMBRE DE PULSATION
!~              NHASKIND=N/2
!~        ELSE
!~              NHASKIND=N
!~        ENDIF
!~        PRINT*, 'N =',N,'NHASKIND =',NHASKIND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
    !!!!!               IMPORTATION DES FONCTIONS DE GREEN             !!!!!!
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
      OPEN(UNIT=44,FILE=ID(1:lID)//'/QTF/GRIN.QAT',FORM='UNFORMATTED',
     # STATUS='OLD')
      READ(44)IR,JZ,(XR(I),I=1,IR),(XZ(J),J=1,JZ)
      DO J=1,JZ
        READ(44)(APD1X(I,J),I=1,IR),(APD1Z(I,J),I=1,IR),
     #(APD2X(I,J),I=1,IR),(APD2Z(I,J),I=1,IR)
      END DO
      CLOSE(UNIT=44)
    !!!!!!!!!!!!!! FIN D'IMPORTATION DES FONCTIONS DE GREEN    !!!!!!!!!!!!!!
      NOMH0=ID(1:lID)//'/QTF/H0.RES'
      NR8=(2*IMX*NJ+3)*2
      OPEN(30,FILE=NOMH0,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR8)
      NR5=(NL*2*NJ+3)
      NR3=IXX*2*NJ
      NOMH5=ID(1:lID)//'/QTF/H5.RES'
      OPEN(25,FILE=NOMH5,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR8)
    !!!!IMPORTATION DES DENSITES DE SOURCES!

    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
    !!!!!               MAILLAGE DE LA SURFACE LIBRE                   !!!!!!
    !!!!!               MAILLAGE RADIAL COMPRIS ENTRE 10 ET 40         !!!!!!
    !!!!!               NOMBRE DE FACES:NFACSL                         !!!!!!
    !!!!!               POINTS CENTRAUX:XMSL(KKK),YMSL(KKK)              !!!!!!
    !!!!!               AIRE DE LA FACE: AIRESL(KK)                    !!!!!!
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
      OPEN(10,FILE=ID(1:lID)//'/mesh/SF_L12_2.dat',STATUS='OLD')
      READ(10,*) NSYMY,NOMBPO,NFACSL,RCEREXT,NPASR
      WRITE(*,*) "NSYMY =",NSYMY,"NFACSL =",NFACSL,"NOMBPO =",NOMBPO
      DO IOP=1,NOMBPO
        READ(10,*) ITERTEMP,XSL(IOP),YSL(IOP),ZSL
      ENDDO
      READ(10,*)
      DO IOP=1,NFACSL
        READ(10,*) M1SL(IOP),M2SL(IOP),M3SL(IOP),M4SL(IOP)
        X1=XSL(M1SL(IOP))
        Y1=YSL(M1SL(IOP))
        X2=XSL(M2SL(IOP))
        Y2=YSL(M2SL(IOP))
        X3=XSL(M3SL(IOP))
        Y3=YSL(M3SL(IOP))
        X4=XSL(M4SL(IOP))
        Y4=YSL(M4SL(IOP))
        XMSL(IOP)=0.25*(X1+X2+X3+X4)
        YMSL(IOP)=0.25*(Y1+Y2+Y3+Y4)
        H1=sqrt((X3-X2)**2+(Y3-Y2)**2)
        G1=sqrt((X2-X1)**2+(Y2-Y1)**2)
        H2=sqrt((X4-X1)**2+(Y4-Y1)**2)
        G2=sqrt((X4-X3)**2+(Y4-Y3)**2)
        AIRESL(IOP)=0.5*(H1+H2)*0.5*(G1+G2)
      !  WRITE(*,*) 'X1:X4=',X1,X2,X3,X4
      !  WRITE(*,*) 'Y1:Y4=',Y1,Y2,Y3,Y4
      !  WRITE(*,*) '(',IOP,') XM,YM=',XMSL(IOP),YMSL(IOP)
      !  CNSL(1,IOP)=0
      !  CNSL(2,IOP)=0
      ENDDO
      CLOSE(10) 
        NTOT=NFACSL
        NCONT=0
        NLIGNE=0
        WRITE(*,*) 'RCEREXT=', RCEREXT
    
      NR9=2*4*2*2*NFASL !MAX((NFACSL+NCONT),NOMBPO)
      OPEN(11,FILE='SF_L12.RES',STATUS='UNKNOWN',
     &ACCESS='DIRECT',RECL=NR9)
      WRITE(11,REC=1)NSYMY,NFACSL,NCONT,NOMBPO,NLIGNE,RCEREXT,NPASR
      PRINT *,RCEREXT
      WRITE(11,REC=2)(XSL(I),YSL(I),I=1,NOMBPO)
      WRITE(11,REC=3)(XMSL(I),YMSL(I),I=1,NFACSL+NCONT+NLIGNE)
      WRITE(11,REC=4)(AIRESL(I),I=1,NFACSL+NCONT)
      WRITE(11,REC=5)(CNSL(1,I),CNSL(2,I),I=NFACSL+1,NFACSL+NCONT)
      WRITE(11,REC=6)(M1SL(I),M2SL(I),I=1,NFACSL+NCONT)
      WRITE(11,REC=7)(M3SL(I),M4SL(I),I=1,NFACSL)
      CLOSE(11)


    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
    !!!!!               FIN MAILLAGE                                    !!!!!
    !!!!!!!!!!!!!                                            !!!!!!!!!!!!!!!!!
    
    
    
    
    
    
    
    !!!!!!!!!!!!!                                            !!!!!!!!!!!!!!!!!
    !!!!!!!         CALCUL, POUR CHAQUE PULSATION                    !!!!!!!!!
    !!!!!!!         PULSATION: N1, DOF: DOF,                         !!!!!!!!!
    !!!!!!!         FACE DE SL: IFACSL                               !!!!!!!!!
    !!!!!!!         POPERSL(N1,IFACSL):POTENTIEL PERTURBE            !!!!!!!!!
    !!!!!!!         PORADSL(IJ,N1,IFACSL):POTENTIEL RADIE            !!!!!!!!!
    !!!!!!!         VPERSLX,Y,Z(N1,IFACSL):VITESSE PERTURBE          !!!!!!!!!
    !!!!!!!         VRADSLX,Y,Z(IJ,N1,LL):VITESSE RADIE,DU AU MVT IJ !!!!!!!!!
    !!!!!!!!!!!!!                                            !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! version matrice carre pour les QTF+!!!!!!!!!!!!!!!!!!!!!!
c~       ALLOCATE(POTPERSLB(NHASKIND,NTOT))
c~       ALLOCATE(POTPERSLS(NHASKIND,NTOT))
c~       ALLOCATE(POTINCSLB(NHASKIND,NTOT))
c~       ALLOCATE(POTINCSLS(NHASKIND,NTOT))
c~       ALLOCATE(POTTOTSLB(NHASKIND,NTOT))
c~       ALLOCATE(POTTOTSLS(NHASKIND,NTOT))
c~       ALLOCATE(POTRADSLB(6,N,NTOT))       ! LES POTENTIELS DE RADIATION SONT PRIS AUX PULSATION W1+-W2
c~       ALLOCATE(POTRADSLS(6,N,NTOT))       ! LES POTENTIELS DE RADIATION SONT PRIS AUX PULSATION W1+-W2
c~       ALLOCATE(VPERSLB(3,NHASKIND,NTOT))
c~       ALLOCATE(VINCSLB(3,NHASKIND,NTOT))
c~       ALLOCATE(VINCSLS(3,NHASKIND,NTOT))
c~       ALLOCATE(VTOTSLB(3,NHASKIND,NTOT))
c~       ALLOCATE(VPERSLS(3,NHASKIND,NTOT))
c~       ALLOCATE(VTOTSLS(3,NHASKIND,NTOT))
c~       ALLOCATE(VRADSLB(3,6,N,NTOT))       ! LES VITESSES DE RADIATION SONT PRISES AUX PULSATION W1+-W2
c~       ALLOCATE(VRADSLS(3,6,N,NTOT))       ! LES VITESSES DE RADIATION SONT PRISES AUX PULSATION W1+-W2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs) !!!!
      ALLOCATE(POTPERSLB(N,NTOT))
      ALLOCATE(POTPERSLS(N,NTOT))
      ALLOCATE(POTINCSLB(N,NTOT))
      ALLOCATE(POTINCSLS(N,NTOT))
      ALLOCATE(POTTOTSLB(N,NTOT))
      ALLOCATE(POTTOTSLS(N,NTOT))
      ALLOCATE(POTRADSLB(6,N,NTOT))       ! LES POTENTIELS DE RADIATION SONT PRIS AUX PULSATION W1+-W2
      ALLOCATE(POTRADSLS(6,N,NTOT))       ! LES POTENTIELS DE RADIATION SONT PRIS AUX PULSATION W1+-W2
      ALLOCATE(VPERSLB(3,N,NTOT))
      ALLOCATE(VINCSLB(3,N,NTOT))
      ALLOCATE(VINCSLS(3,N,NTOT))
      ALLOCATE(VTOTSLB(3,N,NTOT))
      ALLOCATE(VPERSLS(3,N,NTOT))
      ALLOCATE(VTOTSLS(3,N,NTOT))
      ALLOCATE(VRADSLB(3,6,N,NTOT))       ! LES VITESSES DE RADIATION SONT PRISES AUX PULSATION W1+-W2
      ALLOCATE(VRADSLS(3,6,N,NTOT))       ! LES VITESSES DE RADIATION SONT PRISES AUX PULSATION W1+-W2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c~       K=-1
      file8=ID(1:lID)//'/QTF/SLIINF.wat'
      ! check if dir exists
      inquire( file=file8, exist=dir_e )
      if ( dir_e ) then
        CALL SYSTEM('rm '//file8)
      endif
      NB8=8*IMX*2
      OPEN(UNIT=8,FILE=file8,FORM='UNFORMATTED',STATUS='NEW',
     & ACCESS='DIRECT',RECL=4*NB8)
     
      IF(H.LE.0..OR.H.GT.0.9E20) THEN
        CALL VAV(1,H,ZER)
        H=1.E20
      ELSE
        CALL VAV(2,H,ZER)
      ENDIF
     
!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(POTPERSLB,POTPERSLS,POTINCSLB,
!$OMP& POTINCSLS,POTTOTSLB,POTTOTSLS,POTRADSLB,POTRADSLS,VPERSLB,
!$OMP& VINCSLB,VINCSLS,VTOTSLB,VPERSLS,VTOTSLS,VRADSLB,VRADSLS,
!$OMP& H,IMIN,RHO,BB,BC,BD,ZER,IMX,NTOT,NHASKIND,NJ,DPI,XMSL,XEFF,
!$OMP& YMSL,YEFF,NSYMY,lID,ID,NT)
      DO IJK=1,N
        UNIT9=IJK+4000
        WRITE(FMT9,'(I0.5)') IJK
        file9=ID(1:lID)//'/QTF/sliper'//FMT9//'.wat'
        ! check if dir exists
        inquire( file=file9, exist=dir_e )
        if ( dir_e ) then
          CALL SYSTEM('rm '//file9)
        endif
        OPEN(UNIT=UNIT9,FILE=file9,FORM='UNFORMATTED',STATUS='UNKNOWN')
        PRINT *,'PULSATION=',IJK
        N1=NT(IJK)
        READ(21,REC=N1)TR1,BETA1,AM1
        IF(H.GT.0.9E20) THEN
          CALL VNS(TR1,BETA1,VA,H,IMIN,IZ,RHO,BB,BC,BD,ZER,UNIT9)
        ELSE
          CALL VNSF(TR1,BETA1,VA,H,IMIN,IZ,RHO,BB,BC,BD,ZER,UNIT9)
        ENDIF
  
        DO IJ=1,6
          NIJ=(IJK-1)*6+IJ
          READ(30,REC=NIJ)TR,BETA,AK0,
     & ((DNSRADR(JJ,J),J=1,IMX),JJ=1,NJ),
     & ((DNSRADM(JJ,J),J=1,IMX),JJ=1,NJ)
          DO J=1,IMX
            ZSB(J)=CMPLX(DNSRADR(1,J),DNSRADM(1,J))
            ZSS(J)=CMPLX(DNSRADR(2,J),DNSRADM(2,J))
          END DO
          REWIND UNIT9
          DO IK=1,NTOT
            READ(UNIT9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX),
     # (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
            READ(UNIT9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX),
     # (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
            READ(UNIT9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX),
     # (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)
            READ(UNIT9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX),
     # (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
            ZPB=(0.,0.)
            ZPS=(0.,0.)
            ZVXB=(0.,0.)
            ZVXS=(0.,0.)
            ZVYB=(0.,0.)
            ZVYS=(0.,0.)
            ZVZB=(0.,0.)
            ZVZS=(0.,0.)
            DO J=1,IMX
              ZPB=ZPB+0.5*(ZSB(J)*CMPLX(SP1(J)+SM1(J),
     # SP2(J)+SM2(J))+ZSS(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
              ZPS=ZPS+0.5*(ZSS(J)*CMPLX(SP1(J)+SM1(J),
     # SP2(J)+SM2(J))+ZSB(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
              ZVXB=ZVXB+0.5*(ZSB(J)*CMPLX(VSXP1(J)+
     # VSXM1(J),VSXP2(J)+VSXM2(J))+ZSS(J)*CMPLX(VSXP1(J)-VSXM1(J),
     # VSXP2(J)-VSXM2(J)))
              ZVXS=ZVXS+0.5*(ZSS(J)*CMPLX(VSXP1(J)+
     # VSXM1(J),VSXP2(J)+VSXM2(J))+ZSB(J)*CMPLX(VSXP1(J)-VSXM1(J),
     # VSXP2(J)-VSXM2(J)))
              ZVYB=ZVYB+0.5*(ZSB(J)*CMPLX(VSYP1(J)+
     # VSYM1(J),VSYP2(J)+VSYM2(J))+ZSS(J)*CMPLX(VSYP1(J)-VSYM1(J),
     # VSYP2(J)-VSYM2(J)))
              ZVYS=ZVYS-0.5*(ZSS(J)*CMPLX(VSYP1(J)+
     # VSYM1(J),VSYP2(J)+VSYM2(J))+ZSB(J)*CMPLX(VSYP1(J)-VSYM1(J),
     # VSYP2(J)-VSYM2(J)))
              ZVZB=ZVZB+0.5*(ZSB(J)*CMPLX(VSZP1(J)+
     # VSZM1(J),VSZP2(J)+VSZM2(J))+ZSS(J)*CMPLX(VSZP1(J)-VSZM1(J),
     # VSZP2(J)-VSZM2(J)))
              ZVZS=ZVZS+0.5*(ZSS(J)*CMPLX(VSZP1(J)+
     # VSZM1(J),VSZP2(J)+VSZM2(J))+ZSB(J)*CMPLX(VSZP1(J)-VSZM1(J),
     # VSZP2(J)-VSZM2(J)))
            END DO
            POTRADSLB(IJ,IJK,IK)=ZPB
            POTRADSLS(IJ,IJK,IK)=ZPS
            VRADSLB(1,IJ,IJK,IK)=ZVXB
            VRADSLS(1,IJ,IJK,IK)=ZVXS
            VRADSLB(2,IJ,IJK,IK)=ZVYB
            VRADSLS(2,IJ,IJK,IK)=ZVYS
            VRADSLB(3,IJ,IJK,IK)=ZVZB
            VRADSLS(3,IJ,IJK,IK)=ZVZS
          END DO
        END DO
  
        REWIND UNIT9
c~         IF (IJK .LE. NHASKIND) THEN                                  ! version matrice carre pour les QTF+
        IF (IJK .LE. N) THEN                                            ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
        READ(25,REC=IJK)TR,BETA,AK0,
     & ((DNSPERR(JJ,J),J=1,IMX),JJ=1,NJ),
     & ((DNSPERM(JJ,J),J=1,IMX),JJ=1,NJ)
        DO J=1,IMX
          ZSB(J)=CMPLX(DNSPERR(1,J),DNSPERM(1,J))
          ZSS(J)=CMPLX(DNSPERR(2,J),DNSPERM(2,J))
        END DO
  
        DO IK=1,NTOT
          READ(UNIT9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX),
     # (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
          READ(UNIT9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX),
     # (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
          READ(UNIT9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX),
     # (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)
          READ(UNIT9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX),
     # (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
          ZPB=(0.,0.)
          ZPS=(0.,0.)
          ZVXB=(0.,0.)
          ZVXS=(0.,0.)
          ZVYB=(0.,0.)
          ZVYS=(0.,0.)
          ZVZB=(0.,0.)
          ZVZS=(0.,0.)
          DO J=1,IMX
            ZPB=ZPB+0.5*(ZSB(J)*CMPLX(SP1(J)+SM1(J),
     # SP2(J)+SM2(J))+ZSS(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
            ZPS=ZPS+0.5*(ZSS(J)*CMPLX(SP1(J)+SM1(J),
     # SP2(J)+SM2(J))+ZSB(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
            ZVXB=ZVXB+0.5*(ZSB(J)*CMPLX(VSXP1(J)+
     # VSXM1(J),VSXP2(J)+VSXM2(J))+ZSS(J)*CMPLX(VSXP1(J)-VSXM1(J),
     # VSXP2(J)-VSXM2(J)))
            ZVXS=ZVXS+0.5*(ZSS(J)*CMPLX(VSXP1(J)+
     # VSXM1(J),VSXP2(J)+VSXM2(J))+ZSB(J)*CMPLX(VSXP1(J)-VSXM1(J),
     # VSXP2(J)-VSXM2(J)))
            ZVYB=ZVYB+0.5*(ZSB(J)*CMPLX(VSYP1(J)+
     # VSYM1(J),VSYP2(J)+VSYM2(J))+ZSS(J)*CMPLX(VSYP1(J)-VSYM1(J),
     # VSYP2(J)-VSYM2(J)))
            ZVYS=ZVYS-0.5*(ZSS(J)*CMPLX(VSYP1(J)+
     # VSYM1(J),VSYP2(J)+VSYM2(J))+ZSB(J)*CMPLX(VSYP1(J)-VSYM1(J),
     # VSYP2(J)-VSYM2(J)))
            ZVZB=ZVZB+0.5*(ZSB(J)*CMPLX(VSZP1(J)+
     # VSZM1(J),VSZP2(J)+VSZM2(J))+ZSS(J)*CMPLX(VSZP1(J)-VSZM1(J),
     # VSZP2(J)-VSZM2(J)))
            ZVZS=ZVZS+0.5*(ZSS(J)*CMPLX(VSZP1(J)+
     # VSZM1(J),VSZP2(J)+VSZM2(J))+ZSB(J)*CMPLX(VSZP1(J)-VSZM1(J),
     # VSZP2(J)-VSZM2(J)))
          END DO
          POTPERSLB(IJK,IK)=ZPB
          POTPERSLS(IJK,IK)=ZPS
          VPERSLB(1,IJK,IK)=ZVXB
          VPERSLS(1,IJK,IK)=ZVXS
          VPERSLB(2,IJK,IK)=ZVYB
          VPERSLS(2,IJK,IK)=ZVYS
          VPERSLB(3,IJK,IK)=ZVZB
          VPERSLS(3,IJK,IK)=ZVZS
        END DO
        REWIND UNIT9
               
        WH=DPI/TR1
        AK0=WH**2/G
        IF(H.LT.1.E10) THEN
          AM0=X0(AK0*H)/H
        ELSE
          AM0=AK0
        ENDIF
        AKAD=-G*AM0/WH
        C1=-G/WH
        CB=COS(BETA1)
        SB=SIN(BETA1)
        DO IK=1,NTOT
          ZPO=0.
          Q1=CIH(AM0,ZPO,H)*AKAD
          Q2=SIH(AM0,ZPO,H)*AKAD
          Q3=AM0*((XMSL(IK)-XEFF)*CB+(YMSL(IK)-YEFF)*SB)
          VXB1=-Q1*SIN(Q3)*CB
          VYB1=-Q1*SIN(Q3)*SB
          VZB1=+Q2*COS(Q3)
          VXB2=+Q1*COS(Q3)*CB
          VYB2=+Q1*COS(Q3)*SB
          VZB2=+Q2*SIN(Q3)
          Q5=C1*CIH(AM0,ZPO,H)
          POTB1=+Q5*COS(Q3)
          POTB2=+Q5*SIN(Q3)
          POTINCSLB(IJK,IK)=ZI*CMPLX(POTB1,POTB2) ! comments RK: added ZI to be same convention as in the NEMOH First order
          VINCSLB(1,IJK,IK)=ZI*CMPLX(VXB1,VXB2)
          VINCSLB(2,IJK,IK)=ZI*CMPLX(VYB1,VYB2)
          VINCSLB(3,IJK,IK)=ZI*CMPLX(VZB1,VZB2)
          IF (NSYMY.EQ.1) THEN ! SI NSYMY=1, YEFF=0., SIGNE INDIFFERENT
            Q3=AM0*((XMSL(IK)-XEFF)*CB-(YMSL(IK)+YEFF)*SB)
            VXS1=-Q1*SIN(Q3)*CB
            VYS1=+Q1*SIN(Q3)*SB
            VZS1=+Q2*COS(Q3)
            VXS2=+Q1*COS(Q3)*CB
            VYS2=-Q1*COS(Q3)*SB
            VZS2=+Q2*SIN(Q3)
            Q5=C1*CIH(AM0,ZPO,H)
            POTS1=+Q5*COS(Q3)
            POTS2=+Q5*SIN(Q3)
            POTINCSLS(IJK,IK)=ZI*CMPLX(POTS1,POTS2) ! added ZI by RK
            VINCSLS(1,IJK,IK)=ZI*CMPLX(VXS1,VXS2)
            VINCSLS(2,IJK,IK)=ZI*CMPLX(VYS1,VYS2)
            VINCSLS(3,IJK,IK)=ZI*CMPLX(VZS1,VZS2)
          ENDIF
        END DO
        POTTOTSLB=POTINCSLB+POTPERSLB
        VTOTSLB=VINCSLB+VPERSLB
        POTTOTSLS=POTINCSLS+POTPERSLS
        VTOTSLS=VINCSLS+VPERSLS
        ENDIF
        
      ! check if dir exists
      inquire( file=file9, exist=dir_e )
      if ( dir_e ) then
        CALL SYSTEM('rm '//file9)
      endif
        
      END DO
!$OMP END PARALLEL DO
      ! check if dir exists
      inquire( file=file8, exist=dir_e )
      if ( dir_e ) then
        CALL SYSTEM('rm '//file8)
      endif
      NR10=NTOT*2*2*4
      OPEN(12,FILE='POTENTIELS_SL.RES',ACCESS='DIRECT',RECL=2*NR10)
      DO IJK=1,N
        IREC=(IJK-1)*36
c~         IF (IJK .LE. NHASKIND) THEN                                  ! version matrice carre pour les QTF+
        IF (IJK .LE. N) THEN                                            ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
          WRITE(12,REC=IREC+1)(POTTOTSLB(IJK,IFACSL),IFACSL=1,NTOT),
     & (POTTOTSLS(IJK,IFACSL),IFACSL=1,NTOT)
          WRITE(12,REC=IREC+2)(POTINCSLB(IJK,IFACSL),IFACSL=1,NTOT),
     & (POTINCSLS(IJK,IFACSL),IFACSL=1,NTOT)
          WRITE(12,REC=IREC+3)(POTPERSLB(IJK,IFACSL),IFACSL=1,NTOT),
     & (POTPERSLS(IJK,IFACSL),IFACSL=1,NTOT)
        ENDIF
        DO IK=1,6
          WRITE(12,REC=IREC+3+IK)(POTRADSLB(IK,IJK,IFACSL),IFACSL=1,
     & NTOT),(POTRADSLS(IK,IJK,IFACSL),IFACSL=1,NTOT)
        END DO
        DO IDOF=1,3
c~           IF (IJK .LE. NHASKIND) THEN                                ! version matrice carre pour les QTF+
          IF (IJK .LE. N) THEN                                          ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
            WRITE(12,REC=IREC+10+9*(IDOF-1))
     & (VTOTSLB(IDOF,IJK,IFACSL),IFACSL=1,NTOT),
     & (VTOTSLS(IDOF,IJK,IFACSL),IFACSL=1,NTOT)

            WRITE(12,REC=IREC+11+9*(IDOF-1))
     & (VINCSLB(IDOF,IJK,IFACSL),IFACSL=1,NTOT),
     & (VINCSLS(IDOF,IJK,IFACSL),IFACSL=1,NTOT)

            WRITE(12,REC=IREC+12+9*(IDOF-1))
     & (VPERSLB(IDOF,IJK,IFACSL),IFACSL=1,NTOT),
     & (VPERSLS(IDOF,IJK,IFACSL),IFACSL=1,NTOT)
          ENDIF
          DO IK=1,6
            WRITE(12,REC=IREC+12+9*(IDOF-1)+IK)
     & (VRADSLB(IDOF,IK,IJK,IFACSL),IFACSL=1,NTOT),
     & (VRADSLS(IDOF,IK,IJK,IFACSL),IFACSL=1,NTOT)
          END DO
        END DO
      END DO
      CLOSE(12)
125   FORMAT('PULSATION:',2X,1I2.1,1X,'/',1X,1I2.1,4X,
     & 'POURCENTAGE:',1X,1I3.1,1X,'%')
        
      DEALLOCATE(POTPERSLB)
      DEALLOCATE(POTPERSLS)
      DEALLOCATE(POTINCSLB)
      DEALLOCATE(POTINCSLS)
      DEALLOCATE(POTTOTSLB)
      DEALLOCATE(POTTOTSLS)
      DEALLOCATE(POTRADSLB)
      DEALLOCATE(POTRADSLS)
      DEALLOCATE(VPERSLB)
      DEALLOCATE(VPERSLS)
      DEALLOCATE(VINCSLB)
      DEALLOCATE(VINCSLS)
      DEALLOCATE(VTOTSLB)
      DEALLOCATE(VTOTSLS)
      DEALLOCATE(VRADSLB)
      DEALLOCATE(VRADSLS)
      CLOSE(25)
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!
    !!!!!           FIN IMPORTATION POTENTIELS/VITESSES                 !!!!!
    !!!!!!!!!!!!!                                           !!!!!!!!!!!!!!!!!

      STOP
      END
