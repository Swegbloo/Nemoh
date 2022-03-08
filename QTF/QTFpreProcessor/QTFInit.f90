!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991 PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------

MODULE MQTFInit

!    .ECRITURE DES QUANTITES DU PREMIER ORDRE 
!    .CALCUL ET ECRITURE DE CES QUANTITES AU NIVEAU DE LA FLOTTAISON
!    POUR LE CALCUL DES EFFORTS DE 2ND ORDRE
!
!    .CE PROGRAMME UTILISE LES MOUVEMENTS
!    .CE PROGRAMME UTILISE LE FICHIER GRIN.QAT
!
!    .WRITING THE QUANTITIES OF THE FIRST ORDER
!     CALCULATION AND WRITING OF THESE QUANTITIES AT FLOTATION LEVEL
!     FOR CALCULATING 2ND ORDER EFFORTS
!
!     .THIS PROGRAM USES MOVEMENTS
!     .THIS PROGRAM USES THE GRIN.QAT FILE
!

CONTAINS

  SUBROUTINE QTFInit(ID,T,BETA,IMIN,H)
      
      
    !BETA in rad  
      
    USE QTFCOM_VAR
    USE MIdentification
    USE MQTFBaseFunctions
    
    IMPLICIT NONE
    
    TYPE(TID) :: ID
    
    REAL :: VXR(2,NFA),VXM(2,NFA),VYR(2,NFA),VYM(2,NFA) ! velocities (real+imag)
    REAL :: VZR(2,NFA),VZM(2,NFA),VGR(2,NFA),VGM(2,NFA) ! velocities (real+imag)
    REAL :: POR(2,NFA),POM(2,NFA)
    REAL :: DNSR(2,NFA),DNSM(2,NFA)
    COMPLEX :: ZPB(NFA),ZPS(NFA),ZVXB(NFA),ZVXS(NFA),ZTGB(NFA),ZTGS(NFA)
    COMPLEX :: ZVYB(NFA),ZVYS(NFA),ZVZB(NFA),ZVZS(NFA)
    CHARACTER*16 NOMF1,NOMF2,NOMF3,NOMF4,NOMF5,NOMF6
    CHARACTER*16 NOMB1,NOMB2,NOMB3,NOMB4
    CHARACTER*16 NOMH0,NOMH1,NOMH2,NOMH3,NOMH4,NOMH5
    CHARACTER*5 :: strIMIN

    COMPLEX :: ZIGB(NFA,7), ZIGS(NFA,7)
    
    ! maybe still to clean up
    REAL :: DRAD, T, BETA,VA,WH,AK0,WR,TR,ZMAX,ZMIN,TIRAN,BETAD,C1,AKAD,CB,SB,HH,CR
    REAL :: H,AM0,AKK,ZER,CCIH,CP,COEFB,CIB,SIB,COEFS,CIS,SIS
    COMPLEX :: ZI ! the imaginary number i*i=-1
    INTEGER :: JQ,I1,IINC, IMIN, NR,NJJ,N0,NM,NSM,KS,KKK,I,IQ,K,L,NR8,JJ,J,NR1,NR3,II,IJ,NIJ,NR5,NR2,N00
    
    ! -----> INITIALISATION
    ZI=(0.,1.)
    DRAD=PI/180.
    VA=0
    I1=IMX+1
    IINC=1
    
    NR=IMIN
    WH=DPI/T
    AK0=WH**2/G
    
    IF(H.LT.0.9E20)THEN
	AM0=X0(AK0*H)/H     
    ELSE
	AM0=AK0     
    ENDIF
    
    AKK=(8.*ATAN(1.)/T)**2/G            !wave number

    ! PULSATION DE RENCONTRE
    WR=WH-VA*AM0*COS(BETA)
    TR=DPI/WR
    ! SYMETRIE
    NJJ=NSYMY+1                         !NJJ=2 for symmetric, else NJJ=1
    
    ! CALCUL DES POINTS DE GAUSS ET POIDS SUR LES FACETTES
    !CALCULATION OF GAUSS POINTS AND WEIGHT ON THE FACETS
    NG=1
    DO I=1,IMX
	CALL SOMGO(X,Y,Z,M1,M2,M3,M4,I,XJAC,XGA,YGA,ZGA,NG)
	IF(NG.EQ.1)THEN
	    XGA(1,I)=XM(I)
	    YGA(1,I)=YM(I)
	    ZGA(1,I)=ZM(I)
	ENDIF
	DO L=1,NG
	    XJAC(L,I)=XJAC(L,I)/AIRE(I)
	END DO
    END DO
    
    ! CALCUL DU TIRANT DEAU ET DU ZERO DU CALCUL
    ! CALCULATION OF THE DRAFT AND ZERO OF THE CALCULATION
    ZMAX=0.
    DO I=1,IMX
	ZMAX=MIN(ZMAX,ZM(I))
    END DO
    ZMAX=ABS(ZMAX)
    ZMIN=-ZMAX
    ZMIN=ABS(ZMIN)
    TIRAN=0.
    DO I=1,IMX
	IF(ZM(I).NE.0)THEN
	    TIRAN=MAX(TIRAN,ABS(ZM(I))) ! draft
	ENDIF
    END DO
    PRINT *,'TIRAN = ',TIRAN
    ZER=-0.001*TIRAN
    WRITE(*,'(E14.5)') ZER 
    ! OUVERTURE DES FICHIERS DES COEFFICIENTS D'INFLUENCE POTENTIEL ET VITESSES
    OPEN(UNIT=8,FILE=ID%ID(1:ID%lID)//'/QTF/QTFinf.wat', FORM='UNFORMATTED',STATUS='UNKNOWN')
    OPEN(UNIT=9,FILE=ID%ID(1:ID%lID)//'/QTF/QTFper.wat', FORM='UNFORMATTED',STATUS='UNKNOWN')
    ! WRITE(strIMIN,'(I5)') IMIN
    ! OPEN(99,FILE=ID%ID(1:ID%lID)//'/QTF/QTFper'//strIMIN//'.dat') ! checking the coef compared with in NEMOH1
    ! CALCUL DES COEFFICIENTS D'INFLUENCE POTENTIEL ET VITESSES
    ! CALCULATION OF POTENTIAL INFLUENCE COEFFICIENTS AND SPEEDS    
    IF(H.LE.0.OR.H.GT.0.9E20)THEN       ! deep water case
	H=0.
	CALL VAV(1,0.,ZER)              !S^0_ij, K^0_ij
	WRITE(LE,7800)T
	7800 FORMAT(/5X,' T = ',F15.5,' SECONDES   PROFONDEUR INFINIE'/)
	CALL VNS(T,BETA,VA,ZER)         !S_ij=S^0_ij+S^1_ij+S^2_ij, K_ij=K^0_ij+K^1_ij+K^2_ij
    ELSE                                !Finite depth case
	CALL VAV(2,H,ZER)               
	WRITE(LE,7801)T,H                                                         
	7801 FORMAT(/5X,' T = ',F15.5,' SECONDES   PROFONDEUR = ',F15.5,' METRES'/)                                                               
	CALL VNSF(T,BETA,VA,H,ZER)
    ENDIF                                                                  
    ! UNE SEULE INCIDENCE AVEC CETTE VERSION  
    ! on direction incidence wave?
    BETAD=BETA*180./PI
    WRITE(LE,1112)IINC,BETAD
    1112 FORMAT(1X,'INCIDENCE N0 ',I3,' : ANGLE DE LA HOULE AVEC OX = ', F8.3,' DEGRES')
    
    N0=IMIN
    IF(N0.LE.0)STOP
    WRITE(LE,990)N0,T,WR,BETAD,IMX,NP,NSYMY,VA
    990  FORMAT(//1X,'NUMERO D''ENREGISTREMENT N0 = ',I4/         1X,'           PERIODE DE HOULE         T = ',F10.6,' SECONDES'/1X,'           PULSATION (DE RENCONTRE) W = ',F10.6, ' HZ'/ 1X, '           DIRECTION DE LA HOULE BETA = ',F10.6,' DEGRES'/11X,I4 ,' FACETTES',4X,I4,' POINTS',4X,I1,' SYMETRIE'/,12X,'VITESSE D''AVANCE VA=',F10.6,' M/S'//)
    
    ! LECTURE DES SINGULARITES TOTALES DE DIFFRACTION ET DE RADIATION
    NM=6*NC                             !NC 
    NSM=NM+NIN
    OPEN(UNIT=10,FILE=ID%ID(1:ID%lID)//'/QTF/sing.wat',ACCESS='DIRECT', STATUS='UNKNOWN',RECL=4*4*NFA)
    DO KS=1,NSM
	READ(10,REC=KS)(ZIGB(I,KS),I=1,IMX),(ZIGS(I,KS),I=1,IMX)
    END DO
    ! ZIG(I,KS)	: SOURCE DE COURANT DE LA Ieme FACETTE CORRESPONDANT            Current source of the ith panel
    !		   KS=1,6*NC : AU PROBLEME DE RADIATION SUR LE IQeme DOF        Radiation problem on the IQeme DOF
    !			       ET LE Keme CORPS                                 and the Keme body
    !		   KS=NM+1,NM+NIN : AU KKKeme PROBLEME DE RADIATION             at KKKeme radiation problem
    !				    (NIN DIRECTIONS DE HOULE)                   NIN wave direction
    ! ZA	: RAO (DOF,NUMERO DU CORPS)                                     RAO(DOF,Nuber of body)
    ! ZTG(I) 	: SOURCE DE COURANT TOTALE A LA FACETTE I                       ZTG(Total source on the panel I)
    ! B/S --> body / symetrique                                                 body/symmetric
    DO KKK=1,NIN
    DO I=1,IMX
	ZTGB(I)=ZIGB(I,NM+KKK)
	ZTGS(I)=ZIGS(I,NM+KKK)
	DO K=1,NC                               !NC is number of body
        DO IQ=1,6                               !ID DoF
	    KS=IQ+(K-1)*6
	    ZTGB(I)=ZTGB(I)-ZI*WR*ZIGB(I,KS)*ZA(IQ,K)
	    ZTGS(I)=ZTGS(I)-ZI*WR*ZIGS(I,KS)*ZA(IQ,K)
           ! write(*,'(E15.4,E15.4)') ZTGB(I),ZTGS(I)
	END DO
	END DO
    END DO
    END DO

    CLOSE(10)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CC
    !                                                                     !
    ! STOCKAGE DES QUANTITES DU 1ER ORDRE                                 !
    ! POUR LE CALCUL DE L'INTEGRALE D'HASKIND SUR LA SURFACE LIBRE        !
    !                                                                     ! 
    ! STORAGE OF 1ST ORDER QUANTITIES                                     !
    ! FOR THE CALCULATION OF THE INTEGRAL OF HASKIND ON THE FREE SURFACE  !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CC
    ! -----> ECRITURE DES DENSITES TOTALES DE SOURCE
    ! ! ! ! H5 VERIFIE
    DO I=1,IMX
	DNSR(1,I)=REAL(ZTGB(I))
	DNSM(1,I)=AIMAG(ZTGB(I))
	DNSR(2,I)=REAL(ZTGS(I))
	DNSM(2,I)=AIMAG(ZTGS(I))
    END DO
    NR8=(2*IMX*NJJ+3)*2
    NOMH5=ID%ID(1:ID%lID)//'/QTF/H5.RES'
    OPEN(25,FILE=NOMH5,ACCESS='DIRECT', STATUS='UNKNOWN',RECL=4*NR8)
    WRITE(25,REC=N0)T,BETA,AK0, ((DNSR(JJ,J),J=1,IMX),JJ=1,NJJ),((DNSM(JJ,J),J=1,IMX),JJ=1,NJJ)
    CLOSE(25)
    !
    !
    ! CALCUL DES POTENTIELS ET VITESSES SUR LE CORPS ET A LA FLOTTAISON
    ! CALCULATION OF POTENTIALS AND SPEEDS ON THE BODY AND AT FLOTATION
    C1=-G/WH
    AKAD=C1*AM0
    KKK=NIN
    BETA=AIND(KKK)*DRAD
    CB=COS(BETA)
    SB=SIN(BETA)
    HH=H
    IF(H.EQ.0.)HH=1.E20
    REWIND 9
    DO I=1,IXX
	ZPB(I)=(0.,0.)
	ZPS(I)=(0.,0.)
	ZVXB(I)=(0.,0.)
	ZVYB(I)=(0.,0.)
	ZVZB(I)=(0.,0.)
	ZVXS(I)=(0.,0.)
	ZVYS(I)=(0.,0.)
	ZVZS(I)=(0.,0.)
	READ(9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX), (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
	READ(9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX), (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
	READ(9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX), (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)
	READ(9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX), (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
	IF(ZM(I).LT.ZER.OR.I.GT.IMX)THEN
	    CR=SH(AM0,ZM(I),HH)*AKAD
	    CCIH=CH(AM0,ZM(I),HH)
	    CP=CCIH*AKAD
	    COEFB=AM0*((XM(I)-XEFF)*CB+(YM(I)-YEFF)*SB)
	    CIB=COS(COEFB)
	    SIB=SIN(COEFB)
	    ZVXB(I)=CP*CB*CMPLX(-SIB,CIB)
	    ZVYB(I)=CP*SB*CMPLX(-SIB,CIB)
	    ZVZB(I)=CR*CMPLX(CIB,SIB)
	    ZPB(I)=C1*CCIH*CMPLX(CIB,SIB)
	    IF(NSYMY.EQ.1)THEN
		COEFS=AM0*((XM(I)-XEFF)*CB-(YM(I)-YEFF)*SB)
		CIS=COS(COEFS)
		SIS=SIN(COEFS)
		ZVXS(I)=CP*CB*CMPLX(-SIS,CIS)
		ZVYS(I)=CP*SB*CMPLX(-SIS,CIS)
		ZVZS(I)=CR*CMPLX(CIS,SIS)
		ZPS(I)=C1*CCIH*CMPLX(CIS,SIS)
	    ENDIF
	    DO J=1,IMX
		! ATTENTION, 1 ET 2 NE FONT PAS REFERENCE AUX PULSATION 1 ET 2
		ZPB(I)=ZPB(I)+0.5*(ZTGB(J)*CMPLX(SP1(J)+SM1(J),SP2(J)+SM2(J))+ZTGS(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
		ZPS(I)=ZPS(I)+0.5*(ZTGS(J)*CMPLX(SP1(J)+SM1(J),SP2(J)+SM2(J))+ZTGB(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
             !  write(*,'(E15.4,E15.4,E15.4,E15.4)') SP1(J),SM1(J),SP2(J),SM2(J)
             !   write(*,'(E15.4,E15.4,E15.4,E15.4)') VSXP1(J),VSXM1(J),VSXP2(J),VSXM2(J)
             !   write(*,'(E15.4,E15.4,E15.4,E15.4)') VSYP1(J),VSYM1(J),VSYP2(J),VSYM2(J)
		ZVXB(I)=ZVXB(I)+0.5*(ZTGB(J)*CMPLX(VSXP1(J)+VSXM1(J), VSXP2(J)+VSXM2(J))+ZTGS(J)*CMPLX(VSXP1(J)-VSXM1(J), VSXP2(J)-VSXM2(J)))
		ZVXS(I)=ZVXS(I)+0.5*(ZTGS(J)*CMPLX(VSXP1(J)+VSXM1(J), VSXP2(J)+VSXM2(J))+ZTGB(J)*CMPLX(VSXP1(J)-VSXM1(J), VSXP2(J)-VSXM2(J)))
		ZVYB(I)=ZVYB(I)+0.5*(ZTGB(J)*CMPLX(VSYP1(J)+VSYM1(J), VSYP2(J)+VSYM2(J))+ZTGS(J)*CMPLX(VSYP1(J)-VSYM1(J), VSYP2(J)-VSYM2(J)))
		ZVYS(I)=ZVYS(I)-0.5*(ZTGS(J)*CMPLX(VSYP1(J)+VSYM1(J), VSYP2(J)+VSYM2(J))+ZTGB(J)*CMPLX(VSYP1(J)-VSYM1(J), VSYP2(J)-VSYM2(J)))
		ZVZB(I)=ZVZB(I)+0.5*(ZTGB(J)*CMPLX(VSZP1(J)+VSZM1(J), VSZP2(J)+VSZM2(J))+ZTGS(J)*CMPLX(VSZP1(J)-VSZM1(J), VSZP2(J)-VSZM2(J)))
		ZVZS(I)=ZVZS(I)+0.5*(ZTGS(J)*CMPLX(VSZP1(J)+VSZM1(J), VSZP2(J)+VSZM2(J))+ZTGB(J)*CMPLX(VSZP1(J)-VSZM1(J), VSZP2(J)-VSZM2(J)))
	    END DO
	ELSE
	    ZVXB(I)=(0.,0.)
	    ZVYB(I)=(0.,0.)
	    ZVZB(I)=(0.,0.)
	    ZPB(I)=(0.,0.)
	    ZVXS(I)=(0.,0.)
	    ZVYS(I)=(0.,0.)
	    ZVZS(I)=(0.,0.)
	    ZPS(I)=(0.,0.)
	ENDIF
    END DO
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CC
    !                                                                     !
    ! STOCKAGE DES QUANTITES DU 1ER ORDRE                                 !
    ! POUR LE CALCUL DES EFFORTS BICHROMATIQUES (BICK)                    !
    !                                                                     !
    ! STORAGE OF 1ST ORDER QUANTITIES                                     !
    ! FOR THE CALCULATION OF BICHROMATIC EFFORTS (BICK)                   ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CC
    !
    ! -----> CALCUL DES VITESSES SUR LE CORPS ET A LA FLOTTAISON
    ! -----> CALCULATION OF SPEEDS ON THE BODY AND AT FLOTATION 
    !
    ! CALCUL DES VITESSES SUR LES CORPS 
    ! CALCULATION OF SPEEDS ON THE BODIES  
    ! --------------------------------
    DO I=1,IXX
	VXR(1,I)=REAL(ZVXB(I))
	VXR(2,I)=REAL(ZVXS(I))
	VXM(1,I)=AIMAG(ZVXB(I))
	VXM(2,I)=AIMAG(ZVXS(I))
	VYR(1,I)=REAL(ZVYB(I))
	VYR(2,I)=REAL(ZVYS(I))
	VYM(1,I)=AIMAG(ZVYB(I))
	VYM(2,I)=AIMAG(ZVYS(I))
	VZR(1,I)=REAL(ZVZB(I))
	VZR(2,I)=REAL(ZVZS(I))
	VZM(1,I)=AIMAG(ZVZB(I))
	VZM(2,I)=AIMAG(ZVZS(I))
	!       VGR(1,I)=REAL(ZPB(I))
	!       VGR(2,I)=REAL(ZPS(I))
	!       VGM(1,I)=AIMAG(ZPB(I))
	!       VGM(2,I)=AIMAG(ZPS(I))
    END DO
    !
    ! CALCUL DES HAUTEURS A LA FLOTTAISON
    ! CALCULATION On WATERLINE
    ! -----------------------------------
    DO I=1,NFFL
	VGR(1,I)=REAL(ZPB(I+IMX))
	VGR(2,I)=REAL(ZPS(I+IMX))
	VGM(1,I)=AIMAG(ZPB(I+IMX))
	VGM(2,I)=AIMAG(ZPS(I+IMX))
    END DO
    NR1=3+24*NC+72*NC**2+NFFL*2*NJJ
    NR3=IXX*2*NJJ
    NOMB1=ID%ID(1:ID%lID)//'/QTF/B1.RES'
    NOMB2=ID%ID(1:ID%lID)//'/QTF/B2.RES'
    NOMB3=ID%ID(1:ID%lID)//'/QTF/B3.RES'
    NOMB4=ID%ID(1:ID%lID)//'/QTF/B4.RES'
    OPEN(21,FILE=NOMB1,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR1)
    OPEN(22,FILE=NOMB2,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
    OPEN(23,FILE=NOMB3,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
    OPEN(24,FILE=NOMB4,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
    print *,'n0 = ', N0
    WRITE(21,REC=N0)TR,BETA,AM0,((ZF1(L,II,IINC),L=1,6),II=1,NC),((ZA(IQ,J),IQ=1,6),J=1,NC),((((CM(I,J,K,L),J=1,6),L=1,NC), ((CA(I,J,K,L),J=1,6),L=1,NC),I=1,6),K=1,NC),((VGR(JJ,J),J=1,NFFL),JJ=1,NJJ),((VGM(JJ,J),J=1,NFFL),JJ=1,NJJ)
    WRITE(22,REC=N0) ((VXR(JJ,J),J=1,IXX),JJ=1,NJJ),((VXM(JJ,J),J=1,IXX),JJ=1,NJJ)
    WRITE(23,REC=N0) ((VYR(JJ,J),J=1,IXX),JJ=1,NJJ),((VYM(JJ,J),J=1,IXX),JJ=1,NJJ)
    WRITE(24,REC=N0) ((VZR(JJ,J),J=1,IXX),JJ=1,NJJ),((VZM(JJ,J),J=1,IXX),JJ=1,NJJ)
    CLOSE(UNIT=21)
    CLOSE(UNIT=22)
    CLOSE(UNIT=23)
    CLOSE(UNIT=24)   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CC
    !                                                                     !
    ! STOCKAGE DES QUANTITES DU 1ER ORDRE                                 !
    ! POUR LE CALCUL DE L'INTEGRALE D'HASKIND SUR LE CORPS                !
    !                                                                     !
    !                                                                     !
    ! STORAGE OF 1ST ORDER QUANTITIES                                     !
    ! FOR THE CALCULATION OF THE INTEGRAL OF HASKIND ON THE BODY          !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!CC
    !
    ! --------------------------------------------------------------
    ! FICHIERS RELATIFS AUX CALCULS D'INTEGRALES D'HASKIND SUR LE CORPS :
    !       POTENTIELS ADDITIONNELS DE RADIATION ET GRADIENTS
    ! FILES RELATING TO THE CALCULATIONS OF INTEGRALS OF HASKIND ON THE BODY:
    ! ADDITIONAL POTENTIAL FOR RADIATION AND GRADIENTS
    !---------------------------------------------------------------                                                               
    !
    NOMH1=ID%ID(1:ID%lID)//'/QTF/H1.RES'
    NOMH2=ID%ID(1:ID%lID)//'/QTF/H2.RES'
    NOMH3=ID%ID(1:ID%lID)//'/QTF/H3.RES'
    NOMH4=ID%ID(1:ID%lID)//'/QTF/H4.RES'
    !
    ! -----> LECTURE DES DENSITES DES SOURCES DE RADIATION COURANTE IJ
    ! -----> ECRITURE DES DENSITES DES SOURCES DE RADIATION DE LA PERIODE NIJ
    !
    ! ! ! H0 VERIFIE
    OPEN(UNIT=10,FILE=ID%ID(1:ID%lID)//'/QTF/sing.wat',ACCESS='DIRECT', STATUS='UNKNOWN',RECL=4*4*NFA)
    NR8=(2*IMX*NJJ+3)*2
    NOMH0=ID%ID(1:ID%lID)//'/QTF/H0.RES'
    OPEN(30,FILE=NOMH0,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR8)
    DO IJ=1,6
	NIJ=(N0-1)*6+IJ
	READ(10,REC=IJ)(DNSR(1,J),DNSM(1,J),J=1,IMX), (DNSR(2,J),DNSM(2,J),J=1,IMX)
	WRITE(30,REC=NIJ)T,BETA,AK0, ((DNSR(JJ,J),J=1,IMX),JJ=1,NJJ),((DNSM(JJ,J),J=1,IMX),JJ=1,NJJ)
    END DO
    CLOSE(UNIT=30)
    !
    NR5=IXX*2*NJJ+3
    OPEN(26,FILE=NOMH1,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR5)
    NR2=IMX*2*NJJ
    OPEN(27,FILE=NOMH2,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)
    OPEN(28,FILE=NOMH3,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)
    OPEN(29,FILE=NOMH4,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)
    REWIND 9
    DO IJ=1,6
	READ(10,REC=IJ)(S1B(I),S2B(I),I=1,IMX), (S1S(I),S2S(I),I=1,IMX)
	! ---- LECTURE INTENSITE SOURCES DE RADIATION DE LA PERIODE
	DO I=1,IMX
	    ZTGB(I)=S1B(I)+ZI*S2B(I)
	    ZTGS(I)=S1S(I)+ZI*S2S(I)
	END DO
	DO I=1,IXX
	    ZPB(I)=(0.,0.)
	    ZVXB(I)=(0.,0.)
	    ZVYB(I)=(0.,0.)
	    ZVZB(I)=(0.,0.)
	    ZPS(I)=(0.,0.)
	    ZVXS(I)=(0.,0.)
	    ZVYS(I)=(0.,0.)
	    ZVZS(I)=(0.,0.)
	    READ(9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX), (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
	    READ(9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX), (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
	    READ(9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX), (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)
	    READ(9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX), (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
	    IF(ZM(I).LT.ZER.OR.I.GT.IMX)THEN
		DO J=1,IMX
		    ZPB(I)=ZPB(I)+0.5*(ZTGB(J)*CMPLX(SP1(J)+SM1(J),SP2(J)+SM2(J))+ZTGS(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
		    ZPS(I)=ZPS(I)+0.5*(ZTGS(J)*CMPLX(SP1(J)+SM1(J),SP2(J)+SM2(J))+ZTGB(J)*CMPLX(SP1(J)-SM1(J),SP2(J)-SM2(J)))
		    ZVXB(I)=ZVXB(I)+0.5*(ZTGB(J)*CMPLX(VSXP1(J)+VSXM1(J), VSXP2(J)+VSXM2(J))+ZTGS(J)*CMPLX(VSXP1(J)-VSXM1(J), VSXP2(J)-VSXM2(J)))
		    ZVXS(I)=ZVXS(I)+0.5*(ZTGS(J)*CMPLX(VSXP1(J)+VSXM1(J), VSXP2(J)+VSXM2(J))+ZTGB(J)*CMPLX(VSXP1(J)-VSXM1(J), VSXP2(J)-VSXM2(J)))
		    ZVYB(I)=ZVYB(I)+0.5*(ZTGB(J)*CMPLX(VSYP1(J)+VSYM1(J), VSYP2(J)+VSYM2(J))+ZTGS(J)*CMPLX(VSYP1(J)-VSYM1(J), VSYP2(J)-VSYM2(J)))
		    ZVYS(I)=ZVYS(I)-0.5*(ZTGS(J)*CMPLX(VSYP1(J)+VSYM1(J), VSYP2(J)+VSYM2(J))+ZTGB(J)*CMPLX(VSYP1(J)-VSYM1(J), VSYP2(J)-VSYM2(J)))
		    ZVZB(I)=ZVZB(I)+0.5*(ZTGB(J)*CMPLX(VSZP1(J)+VSZM1(J), VSZP2(J)+VSZM2(J))+ZTGS(J)*CMPLX(VSZP1(J)-VSZM1(J), VSZP2(J)-VSZM2(J)))
		    ZVZS(I)=ZVZS(I)+0.5*(ZTGS(J)*CMPLX(VSZP1(J)+VSZM1(J), VSZP2(J)+VSZM2(J))+ZTGB(J)*CMPLX(VSZP1(J)-VSZM1(J), VSZP2(J)-VSZM2(J)))
		END DO
	    ELSE
		ZPB(I)=(0.,0.)
		ZPS(I)=(0.,0.)
		ZVXB(I)=(0.,0.)
		ZVYB(I)=(0.,0.)
		ZVZB(I)=(0.,0.)
		ZVXS(I)=(0.,0.)
		ZVYS(I)=(0.,0.)
		ZVZS(I)=(0.,0.)
	    ENDIF
	END DO
	DO I=1,IXX
	    POR(1,I)=REAL(ZPB(I))
	    POR(2,I)=REAL(ZPS(I))
	    POM(1,I)=AIMAG(ZPB(I))
	    POM(2,I)=AIMAG(ZPS(I))
	    IF(I.LE.IMX)THEN
		VXR(1,I)=REAL(ZVXB(I))
		VXR(2,I)=REAL(ZVXS(I))
		VXM(1,I)=AIMAG(ZVXB(I))
		VXM(2,I)=AIMAG(ZVXS(I))
		VYR(1,I)=REAL(ZVYB(I))
		VYR(2,I)=REAL(ZVYS(I))
		VYM(1,I)=AIMAG(ZVYB(I))
		VYM(2,I)=AIMAG(ZVYS(I))
		VZR(1,I)=REAL(ZVZB(I))
		VZR(2,I)=REAL(ZVZS(I))
		VZM(1,I)=AIMAG(ZVZB(I))
		VZM(2,I)=AIMAG(ZVZS(I))
	    ENDIF
	END DO
	REWIND 9
	! -----> ECRITURE
	N00=(N0-1)*6+IJ
	WRITE(26,REC=N00)TR,BETA,AK0, ((POR(JJ,J),J=1,IXX),JJ=1,NJJ),((POM(JJ,J),J=1,IXX),JJ=1,NJJ)
	WRITE(27,REC=N00)((VXR(JJ,J),J=1,IMX),JJ=1,NJJ),                    ((VXM(JJ,J),J=1,IMX),JJ=1,NJJ)
	WRITE(28,REC=N00)((VYR(JJ,J),J=1,IMX),JJ=1,NJJ),                    ((VYM(JJ,J),J=1,IMX),JJ=1,NJJ)
	WRITE(29,REC=N00)((VZR(JJ,J),J=1,IMX),JJ=1,NJJ),                    ((VZM(JJ,J),J=1,IMX),JJ=1,NJJ)
    END DO
    CLOSE(UNIT=8)
    CLOSE(UNIT=9)
    CLOSE(UNIT=10)
    CLOSE(UNIT=26)
    CLOSE(UNIT=27)
    CLOSE(UNIT=28)
    CLOSE(UNIT=29)
   ! CLOSE(99)
  END SUBROUTINE QTFInit
  
END MODULE MQTFInit

