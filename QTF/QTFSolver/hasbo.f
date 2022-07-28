!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!    - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)
!    - F.VILLEGER 
!    - Fabien Robaux (EDF/INNOSEA)
!    - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------
      PROGRAM HASBO
      
      USE PARA      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                            PROGRAM HASBO                               C
C     ****************************************************************** C
C     *                                                                * C
C     *      PROGRAMME    HASBO         ( 1 CORPS, 0 OU 1 SYMETRIE)    * C
C     *                   *****                                        * C
C     *                                                                * C
C     *        LHN NANTES, THESE DE CHEN XIAO-BO  1988                 * C
C     *        SIREHNA       1989                                      * C
C     *                   VERSION DU 10.12.89                          * C
C     ****************************************************************** C
C                                                                        C
C          EFFORTS DU 2D ORDRE PROVENANT DU POTENTIEL DU 2D ORDRE        C
C                     METHODE PROPOSEE PAR B. MOLIN                      C
C                     FROUDE KRYLOV ET HASKIND CARENE                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      COMMON NSYM,NFAC,NFFL,RHO,
     1 XEFF,YEFF,ZEFF,XCDG,YCDG,ZCDG,TPER,BETA,NDEB,NFIN
      COMMON X(NPT),Y(NPT),Z(NPT)
      COMMON M1(NFA2),M2(NFA2),M3(NFA2),M4(NFA2)
      COMMON XM(NFA2),YM(NFA2),ZM(NFA2),AIRE(NFA2),ALF(NFA2)
      COMMON P(NFA2),Q(NFA2),R(NFA2)
      DIMENSION VXR1(4,NFA),VXM1(4,NFA),VYR1(4,NFA),VYM1(4,NFA)
      DIMENSION VZR1(4,NFA),VZM1(4,NFA),VYR2(4,NFA),VYM2(4,NFA)
      DIMENSION VXR2(4,NFA),VXM2(4,NFA),VZR2(4,NFA),VZM2(4,NFA)
      DIMENSION VGR1(4,NFA),VGM1(4,NFA),VGR2(4,NFA),VGM2(4,NFA)
      DIMENSION PRM(4,NFA),PMM(4,NFA),XRM(4,NFA),XMM(4,NFA)
      DIMENSION YRM(4,NFA),YMM(4,NFA),ZRM(4,NFA),ZMM(4,NFA)
      DIMENSION PRP(4,NFA),PMP(4,NFA),XRP(4,NFA),XMP(4,NFA),TABWR2(NPIN)
      DIMENSION YRP(4,NFA),YMP(4,NFA),ZRP(4,NFA),ZMP(4,NFA),TABWR1(NPIN)
      DIMENSION AMSS(6,6),AMRT(6,6),AHYD(6,6),CN(NFA2,6)
      DIMENSION PRO1(6,4),PRO2(6,4),FP(2),FM(2)
      DIMENSION EFWPS(8,12),EFWMN(8,12)
      DIMENSION AMODP(6),AMODM(6),PHASP(6),PHASM(6)
      DIMENSION SOLP(12),SOLM(12),NT(501)
      DIMENSION PRO3(6,4),PRO4(6,4),A1(12),A2(12)
      DIMENSION DNSR(4,NFA),DNSM(4,NFA)
      DIMENSION IPOS(LN),IMXC(LN),LK(5)
      DIMENSION AIND(NI),IND(2*NFA),DIST(NFA),TDIS(NFA)
      DIMENSION AIN(6,6,LN),ASH(6,6,LN),AMOR(6*LN,6*LN),
     $ RAID(6*LN,6*LN)
      DIMENSION CM(6,6,LN,LN),CA(6,6,LN,LN)
      COMPLEX FOR1(6,LN),AMOU(6,LN),AC1(6),AC2(6),VN1,VN2,VX1,VX2,VY1
     $ ,VY2,VZ1,VZ2,XP1,XP2,YP1,YP2,ZP1,ZP2,XMC,YMC,ZMC,GPM1,GPM2,GPM2C,
     $ FMC,FPC,EFWMNC(8,12)
      CHARACTER*1 QPM(2)
      CHARACTER*3 QRI(2)
      CHARACTER*10 ID,IDENB,CCONT(8)
      CHARACTER*60 FCHA
      CHARACTER*16 NOMF1,NOMF2,NOMF3,NOMF4,NOMF5,NOMF6
      CHARACTER*16 NOMB1,NOMB2,NOMB3,NOMB4,NOMB5,NOMB6
      CHARACTER*16 NOMH1,NOMH2,NOMH3,NOMH4,NOMH5
      LOGICAL FICHEY
      INTEGER FICH,FICH2,lID,IUNI,JUNI(2),DOF
      REAL A,B,NAN
      REAL,DIMENSION(1:12,1:NPIN,1:NPIN)::QTFM,QTFP
      REAL,DIMENSION(1:7,1:12,1:NPIN,1:NPIN,1:7)::CONTRIB
      INTEGER :: NC,NCO,NSYM,NP,NFAC,IXX,IMIN,ILIN,NIN,NDEB,NFIN
      REAL :: XEFF,YEFF,ZEFF,T,H,AMZ,BP,RHO,VA,BETA,WR,AM0
      REAL :: XCDG,YCDG,ZCDG
      INTEGER :: NR1,NR3,NJ,KNC,N,ITEST,NN2,NM,N1,N2,N3,N4,IJK,I1,I,J,K
      INTEGER :: NHASKIND,NPL,Louthasbo,OUT1
      REAL :: CB,SB,XL,AD
      REAL ::dnPHIM_R,dnPHIM_I,PHIM_R,PHIM_I
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      NAN=0.
      NAN=NAN/NAN
 
      OPEN(7,FILE="ID.dat")
      READ(7,*) lID
      READ(7,*) ID
      CLOSE(7)
 
 
      NOMF4=ID(1:lID)//'/QTF/FA.RES'
      INQUIRE(FILE=NOMF4,EXIST=FICHEY)
      INQUIRE(FILE=NOMF4,EXIST=FICHEY)
      OPEN(UNIT=12,FILE=NOMF4,ACCESS='DIRECT',
     #STATUS='UNKNOWN',RECL=4*4*NFA)
      READ(12,REC=1)NC,NCO,NSYM,NP,NFAC,IXX,XEFF,YEFF,ZEFF,
     1(IMXC(I),I=1,NCO),(IPOS(I),I=1,NC),BP,T,IMIN,H,AMZ,ILIN,
     # RHO,VA,BETA,WR,AM0,NIN,(AIND(I),I=1,NIN)
      READ(12,REC=2)(X(I),I=1,NP)
      READ(12,REC=3)(Y(I),I=1,NP)
      READ(12,REC=4)(Z(I),I=1,NP)
      READ(12,REC=5)(M1(I),I=1,NFAC),(M2(I),I=1,NFAC),
     # (M3(I),I=1,NFAC),(M4(I),I=1,NFAC)
      READ(12,REC=6)(P(I),I=1,NFAC),(Q(I),I=1,NFAC),(R(I),I=1,NFAC)
      READ(12,REC=7)(XM(I),I=1,NFAC),(YM(I),I=1,NFAC),(ZM(I),I=1,NFAC)
      READ(12,REC=8)(AIRE(I),I=1,NFAC),(TDIS(I),I=1,NFAC),
     # (DIST(I),I=1,NFAC)
      I1=NFAC+1
      READ(12,REC=9)(XM(I),I=I1,IXX),(YM(I),I=I1,IXX),
     #(ALF(I),I=1,IXX-NFAC),(IND(I),I=1,IXX-NFAC) ! ALF edited by RK 
      JQ=6*NC
      READ(12,REC=10)(((AIN(I,J,K),I=1,6),J=1,6),
     # ((ASH(I,J,K),I=1,6),J=1,6),K=1,NC),((RAID(I,J),I=1,JQ),J=1,JQ),
     # ((AMOR(I,J),I=1,JQ),J=1,JQ)
      CLOSE(12)
      NFFL=IXX-NFAC
      WRITE(*,*)'NFAC=',NFAC,'NFFL=',NFFL,'IXX=',IXX
      NL=NFAC+NFFL
      XCDG=XEFF
      YCDG=YEFF
      ZCDG=ZEFF
      NDEB=0
      NFIN=NFAC+NFFL

     
      DO I=1,NFFL   ! DEFINITION DES GRANDEUR A LA FLOTTAISON (N0 HORIZONTALE --> R=0)
        LK(1)=M1(IND(I))
        LK(2)=M2(IND(I))
        LK(3)=M3(IND(I))
        LK(4)=M4(IND(I))
        LK(5)=LK(1)
        P(NFAC+I)=P(IND(I))
        Q(NFAC+I)=Q(IND(I))
        R(NFAC+I)=0.            !edited by RK
C LES QUANTITES N'ONT PAS ETE CALCULEES EN Z=0 A LA FLOTTAISON
C CEPANDANT ECRIRE ZSL=0 PERMET D'ETRE EXACT DANS LA GEOMETRIE (LES CN)
        ZM(IMX+I)=0.
        DO J=1,4
          TEST=(X(LK(J))-X(LK(J+1)))**2+(Y(LK(J))-Y(LK(J+1)))**2
          IF(ABS(Z(LK(J)))+ABS(Z(LK(J+1))).LT.
     1 1.E-05 .AND.TEST.GT.1.E-5) THEN
            AIRE(NFAC+I)=SQRT(TEST)   !  AIRE(I),I>NFAC, EST LA LONGUEUR DES COTES A LA FLOTTAISON
          ENDIF
        END DO
      END DO

      DO I=1,NL     ! VECTEUR NORMALE GENERALISE
        IF(I.LE.NFAC) THEN  ! since P,Q,R are already unit vector so no effect by normalizing again
                  CN(I,1)=P(I)/SQRT(P(I)**2+Q(I)**2+R(I)**2)
                  CN(I,2)=Q(I)/SQRT(P(I)**2+Q(I)**2+R(I)**2)
                  CN(I,3)=R(I)/SQRT(P(I)**2+Q(I)**2+R(I)**2)
                  CN(I,4)=(YM(I)-YEFF)*CN(I,3)-(ZM(I)-ZEFF)*CN(I,2)
                  CN(I,5)=(ZM(I)-ZEFF)*CN(I,1)-(XM(I)-XEFF)*CN(I,3)
                  CN(I,6)=(XM(I)-XEFF)*CN(I,2)-(YM(I)-YEFF)*CN(I,1)
            ELSE  ! FACETTE HORIZONTALE A LA FLOTTAISON
                  CN(I,1)=P(I)/SQRT(P(I)**2+Q(I)**2)
                  CN(I,2)=Q(I)/SQRT(P(I)**2+Q(I)**2)
                  CN(I,3)=0.
                  CN(I,4)=(YM(I)-YEFF)*CN(I,3)-(ZM(I)-ZEFF)*CN(I,2)
                  CN(I,5)=(ZM(I)-ZEFF)*CN(I,1)-(XM(I)-XEFF)*CN(I,3)
                  CN(I,6)=(XM(I)-XEFF)*CN(I,2)-(YM(I)-YEFF)*CN(I,1)
            ENDIF
      END DO

      BETA=BETA/RADD
      CB=COS(BETA)
      SB=SIN(BETA)
      NJ=NSYM+1

C ---- FICHIERS DES QUANTITES DU PREMIER ORDRE  -----
      NOMB1=ID(1:lID)//'/QTF/B1.RES'
      NOMB2=ID(1:lID)//'/QTF/B2.RES'
      NOMB3=ID(1:lID)//'/QTF/B3.RES'
      NOMB4=ID(1:lID)//'/QTF/B4.RES'
      
C ---- FICHIERS PROBLEMES ADDITIONNELS DE RADIATION  -----
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
      NOMB6=ID(1:lID)//'/QTF/AHP.RES'
      OPEN(UNIT=40,FILE=NOMB5,STATUS='UNKNOWN')
      OPEN(UNIT=41,FILE=NOMB6,STATUS='UNKNOWN')

      XL=1        ! LONGUEUR D''ADIMENSIONNALISATION
      AD=RHO*G*XL
      
      
C SOMMATION SUR 1 SEUL CORPS
      KNC=1
C PROCEDURE VALABLE SI W(NT(I)-NT(I-1))=W(I)-W(I-1)
      READ(LL,*)N
            
c~       WRITE(LE,*)' ENTRER LES NUMEROS D''ENREGISTREMENTS DANS L''ORDRE'
c~       WRITE(LE,*)' DES PULSATIONS CROISSANTES (PERIODES DECROISSANTES)'
      READ(LL,*)(NT(I),I=1,N)
      READ(LL,*) LQTFP
      READ(LL,*) 
      READ(LL,*) 
      READ(LL,*) Louthasbo
!!!!!!!!!!!!!! version matrice carre pour les QTF+!!!!!!!!!!!!!!!!!!!!!!
c~       IF (LQTFP==1) THEN 
c SI LES QTF+ SONT CALCULES, IL FAUT REDUIRE LE NOMBRE DE PULSATION                              
c~             NHASKIND=N/2
c~       ELSE
c~             NHASKIND=N
c~       ENDIF
c~       PRINT*, 'N =',N,'NHASKIND =',NHASKIND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c~       PRINT*,'CALCUL AUTOMATIQUE:1, CAS TEST:#1'
c~       READ(*,*)ITEST
      ITEST=1
      IF(ITEST.EQ.1)THEN
C CALCUL AUTOMATIQUE
        NN2=0
      ELSE
C CAS TEST
C ON NE PREND QUE LES NM DERNIERE FREQUENCES DE LA IJKeme SURDIAGONALE
        PRINT*,'NOMBRE DE PERIODES COMBINEES ?'
        READ(*,*)NM
c~         NN2=NHASKIND-NM                                              ! version matrice carre pour les QTF+
        NN2=N-NM                                                        ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C IJK= DIFFERENCE ENTRE N1 ET N2
C BOUCLE PRINCIPALE QUI FINIT TOUT AU BOUT, ON BOUCLE SUR LES DEUX PULSATIONS
c~       DO IJK=0,NHASKIND-1                                            ! version matrice carre pour les QTF+
c~       DO I1=IJK+1+NN2,NHASKIND                                       ! version matrice carre pour les QTF+
      DO IJK=0,N-1                                                      ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
      DO I1=IJK+1+NN2,N                                                 ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
        N1=NT(I1)         !NUMERO D'ENREGISTREMENT PULSATION I1
        N2=NT(I1-IJK)     !NUMERO D'ENREGISTREMENT PULSATION I2
        IF (IJK>0) THEN
          N3=NT(IJK)      !NUMERO D'ENREGISTREMENT PULSATION DELTA I
        ELSE
          N3=0
        ENDIF
        IF (2*I1-IJK.LE.N .AND. LQTFP==1) THEN
          N4=NT(2*I1-IJK) !NUMERO D'ENREGISTREMENT PULSATION SOMME
        ELSE
          N4=0
        ENDIF
        
        DO I=1,8          ! INITIALISATION DES EFFORTS (0)
        DO J=1,12
          EFWPS(I,J)=0.   !ON MET A 0 LES EFFORTS DE SOMMES
          EFWMN(I,J)=0.   !ON MET A 0 LES EFFORTS DE DIFFERENCE
          EFWMNC(I,J)=0.  !ON MET A 0 LES EFFORTS DE DIFFERENCE COMPLEXE
        END DO
        END DO

        IF(N1.NE.0)THEN   ! IMPORTATION DES QUANTITE DU PREMIER ORDRE POUR LA PULSATION 1
        ! TPE       : PERIODE (PROPRE) DE HOULE
        ! BETA      : ANGLE D'INCIDENCE
        ! AK        : NOMBRE D'ONDE
        ! FOR       : FORCE EXCITATRICE
        ! AMOU      : RAO
        ! CM        : MASSE AJOUTEE
        ! CA        : AMORTISSEMENT AJOUTE
        ! VG        : POTENTIEL(1) A LA FLOTTAISON
        ! VX,VY,VZ  : VITESSES(1) DU FLUIDE
          READ(21,REC=N1)TPE1,BETA,AK1,((FOR1(J,II),J=1,6),II=1,NC),
     1 ((AMOU(J,II),J=1,6),II=1,NC),
     1 ((((CM(I,J,K,L),J=1,6),L=1,NC),
     1 ((CA(I,J,K,L),J=1,6),L=1,NC),I=1,6),K=1,NC),
     1 ((VGR1(JJ,J),J=1,NFFL),JJ=1,NJ),((VGM1(JJ,J),J=1,NFFL),JJ=1,NJ)
          READ(22,REC=N1) ((VXR1(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VXM1(JJ,J),J=1,NL),JJ=1,NJ)
          READ(23,REC=N1) ((VYR1(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VYM1(JJ,J),J=1,NL),JJ=1,NJ)
          READ(24,REC=N1) ((VZR1(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VZM1(JJ,J),J=1,NL),JJ=1,NJ)
          W1=2*PI/TPE1
          DO J=1,6
            A1(J)=REAL(AMOU(J,KNC))
            ! A1(J)=-AIMAG(AMOU(J,KNC))
            A1(J+6)=AIMAG(AMOU(J,KNC))
            ! A1(J+6)=REAL(AMOU(J,KNC))
          END DO
        ENDIF

        IF(N2.NE.0)THEN   ! IMPORTATION DES QUANTITE DU PREMIER ORDRE POUR LA PULSATION 2
          READ(21,REC=N2)TPE2,BETA,AK2,((FOR1(J,II),J=1,6),II=1,NC),
     1((AMOU(J,II),J=1,6),II=1,NC),
     1((((CM(I,J,K,L),J=1,6),L=1,NC),
     1((CA(I,J,K,L),J=1,6),L=1,NC),I=1,6),K=1,NC),
     1((VGR2(JJ,J),J=1,NFFL),JJ=1,NJ),((VGM2(JJ,J),J=1,NFFL),JJ=1,NJ)
          READ(22,REC=N2) ((VXR2(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VXM2(JJ,J),J=1,NL),JJ=1,NJ)
          READ(23,REC=N2) ((VYR2(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VYM2(JJ,J),J=1,NL),JJ=1,NJ)
          READ(24,REC=N2) ((VZR2(JJ,J),J=1,NL),JJ=1,NJ),
     1                    ((VZM2(JJ,J),J=1,NL),JJ=1,NJ)
          W2=2*PI/TPE2
          DO J=1,6
            A2(J)=REAL(AMOU(J,KNC))
            A2(J+6)=AIMAG(AMOU(J,KNC))
          END DO
        ENDIF
        
C FIN D'IMPORTATION DES QUANTITES DU PREMIER ORDRE

        AKMM=AK1-AK2      ! NOMBRE D'ONDE - DE L'ONDE DIFFRACTEE
        AKPP=AK1+AK2      ! NOMBRE D'ONDE + DE L'ONDE DIFFRACTEE

c~         NR8=(2*NFAC*NJ+3)*2
c~         NR5=(NL*2*NJ+3)*2
c~         NR2=NFAC*2*NJ*2
        NR8=(2*NFAC*NJ+3)
        NR5=(NL*2*NJ+3)
        NR2=NFAC*2*NJ
        OPEN(31,FILE=NOMH1,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR5)
        OPEN(32,FILE=NOMH2,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)
        OPEN(33,FILE=NOMH3,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)
        OPEN(34,FILE=NOMH4,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO IJ=1,6         ! ITERATION SUR LES 6 MOUVEMENTS
          IJ6=IJ+6        ! POUR L'ENREGISTREMENT DES PARTIES IMAGINAIRES

C ON IMPORTE LES DONNEES POUR LA PULSATION DELTA I MOUVEMENT IJ
          IF(N3.NE.0)THEN ! QTF-
            NM=(N3-1)*6+IJ  !POUR LA SORTIE DE L'ENREGISTREMENT CORRESPONDANT A LA PULSATION N3 ET LE DOF IJ.
                            
            ! CF QTFinit.f90 
            ! TPE     : PERIODE
            ! BETA    : DIRECTION 
            ! AK      : NOMBRE D'ONDE
            ! P       : POTENTIEL(1) DE RADIATION SUR LE CORPS INDUIT PAR LE DOF IJ A LA PULSATION ~N3
            ! X,Y,Z   : VITESSES(1) DU FLUIDE SUR LE CORPS
            ! 
            ! P ET X,Y,Z : PSI ET SES DERIVEES SPATIALES
            READ(31,REC=NM)TPEM,BETA,AKM,((PRM(IC,J),J=1,NL),IC=1,NJ),
     1                                 ((PMM(IC,J),J=1,NL),IC=1,NJ)
            READ(32,REC=NM)((XRM(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((XMM(IC,J),J=1,NFAC),IC=1,NJ)
            READ(33,REC=NM)((YRM(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((YMM(IC,J),J=1,NFAC),IC=1,NJ)
            READ(34,REC=NM)((ZRM(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((ZMM(IC,J),J=1,NFAC),IC=1,NJ)
            WM=2*PI/TPEM

            IF(IJ.EQ.1) WRITE(LE,999,ADVANCE="NO")WM,W1,W2,CHAR(13)
            IF(ABS(WM-(W1-W2)).GT.5.E-03)THEN
              PRINT*,'WM=',WM,'  <>(W1-W2)=',W1-W2
c~               STOP POURQUOI STOP ?
            ENDIF
          ELSE
            TPEM=NAN
            AKM=NAN
            PRM=NAN
            PMM=NAN
            XRM=NAN
            XMM=NAN
            YRM=NAN
            YMM=NAN
            ZRM=NAN
            ZMM=NAN
          ENDIF

          IF(N4.NE.0)THEN ! QTF+ 
            NPL=(N4-1)*6+IJ
            READ(31,REC=NPL)TPEP,BETA,AKP,((PRP(IC,J),J=1,NL),IC=1,NJ),
     1                                 ((PMP(IC,J),J=1,NL),IC=1,NJ)
            READ(32,REC=NPL)((XRP(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((XMP(IC,J),J=1,NFAC),IC=1,NJ)
            READ(33,REC=NPL)((YRP(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((YMP(IC,J),J=1,NFAC),IC=1,NJ)
            READ(34,REC=NPL)((ZRP(IC,J),J=1,NFAC),IC=1,NJ),
     1                   ((ZMP(IC,J),J=1,NFAC),IC=1,NJ)
            WP=2*PI/TPEP
            IF(ABS(WP-(W1+W2)).GT.1.E-04)THEN
              PRINT*,'WP=',WP,'  <>(W1+W2)=',W1+W2
            ENDIF
          ELSE
            TPEP=NAN
            AKP=NAN
            PRP=NAN
            PMP=NAN
            XRP=NAN
            XMP=NAN
            YRP=NAN
            YMP=NAN
            ZRP=NAN
            ZMP=NAN
          ENDIF
C
C ** INTEGRALE D'HASKIND SUR LE CORPS **
C                             TERME:   (GRAD(PHI).GRAD(PSI))*(MOM.N)
C HVPSI
C --------------------------------------------------------------
          DO JJ=1,NJ
          DO II=1,NFAC
            B=-1    !CELA NE CHANGE ABSOLUMENT RIEN SUR LA SPHERE LIBRE
            PRO1(1,1)=(VXR1(JJ,II)*XRM(JJ,II)-VXM1(JJ,II)*XMM(JJ,II)+   !PARTIE REELLE DE \NABLA PHI_1 * \NABLA PSI-
     1             VYR1(JJ,II)*YRM(JJ,II)-VYM1(JJ,II)*YMM(JJ,II)+
     1             VZR1(JJ,II)*ZRM(JJ,II)-VZM1(JJ,II)*ZMM(JJ,II))
            PRO1(1,2)=(VXR1(JJ,II)*XMM(JJ,II)+VXM1(JJ,II)*XRM(JJ,II)+   !PARTIE IMAGINAIRE
     1             VYR1(JJ,II)*YMM(JJ,II)+VYM1(JJ,II)*YRM(JJ,II)+
     1             VZR1(JJ,II)*ZMM(JJ,II)+VZM1(JJ,II)*ZRM(JJ,II))
            PRO1(1,3)=(VXR2(JJ,II)*XRM(JJ,II)-B*VXM2(JJ,II)*XMM(JJ,II)+ !PARTIE REELLE DE \NABLA PHI_2 * \NABLA PSI- 
     1             VYR2(JJ,II)*YRM(JJ,II)-B*VYM2(JJ,II)*YMM(JJ,II)+     !
     1             VZR2(JJ,II)*ZRM(JJ,II)-B*VZM2(JJ,II)*ZMM(JJ,II))     !TEST DE MODIFICATION DE SIGNE POUR AVOIR LES BONS
            PRO1(1,4)=(B*VXR2(JJ,II)*XMM(JJ,II)+VXM2(JJ,II)*XRM(JJ,II)+ !PARTIE IMAGINAIRE  !JE LA GARDE CAR LA PENSE OK
     1             B*VYR2(JJ,II)*YMM(JJ,II)+VYM2(JJ,II)*YRM(JJ,II)+     !it is correct because it will be conjugated in PRODTN such that
     1             B*VZR2(JJ,II)*ZMM(JJ,II)+VZM2(JJ,II)*ZRM(JJ,II))     !!VXM2/VYM2/VZM2 become negatif and XMM/YMM/ZMM terms become positive

            PRO3(1,1)=(VXR1(JJ,II)*XRP(JJ,II)-VXM1(JJ,II)*XMP(JJ,II)+
     1             VYR1(JJ,II)*YRP(JJ,II)-VYM1(JJ,II)*YMP(JJ,II)+
     1             VZR1(JJ,II)*ZRP(JJ,II)-VZM1(JJ,II)*ZMP(JJ,II))
            PRO3(1,2)=(VXR1(JJ,II)*XMP(JJ,II)+VXM1(JJ,II)*XRP(JJ,II)+
     1             VYR1(JJ,II)*YMP(JJ,II)+VYM1(JJ,II)*YRP(JJ,II)+
     1             VZR1(JJ,II)*ZMP(JJ,II)+VZM1(JJ,II)*ZRP(JJ,II))
            PRO3(1,3)=(VXR2(JJ,II)*XRP(JJ,II)-VXM2(JJ,II)*XMP(JJ,II)+
     1             VYR2(JJ,II)*YRP(JJ,II)-VYM2(JJ,II)*YMP(JJ,II)+
     1             VZR2(JJ,II)*ZRP(JJ,II)-VZM2(JJ,II)*ZMP(JJ,II))
            PRO3(1,4)=(VXR2(JJ,II)*XMP(JJ,II)+VXM2(JJ,II)*XRP(JJ,II)+
     1             VYR2(JJ,II)*YMP(JJ,II)+VYM2(JJ,II)*YRP(JJ,II)+
     1             VZR2(JJ,II)*ZMP(JJ,II)+VZM2(JJ,II)*ZRP(JJ,II))

            XOII=XM(II)-XEFF
            YOII=YM(II)*(-1.)**(JJ+1)-YEFF
            ZOII=ZM(II)-ZEFF
            
            CN1=CN(II,1)
            CN2=CN(II,2)*(-1.)**(JJ+1)
            CN3=CN(II,3)
            
            XPR1=A1(1)+A1(5)*ZOII-A1(6)*YOII      !AMPLITUDE DU MOUVEMENT /POINT DE CALCUL DES EFFORTS MOM
            YPR1=A1(2)+A1(6)*XOII-A1(4)*ZOII
            ZPR1=A1(3)+A1(4)*YOII-A1(5)*XOII
            XPM1=A1(7)+A1(11)*ZOII-A1(12)*YOII    !PARTIE IMAGINAIRE
            YPM1=A1(8)+A1(12)*XOII-A1(10)*ZOII
            ZPM1=A1(9)+A1(10)*YOII-A1(11)*XOII

            XPR2=A2(1)+A2(5)*ZOII-A2(6)*YOII      !PAREIL POUR LE POINT 2
            YPR2=A2(2)+A2(6)*XOII-A2(4)*ZOII
            ZPR2=A2(3)+A2(4)*YOII-A2(5)*XOII
            XPM2=A2(7)+A2(11)*ZOII-A2(12)*YOII
            YPM2=A2(8)+A2(12)*XOII-A2(10)*ZOII
            ZPM2=A2(9)+A2(10)*YOII-A2(11)*XOII

            PRO2(1,1)=XPR1*CN1+YPR1*CN2+ZPR1*CN3  !ON CALCULE RE(MOM_1.N)
            PRO2(1,2)=XPM1*CN1+YPM1*CN2+ZPM1*CN3  !ON CALCULE IM(MOM_1.N)
            PRO2(1,3)=XPR2*CN1+YPR2*CN2+ZPR2*CN3  !ON CALCULE RE(MOM_2.N)
            PRO2(1,4)=XPM2*CN1+YPM2*CN2+ZPM2*CN3  !
            DO IP=1,4
              PRO4(1,IP)=PRO2(1,IP)
            END DO
            CALL PRODTN(PRO3,PRO4,PRO1,PRO2,FP,FM,1) ! not this products include conjugation for QTF- , RK 
            EFWPS(1,IJ)=EFWPS(1,IJ)+FP(1)*AIRE(II)
            EFWPS(1,IJ6)=EFWPS(1,IJ6)+FP(2)*AIRE(II)
            EFWMN(1,IJ)=EFWMN(1,IJ)+FM(1)*AIRE(II)  !ON CALCULE LA CONTRIBUTION A L'EFFORT SUIVANT X (RE) EN AJOUTANT POUR CHAQUE FACE CE QUE L'ON A CALCULE
            EFWMN(1,IJ6)=EFWMN(1,IJ6)+FM(2)*AIRE(II)
           ! IF (IJ==1.AND. JJ==1) THEN
           !     PRINT*,II,FM(1),FM(2),FP(1),FP(2)
           ! ENDIF

          END DO
          END DO
          EFP=EFWPS(1,IJ)
          EFM=EFWMN(1,IJ)                            !ON SAUVEGARDE AVANT DE MULTIPLIER PAR I
          EFWPS(1,IJ)=-EFWPS(1,IJ6)*RHO*WP            !1i*w*
          EFWPS(1,IJ6)=EFP*RHO*WP
          EFWMN(1,IJ)=-EFWMN(1,IJ6)*RHO*WM             
          EFWMN(1,IJ6)=EFM*RHO*WM

C ** INTEGRALE D'HASKIND SUR LE CORPS **
C                             TERME:   -(GRAD(PHI).N) (GRAD(PSI).MOM)
C HVN by RK
C --------------------------------------------------------------
          
          DO JJ=1,NJ
            DO II=1,NFAC
            XOII=XM(II)-XEFF
            YOII=YM(II)*(-1.)**(JJ+1)-YEFF
            ZOII=ZM(II)-ZEFF
            
            CN1=CN(II,1)
            CN2=CN(II,2)*(-1.)**(JJ+1)
            CN3=CN(II,3)
            
            XPR1=A1(1)+A1(5)*ZOII-A1(6)*YOII      !AMPLITUDE DU MOUVEMENT /POINT DE CALCUL DES EFFORTS MOM
            YPR1=A1(2)+A1(6)*XOII-A1(4)*ZOII
            ZPR1=A1(3)+A1(4)*YOII-A1(5)*XOII
            XPM1=A1(7)+A1(11)*ZOII-A1(12)*YOII    !PARTIE IMAGINAIRE
            YPM1=A1(8)+A1(12)*XOII-A1(10)*ZOII
            ZPM1=A1(9)+A1(10)*YOII-A1(11)*XOII

            XPR2=A2(1)+A2(5)*ZOII-A2(6)*YOII      !PAREIL POUR LE POINT 2
            YPR2=A2(2)+A2(6)*XOII-A2(4)*ZOII
            ZPR2=A2(3)+A2(4)*YOII-A2(5)*XOII
            XPM2=A2(7)+A2(11)*ZOII-A2(12)*YOII
            YPM2=A2(8)+A2(12)*XOII-A2(10)*ZOII
            ZPM2=A2(9)+A2(10)*YOII-A2(11)*XOII

            PRO1(1,1)=VXR1(JJ,II)*CN1+VYR1(JJ,II)*CN2+VZR1(JJ,II)*CN3 !PARTIE REELLE DE \NABLA PHI_1.N
            PRO1(1,2)=VXM1(JJ,II)*CN1+VYM1(JJ,II)*CN2+VZM1(JJ,II)*CN3 !Imaginer part
            PRO1(1,3)=VXR2(JJ,II)*CN1+VYR2(JJ,II)*CN2+VZR2(JJ,II)*CN3 !PARTIE REELLE DE \NABLA conj(PHI_2).N
            PRO1(1,4)=VXM2(JJ,II)*CN1+VYM2(JJ,II)*CN2+VZM2(JJ,II)*CN3 !Imaginer part
            DO IP=1,4
              PRO3(1,IP)=PRO1(1,IP)
            END DO
            B=-1    !CELA NE CHANGE ABSOLUMENT RIEN SUR LA SPHERE LIBRE

            PRO2(1,1)=(XRM(JJ,II)*XPR1-XMM(JJ,II)*XPM1+   !PARTIE REELLE DE (GRADPSI.MOM_1)
     1             YRM(JJ,II)*YPR1-YMM(JJ,II)*YPM1+
     1             ZRM(JJ,II)*ZPR1-ZMM(JJ,II)*ZPM1)
            PRO2(1,2)=(XMM(JJ,II)*XPR1+XRM(JJ,II)*XPM1+   !PARTIE IMAGINAIRE
     1             YMM(JJ,II)*YPR1+YRM(JJ,II)*YPM1+
     1             ZMM(JJ,II)*YPR1+ZRM(JJ,II)*ZPM1)
            PRO2(1,3)=(XRM(JJ,II)*XPR2-B*XMM(JJ,II)*XPM2+   !PARTIE REELLE DE (GRADPSI.Conj(MOM_2))
     1             YRM(JJ,II)*YPR2-B*YMM(JJ,II)*YPM2+
     1             ZRM(JJ,II)*ZPR2-B*ZMM(JJ,II)*ZPM2)
            PRO2(1,4)=(B*XMM(JJ,II)*XPR2+XRM(JJ,II)*XPM2+   !PARTIE IMAGINAIRE
     1             B*YMM(JJ,II)*YPR2+YRM(JJ,II)*YPM2+       !!it is correct because it will be conjugated in PRODTN
     1             B*ZMM(JJ,II)*YPR2+ZRM(JJ,II)*ZPM2)

            PRO4(1,1)=(XRP(JJ,II)*XPR1-XMP(JJ,II)*XPM1+   !PARTIE REELLE DE (GRADPSI.MOM_1) for QTF+
     1             YRP(JJ,II)*YPR1-YMP(JJ,II)*YPM1+
     1             ZRP(JJ,II)*ZPR1-ZMP(JJ,II)*ZPM1)
            PRO4(1,2)=(XMP(JJ,II)*XPR1+XRP(JJ,II)*XPM1+   !PARTIE IMAGINAIRE
     1             YMP(JJ,II)*YPR1+YRP(JJ,II)*YPM1+
     1             ZMP(JJ,II)*YPR1+ZRP(JJ,II)*ZPM1)
            PRO4(1,3)=(XRP(JJ,II)*XPR2-XMP(JJ,II)*XPM2+   !PARTIE REELLE DE (GRADPSI.(MOM_2))
     1             YRP(JJ,II)*YPR2-YMP(JJ,II)*YPM2+
     1             ZRP(JJ,II)*ZPR2-ZMP(JJ,II)*ZPM2)
            PRO4(1,4)=(XMP(JJ,II)*XPR2+XRP(JJ,II)*XPM2+   !PARTIE IMAGINAIRE
     1             YMP(JJ,II)*YPR2+YRP(JJ,II)*YPM2+
     1             ZMP(JJ,II)*YPR2+ZRP(JJ,II)*ZPM2)

     
            CALL PRODTN(PRO3,PRO4,PRO1,PRO2,FP,FM,1)    ! note this products includes conjugation for QTF-
            EFWPS(2,IJ)=EFWPS(2,IJ)+FP(1)*AIRE(II)
            EFWPS(2,IJ6)=EFWPS(2,IJ6)+FP(2)*AIRE(II)
            EFWMN(2,IJ)=EFWMN(2,IJ)+FM(1)*AIRE(II)  !ON CALCULE LA CONTRIBUTION A L'EFFORT SUIVANT X (RE) EN AJOUTANT POUR CHAQUE FACE CE QUE L'ON A CALCULE
            EFWMN(2,IJ6)=EFWMN(2,IJ6)+FM(2)*AIRE(II)
           ! IF (IJ==1.AND. II.LE.100) THEN
           !     !PRINT*,II,AIRE(II)
           !     PRINT*,JJ,II,FM(1),FM(2),FP(1),FP(2)
           !     PRINT*,II,VXR1(JJ,II),VYR1(JJ,II),VZR1(JJ,II)
           !     PRINT*,II,VXM1(JJ,II),VYM1(JJ,II),VZM1(JJ,II)
           !    ! PRINT*,JJ,II,EFWMN(2,IJ),EFWMN(2,IJ6)
           ! ENDIF

          END DO
          END DO
          EFP=EFWPS(2,IJ)
          EFM=EFWMN(2,IJ)                            !ON SAUVEGARDE AVANT DE MULTIPLIER PAR -I
          EFWPS(2,IJ)=EFWPS(2,IJ6)*RHO*WP            !-1i*w*
          EFWPS(2,IJ6)=-EFP*RHO*WP
          EFWMN(2,IJ)=EFWMN(2,IJ6)*RHO*WM             
          EFWMN(2,IJ6)=-EFM*RHO*WM
                
       !  print*,EFWMN(2,IJ),EFWMN(2,IJ6)


CC ** INTEGRALE D'HASKIND SUR LE CORPS **
CC                             TERME:   -GRAD(PHI).N*GRAD(PSI).MOM
CC HVN
CC ---------------------------------------------------------------
CC CODE EN COMPLEXE PAR FABIEN ROBAUX, POUR INNOSEA/EDF
C
C          DO JJ=1,NJ
C          DO II=1,NFAC
C            CN1=CN(II,1)
C            CN2=CN(II,2)*(-1.)**(JJ+1)
C            CN3=CN(II,3)
C
C            VX1=CMPLX(VXR1(JJ,II),VXM1(JJ,II))
C            VY1=CMPLX(VYR1(JJ,II),VYM1(JJ,II))
C            VZ1=CMPLX(VZR1(JJ,II),VZM1(JJ,II))
C            VN1=VX1*CN1+VY1*CN2+VZ1*CN3    !CALCUL DE \NABLA{PHI}.N
C            
C            VX2=CMPLX(VXR2(JJ,II),VXM2(JJ,II))
C            VY2=CMPLX(VYR2(JJ,II),VYM2(JJ,II))
C            VZ2=CMPLX(VZR2(JJ,II),VZM2(JJ,II))
C            VN2=VX2*CN1+VY2*CN2+VZ2*CN3    !CALCUL DE \NABLA{PHI}.N
C
C            DO I=1,6
C                AC1(I)=CMPLX(A1(I),A1(I+6))    ! DEPLACEMENT 1 COMPLEXE
C                AC2(I)=CMPLX(A2(I),A2(I+6))    ! DEPLACEMENT 2 COMPLEXE
C            END DO
C            XOII=XM(II)-XEFF                   ! added by RK
C            YOII=YM(II)*(-1.)**(JJ+1)-YEFF
C            ZOII=ZM(II)-ZEFF
C
C            XP1=AC1(1)+AC1(5)*ZOII-AC1(6)*YOII     !
C            YP1=AC1(2)+AC1(6)*XOII-AC1(4)*ZOII     !CALCUL DE MOM_1
C            ZP1=AC1(3)+AC1(4)*YOII-AC1(5)*XOII     !
C            
C            XP2=AC2(1)+AC2(5)*ZOII-AC2(6)*YOII     !CALCUL DE MOM_2
C            YP2=AC2(2)+AC2(6)*XOII-AC2(4)*ZOII     !
C            ZP2=AC2(3)+AC2(4)*YOII-AC2(5)*XOII     !
C            
C            XMC=CMPLX(XRM(JJ,II),XMM(JJ,II))         !
C            YMC=CMPLX(YRM(JJ,II),YMM(JJ,II))         ! NABLA PSI-
C            ZMC=CMPLX(ZRM(JJ,II),ZMM(JJ,II))         !
C            
C            GPM1=XMC*XP1+YMC*YP1+ZMC*ZP1                             !CALCUL GRAD(PSI)(MOM_1)
C            GPM2C=XMC*CONJG(XP2)+YMC*CONJG(YP2)+ZMC*CONJG(ZP2)       !CALCUL GRAD(PSI)(MOM_2*)
C            GPM2=XMC*XP2+YMC*YP2+ZMC*ZP2                             !CALCUL GRAD(PSI)(MOM_2)
C
C
C
C
C            FMC=1./2.*(GPM2C*VN1+CONJG(VN2)*GPM1)
C            EFWMNC(2,IJ)=EFWMNC(2,IJ)+FMC*AIRE(II)
C          END DO
C          END DO
C
C          EFWMNC(2,IJ)=ZI*EFWMNC(2,IJ)*RHO*WM
C          
C          EFWMN(2,IJ)=REAL(EFWMNC(2,IJ))
C          EFWMN(2,IJ6)=IMAG(EFWMNC(2,IJ))
C


C
C ** INTEGRALE D'HASKIND SUR LE CORPS **
C                  TERME:  INTEGRALE SUR LA LIGNE DE FLOTTAISON
C HFLO
C -------------------------------------------------------------
          DO JJ=1,NJ
          DO II=NFAC+1,NL
            II1=II-NFAC
            XOII=XM(II)-XEFF
            YOII=YM(II)*(-1.)**(JJ+1)-YEFF
            ZOII=ZM(II)-ZEFF
            GAMMA0=ALF(II1)  !added by RK, LONGUEUR DE LA LIGNE DE FLOTTAISON --> OUI! cf QTFGeom l.99
C ATTENTION A LA DEFINITION DE DL SUR LA PARTIE SYMETRIQUE DU CORPS
C LE CONTOUR EST ORIENTE VERS L'EXTERIEUR DU CORPS (INVERSE DE Z!)
            CL1=Q(II)*(-1.)**(JJ+1)!!DL_X/|DL| ny/|n|
            
            CL2=-P(II)             !!DL_Y/|DL|  -nx/|n|

            XPR1=A1(1)+A1(5)*ZOII-A1(6)*YOII
            YPR1=A1(2)+A1(6)*XOII-A1(4)*ZOII
            ZPR1=A1(3)+A1(4)*YOII-A1(5)*XOII
            XPM1=A1(7)+A1(11)*ZOII-A1(12)*YOII
            YPM1=A1(8)+A1(12)*XOII-A1(10)*ZOII
            ZPM1=A1(9)+A1(10)*YOII-A1(11)*XOII

            XPR2=A2(1)+A2(5)*ZOII-A2(6)*YOII        !
            YPR2=A2(2)+A2(6)*XOII-A2(4)*ZOII        !
            ZPR2=A2(3)+A2(4)*YOII-A2(5)*XOII        !ON CALCULE LES PARTIES REELLES ET IMAGINAIRES DES DEPLACEMENTS/POINT DE CALCUL DES EFFORTS
            XPM2=A2(7)+A2(11)*ZOII-A2(12)*YOII      !
            YPM2=A2(8)+A2(12)*XOII-A2(10)*ZOII      !
            ZPM2=A2(9)+A2(10)*YOII-A2(11)*XOII      !

            PRO1(1,1)=ZPR1*CL1                      !ON FAIT MOM1_Z*DL_X
            PRO1(1,2)=ZPM1*CL1
            PRO1(1,3)=ZPR2*CL1                      !ON FAIT MOM2_Z*DL_X
            PRO1(1,4)=ZPM2*CL1

            PRO1(2,1)=-YPR1*CL1                     !ON FAIT -MOM1_Y*DL_X
            PRO1(2,2)=-YPM1*CL1
            PRO1(2,3)=-YPR2*CL1
            PRO1(2,4)=-YPM2*CL1

            PRO1(3,1)=XPR1*CL2                       !ON FAIT MOM1_X*DL_Y
            PRO1(3,2)=XPM1*CL2
            PRO1(3,3)=XPR2*CL2
            PRO1(3,4)=XPM2*CL2

            PRO1(4,1)=-ZPR1*CL2                      !ON FAIT -MOM1_Z*DL_Y
            PRO1(4,2)=-ZPM1*CL2
            PRO1(4,3)=-ZPR2*CL2
            PRO1(4,4)=-ZPM2*CL2

            PRO2(1,1)=VYR1(JJ,II)                   !\NABLA PHI_Y
            PRO2(1,2)=VYM1(JJ,II)
            PRO2(1,3)=VYR2(JJ,II)
            PRO2(1,4)=VYM2(JJ,II)

            PRO2(2,1)=VZR1(JJ,II)                   !\NABLA PHI_Z
            PRO2(2,2)=VZM1(JJ,II)
            PRO2(2,3)=VZR2(JJ,II)
            PRO2(2,4)=VZM2(JJ,II)

            PRO2(3,1)=VZR1(JJ,II)                   !\NABLA PHI_Z
            PRO2(3,2)=VZM1(JJ,II)
            PRO2(3,3)=VZR2(JJ,II)
            PRO2(3,4)=VZM2(JJ,II)

            PRO2(4,1)=VXR1(JJ,II)                   !\NABLA PHI_X
            PRO2(4,2)=VXM1(JJ,II)
            PRO2(4,3)=VXR2(JJ,II)
            PRO2(4,4)=VXM2(JJ,II)
           ! curl(gradPhi,MOM) 
            CALL PRODTN(PRO1,PRO2,PRO1,PRO2,FP,FM,4)!NORMALEMENT VERIFIE, C'EST OK ON OBTIENT BIEN LES BONS SIGNES:(POUR FM) +1/2 DL/|DL|.(\NABLA PHI_1 ^ MOM_2*   + \NABLA PHI_2* ^ MOM_1)
            AIR=AIRE(II)
           ! print*,GAMMA0,AIR  !GAMMA0=AIR for II>NFAC 
            EFWPS(3,IJ)=EFWPS(3,IJ)+GAMMA0*(FP(1)*PRP(JJ,II) !AIR replaced with GAMMA0 by RK
     1 -FP(2)*PMP(JJ,II))                            !! PRP ET PMP SONT LES PARTIES R ET I DE PSI+ POUR LA PULSATION W+ ET DE DOF K
            EFWPS(3,IJ6)=EFWPS(3,IJ6)+GAMMA0*(FP(1)*PMP(JJ,II)
     1 +FP(2)*PRP(JJ,II))
            EFWMN(3,IJ)=EFWMN(3,IJ)+GAMMA0*(FM(1)*PRM(JJ,II)
     1 -FM(2)*PMM(JJ,II))                            ! ON FAIT TERME 4LIGNES PLUS HAUT *PSI_K ET ON INTEGRE
            EFWMN(3,IJ6)=EFWMN(3,IJ6)+GAMMA0*(FM(1)*PMM(JJ,II)
     1 +FM(2)*PRM(JJ,II))
            IF (IJ==1 .AND. JJ==2) THEN
!                PRINT*,II-NFAC,GAMMA0*(FM(1)*PRM(JJ,II)
!     1 -FM(2)*PMM(JJ,II)), GAMMA0*(FM(1)*PMM(JJ,II)
!     1 +FM(2)*PRM(JJ,II)) 
!                PRINT*,II-NFAC,RHO*WM*GAMMA0*(FM(1)*PMM(JJ,II)
!     1 +FM(2)*PRM(JJ,II)), 
!     1   -RHO*WM*GAMMA0*(FM(1)*PRM(JJ,II)
!     1 -FM(2)*PMM(JJ,II)) 
!               PRINT*,II-NFAC,FM(1)*GAMMA0,FM(2)*GAMMA0
            ENDIF

          END DO
          END DO
          EFP=EFWPS(3,IJ)
          EFM=EFWMN(3,IJ)
          EFWPS(3,IJ) =-EFWPS(3,IJ6)*RHO*WP
          EFWPS(3,IJ6)=EFP*RHO*WP
          EFWMN(3,IJ) =-EFWMN(3,IJ6)*RHO*WM            !ON MULTIPLIE PAR IWM RHO
          EFWMN(3,IJ6)=EFM*RHO*WM                 
C
C
C ** INTEGRALE D'HASKIND SUR LE CORPS **
C                             TERME:  PSI*(2.GRAD(PHI)-VE).R.N
C HVRN
C --------------------------------------------------------------
          DO JJ=1,NJ
          DO II=1,NFAC
            CN1=CN(II,1)
            CN2=CN(II,2)*(-1.)**(JJ+1)
            CN3=CN(II,3)

            XOII=XM(II)-XEFF
            YOII=YM(II)*(-1.)**(JJ+1)-YEFF
            ZOII=ZM(II)-ZEFF

            XPR1=A1(1)+A1(5)*ZOII-A1(6)*YOII
            YPR1=A1(2)+A1(6)*XOII-A1(4)*ZOII
            ZPR1=A1(3)+A1(4)*YOII-A1(5)*XOII
            XPM1=A1(7)+A1(11)*ZOII-A1(12)*YOII
            YPM1=A1(8)+A1(12)*XOII-A1(10)*ZOII
            ZPM1=A1(9)+A1(10)*YOII-A1(11)*XOII

            XPR2=A2(1)+A2(5)*ZOII-A2(6)*YOII
            YPR2=A2(2)+A2(6)*XOII-A2(4)*ZOII
            ZPR2=A2(3)+A2(4)*YOII-A2(5)*XOII
            XPM2=A2(7)+A2(11)*ZOII-A2(12)*YOII
            YPM2=A2(8)+A2(12)*XOII-A2(10)*ZOII
            ZPM2=A2(9)+A2(10)*YOII-A2(11)*XOII

            A=1                               !!ON FIXE CE PARAMETRE A 1 CAR LE DEVEPNT EN SERIE DU PB DE BASE (GRAD(PHI).N=VE.N)
            PRO1(1,1)=(YPM1*W1-A*VYR1(JJ,II))           ! dtMoM=-iw*MoM
            PRO1(1,2)=(-YPR1*W1-A*VYM1(JJ,II))
            PRO1(1,3)=(YPM2*W2-A*VYR2(JJ,II))
            PRO1(1,4)=(-YPR2*W2-A*VYM2(JJ,II))

            PRO1(2,1)=(ZPM1*W1-A*VZR1(JJ,II))
            PRO1(2,2)=(-ZPR1*W1-A*VZM1(JJ,II))
            PRO1(2,3)=(ZPM2*W2-A*VZR2(JJ,II))
            PRO1(2,4)=(-ZPR2*W2-A*VZM2(JJ,II))

            PRO1(3,1)=(XPM1*W1-A*VXR1(JJ,II))
            PRO1(3,2)=(-XPR1*W1-A*VXM1(JJ,II))
            PRO1(3,3)=(XPM2*W2-A*VXR2(JJ,II))
            PRO1(3,4)=(-XPR2*W2-A*VXM2(JJ,II))
            DO IP=1,4
              PRO1(4,IP)=PRO1(2,IP)
              PRO1(5,IP)=PRO1(3,IP)
              PRO1(6,IP)=PRO1(1,IP)
            END DO
            
            PRO2(1,1)=-A1(6)*CN1                        !! R.N
            PRO2(1,2)=-A1(12)*CN1
            PRO2(1,3)=-A2(6)*CN1
            PRO2(1,4)=-A2(12)*CN1

            PRO2(2,1)=A1(5)*CN1
            PRO2(2,2)=A1(11)*CN1
            PRO2(2,3)=A2(5)*CN1
            PRO2(2,4)=A2(11)*CN1

            PRO2(3,1)=A1(6)*CN2
            PRO2(3,2)=A1(12)*CN2
            PRO2(3,3)=A2(6)*CN2
            PRO2(3,4)=A2(12)*CN2

            PRO2(4,1)=-A1(4)*CN2
            PRO2(4,2)=-A1(10)*CN2
            PRO2(4,3)=-A2(4)*CN2
            PRO2(4,4)=-A2(10)*CN2

            PRO2(5,1)=-A1(5)*CN3
            PRO2(5,2)=-A1(11)*CN3
            PRO2(5,3)=-A2(5)*CN3
            PRO2(5,4)=-A2(11)*CN3

            PRO2(6,1)=A1(4)*CN3
            PRO2(6,2)=A1(10)*CN3
            PRO2(6,3)=A2(4)*CN3
            PRO2(6,4)=A2(10)*CN3
            !edited by RK 220712 to be same as in the eqs
            CALL PRODTN(PRO1,PRO2,PRO1,-PRO2,FP,FM,6)  !-PRO2 to be consistent (RK) 
            AIR=AIRE(II)
            EFWPS(4,IJ)=EFWPS(4,IJ)+AIR*(FP(1)*PRP(JJ,II)
     1 -FP(2)*PMP(JJ,II))                                   !MULTIPLICATION PAR PSI
            EFWPS(4,IJ6)=EFWPS(4,IJ6)+AIR*(FP(1)*PMP(JJ,II)
     1 +FP(2)*PRP(JJ,II))
            EFWMN(4,IJ)=EFWMN(4,IJ)+AIR*(FM(1)*PRM(JJ,II)
     1 -FM(2)*PMM(JJ,II))
            EFWMN(4,IJ6)=EFWMN(4,IJ6)+AIR*(FM(1)*PMM(JJ,II)
     1 +FM(2)*PRM(JJ,II))
            !IF (IJ==1) THEN
            !    PRINT*,II,FM(1),FM(2)
            !    PRINT*,II,FP(1),FP(2)
            !ENDIF
          END DO
          END DO
          EFP=EFWPS(4,IJ)
          EFM=EFWMN(4,IJ)
          EFWPS(4,IJ)=EFWPS(4,IJ6)*RHO*WP
          EFWPS(4,IJ6)=-EFP*RHO*WP
          EFWMN(4,IJ)=EFWMN(4,IJ6)*RHO*WM                   !MULTIPLICATION PAR -I*WM*RHO
          EFWMN(4,IJ6)=-EFM*RHO*WM
C
C
C ***  CONTRIBUTION DE LA DERIVE NORMALE DU
C                                POTENTIEL INCIDENT DU 2E ORDRE ***
C HIDN
C--------------------------------------------------------------------
          IF(H.GT.0)THEN
            HII=H
          ELSE
            HII=1.E20
          ENDIF
C          IF(W1.EQ.W2)THEN ! closed by RK
C            PHIM=0.
C            PHIP=0.
C          ELSE
            IF(HII.GE.1.E19)THEN
              PHIM=AMAX1(W1,W2)
              PHIP=0
            ELSE
              AKH1=AK1*HII
              AKH2=AK2*HII
              AKHMM=AKMM*HII
              AKHPP=AKPP*HII
              PHIM=WM*AK1*AK2*(1+TANH(AKH1)*TANH(AKH2))/(W1*W2)+
     1    (AK1*AK1/(W1*COSH(AKH1)**2)-AK2*AK2/(W2*COSH(AKH2)**2))/2
              PHIM=PHIM*G*G/(G*AKMM*TANH(AKHMM)-WM**2)
              PHIP=WP*AK1*AK2*(TANH(AKH1)*TANH(AKH2)-1)/(W1*W2)-
     1    (AK1*AK1/(W1*COSH(AKH1)**2)+AK2*AK2/(W2*COSH(AKH2)**2))/2
              PHIP=-PHIP*G*G/(G*AKPP*TANH(AKHPP)-WP**2)
C due to convention of the first order incoming potential in NEMOH
C PHIP is multiplied by -1, by RK  
            ENDIF
C         ENDIF

          DO JJ=1,NJ
          DO II=1,NFAC
            CN1=CN(II,1)
            CN2=CN(II,2)*(-1.)**(JJ+1)
            CN3=CN(II,3)
            COEFM=AKMM*((XM(II)-XEFF)*CB+(YM(II)*(-1.)**(JJ+1)-YEFF)*SB)
            COEFP=AKPP*((XM(II)-XEFF)*CB+(YM(II)*(-1.)**(JJ+1)-YEFF)*SB)
            AZCHM=CH(AKMM,ZM(II),HII)
            AZCHP=CH(AKPP,ZM(II),HII)
            AZSHM=SH(AKMM,ZM(II),HII)
            AZSHP=SH(AKPP,ZM(II),HII)
            PHIM_R=-PHIM*AZCHM*SIN(COEFM)
            PHIM_I=PHIM*AZCHM*COS(COEFM)
            dnPHIM_R=(-AKMM*CB*PHIM_I)*CN1+(-AKMM*SB*PHIM_I)*CN2+
     1               AKMM*TANH(AKMM*(HII+ZM(II)))*PHIM_R*CN3 
            dnPHIM_I=(AKMM*CB*PHIM_R)*CN1+(AKMM*SB*PHIM_R)*CN2+
     1               AKMM*TANH(AKMM*(HII+ZM(II)))*PHIM_I*CN3  



            ACQM=AIRE(II)*PHIM*AKMM
            DFRM=-AZCHM*(CN1*CB+CN2*SB)*COS(COEFM)-AZSHM*CN3*SIN(COEFM)
            DFMM=-AZCHM*(CN1*CB+CN2*SB)*SIN(COEFM)+AZSHM*CN3*COS(COEFM)
            ACQP=AIRE(II)*PHIP*AKPP
            DFRP=-AZCHP*(CN1*CB+CN2*SB)*COS(COEFP)-AZSHP*CN3*SIN(COEFP)
            DFMP=-AZCHP*(CN1*CB+CN2*SB)*SIN(COEFP)+AZSHP*CN3*COS(COEFP)
            IF (ZM(II).LT.0.) THEN   ! not calculating on lid panels for irreg. freq removal ,added by RK
            EFWPS(5,IJ)=EFWPS(5,IJ)+ACQP*(DFRP*PRP(JJ,II)
     1 -DFMP*PMP(JJ,II))
            EFWPS(5,IJ6)=EFWPS(5,IJ6)+ACQP*(DFRP*PMP(JJ,II)
     1 +DFMP*PRP(JJ,II))
            EFWMN(5,IJ)=EFWMN(5,IJ)+ACQM*(DFRM*PRM(JJ,II)
     1 -DFMM*PMM(JJ,II))
            EFWMN(5,IJ6)=EFWMN(5,IJ6)+ACQM*(DFRM*PMM(JJ,II)
     1 +DFMM*PRM(JJ,II))
C             IF (IJ==1) THEN
C             print*,II,(-AKMM*CB*PHIM_I),(AKMM*CB*PHIM_R)
C             print*,II,(-AKMM*SB*PHIM_I),(AKMM*SB*PHIM_R)
C             print*,II,AKMM*TANH(AKMM*(HII+ZM(II)))*PHIM_R,
C     1                 AKMM*TANH(AKMM*(HII+ZM(II)))*PHIM_I

C              print*,II,dnPHIM_R,dnPHIM_I
C              print*,II,PRM(JJ,II),PMM(JJ,II)
C               print*,II,ACQM*(DFRM*PRM(JJ,II)
C     1 -DFMM*PMM(JJ,II)),ACQM*(DFRM*PMM(JJ,II)
C     1 +DFMM*PRM(JJ,II))
C              ENDIF

            ENDIF
          END DO
          END DO
          EFP=EFWPS(5,IJ)
          EFM=EFWMN(5,IJ)
          EFWPS(5,IJ)=-EFWPS(5,IJ6)*RHO*WP              !i omega
          EFWPS(5,IJ6)=EFP*RHO*WP
          EFWMN(5,IJ)=-EFWMN(5,IJ6)*RHO*WM
          EFWMN(5,IJ6)=EFM*RHO*WM
             
C FIN DE LA SOMMATION SUR LES 6 MOUVEMENTS
        END DO
        
        
        
        
        
C
C
C ****    CONTRIBUTION DE LA HOULE INCIDENTE DU 2D ORDRE   ****
C                                EFFORTS DE FROUDE KRYLOV
C HFRK
C--------------------------------------------------------------
C  PROFONDEUR FINIE , 2 HOULES DE MEME DIRECTION
C                      E = A*SIN(AK*(X.COSB+Y.SINB)-WT)
C  Replace convention by RK to be same as in first order NEMOH
C                     Phi=-iAg/w f(z) e^(ik.x)
C  difference frequency part is same but not in the sum freq
C So PHIP is multiply by - sign
C          IF(W1.EQ.W2)THEN ! closed by RK
C            PHIM=0.
C            PHIP=0.
C          ELSE
        IF(HII.GE.1.E19)THEN
            PHIM=AMAX1(W1,W2)
            PHIP=0
          ELSE
            AKH1=AK1*HII
            AKH2=AK2*HII
            AKHMM=AKMM*HII
            AKHPP=AKPP*HII
            PHIM=WM*AK1*AK2*(1+TANH(AKH1)*TANH(AKH2))/(W1*W2)+
     1    (AK1*AK1/(W1*COSH(AKH1)**2)-AK2*AK2/(W2*COSH(AKH2)**2))/2
            PHIM=PHIM*G*G/(G*AKMM*TANH(AKHMM)-WM**2)
            PHIP=WP*AK1*AK2*(TANH(AKH1)*TANH(AKH2)-1)/(W1*W2)-
     1    (AK1*AK1/(W1*COSH(AKH1)**2)+AK2*AK2/(W2*COSH(AKH2)**2))/2
            PHIP=-PHIP*G*G/(G*AKPP*TANH(AKHPP)-WP**2)
C due to convention of the first order incoming potential in NEMOH
C PHIP is multiplied by -1, by RK  
          ENDIF
C       ENDIF

        DO JJ=1,NJ
        DO II=1,NFAC
          COEFM=AKMM*((XM(II)-XEFF)*CB+(YM(II)*(-1.)**(JJ+1)-YEFF)*SB)
          COEFP=AKPP*((XM(II)-XEFF)*CB+(YM(II)*(-1.)**(JJ+1)-YEFF)*SB)
          AZCHM=CH(AKMM,ZM(II),HII)
          AZCHP=CH(AKPP,ZM(II),HII)
          DO IJ=1,6
            IJ6=IJ+6
            ACQM=AIRE(II)*RHO*WM*PHIM*AZCHM*((-1.)**(IJ+1))**(JJ+1)
     1 *CN(II,IJ)
            ACQP=AIRE(II)*RHO*WP*PHIP*AZCHP*((-1.)**(IJ+1))**(JJ+1)
     1 *CN(II,IJ)
            IF (ZM(II).LT.0.) THEN   ! not calculating on lid panels for irreg. freq removal
            EFWPS(6,IJ)=EFWPS(6,IJ)+ACQP*COS(COEFP)
            EFWPS(6,IJ6)=EFWPS(6,IJ6)+ACQP*SIN(COEFP)
            EFWMN(6,IJ)=EFWMN(6,IJ)+ACQM*COS(COEFM)
            EFWMN(6,IJ6)=EFWMN(6,IJ6)+ACQM*SIN(COEFM)
            ENDIF
          END DO
        END DO
        END DO
C *** FONCTION DE TRANSTERT H2 (DIVISION PAR 2) ***
C -------------------------------------------------
        DO J2=1,12
        DO J1=1,6
          EFWPS(J1,J2)=EFWPS(J1,J2)/2
          EFWMN(J1,J2)=EFWMN(J1,J2)/2
        END DO
        END DO
C
C*******************************************************************
C    ----    SOMMATION  DES  DIFFERENTES  CONTRIBUTIONS    ----
C*******************************************************************
        DO J2=1,12  !SOMMATION SUR LES DOF REELS PUIS IMAGINAIRES
        DO J1=1,6   !SOMMATION SUR LES CONTRIBUTIONS
          EFWPS(7,J2)=EFWPS(7,J2)+EFWPS(J1,J2)
          EFWMN(7,J2)=EFWMN(7,J2)+EFWMN(J1,J2)
        END DO
        END DO
C
C********************************************************************
C -----   ECRITURE   ----
C ON ECRIT EN FAIT W1,W2,EFXREEL/ADIMENSION,EFXIM/ADIMENSION,EFYREEL/ADIMENSION,EFYIM/ADIMENSION
C ...MOMENTXREEL/ADIM..ETC DANS LE FICHER AH.RES ET AHP.RES
C********************************************************************
        WRITE(40,798)W1,W2,(EFWMN(7,J)/AD,EFWMN(7,J+6)/AD,J=1,3),
     1(EFWMN(7,J)/AD/XL,EFWMN(7,J+6)/AD/XL,J=4,6)
        WRITE(41,798)W1,W2,(EFWPS(7,J)/AD,EFWPS(7,J+6)/AD,J=1,3),
     1(EFWPS(7,J)/AD/XL,EFWPS(7,J+6)/AD/XL,J=4,6)
 798    FORMAT(14E11.3)
 
 !REMPLISSAGE DE CE QUI NOUS SERA UTILE POUR CREER LES FICHER QTFS
        TABWR2(N2)=W2 !!!!PROBLEME , LA DERNIERE VALEUR NE SERA PAS REMPLIE, ON SUPPOSERA SYMETRIQUE.
        TABWR1(N1)=W1

        DO J=1,6
          QTFM(J,N1,N2)=EFWMN(7,J)        !PARTIE REELE DE LA QTF MOINS SUIVANT LE DOF J
          QTFM(J+6,N1,N2)=EFWMN(7,J+6)    !PARTIE IMAGINAIRE DE LA QTF MOINS SUIVANT LE DOF J
          QTFP(J,N1,N2)=EFWPS(7,J)        !PARTIE REELE DE LA QTF PLUS SUIVANT LE DOF J
          QTFP(J+6,N1,N2)=EFWPS(7,J+6)    !PARTIE IMAGINAIRE DE LA QTF PLUS SUIVANT LE DOF J
          DO K=1,7
            CONTRIB(K,J,N1,N2,1)=EFWMN(K,J)
            CONTRIB(K,J+6,N1,N2,1)=EFWMN(K,J+6)
            CONTRIB(K,J,N1,N2,2)=EFWPS(K,J)
            CONTRIB(K,J+6,N1,N2,2)=EFWPS(K,J+6)
          END DO 
        END DO

        J11=7
        J12=7
        
        IF(J11.NE.0)THEN
          DO J1=J11,J12                 ! AU VU DE LA BOUCLE ON LA FAIT UNE FOIS POUR J1=7 (UNIQUEMENT LA SOMME DES CONTRIBS)
            DO J2=1,12
              SOLP(J2)=EFWPS(J1,J2)     ! ON ECRIT DANS SOLP(J2) EFWPS(7,J2)
              SOLM(J2)=EFWMN(J1,J2)
            END DO
            CALL AFIN(12,SOLP,1.E-06)   ! ROUTINE QUI SUPPRIME LES VALEURS DE SOLP TQ SOLP(I)<E-06*MAX(SOLP)
            CALL AFIN(12,SOLM,1.E-06)
            DO J2=1,12
              EFWPS(J1,J2)=SOLP(J2)     ! ON REECRIT APRES AVOIR SUPPRIME LES TROP PETITES VALEURS DANS LE TABLEAU
              EFWMN(J1,J2)=SOLM(J2)
            END DO

            DO J2=1,6
              AMODP(J2)=SQRT(EFWPS(J1,J2)**2+EFWPS(J1,J2+6)**2)         ! CALCUL MODULE DE EFWPS POUR LE DOF J2
              IF(ABS(EFWPS(J1,J2+6)).LE.1.E-20)EFWPS(J1,J2+6)=
     $ SIGN(1.E-20,EFWPS(J1,J2+6))
              PHASP(J2)=ATAN2(EFWPS(J1,J2),EFWPS(J1,J2+6))*RADD+180.
              PHASP(J2)=AMOD(PHASP(J2),360.)
              IF(ABS(EFWPS(J1,J2+6)).LE.1.E-20)EFWPS(J1,J2+6)=0.000000
              AMODM(J2)=SQRT(EFWMN(J1,J2)**2+EFWMN(J1,J2+6)**2)
              IF(ABS(EFWMN(J1,J2+6)).LE.1.E-20)EFWMN(J1,J2+6)=
     $ SIGN(1.E-20,EFWMN(J1,J2+6))
              PHASM(J2)=ATAN2(EFWMN(J1,J2),EFWMN(J1,J2+6))*RADD+180.
              PHASM(J2)=AMOD(PHASM(J2),360.)
              IF(ABS(EFWMN(J1,J2+6)).LE.1.E-20)EFWMN(J1,J2+6)=0.0000000
            END DO
          END DO
        ENDIF
      END DO 
      END DO
      
      WRITE(LE,*) 
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
      
      WRITE(*,*) 'QTF- CALCULEES POUR W- = 0 --> ',TABWR1(N-1)
      IF (LQTFP==1) THEN
        WRITE(*,*) 'QTF+ CALCULEES POUR W+ = ',TABWR1(2),' --> ',
     1 TABWR1(N)
      ENDIF
C********************************************************************
C -----   ECRITURE   ----
C ON ECRIT DANS 6 LES CONTRIBS
C********************************************************************
      QPM(1)='M'
      QPM(2)='P'
      QRI(1)='_PR'
      QRI(2)='_PI'
      
      CCONT(1)='HVPSI'
      CCONT(2)='HVN'
      CCONT(3)='HFLO'
      CCONT(4)='HVRN'
      CCONT(5)='HIDN'
      CCONT(6)='HFRK'
      CCONT(7)='HASBO'
      IF (LQTFP .EQ. 1) THEN
            NQTFP=2
      ELSE
            NQTFP=1
      ENDIF
      IF (Louthasbo==1) THEN
	    OUT1=1	! ON ECRIT TOUTES LES CONTRIBUTIONS
      ELSE
	    OUT1=7	! ON ECRIT SUELEMENT HASBO
      ENDIF
      DO DOF=1,6      ! ECRITURE DES RESULTATS
      DO L=1,2        ! PR PI (L=1 : PR / L=2 : PI)
      DO M=OUT1,7     ! CONTRIB
      DO K=1,NQTFP    ! QTF+- (K=1 : QTF- / K=2 : QTF+)
        IF (M<7) THEN          
          FCHA=ID(1:lID)//'/results/QTF/Appendix/QTF'//QPM(K)//'_'
     1 //TRIM(CCONT(M))//'_DOF_'//CHAR(DOF+48)//QRI(L)//'.dat'
        ELSE
          FCHA=ID(1:lID)//'/results/QTF/QTF'//QPM(K)//'_'
     1 //TRIM(CCONT(M))//'_DOF_'//CHAR(DOF+48)//QRI(L)//'.dat'
        ENDIF
        IUNI=60+L+2*(DOF-1)+12*(K-1)+24*(M-1)
        OPEN(IUNI,FILE=TRIM(FCHA),STATUS='UNKNOWN')
        
c~         WRITE(IUNI,1987)TABWR2(1),(TABWR1(I),I=2,NHASKIND)         ! version matrice carre pour les QTF+
        WRITE(IUNI,1987)TABWR2(1),(TABWR1(I),I=2,N)                   ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)

!!!! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs) !!!!
!JND : attention, critere a changer si le tri des pulsation n'est pas croissant
        DO I1=1,N                                                       
         DO I2=1,I1
                IF (K==2 .AND. I1+I2>N) THEN                      
                   CONTRIB(M,DOF+6*(L-1),I1,I2,K)=0.
                   CONTRIB(M,DOF+6*(L-1),I2,I1,K)=0.
                ENDIF
          ENDDO
          IF (K==1) THEN
                CONTRIB(M,DOF+6*(L-1),I1,I1,K)=0.
          ENDIF
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c~         DO I1=1,NHASKIND                                             ! version matrice carre pour les QTF+
        DO I1=1,N                                                       ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
          IF (L==2 .AND. K==1) THEN
            !!!  QTFij-=QTFji-* !!!
            ! ON PREND LE CONJUGE DU TRIANGLE SUP
            WRITE(IUNI,1991)TABWR1(I1)
     1 ,(CONTRIB(M,DOF+6*(L-1),I1,I2,K),I2=1,I1)
     1 ,(-CONTRIB(M,DOF+6*(L-1),I2,I1,K),I2=I1+1,N)                     ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
c~      1 ,(-CONTRIB(M,DOF+6*(L-1),I2,I1,K),I2=I1+1,NHASKIND)           ! version matrice carre pour les QTF+
          ELSE    
            WRITE(IUNI,1991)TABWR1(I1)
     1 ,(CONTRIB(M,DOF+6*(L-1),I1,I2,K),I2=1,I1)
     1 ,(CONTRIB(M,DOF+6*(L-1),I2,I1,K),I2=I1+1,N)                      ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
c~      1 ,(CONTRIB(M,DOF+6*(L-1),I2,I1,K),I2=I1+1,NHASKIND)            ! version matrice carre pour les QTF+
          ENDIF
        END DO                 
        WRITE(IUNI,*) "# FORMAT:",CHAR(10),
     1 "# 0.    W2[1]    --- W2[j]    --- W2[N]",CHAR(10),
     1 "# W1[1] QTF[1,1] --- QTF[1,j] --- QTF[1,N]",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# W1[i] QTF[i,1] --- QTF[i,j] --- QTF[i,N]",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# .      .            .            .",CHAR(10),
     1 "# W1[N] QTF[N,1] --- QTF[N,j] --- QTF[N,N]",CHAR(10),CHAR(10)
!!!! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs) !!!!
        IF (K==2) THEN
          WRITE(IUNI,*) CHAR(10),CHAR(10),"# NB: Les valeurs des QTF+",
     1 "pour w1+w2>n*dw sont definies a zeros. ", CHAR(10),
     1 "# en effet, un calcul de ces valeurs ",
     1 "#necessiterat un calcul d'ordre 1 sur 2 fois plus de pulsation."
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CLOSE(IUNI)
      END DO
      END DO
      END DO
      END DO






1991  FORMAT(1E13.6,4000(3X,1E13.6))
1987  FORMAT('0.',9X,4000(3X,1E13.6))
c~ 1988  FORMAT('QTF+',9X,4000(3X,1E13.6))

C ---- FORMATS -----
 1000 FORMAT(1X,'DEFINITION',4X,'PARTIE REAL',5X,'PARTIE IMAG.',
     1 6X,'MODE',13X,'PHASE'/
     1 1X,'FORCE DIR X :',2(E13.6,3X),2X,E13.6,3X,F8.2/
     1 1X,'FORCE DIR Y :',2(E13.6,3X),2X,E13.6,3X,F8.2/
     2 1X,'FORCE DIR Z :',2(E13.6,3X),2X,E13.6,3X,F8.2//
     3 1X,'MOMENT AUT X:',2(E13.6,3X),2X,E13.6,3X,F8.2/
     4 1X,'MOMENT AUT Y:',2(E13.6,3X),2X,E13.6,3X,F8.2/
     5 1X,'MOMENT AUT Z:',2(E13.6,3X),2X,E13.6,3X,F8.2/)
  999 FORMAT('**HASBO**  WM=',F10.4,'  W1=',F10.4,'   W2=',F10.4,1A)
 1013	FORMAT(//1X,I1,'EME CONTRIBUTION: EN MODE SOMME      (W1+W2)')
 1014	FORMAT(1X,   '                  EN MODE DIFFERENCE (W1-W2)')
 1015	FORMAT(' ')
 
      STOP
      END 
      
      FUNCTION CH(AK,Z,H)
        IF(AK*H.LE.0.OR.AK*H.GT.20)THEN
          CH=EXP(AK*Z)
        ELSE
          CH=COSH(AK*(Z+H))/COSH(AK*H)
        ENDIF
        RETURN
      END
      
      FUNCTION SH(AK,Z,H)
        IF(AK*H.LE.0.OR.AK*H.GT.20)THEN
            SH=EXP(AK*Z)
        ELSE
            SH=SINH(AK*(Z+H))/COSH(AK*H)
        ENDIF
        RETURN
      END

      SUBROUTINE PRODTN(A,B,C,D,E,F,N)
        DIMENSION A(6,4),B(6,4),C(6,4),D(6,4),E(2),F(2)
        E(1)=0.
        E(2)=0.
        F(1)=0.
        F(2)=0.
        DO I=1,N
          E(1)=E(1)+.5*(A(I,1)*B(I,3)-A(I,2)*B(I,4)+
     1                A(I,3)*B(I,1)-A(I,4)*B(I,2))
          E(2)=E(2)+.5*(A(I,1)*B(I,4)+A(I,2)*B(I,3)+
     1                A(I,3)*B(I,2)+A(I,4)*B(I,1))
          F(1)=F(1)+.5*(C(I,1)*D(I,3)+C(I,2)*D(I,4)+
     1                C(I,3)*D(I,1)+C(I,4)*D(I,2))
          F(2)=F(2)+.5*(C(I,2)*D(I,3)-C(I,1)*D(I,4)+
     1                C(I,3)*D(I,2)-C(I,4)*D(I,1))
        END DO
        RETURN
      END

      SUBROUTINE AFIN(N,SOL,ERRS)
        DIMENSION SOL(N)
        SMAX=ABS(SOL(1))
        DO I=2,N
          IF(ABS(SOL(I)).GT.SMAX)SMAX=ABS(SOL(I))
        END DO
        IF(SMAX.EQ.0)RETURN
        DO I=1,N
          IF(ABS(SOL(I))/SMAX.LT.ERRS)SOL(I)=0.
        END DO
        RETURN
      END


      FUNCTION IL(NOM)
        CHARACTER*10 NOM
        IL=INDEX(NOM,' ')-1
        IF(IL.LT.0)IL=LEN(NOM)
        IF(IL.EQ.0)THEN
          NOM='DEFAUT'
          IL=6
        ENDIF
        RETURN
      END
