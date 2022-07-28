!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!    - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)
!    - F.VILLEGER 		
!    - Fabien Robaux (EDF/INNOSEA)
!    - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------


      PROGRAM DUOK
      USE PARA
C***********************************************************************
C                                                                      *
C                      PROGRAMME  DUOK                                 *
C                      ---------------                                 *
C                                                                      *
C     CALCUL DE LA PARTIE QUADRATIQUE DES EFFORTS DU SECOND ORDRE      *
C     POUR UNE HOULE BICHROMATIQUE.                                    *
C     - 2 HOULES DIRECTIONNELLES (ANGLES EVENTUELLEMENT DIFFERENTS)    *
C     - UN OU PLUSIEURS CORPS                                          *
C     - 0 OU 1 SYMETRIE                                                *
C     - PROFONDEUR FINIE OU INFINIE                                    *
C     - FAIBLE VITESSE D'AVANCE DU CORPS                               *
C     CALCULATION OF THE QUADRATIC PART OF SECOND ORDER EFFORTS        *
C     FOR A BICHROMATIC SWELL.                                         *
C     - 2 DIRECTIONAL SWELLS (POSSIBLE DIFFERENT ANGLES)               *
C     - ONE OR MORE BODIES                                             *
C     - 0 OR 1 SYMMETRY                                                *
C     - FINITE OR INFINITE DEPTH                                       *
C     - LOW BODY Forward SPEED                                                                              *
C     G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)                                 *
C     LABORATOIRE D'HYDRODYNAMIQUE NAVALE .                            *
C     E.N.S.M. 1 rue de la Noe , 44072 Nantes CEDEX .                  *
C     Tel : 40-37-16-00                                                *
C                                                                      *
C     F.VILLEGER                                                       *
C     SIREHNA                                                          *
C     2 quai de la Joneliere , 44300 Nantes .                          *
C     Tel : 40-74-61-71                                                *
C                                                                      *
C                                                                      *
C     Version 1.1 ( 04/91 )                                            *
C                                                                      *
C***********************************************************************
C   LA CONTRIBUTION DES EFFORTS DU SECOND ORDRE CALCULE PAR DUOG       *
C CONCERNE LES PRODUITS DE TERMES DU PREMIER ORDRE POUR UNE HOULE      *
C BICHROMATIQUE (2 PULSATIONS WI ET WJ) DE DIRECTION  QUELCONQUE       *
C (LES ANGLES DES DEUX HOULES PEUVENT ETRE DIFFERENTS(PAS SUR)). LES EFFORTS    *
C CALCULES SONT DES EFFORTS DE BASSE FREQUENCE (PULSATION WI-WJ) ET    *
C DES EFFORTS DE HAUTE FREQUENCE (PULSATION WI+WJ). CE MODULE          *
C PERMET DE PRENDRE EN COMPTE DES MULTICORPS (<5), D'EXPLOITER UNE     *
C SYMETRIE, CECI POUR UNE PROFONDEUR D'EAU FINIE OU NON.IL EST POSSIBLE*
C DE TENIR COMPTE D'UNE FAIBLE VITESSE D'AVANCE AVEC LES MEMES         *
C HYPOTHESES QUE DANS AQUA+ (SELON LA METHODE UTILISEE PAR G.DELHOMMEAU*
C ET J.M.KOBUS).                                                       *
C                                                                      *
C L'ELEVATION DE LA HOULE AU POINT (XEFF,YEFF,0) SUR LA SURFACE LIBRE  *
C EST DONNEE PAR:                                                      *
C                        H=-A*SIN(W*T)                                 *
C LES CALCULS SONT EFFECTUES AVEC LES GRANDEURS DIMENSIONNELLES        *
C***********************************************************************
C THE CONTRIBUTION OF SECOND ORDER EFFORTS CALCULATED BY DUOG          *
C CONCERNS PRODUCTS OF TERMS OF THE FIRST ORDER FOR A BICHROMATIC WAVE *
C BICHROMATIC (2 PULSATIONS WI AND WJ) WITH ANY DIRECTION              *
C (THE ANGLES OF THE TWO SWELLS MAY BE DIFFERENT (NOT SURE)). THE EFFORTS       *
C CALCULATIONS ARE LOW FREQUENCY EFFORTS (WI-WJ PULSATION) AND         *
C HIGH FREQUENCY EFFORTS (WI + WJ PULSATION). THIS MODULE              *
C ALLOWS TO TAKE INTO ACCOUNT MULTIBODY (<5), TO OPERATE A             *
C SYMMETRY, THIS FOR A FINITE OR INFINITE WATER DEPTH. IT IS POSSIBLE  *
C TAKE INTO ACCOUNT A LOW FEED SPEED WITH THE SAMES                    *
C ASSUMPTIONS THAT IN AQUA + (ACCORDING TO THE METHOD USED BY G. DELHOMMEAU *
C AND J.M. KOBUS).                                                     *                                                                      *
C SWELL LEVEL AT POINT (XEFF, YEFF, 0) OVER THE OPEN SURFACE           *
C IT IS GIVEN BY:                                                      *
C H = -A * SIN (W * T)                                                 *
C THE CALCULATIONS ARE CARRIED OUT WITH THE DIMENSIONAL QUANTITIES      
C***********************************************************************
C   NPT   : NOMBRE MAXIMUM DE SOMMETS DU MAILLAGE 
C         : MAXIMUM NUMBER OF MESH TOPS
C   NFA   : NOMBRE MAXIMUM DE FACETTES
C   NFB   : NOMBRE MAXIMUM DE FACETTES DU MAILLAGE ET DE
C                                SEGMENTS A LA FLOTTAISON
C         : MAXIMUM NUMBER OF MESH FACETS AND
C                               SEGMENTS AT THE FLOTATION 
C   LN    : NOMBRE MAXIMUM DE CORPS
C           Maximum number of bodies      
C   NI    : NOMBRE MAXIMUM D'ANGLES DE HOULE
C           MAXIMUM NUMBER OF wave ANGLES
C   NC    : NOMBRE DE CORPS POUVANT OSCILLER INDEPENDAMMENT
C           NUMBER OF BODIES THAT CAN OSCILLATE INDEPENDENTLY 
C   NCO   : NOMBRE DE CORPS EFFECTIVEMENT REPRESENTES
C           NUMBER OF BODIES ACTUALLY REPRESENTED
C   NSYMY : NOMBRE DE SYMETRIE PAR RAPPORT AU PLAN XOZ (0 OU 1)
C           NUMBER OF SYMMETRY IN RELATION TO THE XOZ PLAN (0 OR 1)  
C   NP    : NOMBRE DE SOMMETS DE FACETTES
C           NUMBER OF SUMMITS OF FACETS      
C   IMX   : NOMBRE DE FACETTES
C   NFFL  : NOMBRE DE POINTS DE CONTROLE A LA FLOTTAISON (1 PAR FACETTE
C                                              POUR COMPATIBILITE AQUA+)
C  NUMBER OF FLOTATION CONTROL POINTS (1 PER FACET)
C                                              FOR AQUA + COMPATIBILITY)  
C   IXX   : NOMBRE DE POINTS DE CONTROLE SUR LE CORPS ET A LA FLOTTAISON
C         :NUMBER OF CONTROL POINTS ON THE BODY AND ON THE FLOTATION
C   XEFF,YEFF,ZEFF : COORDONNEES DU POINT DE CALCUL DES EFFORTS
C                  COORDINATES OF THE CALCULATION POINT OF Forces
C   IDEN  : IDENTIFICATEUR DU CALCUL
C           CALCULATION IDENTIFIER
C  IMXC(I): NOMBRE DE FACETTES SUR LE CORPS I
C           NUMBER OF FACETS ON THE BODY I
C  IPOS(I): 0 SI LA MOITIE DU CORPS I EST DECRITE
C           1 S'IL EXISTE UN CORPS SYMETRIQUE DU CORPS I DONT LES
C             OSCILLATIONS SONT INDEPENDANTES DES OSCILLATIONS DU CORPS I
C             (PAS DE LIAISON MECANIQUE RIGIDE ENTRE LES DEUX CORPS)
C           0 IF HALF OF BODY I IS DESCRIBED
C           1 IF THERE IS A SYMMETRICAL BODY OF BODY I WHOSE
C             OSCILLATIONS ARE INDEPENDENT OF BODY OSCILLATIONS I
C             (NO RIGID MECHANICAL LINK BETWEEN THE TWO BODIES)
C   BP    : ORDRE DE L'INTERPOLATION DANS LES FICHIERS
C           BP: ORDER OF INTERPOLATION IN THE FILES
C   TR    : PERIODE (DE RENCONTRE SI LA VITESSE D'AVANCE EST NON NULLE)
C           PERIOD (ENCOUNTER IF THE FEED SPEED IS NON ZERO)
C   IMIN  : COMPTEUR DE BOUCLE
C           LOOP COUNTER           
C   H     : PROFONDEUR DU FOND EN METRES
C           BOTTOM DEPTH IN METERS
C   AMZ   : ?NOMBRE D'ONDE?
C            NUMBER OF WAVE?
C   ILIN  : 1 POUR COEFFICIENTS HYDROSTATIQUES LINEARISES
C           1 FOR LINEARIZED HYDROSTATIC COEFFICIENTS
C   RHO   : DENSITE DU FLUIDE
C   VA    : VITESSE D'AVANCE
C           Forward SPEED
C   BETAI : ANGLE DE LA HOULE I (I=1 OU 2) AVEC L'AXE OX
C           WAVE ANGLE I (I = 1 OR 2) WITH OX AXIS
C   WRI   : PULSATION (DE RENCONTRE SI VA#0) DE LA HOULE I (I=1 OU2)
C           frequency (ENCOUNTER IF VA not 0) OR WAVE I (I = 1 OR 2)
C   WI    : PULSATION DANS LE REPERE ABSOLU DE LA HOULE I
C           FREQUENCY IN THE ABSOLUTE WAVE I
C   AM0   : NOMBRE D'ONDE
C           NUMBER OF WAVE      
C   NIN   : NOMBRE D'ANGLES DE HOULE
C           NUMBER OF WAVE DIRECTION
C   AIND  : VALEURS DES ANGLES DE LA HOULE
c           SWELL angle value
C   X,Y,Z : COORDONNEES DES SOMMETS DES FACETTES
C           COORDINATES OF THE SUMMITS OF THE FACETS
C   XM,YM,ZM: COORDONNEES DES CENTRES DES FACETTES ET COORDONNEES DES
C             CENTRES DES FACETTES SITUEES A LA FLOTTAISON
C             COORDINATES OF FACET CENTERS AND COORDINATES OF 
C             FACET CENTERS LOCATED AT FLOTATION
C   AIRE  : AIRE DES FACETTES ET LONGUEUR DES SEGMENTS SITUES A LA
C                                                        FLOTTAISON
C           AREA OF THE FACETS AND LENGTH OF THE SEGMENTS LOCATED AT THE
C                                                        FLOTATION
C  (P,Q,R): NORMALES DES FACETTES ET DES FACETTES DE LA FLOTTAISON
C           NORMALS OF FLOTATION FACETS AND FACETS
C      
      COMMON NC,NCO,NSYMY,NP,IMX,NFFL,IXX,XEFF,YEFF,ZEFF,
     #RHO,TPER,BETA1,BETA2
      COMMON X(NPT),Y(NPT),Z(NPT),M1(NFA),M2(NFA)
      COMMON M3(NFA),M4(NFA),P(NFA),Q(NFA),R(NFA)
      COMMON XM(2*NFA),YM(2*NFA),ZM(2*NFA),AIRE(NFA),ALF(NFA)
      COMMON IPOS(LN),IMXC(LN)
      DIMENSION VXR1(4,NFB),VXM1(4,NFB),VYR1(4,NFB),VYM1(4,NFB)
      DIMENSION VZR1(4,NFB),VZM1(4,NFB),VYR2(4,NFB),VYM2(4,NFB)
      DIMENSION VXR2(4,NFB),VXM2(4,NFB),VZR2(4,NFB),VZM2(4,NFB),FKI(6)
      DIMENSION VGR1(4,NFA),VGM1(4,NFA),VGR2(4,NFA),VGM2(4,NFA),FKR(6)
      DIMENSION CM1(6,6,LN,LN),CA1(6,6,LN,LN)
      DIMENSION CM2(6,6,LN,LN),CA2(6,6,LN,LN)
      DIMENSION TABWR1(NFA),TABWR2(NFA)
      DIMENSION AINTV(6,6),AHYD(6,6),CN(NFB,6)
      DIMENSION FINT1(12),FINT2(12),PRO1(6,4),PRO2(6,4),FP(2),FM(2)
      REAL,DIMENSION(1:6,1:6):: AINTR1,AINTM1,AINTR2,AINTM2
      REAL,DIMENSION(1:6,1:6):: AINTINERTIER,AINTINERTIEM
      REAL,DIMENSION(1:12):: FINT11,FINT21,FINT12,FINT22
      REAL,DIMENSION(1:12):: FHY1,FHY2
      DIMENSION EFWPS(10,12,LN),EFWMN(10,12,LN),K1(6),K2(6),K3(6),K4(6)
      DIMENSION K5(4),K6(4),K7(4),X8(4),PRPR(6),PRPI(6)
      DIMENSION AMODP(6),AMODM(6),PHASP(6),PHASM(6)
      DIMENSION SOLP(12),SOLM(12),A1(12),A2(12),NT(NPIN),LK(5)
      DIMENSION AIND(NI),IND(NFA),DIST(NFA),TDIS(NFA)
      DIMENSION AIN(6,6,LN),ASH(6,6,LN),AMOR(6*LN,6*LN),RAID(6*LN,6*LN)
      REAL,DIMENSION(1:6,1:2,1:NPIN,1:NPIN,1:2,1:7):: FCONT
      DIMENSION EFFP(12,NPIN,NPIN,3)
      DIMENSION FH(NPIN),WWX(NPIN)
      COMPLEX FOR1(6,LN),FOR2(6,LN)
      COMPLEX ZA1(6,LN),ZA2(6,LN),ZAX(6,LN) ! ,ZET,ZZZ,ZF2 INNOSEA UNUSED !
      COMPLEX ZZA(6),ZZX(6,LN,NPIN)
      CHARACTER*1 QPM(2)
      CHARACTER*3 QRI(2)
      CHARACTER*10 ID,IDENB,CCONT(8)
      CHARACTER*16 NOMF4,NOMB1,NOMB2,NOMB3,NOMB4,NOMB5,NOMB6
      CHARACTER*60 FCHA
      LOGICAL FICHEX,FICHEY
      DATA K1/2,3,1,5,6,4/
      DATA K2/3,1,2,6,4,5/
      DATA K3/6,4,5,6,4,5/
      DATA K4/5,6,4,5,6,4/
      DATA K5/4,5,5,4/
      DATA K6/4,5,6,6/
      DATA K7/3,3,4,5/
      DATA X8/1,1,-1,1/
      INTEGER lID,FICH,CHOIXAPPROX,Loutduok !CHOIX D'APPROXIMATION DE PINKSTER : 1POUR FROUDE-KRYLOV   2:POUR FEX
      CHARACTER*1 FMT,FMT1
      REAL,DIMENSION(1:6,1:NFA,1:NFA)::QTF1,QTFAP !LES QTF APPROX PINKSTER AVEC DUOK
      REAL,DIMENSION(1:2,1:NFA,1:NFA)::Q1,Q2,Q3,Q4,Q5,Q8 !SEPARES (+)
      REAL,DIMENSION(1:6,1:NFA)::FEX,FEXPH
      REAL,DIMENSION(1:NFA)::PULS
      INTEGER REE,FICH2,DOF,IUNI,NQTFP,LQTFP,LHASK,OUT1,DOFSYM,N1,NdofIt
      REAL,DIMENSION(1:6,1:12,1:NFA,1:NFA)::CONTRIB
      REAL PHIR,PHII,C1H,COST




!#######################################################################
      ISTOK=2           !CHOIX DES CONTRIBUTIONS AUX EFFORTS CALCULEES
      IF(ISTOK.LT.1 .OR.ISTOK.GT.3)ISTOK=1
      IECR1=1           !ECRITURE DES EFFORTS PAR RAPPORT A UN POINT FIXE (1:OUI)
	!					 WRITING OF EFFORTS IN RELATION TO A FIXED POINT (1: YES)
      IECR2=0           !ECRITURE DES DIFFERENTES CONTRIBUTIONS AUX EFFORTS (1:OUI)
	!					 WRITING OF THE VARIOUS CONTRIBUTIONS TO THE EFFORTS (1: YES)
      IMZERO=0          !MISE A ZERO DES COMPOSANTES TROP PETITES (1:OUI)
	!					 ZEROING OF COMPONENTS TOO SMALL (1: YES)
      IECRBF=1          !ECRITURE FONCTION DE TRANSFERT ADIMENSIONNELLE BASSE FREQUENCE
	!					 WRITE LOW FREQUENCY ADIMENSIONAL TRANSFER FUNCTION    
	  IECRHF=0          !ECRITURE FONCTION DE TRANSFERT ADIMENSIONNELLE HAUTE FREQUENCE
	!					WRITING HIGH FREQUENCY ADIMENSIONAL TRANSFER FUNCTION
      OPEN(7,FILE="ID.dat")
      READ(7,*) lID
      READ(7,*) ID
      CLOSE(7)
      
C *** INITIALISATION DES VARIABLES ***
      NOMF4=ID(1:lID)//'/QTF/FA.RES'		 !GEOMETRY data,Hydrostatic stifness matric etc from QTFpreprocessing
      INQUIRE(FILE=NOMF4,EXIST=FICHEY)
      OPEN(UNIT=12,FILE=NOMF4,ACCESS='DIRECT',
     #STATUS='UNKNOWN',RECL=4*4*NFA)
      ! DEFINITION EN PREAMBULE				! PREAMBLE DEFINITION
      READ(12,REC=1)NC,NCO,NSYMY,NP,IMX,IXX,XEFF,YEFF,ZEFF,
     #(IMXC(I),I=1,NCO),(IPOS(I),I=1,NC),BP,TPER,IMIN,H,AMZ,ILIN,
     #RHO,VA,BETA,WR,AM0,NIN,(AIND(I),I=1,NIN)

      NFFL=IXX-IMX			! number of floatation control point

      NJ=NSYMY+1
      DOFSYM=1
      IF (NJ==2 .AND. BETA==0) DOFSYM=2
      XL=1.             !LONGUEUR D'ADIMENSIONNALISATION ! ADIMENSIONALIZATION LENGTH
      AD=RHO*G*XL

C *** LECTURE DE LA GEOMETRIE ***
      ! X,Y,Z           : COORDONNEES DU Ieme POINT
      ! Mi (i=1:4)      : INDICE DU ieme POINT DE LA Ieme FACETTE
      ! P,Q,R           : COORDONNEES DE LA NORMALE A LA Ieme FACETTE
      ! XM,YM,ZM        : COORDONNEES DU CENTRE DE LE Ieme FACETTE
      ! AIRE            : SURFACE DE LA Ieme FACETTE
      ! TDIS            : ?  non utilisee
      ! DIST            : ?  non utilisee
      ! XM,YM           : COORDONNEES DU Ieme SEGMENT A LA FLOTTAISON
      ! ALF             : LONGUEUR DU Ieme SEGMENT A LA FLOTTAISON ?
      ! IND             : INDICE DE LA Ieme FACETTE A LA FLOTTAISON
      ! AIN             : MATRICE D'INERTIE
      ! ASH             : MATRICE DE RAIDEUR HYDROSTATIQUE
      ! RAID            : RAIDEUR D'ANCRAGE (EXT)
      ! AMOR            : AMORTISSEMENT VISQUEUX (EXT)
C *** READING THE GEOMETRY ***
       ! X, Y, Z: COORDINATES OF THE ITH POINT
       ! Mi (i = 1: 4): INDEX OF THE iTH POINT OF THE ITH FACET
       ! P, Q, R: COORDINATES FROM NORMAL TO ITH FACET
       ! XM, YM, ZM: COORDINATES OF THE CENTER OF THE ITH FACET
       ! AIRE: AREA SURFACE OF THE ITH FACET
       ! TDIS :? not used
       ! DIST:? not used
       ! XM, YM: COORDINATES OF THE ITH SEGMENT AT FLOTATION
       ! ALF: LENGTH OF THE ITH SEGMENT AT THE FLOTATION?
       ! IND: INDEX OF THE ITH FACET AT FLOTATION
       ! AIN: INERTIA MATRIX
       ! ASH: HYDROSTATIC STIFFNESS MATRIX
       ! RAID: ANCHOR STIFFNESS (EXT)
       ! AMOR: VISCOUS CUSHIONING (EXT)           
      
      READ(12,REC=2)(X(I),I=1,NP)
      READ(12,REC=3)(Y(I),I=1,NP)
      READ(12,REC=4)(Z(I),I=1,NP)
      READ(12,REC=5)(M1(I),I=1,IMX),(M2(I),I=1,IMX),
     #(M3(I),I=1,IMX),(M4(I),I=1,IMX)
      READ(12,REC=6)(P(I),I=1,IMX),(Q(I),I=1,IMX),(R(I),I=1,IMX)   !is unit normal vector as prepared in QTFpreProcessor
      READ(12,REC=7)(XM(I),I=1,IMX),(YM(I),I=1,IMX),(ZM(I),I=1,IMX)
      READ(12,REC=8)(AIRE(I),I=1,IMX),(TDIS(I),I=1,IMX),
     #(DIST(I),I=1,IMX)
      I1=IMX+1
      READ(12,REC=9)(XM(I),I=I1,IXX),(YM(I),I=I1,IXX),
     #(ALF(I),I=1,IXX-IMX),(IND(I),I=1,IXX-IMX)
      JQ=6*NC
      READ(12,REC=10)(((AIN(I,J,K),I=1,6),J=1,6),
     #((ASH(I,J,K),I=1,6),J=1,6),K=1,NC),((RAID(I,J),I=1,JQ),J=1,JQ),
     #((AMOR(I,J),I=1,JQ),J=1,JQ)
      CLOSE(UNIT=12)
      DO I=1,NFFL
            LK(1)=M1(IND(I))              ! INDICE DU 1ER PT DE LA FCETTE #IND
            LK(2)=M2(IND(I))              ! "         2EME "
            LK(3)=M3(IND(I))              ! "         3EME "
            LK(4)=M4(IND(I))
            LK(5)=LK(1)
      
            P(IMX+I)=P(IND(I))            ! NORMALE A LA FLOTTAISON (NX)
            Q(IMX+I)=Q(IND(I))            ! NORMALE A LA FLOTTAISON (NY)
            R(IMX+I)=0                    ! NORMALE A LA FLOTTAISON (NZ=0 --> FACETTE VERTICALE)
            ZM(IMX+I)=0.
            DO J=1,4
                  TEST=(X(LK(J))-X(LK(J+1)))**2+(Y(LK(J))-Y(LK(J+1)))**2
                  IF(ABS(Z(LK(J)))+ABS(Z(LK(J+1))).LT. 1.E-05
     1 .AND.TEST.GT. 1.E-5) THEN
C REMARQUE: POUR PLUS DE PRECISION DANS LE CALCUL DE LA CONTRIBUTION
C DE L'INTEGRALE DE LIGNE, IL PEUT ETRE INTERESSANT DE PLACER DEUX POINTS
C DE CONTROLE (OU PLUS) SUR CHAQUE SEGMENT A LA FLOTTAISON. IL FAUT ALORS
C UTILISER DES VERSIONS ADAPTEES DE LECK ET DE MOUK(EVENTUELLEMENT). LA
C LONGUEUR DES SEGMENTS A LA FLOTTAISON EST ALORS DIVISE PAR 2.
C NOTE: FOR MORE PRECISION IN CALCULATING THE CONTRIBUTION
C OF THE INTEGRAL LINE, IT MAY BE INTERESTING TO PLACE TWO POINTS
C OF CONTROL (OR MORE) ON EACH SEGMENT AT FLOTATION. THEN YOU MUST
C USE ADAPTED VERSIONS OF LECK AND MOUK (IF POSSIBLE). THE
C LENGTH OF THE SEGMENTS AT THE FLOTATION IS THEN DIVIDED BY 2.
                        AIRE(IMX+I)=SQRT(TEST)
                  ENDIF
            END DO
      END DO
      
      DO I=1,IXX  ! NORMALE GENERALISEE
            IF(I.LE.IMX) THEN  ! since P,Q,R are already unit vector so no effect by normalizing again
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
C *** OUVERTURE DES FICHIERS QUANTITES DU 1ER ORDRE ***
C *** OPENING OF FIRST ORDER QUANTITY FILES ***
      NOMB1=ID(1:lID)//'/QTF/B1.RES' !! TR,beta,AM0, ZF,ZA,CM,CA,VGR,VGM
      NOMB2=ID(1:lID)//'/QTF/B2.RES' !! VXR(1:NJJ,1:IXX) VXM(1:NJJ,1:IXX)
      NOMB3=ID(1:lID)//'/QTF/B3.RES' !! VYR VYM
      NOMB4=ID(1:lID)//'/QTF/B4.RES' !! VZR VZM
C FV NCO AU LIEU DE NC!!??
C FV NCO INSTEAD OF NC !! ??     
      NR1=3+24*NC+72*NC**2+NFFL*2*NJ
      NR3=IXX*2*NJ
      OPEN(21,FILE=NOMB1,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR1)
      OPEN(22,FILE=NOMB2,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
      OPEN(23,FILE=NOMB3,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
      OPEN(24,FILE=NOMB4,STATUS='UNKNOWN',ACCESS='DIRECT',RECL=4*NR3)
C *** OUVERTURE FICHIER ECRITURE FONCTION DE TRANSFERT ***
C *** OPEN FILE WRITE TRANSFER FUNCTION ***
C BASSE FREQUENCE
C LOW FREQUENCY
      IF(IECRBF.EQ.1) THEN
            NOMB5=ID(1:lID)//'/QTF/A1.RES'
            OPEN(UNIT=26,FILE=NOMB5,STATUS='UNKNOWN')
      ENDIF
C HAUTE FREQUENCE
C HIGH FREQUENCY
      IF(IECRHF.EQ.1) THEN
            NOMB6=ID(1:lID)//'/QTF/A2.RES'
            OPEN(UNIT=27,FILE=NOMB6,STATUS='UNKNOWN')
      ENDIF

c~       WRITE(LE,*)' NOMBRE D''ENREGISTREMETS N=?'
      READ(LL,*)N
c~       WRITE(LE,*)' ENTRER LES NUMEROS D''ENREGISTREMENTS DANS L''ORDRE'
c~       WRITE(LE,*)' DES PULSATIONS CROISSANTES (PERIODES DECROISSANTES)'
      READ(LL,*) (NT(I),I=1,N) 
      READ(LL,*) LQTFP
      READ(LL,*) LHASK
      READ(LL,*) Loutduok
c~       PRINT*, LHASK>1,LQTFP==1,LHASK>1 .AND. LQTFP==1
!!!!!!!!!!!!!! version matrice carre pour les QTF+!!!!!!!!!!!!!!!!!!!!!!
!~        IF (LHASK>1 .AND. LQTFP==1) THEN                               
! SI LES QTF+ SONT CALCULES, IL FAUT REDUIRE LE NOMBRE DE PULSATION
!~              NHASKIND=N/2
!~        ELSE
!~              NHASKIND=N
!~        ENDIF
!~        PRINT*, 'N =',N,'NHASKIND =',NHASKIND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      OPEN(99,FILE=ID(1:lID)//'/QTF/WRH.dat',STATUS='UNKNOWN')
c~       WRITE(LE,2003)RHO
c~       WRITE(LE,*)'  VITESSE D''AVANCE :',VA
c~       WRITE(LE,*)'  POINT DE CALCUL DES EFFORTS(XEFF,YEFF,ZEFF):',
c~      #XEFF,YEFF,ZEFF
c~       DO I=1,NHASKIND    ! initialisation des sorties                ! version matrice carre pour les QTF+
					! initialisation output                            ! antitriangular matrix version sup for QTF + (0 elsewhere)
      DO I=1,N    ! initialisation des sorties                          ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
            READ(21,REC=NT(I))TRB1,BETAB1,AMB1
            FH(I)=2*PI/TRB1
            
            DO DOF=1,6
            DO L=1,2          ! PARTIE REEL / IMAGINAIRE
                  DO K=1,2    ! QTF-+
                  DO M=1,7    ! CONTRIB
                        FCONT(DOF,K,1,I+1,L,M)=FH(I)
                  END DO
                  END DO
                  DO K=1,3    ! APPROXs PINKSTER
                        EFFP(DOF+6*(L-1),1,I+1,K)=FH(I)
                  END DO
            END DO
            END DO
            
      END DO
      
      WRITE(*,*) 'CALCUL DES QTF- POUR W- = 0 --> ',FH(N-1)
      IF (LQTFP==1) THEN
            WRITE(*,*) 'CALCUL DES QTF+ POUR W+ = ',FH(2),' --> ',FH(N)
      ENDIF
      
      DO DOF=1,6  
      DO L=1,2
            DO K=1,2    ! QTF-+
            DO M=1,7
                  FCONT(DOF,K,1,1,L,M)=0.
            END DO
            END DO
            DO K=1,3    ! APPROXs PINKSTRER
                  EFFP(DOF+6*(L-1),1,1,K)=0.
            END DO
      END DO
      END DO
      
c~       DO I=1,NHASKIND                                                ! version matrice carre pour les QTF+
      DO I=1,N                                                          ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
            READ(21,REC=NT(I))TR,BETA,AM,((FOR1(J,II),J=1,6),II=1,NC),
     1 ((ZAX(J,II),J=1,6),II=1,NC)
            WWX(I)=2*PI/TR
            DO II=1,NC
            DO J=1,6
                  ZZX(J,II,I)=ZAX(J,II)
            END DO
            END DO
      END DO
   
      ! lecture des forces d'excitation pour pinkster
      OPEN(4444,FILE=ID(1:lID)//'/results/Fe.dat')
      READ(4444,*) 
      READ(4444,*)
      READ(4444,*)
c~       DO I=1,NHASKIND                                                ! version matrice carre pour les QTF+
      DO I=1,N                                                          ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
!   TODO: CHECK PHASES CONVENTION !!!!!!!!!!!  
            READ(4444,*) PULS(I),(FEX(J,I),J=1,6),(FEXPH(J,I),J=1,6)  ! FEXPH expected in degrees  
!!! A.C: WARNING : change to AQ+ phase convention !!!!!
            DO J=1,6
                  FEXPH(J,I) = FEXPH(J,I)/RADD - PI/2.0 
            END DO
      END DO
      CLOSE(4444)      
      
C ______________________________________________
C *** COMBINAISON DES DIFFERENTES PULSATIONS ***
C ----------------------------------------------
C NDD : ECART ENTRE LES PULSATIONS
c~       DO NDD=0,NHASKIND-1                                            ! version matrice carre pour les QTF+
      DO NDD=0,N-1                                                      ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
      WRITE(LE,'(A,I3,A)',ADVANCE="NO")
     & 'ECART ENTRE LES PULSATIONS:', NDD,CHAR(13)

c~       DO I1=NDD+1,NHASKIND                                           ! version matrice carre pour les QTF+
      DO I1=NDD+1,N                                                     ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
      I2=I1-NDD
      ! TRi       : periode de la houle   (s)
      ! BETAi     : angle d'incidence de la houle (rad)
      ! AMi       : nombre d'onde
      ! FORi      : force excitatrice d'ordre 1
      ! ZAi       : RAO complexe
      ! CMi       : masse ajoutee
      ! CAi       : amortissement ajoute
      ! VGci      : potentiel des vitesse sur la Ieme facette
      ! VXci,VYci,VZci : vitesse du fluide sur la Ieme facette
      !
      ! 
      ! i=1,2     : indice de la pulsation du couple etudie
      ! c=R,M     : partie reelle / imaginaire
       ! TRi: period of the swell (s)
       ! BETAi: angle of incidence of the swell (rad)
       ! AMi: wave number
       ! FORi: excitatory force of order 1
       ! ZAi: complex RAO
       ! CMi: added mass
       ! CAi: damping adds
       ! VGci: speed potential on the Ith facet
       ! VXci, VYci, VZci: velocity of the fluid on the Ith facet
       !
       !
       ! i = 1,2: index of the pulsation of the couple studied
       ! c = R, M: real / imaginary part 
      READ(21,REC=NT(I1))TR1,BETA1,AM1,                                 
     1  ((FOR1(J,II),J=1,6),II=1,NC),
     1  ((ZA1(J,II),J=1,6),II=1,NC),
     1  ((((CM1(I,J,K,L),J=1,6),L=1,NC),
     1  ((CA1(I,J,K,L),J=1,6),L=1,NC),I=1,6),K=1,NC),
     1  ((VGR1(IC,J),J=1,NFFL),IC=1,NJ),((VGM1(IC,J),J=1,NFFL),IC=1,NJ)
      READ(22,REC=NT(I1))((VXR1(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VXM1(IC,J),J=1,IXX),IC=1,NJ)
      READ(23,REC=NT(I1))((VYR1(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VYM1(IC,J),J=1,IXX),IC=1,NJ)
      READ(24,REC=NT(I1))((VZR1(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VZM1(IC,J),J=1,IXX),IC=1,NJ)
      WR1=2*PI/TR1
      TABWR1(I1)=WR1
      BETA1D=BETA1*RADD
      READ(21,REC=NT(I2))TR2,BETA2,AM2,
     1  ((FOR2(J,II),J=1,6),II=1,NC),
     1  ((ZA2(J,II),J=1,6),II=1,NC),
     1  ((((CM2(I,J,K,L),J=1,6),L=1,NC),
     1  ((CA2(I,J,K,L),J=1,6),L=1,NC),I=1,6),K=1,NC),
     1  ((VGR2(IC,J),J=1,NFFL),IC=1,NJ),((VGM2(IC,J),J=1,NFFL),IC=1,NJ)
      READ(22,REC=NT(I2))((VXR2(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VXM2(IC,J),J=1,IXX),IC=1,NJ)
      READ(23,REC=NT(I2))((VYR2(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VYM2(IC,J),J=1,IXX),IC=1,NJ)
      READ(24,REC=NT(I2))((VZR2(IC,J),J=1,IXX),IC=1,NJ),
     1                       ((VZM2(IC,J),J=1,IXX),IC=1,NJ)
      WR2=2*PI/TR2
      TABWR1(I2)=WR2
      BETA2D=BETA2*RADD
C *** INITIALISATION DES EFFORTS ***
      DO KNC=1,NC
      DO I=1,10
      DO J=1,12
            EFWPS(I,J,KNC)=0.
            EFWMN(I,J,KNC)=0.
      END DO
      END DO
      END DO
      W1=WR1+VA*AM1*COS(BETA1)
      T1=2*PI/W1
      C1=AM1**2*COS(BETA1)/W1
      W2=WR2+VA*AM2*COS(BETA2)
      C2=AM2**2*COS(BETA2)/W2
      T2=2*PI/W2
C      DO I=1,NFFL
C        print*,I,VGR1(1,I),VGM1(1,I)
C        print*,I,VGR1(2,I),VGM1(2,I)
C      ENDDO
C      DO I=1,IXX
C        print*,I,VXR1(1,I),VXM1(1,I),VXR1(2,I),VXM1(2,I)
C        print*,I,VYR1(1,I),VYM1(1,I),VYR1(2,I),VYM1(2,I)
C        print*,I,VZR1(1,I),VZM1(1,I),VZR1(2,I),VZM1(2,I)
C      ENDDO

C
C *** ITERATION SUR LE NOMBRE DE CORPS INDEPENDANTS DECRITS ***
C -------------------------------------------------------------
      DO KNC=1,NCO
      IF(NC.GT.1) THEN
            WRITE(LE,*)' *** CORPS N0 ',KNC,' ***'
      ENDIF
      IA=1
      IF(KNC.NE.1)IA=IMXC(KNC-1)+1
      IB=IMXC(KNC)
      NMI=1
      NMA=NJ
      IF(IPOS(KNC).GT.0) THEN
            IF(KNC.LE.NCO) THEN
                  NMA=1
            ELSE
C IA=1 A VOIR!! IA=IMXC(KNC-1)+1
C IA = 1 SEE !! IA = IMXC (KNC-1) +1
                  IA=1
C EST IL POSSIBLE QUE IPOS(KNC)>1, NON !
                  IF(IPOS(KNC).NE.1)IA=IMXC(IPOS(KNC)-1)+1
                  IB=IMXC(IPOS(KNC))
                  NMI=2
                  NMA=2
            ENDIF
      ENDIF
C
C *** MOUVEMENTS DU CORPS N0 KNC ***
      DO J=1,6
            A1(J)=REAL(ZA1(J,KNC))
            A1(J+6)=AIMAG(ZA1(J,KNC))
            A2(J)=REAL(ZA2(J,KNC))
            A2(J+6)=AIMAG(ZA2(J,KNC))
      END DO
C
C --------------------------------------------------------
C *** CALCUL DES DIFFERENTES CONTRIBUTIONS AUX EFFORTS ***
C --------------------------------------------------------
C
C *** CONTRIBUTION DUE A LA VITESSE QUADRATIQUE ***
C ADV
C -------------------------------------------------
      DO JJ=NMI,NMA
      DO II=IA,IB
            PRO1(1,1)=VXR1(JJ,II)
            PRO1(2,1)=VYR1(JJ,II)
            PRO1(3,1)=VZR1(JJ,II)
            PRO1(1,2)=VXM1(JJ,II)
            PRO1(2,2)=VYM1(JJ,II)
            PRO1(3,2)=VZM1(JJ,II)
      
            PRO1(1,3)=VXR2(JJ,II)
            PRO1(2,3)=VYR2(JJ,II)
            PRO1(3,3)=VZR2(JJ,II)
            PRO1(1,4)=VXM2(JJ,II)
            PRO1(2,4)=VYM2(JJ,II)
            PRO1(3,4)=VZM2(JJ,II)
            DO J1=1,3
            DO J2=1,4
                  PRO2(J1,J2)=PRO1(J1,J2)
            END DO
            END DO
            CALL PRODT(PRO1,PRO2,FP,FM,3) !CALCUL NOUS DONNE FP=V1.V2 ET FM=1/2(V1*.V2+V1.V2*)
            DO J=1,6
                  ACQ=AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))/2.*RHO 
!! CN(II,J): JEME COORDONNEES DE LA NORMALE A LA IIEME FACETTE, ON BOUCLE EN CHANGEANT LE SIGNE DE CNY DANS LA DEUXIEME BOUCLE SUR LA SYMETRIE
!! EFWPS(NO CONTRIB , DOF(+6 POUR PARTIE IMAG) , NO CORPS)
                  EFWPS(1,J,KNC)=EFWPS(1,J,KNC)+FP(1)*ACQ               !! INTEGRATION CORRECTE
                  EFWPS(1,J+6,KNC)=EFWPS(1,J+6,KNC)+FP(2)*ACQ           !!
                  EFWMN(1,J,KNC)=EFWMN(1,J,KNC)+FM(1)*ACQ               !!
                  EFWMN(1,J+6,KNC)=EFWMN(1,J+6,KNC)+FM(2)*ACQ           !!

              ! IF (J==1) THEN
              ! !print*,II,AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))
              !  print*,II,FM*ACQ,FP*ACQ
              ! ENDIF
            END DO
                   !print*,II,FM,FP
      END DO
      END DO
      Q1(1,I1,I2)=EFWMN(1,3,KNC)/2
      Q1(2,I1,I2)=EFWMN(1,9,KNC)/2
      !print*,NMI,NMA
      !print*,I1,I2,EFWMN(1,1,KNC),EFWMN(1,7,KNC)
C
C *** CONTRIBUTION DUE AU PRODUIT DU GRADIENT DE PRESSION ET
C                                                    DU MOUVEMENT ***
C DPV
C -------------------------------------------------------------------
      DO JJ=NMI,NMA
      DO II=IA,IB
            XMII=XM(II)-XEFF
            YMII=(YM(II)-YEFF)*((-1.)**(JJ+1))
            ZMII=ZM(II)-ZEFF
            PRO1(1,1)=A1(1)+(A1(5)*ZMII-A1(6)*YMII) !MoM=L+\theta \cross r_0
            PRO1(2,1)=A1(2)+(A1(6)*XMII-A1(4)*ZMII)
            PRO1(3,1)=A1(3)+(A1(4)*YMII-A1(5)*XMII)
            PRO1(1,2)=A1(7)+(A1(11)*ZMII-A1(12)*YMII)
            PRO1(2,2)=A1(8)+(A1(12)*XMII-A1(10)*ZMII)
            PRO1(3,2)=A1(9)+(A1(10)*YMII-A1(11)*XMII)
      
            PRO1(1,3)=A2(1)+(A2(5)*ZMII-A2(6)*YMII)
            PRO1(2,3)=A2(2)+(A2(6)*XMII-A2(4)*ZMII)
            PRO1(3,3)=A2(3)+(A2(4)*YMII-A2(5)*XMII)
            PRO1(1,4)=A2(7)+(A2(11)*ZMII-A2(12)*YMII)
            PRO1(2,4)=A2(8)+(A2(12)*XMII-A2(10)*ZMII)
            PRO1(3,4)=A2(9)+(A2(10)*YMII-A2(11)*XMII)
      
            PRO2(1,1)=VXM1(JJ,II)*WR1		!\partial_t \grad\Phi=-i\omega \grad \Phi
            PRO2(2,1)=VYM1(JJ,II)*WR1
            PRO2(3,1)=VZM1(JJ,II)*WR1
            PRO2(1,2)=-VXR1(JJ,II)*WR1
            PRO2(2,2)=-VYR1(JJ,II)*WR1
            PRO2(3,2)=-VZR1(JJ,II)*WR1
      
            PRO2(1,3)=VXM2(JJ,II)*WR2
            PRO2(2,3)=VYM2(JJ,II)*WR2
            PRO2(3,3)=VZM2(JJ,II)*WR2
            PRO2(1,4)=-VXR2(JJ,II)*WR2
            PRO2(2,4)=-VYR2(JJ,II)*WR2
            PRO2(3,4)=-VZR2(JJ,II)*WR2
      
            CALL PRODT(PRO1,PRO2,FP,FM,3)
            DO J=1,6
                  ACQ=AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))*RHO
                  EFWPS(2,J,KNC)=EFWPS(2,J,KNC)+FP(1)*ACQ
                  EFWPS(2,J+6,KNC)=EFWPS(2,J+6,KNC)+FP(2)*ACQ
                  EFWMN(2,J,KNC)=EFWMN(2,J,KNC)+FM(1)*ACQ
                  EFWMN(2,J+6,KNC)=EFWMN(2,J+6,KNC)+FM(2)*ACQ
              ! IF (J==1) THEN
              ! !print*,II,AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))
              ! print*,II,FM*ACQ,FP*ACQ
              ! ENDIF
            END DO
      END DO
      END DO
      Q2(1,I1,I2)=EFWMN(2,3,KNC)/2
      Q2(2,I1,I2)=EFWMN(2,9,KNC)/2

C     print*,I1,I2,EFWMN(2,1,KNC),EFWMN(2,7,KNC)
C
C *** CONTRIBUTION DE LA VITESSE D'AVANCE DANS LE PRODUIT
C                               GRADIENT DE VITESSE ET MOUVEMENT ***
C -------------------------------------------------------------------
C PAS PRIS EN COMPTE POUR LE MOMENT
      IVA=0
      IF(IVA.NE.0) THEN
      IF(ABS(VA).GE..05) THEN
      WRITE(LE,*)'PROFONDEUR H',H
      DO JJ=NMI,NMA
      DO II=IA,IB
            XMII=XM(II)-XEFF
            YMII=(YM(II)-YEFF)*((-1.)**(JJ+1))
            ZMII=ZM(II)-ZEFF
            PRO1(1,1)=A1(1)+(A1(5)*ZMII-A1(6)*YMII)
            PRO1(2,1)=A1(2)+(A1(6)*XMII-A1(4)*ZMII)
            PRO1(3,1)=A1(3)+(A1(4)*YMII-A1(5)*XMII)
            PRO1(1,2)=A1(7)+(A1(11)*ZMII-A1(12)*YMII)
            PRO1(2,2)=A1(8)+(A1(12)*XMII-A1(10)*ZMII)
            PRO1(3,2)=A1(9)+(A1(10)*YMII-A1(11)*XMII)
      
            PRO1(1,3)=A2(1)+(A2(5)*ZMII-A2(6)*YMII)
            PRO1(2,3)=A2(2)+(A2(6)*XMII-A2(4)*ZMII)
            PRO1(3,3)=A2(3)+(A2(4)*YMII-A2(5)*XMII)
            PRO1(1,4)=A2(7)+(A2(11)*ZMII-A2(12)*YMII)
            PRO1(2,4)=A2(8)+(A2(12)*XMII-A2(10)*ZMII)
            PRO1(3,4)=A2(9)+(A2(10)*YMII-A2(11)*XMII)
      
            E1=AM1*(XMII*COS(BETA1)+YMII*SIN(BETA1))
            PRO2(1,1)=-C1*CH(AM1,ZMII,H)*COS(BETA1)*COS(E1)
            PRO2(2,1)=-C1*CH(AM1,ZMII,H)*SIN(BETA1)*COS(E1)
            PRO2(3,1)=-C1*SH(AM1,ZMII,H)*SIN(E1)
            PRO2(1,2)=-C1*CH(AM1,ZMII,H)*COS(BETA1)*SIN(E1)
            PRO2(2,2)=-C1*CH(AM1,ZMII,H)*SIN(BETA1)*SIN(E1)
            PRO2(3,2)=C1*SH(AM1,ZMII,H)*COS(E1)
      
            E2=AM2*(XMII*COS(BETA2)+YMII*SIN(BETA2))
            PRO2(1,3)=-C2*CH(AM2,ZMII,H)*COS(BETA2)*COS(E2)
            PRO2(2,3)=-C2*CH(AM2,ZMII,H)*SIN(BETA2)*COS(E2)
            PRO2(3,3)=-C2*SH(AM2,ZMII,H)*SIN(E2)
            PRO2(1,4)=-C2*CH(AM2,ZMII,H)*COS(BETA2)*SIN(E2)
            PRO2(2,4)=-C2*CH(AM2,ZMII,H)*SIN(BETA2)*SIN(E2)
            PRO2(3,4)=C2*SH(AM2,ZMII,H)*COS(E2)
      
            CALL PRODT(PRO1,PRO2,FP,FM,3)
            DO J=1,6
                  ACQ=AIRE(II)*CN(II,J)*
     1 (((-1.)**(J+1))**(JJ+1))*RHO*G*VA
                  EFWPS(6,J,KNC)=EFWPS(6,J,KNC)+FP(1)*ACQ !edited by RK, replaced index 8 by 6 in the right hand side
                  EFWPS(6,J+6,KNC)=EFWPS(6,J+6,KNC)+FP(2)*ACQ
                  EFWMN(6,J,KNC)=EFWMN(6,J,KNC)+FM(1)*ACQ
                  EFWMN(6,J+6,KNC)=EFWMN(6,J+6,KNC)+FM(2)*ACQ
            END DO
      END DO
      END DO
      ENDIF
      ENDIF
      
C *** CONTRIBUTION DUE AU PRODUIT DES MOUVEMENTS DE ROTATION ET
C                                            DES EFFORTS INERTIES ***
C     Fh = Mẍ + Adx/dt - kh x + Fex
C     CMi : masse ajoutée à la pulsation i
C     ASH : raideur hydrostatique
C     CAi : amortissement ajouté à la pulsation i
C     FORi : efforts d'exitation du premier ordre à la pulsation i
C
C     attention !=Iẍ!!!!!
C     
C     Ai(1,6) : partie reelle
C     Ai(7,12) : partie imaginaire
C
C --> FINT1i(1,6 ) : partie reelle
C --> FINT1i(7,12) : partie imaginaire
C --> FINT2i(1,6 ) : partie imaginaire
C --> FINT2i(7,12) : - partie reelle
C
C RFI
C
C
C -------------------------------------------------------------------
C *** CONTRIBUTION DUE TO THE PROCEEDS OF ROTATION MOVEMENTS AND
C INERTED EFFORTS ***
C Fh = Mẍ + Adx / dt - kh x + Fex
C CMi: mass added to the pulsation i
C ASH: hydrostatic stiffness
C CAi: damping added to the pulsation i
C FORi: exitation forces of the first order at the pulsation i
C
C attention! = Iẍ !!!!!
C
C Ai (1,6): real part
C Ai (7,12): imaginary part
C
C -> FINT1i (1,6): real part
C -> FINT1i (7,12): imaginary part
C -> FINT2i (1,6): imaginary part
C -> FINT2i (7,12): - real part
C
C RFI
C
C
C -------------------------------------------------------------------
      DIM1=WR1*WR1
      DIM2=WR2*WR2
      DO IJ=1,6
      DO IK=1,6
            AINTR1(IJ,IK)=ASH(IJ,IK,KNC)-DIM1*CM1(IJ,IK,KNC,KNC)
            AINTM1(IJ,IK)=-WR1*CA1(IJ,IK,KNC,KNC)
            AINTR2(IJ,IK)=ASH(IJ,IK,KNC)-DIM2*CM2(IJ,IK,KNC,KNC)
            AINTM2(IJ,IK)=-WR2*CA2(IJ,IK,KNC,KNC)
c~             AINTINERTIER(IJ,IK)=DIM1*AIN(IJ,IK,KNC)
c~      1 +RAID(IJ+6*(KNC-1),IK+6*(KNC-1))
c~             AINTINERTIEM(IJ,IK)=-WR1*AMOR(IJ+6*(KNC-1),IK+6*(KNC-1))
      END DO
      END DO
c~       CALL FORINE(AINTINERTIER,A1,FINT1)
c~       CALL FORINE(AINTINERTIEM,A1,FINT2)

      CALL FORINE(AINTR1,A1,FINT11)     ! Matrix multiplication F=[-w^2M+K]_6x6 [X]_6X1         
      CALL FORINE(AINTM1,A1,FINT21)     ! F=[-wB]_6x6 [X]_6X1
      CALL FORINE(AINTR2,A2,FINT12)
      CALL FORINE(AINTM2,A2,FINT22)
      DO J=1,6
            FHY1(J)=REAL(FOR1(J,KNC))+FINT11(J)-FINT21(J+6)
            FHY1(J+6)=IMAG(FOR1(J,KNC))+FINT11(J+6)+FINT21(J)
            FHY2(J)=REAL(FOR2(J,KNC))+FINT12(J)-FINT22(J+6)
            FHY2(J+6)=IMAG(FOR2(J,KNC))+FINT12(J+6)+FINT22(J)
            PRO2(J,1)=FHY1(J)
            PRO2(J,2)=FHY1(J+6)
            PRO2(J,3)=FHY2(J)
            PRO2(J,4)=FHY2(J+6)    
            !print*,A1(J),A1(J+6)
            ! print*,FOR1(J,KNC)        
            !print*,PRO2(J,1),PRO2(J,2),PRO2(J,3),PRO2(J,4) 
            !print*, (AINTR1(J,IK),IK=1,6)
            !WRITE(*,'(6(X,F12.8))') (AINTM1(J,IK),IK=1,6)
      END DO
      
      DO J=1,6
            DO J1=1,6
            DO J2=1,4
                  PRO1(J1,J2)=0.
            END DO
            END DO
            PRO1(K1(J),1)=-A1(K3(J))
            PRO1(K1(J),2)=-A1(K3(J)+6)
            PRO1(K2(J),1)=A1(K4(J))
            PRO1(K2(J),2)=A1(K4(J)+6)
            PRO1(K1(J),3)=-A2(K3(J))
            PRO1(K1(J),4)=-A2(K3(J)+6)
            PRO1(K2(J),3)=A2(K4(J))
            PRO1(K2(J),4)=A2(K4(J)+6)
            CALL PRODT(PRO1,PRO2,FP,FM,6)
            EFWPS(3,J,KNC)=FP(1)
            EFWPS(3,J+6,KNC)=FP(2)
            EFWMN(3,J,KNC)=FM(1)
            EFWMN(3,J+6,KNC)=FM(2)
          !  print*,J 
          !  print*,K1(J),PRO1(K1(J),1),PRO1(K1(J),2)
          !  print*,K2(J),PRO1(K2(J),1),PRO1(K2(J),2)
          !  print*,K1(J),PRO1(K1(J),3),PRO1(K1(J),4)
          !  print*,K2(J),PRO1(K2(J),3),PRO1(K2(J),4)

          !  print*,FM(1),FM(2),FP(1),FP(2)
      END DO
      Q3(1,I1,I2)=EFWMN(3,3,KNC)/2
      Q3(2,I1,I2)=EFWMN(3,9,KNC)/2

C
C *** CONTRIBUTION DE L'ELEVATION DE VAGUE LE LONG DE LA LIGNE
C                                                  DE LA FLOTTAISON ***
C *** CONTRIBUTION OF WAVE ELEVATION ALONG THE LINE
C OF FLOTATION ***
C ---------------------------------------------------------------------
      CTB=RHO*G*0.5
      DO JJ=NMI,NMA     !!
      DO II=IMX+1,IXX   !!D'APRES CE QUE J'AI COMPRIS CELA SERT A
            II1=II-IMX          !!NE GARDER QUE LA LIGNE DE FLOTTAISON
            IF(IND(II1).GE.IA.AND.IND(II1).LE.IB) THEN  
! ON VERIFIE QUE L'INDICE EST DANS [IA,IB] (CORPS CONSIDERE)
                  GAMMA0=ALF(II1)         
!LONGUEUR DE LA LIGNE DE FLOTTAISON --> OUI! cf QTFGeom l.99
                  XMII=XM(II)-XEFF
                  YMII=(YM(II)-YEFF)*((-1.)**(JJ+1))
                  H1W1=(-WR1*VGM1(JJ,II1)+VA*VXR1(JJ,IND(II1)))/G       ! attention non-homogene
                  H2W1=(WR1*VGR1(JJ,II1)+VA*VXM1(JJ,IND(II1)))/G       
                  H1W2=(-WR2*VGM2(JJ,II1)+VA*VXR2(JJ,IND(II1)))/G!CORRECTION DE BUG FROBAUX
                  H2W2=(WR2*VGR2(JJ,II1)+VA*VXM2(JJ,IND(II1)))/G!CORRECTION DE BUG FROBAUX
                  PRO1(1,1)=H1W1-A1(3)-(A1(4)*YMII-A1(5)*XMII)
                  PRO1(1,2)=H2W1-A1(9)-(A1(10)*YMII-A1(11)*XMII)
                  PRO1(1,3)=H1W2-A2(3)-(A2(4)*YMII-A2(5)*XMII)
                  PRO1(1,4)=H2W2-A2(9)-(A2(10)*YMII-A2(11)*XMII)
                  DO J=1,4
                        PRO2(1,J)=PRO1(1,J)
                  END DO
                  CALL PRODT(PRO1,PRO2,FP,FM,1)
                  DO J=1,6
c~                         COST=SQRT(1-CN(II,3)**2)
                        ACQ=CTB*GAMMA0*CN(II,J)*((-1.)**(J+1))**(JJ+1)
c~      1 /COST      ! attention cos(theta)
                        EFWPS(4,J,KNC)=EFWPS(4,J,KNC)-FP(1)*ACQ
                        EFWPS(4,J+6,KNC)=EFWPS(4,J+6,KNC)-FP(2)*ACQ
                        EFWMN(4,J,KNC)=EFWMN(4,J,KNC)-FM(1)*ACQ
                        EFWMN(4,J+6,KNC)=EFWMN(4,J+6,KNC)-FM(2)*ACQ
                    !IF (J==1.AND.JJ==1) THEN
                    !!print*,II,AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))
                    ! print*,II-IMX,-FM*ACQ,-FP*ACQ
                    !!print*,II-IMX,FM(1),FM(2),FP(1),FP(2)
                    !! print*,II-IMX,VGR1(JJ,II1),VGM1(JJ,II1)
                    !  ENDIF
                  END DO
            ENDIF
      END DO 
      END DO
      Q4(1,I1,I2)=EFWMN(4,3,KNC)/2
      Q4(2,I1,I2)=EFWMN(4,9,KNC)/2

!      print*,I1,I2,EFWMN(4,1,KNC),EFWMN(4,7,KNC)
C
C *** CONTRIBUTION COMPLEMENTAIRE HYDROSTATIQUE ***
CC *** ADDITIONAL HYDROSTATIC CONTRIBUTION ***
C Eq 4.51, 4.52 thesis
C -------------------------------------------------
      DO J=1,4
            PRO1(1,1)=A1(K5(J))     !DEPLACEMENT SUIVANT K5(J):K5/4,5,5,4/
            PRO1(1,2)=A1(K5(J)+6)   !   DONC NUL SI RAO EN ANGLES SONT NULLES
            PRO1(1,3)=A2(K5(J))
            PRO1(1,4)=A2(K5(J)+6)
            PRO2(1,1)=A1(K6(J))     !K6/4,5,6,6/
            PRO2(1,2)=A1(K6(J)+6)
            PRO2(1,3)=A2(K6(J))
            PRO2(1,4)=A2(K6(J)+6)
            CALL PRODT(PRO1,PRO2,FP,FM,1)
            ZF=1
            IF(J.LE.2) THEN
                  ZF=-ZEFF
            ENDIF
            DO J1=3,5
                  ACQ=X8(J)/2.*ZF
                  EFWPS(5,J1,KNC)=EFWPS(5,J1,KNC)+
     1 ACQ*ASH(J1,K7(J),KNC)*FP(1)
                  EFWPS(5,J1+6,KNC)=EFWPS(5,J1+6,KNC)+
     1 ACQ*ASH(J1,K7(J),KNC)*FP(2)
                  EFWMN(5,J1,KNC)=EFWMN(5,J1,KNC)+
     1 ACQ*ASH(J1,K7(J),KNC)*FM(1)
                  EFWMN(5,J1+6,KNC)=EFWMN(5,J1+6,KNC)+
     1 ACQ*ASH(J1,K7(J),KNC)*FM(2)       
            END DO 
      END DO
      Q5(1,I1,I2)=EFWMN(5,3,KNC)/2
      Q5(2,I1,I2)=EFWMN(5,9,KNC)/2
C
C *** TERME COMPLEMENTAIRE: EFFORTS / POINT FIXE ***
C FIX
C 
C Moments induit par la résultante d'ordre 1  qui s'exprime en G 
C et non en G au repos --> produit vectoriel {GoG^(1)}{F}
C*** ADDITIONAL TERM: EFFORTS / FIXED POINT ***
C FIX
C
C Moments induced by the resultant of order 1 which is expressed in G
C and not in G at rest -> vector product {GoG ^ (1)} {FInertia}
C L X FI  see Eq 4.46 Xiabo Chen's thesis
C _________________________________________________
      DO J=1,3
            PRO2(J,1)=FHY1(J)
            PRO2(J,2)=FHY1(J+6)
            PRO2(J,3)=FHY2(J)
            PRO2(J,4)=FHY2(J+6)
      END DO

      DO J=1,3
            DO J1=1,6
            DO J2=1,4
                  PRO1(J1,J2)=0.
            END DO
            END DO
            PRO1(K1(J),1)=-A1(K3(J)-3)
            PRO1(K1(J),2)=-A1(K3(J)+3)
            PRO1(K2(J),1)=A1(K4(J)-3)
            PRO1(K2(J),2)=A1(K4(J)+3)
            PRO1(K1(J),3)=-A2(K3(J)-3)
            PRO1(K1(J),4)=-A2(K3(J)+3)
            PRO1(K2(J),3)=A2(K4(J)-3)
            PRO1(K2(J),4)=A2(K4(J)+3)
            CALL PRODT(PRO1,PRO2,FP,FM,3)
            EFWPS(8,J+3,KNC)=FP(1)
            EFWPS(8,J+9,KNC)=FP(2)
            EFWMN(8,J+3,KNC)=FM(1)
            EFWMN(8,J+9,KNC)=FM(2)
            !IF (I1==30.AND.I2==39) THEN
            !IF (J==1)  print*,I1,I2
            !print*,J+3,FM(1),FM(2),FP(1),FP(2)
            !STOP
            !ENDIF
      END DO
      Q8(1,I1,I2)=EFWMN(8,3,KNC)
      Q8(2,I1,I2)=EFWMN(8,7,KNC)
C
C *** SOMMATION DES EFFORTS ***
C -----------------------------
C
      DO J2=1,12
            EFWPS(9,J2,KNC)=0.
            EFWMN(9,J2,KNC)=0.
            EFWPS(10,J2,KNC)=0.
            EFWMN(10,J2,KNC)=0.
! EFFORTS PAR RAPPORT AU POINT DE CALCUL EN MOUVEMENT AVEC LE CORPS
            DO J1=1,5 ! reduit a 3 pour verification (5 normalement)
! ON SOMME LES DIFFERENTES CONTRIBUTIONS CALCULEES DE 1 A 5
!! WE SUM UP THE DIFFERENT CONTRIBUTIONS CALCULATED FROM 1 TO 5   
            IF (J1 /= 0.) THEN ! /= not equal
                  EFWPS(9,J2,KNC)=EFWPS(9,J2,KNC)+EFWPS(J1,J2,KNC)
                  EFWMN(9,J2,KNC)=EFWMN(9,J2,KNC)+EFWMN(J1,J2,KNC)
            ENDIF     
            END DO
! EFFORTS PAR RAPPORT AU POINT DE CALCUL CONSIDERE COMME FIXE
! EFFORTS IN RELATION TO THE CALCULATION POINT CONSIDERED AS FIXED  
            EFWPS(10,J2,KNC)=EFWPS(9,J2,KNC)+EFWPS(8,J2,KNC)
            EFWMN(10,J2,KNC)=EFWMN(9,J2,KNC)+EFWMN(8,J2,KNC)
      END DO

      DO J=1,6
      DO K=1,5
            CONTRIB(K,J,I1,I2)=EFWMN(K,J,KNC)
            CONTRIB(K,J+6,I1,I2)=EFWMN(K,J+6,KNC)
      END DO
      END DO




C
C *** ECRITURES DES DIFFERENTES CONTRIBUTIONS ***
C -----------------------------------------------
      IF(IECR2.EQ.1) THEN
            WRITE(*,3010)
            WRITE(*,3012)
            WRITE(*,3020)
            WRITE(*,3001)(EFWMN(1,I,KNC),I=1,6)
            WRITE(*,3002)(EFWMN(2,I,KNC),I=1,6)
            WRITE(*,3003)(EFWMN(3,I,KNC),I=1,6)
            WRITE(*,3004)(EFWMN(4,I,KNC),I=1,6)
            WRITE(*,3005)(EFWMN(5,I,KNC),I=1,6)
            WRITE(*,3006)(EFWMN(6,I,KNC),I=1,6)
            WRITE(*,3008)(EFWMN(8,I,KNC),I=1,6)
      ENDIF
C
C
C *** ECRITURE DES EFFORTS PAR RAPPORT AU POINT DE CALCUL EN
C                                   MOUVEMENT AVEC LE CORPS CONSIDERE***
C *** WRITING OF FORCES IN RELATION TO THE CALCULATION POINT IN
C MOVEMENT WITH THE BODY CONSIDERED ***
C ______________________________________________________________________
      J1=9
C *** MISE A 0 DES CONTRIBUTIONS TROP PETITES ***
C *** SET TO 0 CONTRIBUTIONS TOO SMALL ***
      IF(IMZERO.EQ.1) THEN
            DO J2=1,12
                  SOLP(J2)=EFWPS(J1,J2,KNC)
                  SOLM(J2)=EFWMN(J1,J2,KNC)
            END DO
            CALL AFIN(12,SOLP,ERRE)
            CALL AFIN(12,SOLM,ERRE)
            DO J2=1,12
                  EFWPS(J1,J2,KNC)=SOLP(J2)
                  EFWMN(J1,J2,KNC)=SOLM(J2)
            END DO
      ENDIF
C *** VALEURS ABSOLUES DES EFFORTS ET PHASES ***
C *** ABSOLUTE VALUES OF EFFORTS AND PHASES ***
      DO J2=1,6
            AMODP(J2)=SQRT(EFWPS(J1,J2,KNC)**2+EFWPS(J1,J2+6,KNC)**2)
            IF(ABS(EFWPS(J1,J2+6,KNC)).LE.1.E-20) THEN
                  EFWPS(J1,J2+6,KNC)=SIGN(1.E-20,EFWPS(J1,J2+6,KNC))
            ENDIF
            PHASP(J2)=ATAN2(EFWPS(J1,J2,KNC),EFWPS(J1,J2+6,KNC))*RADD
     1 +180.
            PHASP(J2)=AMOD(PHASP(J2),360.)
            IF(ABS(EFWPS(J1,J2+6,KNC)).LE.1.E-20) THEN
                  EFWPS(J1,J2+6,KNC)=0.000000
            ENDIF
            AMODM(J2)=SQRT(EFWMN(J1,J2,KNC)**2+EFWMN(J1,J2+6,KNC)**2)
            IF(ABS(EFWMN(J1,J2+6,KNC)).LE.1.E-20) THEN
                  EFWMN(J1,J2+6,KNC)=SIGN(1.E-20,EFWMN(J1,J2+6,KNC))
            ENDIF
            PHASM(J2)=ATAN2(EFWMN(J1,J2,KNC),EFWMN(J1,J2+6,KNC))*RADD
     1 +180.
            PHASM(J2)=AMOD(PHASM(J2),360.)
            IF(ABS(EFWMN(J1,J2+6,KNC)).LE.1.E-20) THEN
                  EFWMN(J1,J2+6,KNC)=0.0000000
            ENDIF
      END DO

C
C ECRITURE DES FONCTIONS DE TRANSFERT BICHROMATIQUES ADIMENSIONNELLES
C -------------------------------------------------------------------
C    FONCTION DE TRANSFERT=EFFORT BICHROMATIQUE/2
C    ECRITES EN FONCTION DES PULSATIONS DANS LE REPERE ABSOLU
C BASSE FREQUENCE                                      ------
C WRITING ADIMENSIONAL BICHROMATIC TRANSFER FUNCTIONS
C  ------------------------------------------------- ------------------
C TRANSFER FUNCTION = BICHROMATIC FORCES / 2
C WRITTEN ACCORDING TO THE PULSATIONS IN THE ABSOLUTE
C LOW FREQUENCY ------
C ----------------------------------------------------------------------      
      IF(IECRBF.EQ.1) THEN
            WRITE(26,799)1/T1,EFWMN(J1,1,KNC)/2
      ENDIF
C HAUTE FREQUENCE
      IF(IECRHF.EQ.1) THEN
            WRITE(27,799)W1,W2,
     1 (EFWPS(J1,J,KNC)/AD/2,EFWPS(J1,J+6,KNC)/AD/2,J=1,3),
     1 (EFWPS(J1,J,KNC)/AD/XL/2,EFWPS(J1,J+6,KNC)/AD/XL/2,J=4,6)
      ENDIF
 799  FORMAT(1P,7E11.3)
C
C *** ECRITURE DES EFFORTS PAR RAPPORT AU POINT DE CALCUL
C                                          CONSIDERE COMME FIXE***
C *** WRITING OF EFFORTS IN RELATION TO THE CALCULATION POINT
C CONSIDERED AS FIXED ***
C ________________________________________________________________
      IF(IECR1.EQ.1) THEN
            J1=10!TOUJOURS DANS CETTE BOUCLE
C *** MISE A 0 DES CONTRIBUTIONS TROP PETITES ***
            IF(IMZERO.EQ.1) THEN
                  DO J2=1,12
                        SOLP(J2)=EFWPS(J1,J2,KNC)
                        SOLM(J2)=EFWMN(J1,J2,KNC)
                  END DO
                  CALL AFIN(12,SOLP,ERRE)
                  CALL AFIN(12,SOLM,ERRE)
                  DO J2=1,12
                        EFWPS(J1,J2,KNC)=SOLP(J2)
                        EFWMN(J1,J2,KNC)=SOLM(J2)
                  END DO
            ENDIF
C *** VALEURS ABSOLUES DES EFFORTS ET PHASES ***
            DO J2=1,6
                  AMODP(J2)=SQRT(EFWPS(J1,J2,KNC)**2
     1 +EFWPS(J1,J2+6,KNC)**2)
                  IF(ABS(EFWPS(J1,J2+6,KNC)).LE.1.E-20) THEN
                        EFWPS(J1,J2+6,KNC)=
     1 SIGN(1.E-20,EFWPS(J1,J2+6,KNC))
                  ENDIF
                  PHASP(J2)=ATAN2(EFWPS(J1,J2,KNC),EFWPS(J1,J2+6,KNC))
     1 *RADD+180.
                  PHASP(J2)=AMOD(PHASP(J2),360.)
                  IF(ABS(EFWPS(J1,J2+6,KNC)).LE.1.E-20) THEN
                        EFWPS(J1,J2+6,KNC)=0.000000
                  ENDIF
                  AMODM(J2)=SQRT(EFWMN(J1,J2,KNC)**2
     1 +EFWMN(J1,J2+6,KNC)**2)
                  IF(ABS(EFWMN(J1,J2+6,KNC)).LE.1.E-20) THEN
                        EFWMN(J1,J2+6,KNC)=
     1 SIGN(1.E-20,EFWMN(J1,J2+6,KNC))
                  ENDIF
                  PHASM(J2)=ATAN2(EFWMN(J1,J2,KNC),EFWMN(J1,J2+6,KNC))
     1 *RADD+180.
                  PHASM(J2)=AMOD(PHASM(J2),360.)
                  IF(ABS(EFWMN(J1,J2+6,KNC)).LE.1.E-20) THEN
                        EFWMN(J1,J2+6,KNC)=0.0000000
                  ENDIF
            END DO
       !!CALCUL DES FORCES TOTALES DUES AU PREMIER ORDRE
      ENDIF      
      !! CALCULATION OF TOTAL FORCES DUE TO FIRST ORDER
      DO DOF=1,6  ! BOUCLE INITIALISER FCONT ET EFFP
      DO L=1,2    ! index real, imajiner 
            DO M=1,7
            DO K=1,2
                FCONT(DOF,1,I2+1,1,L,M)=WR2
            END DO
            END DO
            DO M=1,5 ! force term contribution
                  FCONT(DOF,1,I2+1,I1+1,L,M)=EFWMN(M,DOF+6*(L-1),KNC)/2
                  FCONT(DOF,2,I2+1,I1+1,L,M)=EFWPS(M,DOF+6*(L-1),KNC)/2
            END DO
            FCONT(DOF,1,I2+1,I1+1,L,6)=EFWMN(8,DOF+6*(L-1),KNC)/2
            FCONT(DOF,2,I2+1,I1+1,L,6)=EFWPS(8,DOF+6*(L-1),KNC)/2
            FCONT(DOF,1,I2+1,I1+1,L,7)=EFWMN(10,DOF+6*(L-1),KNC)/2
            FCONT(DOF,2,I2+1,I1+1,L,7)=EFWPS(10,DOF+6*(L-1),KNC)/2
                       
            DO K=1,3
                  EFFP(DOF+6*(L-1),I2+1,1,K)=WR2
            END DO
      
      END DO
      END DO

      
      
      



! ################## APROXIMATION DE PINKSTER #########################
      IF(I2.NE.I1) THEN
      
      C12C=FIJ(AM1,AM2,H)
      CHOIXAPPROX=3
      
      ! APPROXIMATION PAR LES EFFORTS DE FROUDE KRYLOV
      DO J=1,6
            FKR(J)=0
            FKI(J)=0
      END DO
      DO JJ=NMI,NMA!INTEGRATION SUR LES FACETTES DU CORPS
      DO II=IA,IB
            XMII=XM(II)-XEFF
            YMII=(YM(II)-YEFF)*((-1.)**(JJ+1))
            ZMII=ZM(II)-ZEFF
            C1H=-CH(ABS(AM1-AM2),ZM(II),H)*G*RHO
            PHIR=C1H*COS((AM1-AM2)*(XMII*COS(0.)+YMII*SIN(0.)))
            PHII=C1H*SIN((AM1-AM2)*(XMII*COS(0.)+YMII*SIN(0.))) !this is only for heading 0
            DO J=1,6!INTEGRATION POUR CHAQUE DOF
                  ACQ=AIRE(II)*CN(II,J)*(((-1.)**(J+1))**(JJ+1))!!N*DS
                  FKR(J)=FKR(J)+ACQ*PHIR/2                          !!-\IINT(PI.N.DS)
                  FKI(J)=FKI(J)+ACQ*PHII/2              !edited by RK
            END DO
      END DO
      END DO
      
      DO J=1,6
            EFFP(J,I2+1,1,1)=WR2
            EFFP(J,I2+1,I1+1,1)=C12C*FKR(J)
            EFFP(J+6,I2+1,1,1)=WR2
            EFFP(J+6,I2+1,I1+1,1)=C12C*FKI(J)
      END DO



      IF(H.NE.0 .AND.H.LT.0.9E20) THEN
            AK1=AM1*H
            AK2=AM2*H
            AKK=(AK1-AK2)
            AKM=X01(AKK)
            WM=SQRT(AKM*G/H)
      ELSE
            AK1=AM1
            AK2=AM2
            WM=SQRT(G*ABS(AK2-AK1))
      ENDIF
      WM=SQRT((AM1-AM2)*G*TANH((AM1-AM2)*H))!!!EFFORTS D'EXCITATIONS A LA PULSATION K1-K2
c~       DO IPP=2,NHASKIND                                              ! version matrice carre pour les QTF+
      DO IPP=2,N                                                        ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
            IF(WWX(IPP).GT.WM) THEN
                  IPM=IPP
                  GOTO 3886
            ENDIF
      END DO
c~       IPM=NHASKIND                                                   ! version matrice carre pour les QTF+
      IPM=N                                                             ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
 3886 CONTINUE

      DO J=1,6 !ReConstructing the excitation force with linear interpolation for freq WM
            FKR(J)=(FEX(J,IPM-1)*COS((FEXPH(J,IPM-1))))
     1 *(WWX(IPM)-WM)/(WWX(IPM)-WWX(IPM-1))
     1 +(FEX(J,IPM)*COS((FEXPH(J,IPM))))
     1 *(WM-WWX(IPM-1))/(WWX(IPM)-WWX(IPM-1))
            FKI(J)=FEX(J,IPM-1)*SIN((FEXPH(J,IPM-1)))
     1 *(WWX(IPM)-WM)/(WWX(IPM)-WWX(IPM-1))+(FEX(J,IPM)
     1 *SIN((FEXPH(J,IPM))))
     1 *(WM-WWX(IPM-1))/(WWX(IPM)-WWX(IPM-1))
      END DO
      DO J=1,6
            EFFP(J,I2+1,1,2)=WR2
            EFFP(J,I2+1,I1+1,2)=-C12C*FKI(J)
!! FONCTIONNE, MULTIPLICATION PAR i ETRANGE, ON DOIT FAIRE LA MEME AU EFFORTS D'INERTIE (3EME APPROX)
!! WORKS, MULTIPLICATION BY I ETRANGE, WE MUST DO THE SAME FOR INERTIAL EFFORTS (3RD APPROX)            
            EFFP(J+6,I2+1,1,2)=WR2
            EFFP(J+6,I2+1,I1+1,2)=-C12C*FKR(J) ! added by RK so Conjugate(1i*FK) then the results same as HASBO
      END DO

      DO J=1,6 !interpolating RAO
            ZZA(J)=ZZX(J,KNC,IPM-1)+(WM-WWX(IPM-1))
     1 /(WWX(IPM)-WWX(IPM-1))*(ZZX(J,KNC,IPM)-ZZX(J,KNC,IPM-1))
            A1(J)=REAL(ZZA(J))
            A1(J+6)=AIMAG(ZZA(J))
            WRITE(99,*)J,A1(J),A1(J+6)
      END DO
      DO IJ=1,6
      DO IK=1,6
            AINTV(IJ,IK)=AIN(IJ,IK,KNC) !Inertia matrix
      END DO
      END DO
      CALL FORINE(AINTV,A1,FINT1) ! M*X
      DIMM=WM*WM
      C12=G*WM/((W1**2-W2**2)*TANH((W1**2-W2**2)/G*H)-(W1-W2)**2)
      C12B=G*WM/(2*W2*(W1-W2))
      C12C=FIJ(AM1,AM2,H)
      DO IIKL=1,6
            WRITE(99,*)IIKL,WM**2,FINT1(IIKL),FINT1(IIKL+6)
            WRITE(99,*)IIKL,FINT1(IIKL)*WM**2,FINT1(IIKL+6)*WM**2
      END DO
      WRITE(99,*)H,TANH((W1**2-W2**2)/G*H)
      WRITE(99,*)H,((W1**2-W2**2)-(W1-W2)**2),2*W2*(W1-W2)
      WRITE(99,*)'C12 = ',C12,C12B,C12C
      DO J=1,6
             FKR(J)=-FINT1(J+6)*DIMM/2 !-w^2 MX
             FKI(J)=FINT1(J)*DIMM/2
      END DO

      DO J=1,6
            EFFP(J,I2+1,1,3)=WR2
            EFFP(J,I2+1,I1+1,3)=C12C*FKR(J)
            EFFP(J+6,I2+1,1,3)=WR2
            EFFP(J+6,I2+1,I1+1,3)=C12C*FKI(J)
      END DO
      ENDIF






! WE SUM THE VALUES IN QTF FOUND AT THE DEPARTURE BY DUOK, TO WHICH WE ADD THE VALUES FOUND WITH THE APPROXIMATION OF PINKSTER
C *** END OF ITERATIONS ON THE NUMBER OF BODIES ***
!ON SOMME LES VALEURS EN QTF TROUVES AU DEPART PAR DUOK, AUXQUELS ON AJOUTE LES VALEURS TROUVEES AVEC L'APPROXIMATION DE PINKSTER
C *** FIN DES ITERATIONS SUR LE NOMBRE DE CORPS ***
      END DO
C
C
      DO I=1,6
            QTF1(I,I1,I2)=AMODM(I)!MODULE DES QTF DUOK SEUL
      END DO

C
C *** FIN DES ITERATIONS SUR LES PULSATIONS ***
      END DO
      END DO

      QPM(1)='M'
      QPM(2)='P'
      QRI(1)='_PR'
      QRI(2)='_PI'
      
      CCONT(7)='DUOK'
      CCONT(1)='DADV'
      CCONT(2)='DDPV'
      CCONT(3)='DRFI'
      CCONT(4)='DFLO'
      CCONT(5)='DSTA'
      CCONT(6)='DFIX'
      CCONT(8)='APINK'
      
      IF (Loutduok==1) THEN
	    OUT1=1	! ON ECRIT TOUTES LES CONTRIBUTIONS
      ELSE
	    OUT1=7	! ON ECRIT SUELEMENT DUOK ET PINKSTER
      ENDIF
      IF (LQTFP .EQ. 1) THEN
            NQTFP=2
      ELSE
            NQTFP=1
      ENDIF
      
      DO DOF=1,6,DOFSYM        ! ECRITURE DES RESULTATS, if symmetric then writes only 1,3,5
      DO L=1,2                ! PR PI (L=1 : PR / L=2 : PI)
            DO K=1,NQTFP      ! QTF+- (K=1 : QTF- / K=2 : QTF+)
            DO M=OUT1,7       ! CONTRIB
                IF (M<7) THEN
                  FCHA=ID(1:lID)//'/results/QTF/Appendix/QTF'//QPM(K)
     1 //'_'//TRIM(CCONT(M))//'_DOF_'//CHAR(DOF+48)//QRI(L)//'.dat'
                ELSE 
                  FCHA=ID(1:lID)//'/results/QTF/QTF'//QPM(K)//'_'
     1 //TRIM(CCONT(M))//'_DOF_'//CHAR(DOF+48)//QRI(L)//'.dat'
                ENDIF
                IUNI=1116+L+2*DOF+12*K+M*24
                OPEN(IUNI,FILE=TRIM(FCHA),STATUS='UNKNOWN')
                  
                  
c~                 DO I=1,NHASKIND+1                                    ! version matrice carre pour les QTF+
                DO I=1,N+1                                              ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
                  DO J=1,I
                    FCONT(DOF,K,I,J,L,M)=FCONT(DOF,K,J,I,L,M)       !SYMETRISATION
                  ENDDO
                ENDDO
                DO I=2,N+1                                              ! 2 car sinon on change aussi la pulsation 
                  IF (L==2 .AND. K==1) THEN
                     !!!  QTFij-=QTFji-* !!!
                     ! ON PREND LE CONJUGE DU TRIANGLE SUP
                     FCONT(DOF,K,I,I,L,M)=0.
                     DO J=I+1,N+1
                        FCONT(DOF,K,I,J,L,M)=-FCONT(DOF,K,I,J,L,M)
                     ENDDO
                  ENDIF
!!!! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs) !!!!
!JND : attention, critere a changer si le tri des pulsation n'est pas croissant
!!!! antitriangular matrix version sup for QTF + (0 elsewhere) !!!!
! JND: attention, criteria to change if the sorting of the pulsation is not increasing
                  IF (K==2) THEN
                        DO J=1+(N+1)-(I-1),N+1
                           FCONT(DOF,K,I,J,L,M)=0.
                        ENDDO
                  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ENDDO

                N1=N+1
                DO I=1,N1
c~                   WRITE(IUNI,*)(FCONT(DOF,K,I,J,L,M),J=1,NHASKIND+1)         ! version matrice carre pour les QTF+
                 WRITE(IUNI,'(<N1>E15.4)')(FCONT(DOF,K,I,J,L,M),J=1,N1)          ! version matrice antitriangulaire sup pour les QTF+ (0 ailleurs)
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
                   WRITE(IUNI,*) CHAR(10),CHAR(10),
     1 "# NB: Les valeurs des QTF+ pour ", 
     1 "w1+w2>n*dw sont definies a zeros.",CHAR(10),
     1 "# en effet, un calcul de ", 
     1 "#ces valeurs necessiterait un calcul d'ordre 1 sur 2 fois plus ",
     1 "de pulsation."
                ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                CLOSE(IUNI)
            END DO
            END DO
            
            
            M=8   ! Pinkster
            FCHA=ID(1:lID)//'/results/QTF/QTF'//QPM(1)//'_'
     1 //TRIM(CCONT(M))//'_DOF_'//CHAR(DOF+48)//QRI(L)//'.dat'
            IUNI=1116+L+2*DOF+12+M*24
            OPEN(IUNI,FILE=TRIM(FCHA),STATUS='UNKNOWN')
            NdofIt=DOF+6*(L-1) 
c~             DO I=1,NHASKIND+1                                       
           DO I=1,N+1                                                  
                  DO J=1,I
                  DO K=1,3
                        EFFP(DOF,I,J,K)=EFFP(DOF,J,I,K)
                        EFFP(DOF+6,I,J,K)=-EFFP(DOF+6,J,I,K) ! changed by RK such that the imaginer part is skewsymmetric matric as should be!
                  END DO
                  END DO 
c~                   WRITE(IUNI,*) (EFFP(DOF+6*(L-1),I,J,2),J=1,NHASKIND+1)! (2 : # APPROX)     
                  WRITE(IUNI,'(<N1>E15.4)') ( EFFP(NdofIt,I,J,2),J=1,N1)
            END DO

            CLOSE(IUNI)
      END DO
      END DO      
    
C *** ECRITURE REPERE FIN DE FICHIER FONCTION DE TRANSFERT ***
      ZERO=0.
C BASSE FREQUENCE
      IF(IECRBF.EQ.1) THEN
            WRITE(26,799)(ZERO,J=1,14)
            CLOSE(UNIT=26)
      ENDIF
C HAUTE FREQUENCE
      IF(IECRHF.EQ.1) THEN
            WRITE(27,799)(ZERO,J=1,14)
            CLOSE(UNIT=27)
      ENDIF

1991  FORMAT(1E13.6,4000(3X,1E13.6))
1989  FORMAT('QTF-',9X,4000(3X,1E13.6))
      REE=2

C *** FORMAT ***
 2003 FORMAT(/2X,'AMPLITUDE DES EFFORTS 1,2,3 (RESPECTIVEMENT X,Y,',
     1 'Z) EN NEWTONS'/2X,'AMPLITUDE DES EFFORTS 4,5,6 (MOMENTS AUTOUR ',
     1 'DE OX,OY,OZ) EN MN'/ 2X,'DANS UN LIQUIDE DE MASSE VOLUMIQUE ',
     1 F6.0,' KG/M**3'/ 2X,'POUR DEUX HOULES D''AMPLITUDE A = 1 M'/ 2X,
     1 'EFFORTS  : F2(T)=AMPCOS*COS(WR*T)+AMPSIN*SIN(WR*T)'/2X,
     1 '           F2(T)=-MODULE(F2)*SIN(WR*T+PHI) , 0<PHI<360'/2X,
     1 'HOULE AU POINT XEFF,YEFF,0. : H1(T)=-A*SIN(W1*T)',/,
     1 '                                   H2(T)=-A*SIN(W2*T)'/)
 999  FORMAT(///,' EFFORTS DU SECOND ORDRE :',
     1 '  PRODUITS DE TERMES DU 1ER ORDRE',/,
     1 2X,'PULSATIONS DE RENCONTRE : WR1=',F7.4,' RD/S   WR2=',
     1 F7.4,' RD/S',/,
     1 2X,'PULSATIONS REPERE ABSOLU: W1 =',F7.4,' RD/S   W2 =',
     1 F7.4,' RD/S',/,
     1 2X,'PERIODES REPERE ABSOLU  : T1 =',F7.4,' RD/S   T2 =',
     1 F7.4,' RD/S',/,
     1 2X,'ANGLE DE LA HOULE       :BETA1=',F6.3,' DEG   BETA2=',
     1 F6.3,' DEG')
 2004 FORMAT(' EFFORTS PAR RAPPORT AU POINT(XEFF,YEFF,ZEFF)',
     1 'EN MOUVEMENT AVEC LE CORPS CONSIDERE',/)
 1014	FORMAT(' CONTRIBUTION EN MODE DIFFERENCE A LA PULSATION WR1-WR2=',
     1 F6.3,' RD/S')
 1013	FORMAT(' CONTRIBUTION EN MODE SOMME A LA PULSATION WR1+WR2=',F6.3,
     1 ' RD/S')
 1000 FORMAT(1X,'   EFFORT  ',5X,' AMPCOS  ',7X,' AMPSIN  ',
     1 8X,'  MODULE',9X,'PHASE'/
     1 1X,'FORCE DIR X :',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/
     1 1X,'FORCE DIR Y :',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/
     2 1X,'FORCE DIR Z :',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/
     3 1X,'MOMENT AUT X:',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/
     4 1X,'MOMENT AUT Y:',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/
     5 1X,'MOMENT AUT Z:',1P,2(E13.6,3X),2X,E13.6,3X,0PF8.2/)
 2005 FORMAT(' EFFORTS PAR RAPPORT AU POINT(XEFF,YEFF,ZEFF)',
     1 ' CONSIDERE COMME FIXE',/)

 3020 FORMAT(18X,'FX',10X,'FY',10X,'FZ',10X,'MX',10X,'MY',10X,'MZ')
 3010 FORMAT(' * EFFORTS BASSE FREQUENCE *')
 3011 FORMAT(' * EFFORTS HAUTE FREQUENCE *')
 3012 FORMAT(' PARTIE REELLE DES DIFFERENTES CONTRIBUTIONS')
 3013 FORMAT(' PARTIE IMAGINAIRE DES DIFFERENTES CONTRIBUTIONS')
 3001 FORMAT('  SSV**2    ',6(1X,1PE11.4))
 3002 FORMAT(' MOM*DV/DT  ',6(1X,1PE11.4))
 3003 FORMAT(' R(1)*FI(1) ',6(1X,1PE11.4))
 3004 FORMAT('  INT FLO   ',6(1X,1PE11.4))
 3005 FORMAT('   HYDS     ',6(1X,1PE11.4))
 3006 FORMAT(' VA.DGRAD/DX',6(1X,1PE11.4))
 3008 FORMAT(' POINT FIXE ',6(1X,1PE11.4))

      STOP

      END
      

      FUNCTION CH(AK,Z,H)
            IF(AK*H.LE.0.OR.AK*H.GT.20) THEN
                  CH=EXP(AK*Z)
            ELSE
                  CH=COSH(AK*(Z+H))/COSH(AK*H)
            ENDIF
            RETURN
      END

      FUNCTION SH(AK,Z,H)
            IF(AK*H.LE.0 .OR.AK*H.GT.20) THEN
                  SH=EXP(AK*Z)
            ELSE
                  SH=SINH(AK*(Z+H))/COSH(AK*H)
            ENDIF
            RETURN
      END

      SUBROUTINE PRODT(A,B,C,D,N)
            DIMENSION A(6,4),B(6,4),C(2),D(2)
            C(1)=0.
            C(2)=0.
            D(1)=0.
            D(2)=0.
            DO I=1,N
                  C(1)=C(1)+.5*(A(I,1)*B(I,3)-A(I,2)*B(I,4)+
     1                A(I,3)*B(I,1)-A(I,4)*B(I,2))
                  C(2)=C(2)+.5*(A(I,1)*B(I,4)+A(I,2)*B(I,3)+
     1                A(I,3)*B(I,2)+A(I,4)*B(I,1))
                  D(1)=D(1)+.5*(A(I,1)*B(I,3)+A(I,2)*B(I,4)+
     1                A(I,3)*B(I,1)+A(I,4)*B(I,2))
                  D(2)=D(2)+.5*(A(I,2)*B(I,3)-A(I,1)*B(I,4)+
     1                A(I,3)*B(I,2)-A(I,4)*B(I,1))
            END DO
            RETURN
      END

      SUBROUTINE FORINE(AI,A,CF)
      ! CF=-[AI]*A
            DIMENSION AI(6,6),A(12),CF(12)
            DO I=1,6
                  CF(I)=0.
                  CF(I+6)=0.
                  DO J=1,6
                        CF(I)=CF(I)-AI(I,J)*A(J)
                        CF(I+6)=CF(I+6)-AI(I,J)*A(J+6)
                  END DO
 	      END DO
            RETURN
      END

      FUNCTION PHASE(Z)
C	      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
            COMPLEX Z
            PI=3.141592653589793
            DZ=AIMAG(Z)
            RZ=REAL(Z)
            IF(RZ.EQ.0) THEN
                  PHASE=0.
                  IF(DZ.GT.0)PHASE=PI/2.
                  IF(DZ.LT.0)PHASE=3*PI/2.
            ELSE
                  PHASE=ATAN2(DZ,RZ)
                  IF(PHASE.LT.0)PHASE=2*PI+PHASE
            ENDIF
            PHASE=PHASE*180/PI
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

      FUNCTION X0(AK,G,H)
            F(T)=AK-T*G/H*TANH(T)
            EPS=5.E-6
            ITOUR=0
            XI=0.
            XS=XI
            PAS=AMAX1(AK,SQRT(AK))
30          XS=XS+PAS
            ITOUR=ITOUR+1
            IF(ITOUR.LE.1000) THEN
                  IF(F(XS)*F(XI).GT.0)GOTO 30
            ENDIF
            IITER=0
10          XM=(XI+XS)*0.5
            IITER=IITER+1
            IF(IITER.GT.1000.OR.ITOUR.GT.1000) THEN
                  WRITE(*,110)ITOUR,IITER
110               FORMAT(2X,'ERREUR DANS LA RECHERCHE DE LA RACINE',
     1 /2X,'APPROXIMATION =',I5,'   DICHOTOMIE = ',I5)
                  STOP
            ELSE
                  IF(ABS((XS-XI)/XM).GT.EPS) THEN
                        IF(F(XM)*F(XI).LT.0) THEN
                              XS=XM
                        ELSE
                              XI=XM
                        ENDIF
                        GOTO 10
                  ELSE
                        X0=XM!! RENVOIE LA RACINE A EPS PRES DE F(T)
                  ENDIF
            ENDIF
            RETURN
      END


      FUNCTION FIJ(AKI,AKJ,H)
            G=9.81
            WI=SQRT(G*AKI*TANH(AKI*H))
            WJ=SQRT(G*AKJ*TANH(AKJ*H))
            BIJ=AKI**2/(WI*COSH(AKI*H)**2)-AKJ**2/(WJ*COSH(AKJ*H)**2)
            CIJ=2*AKI*AKJ*(WI-WJ)*(1+TANH(AKI*H)*TANH(AKJ*H))/(WI*WJ)
            AIJ=0.5*G**2*(BIJ+CIJ)/((WI-WJ)**2-(AKI-AKJ)*G
     1 *TANH((AKI-AKJ)*H))
            FIJ=AIJ*(WI-WJ)/G
            RETURN
      END


      FUNCTION IL(NOM)
            CHARACTER*10 NOM
            IL=INDEX(NOM,' ')-1
            IF(IL.LT.0)IL=LEN(NOM)
            IF(IL.EQ.0) THEN
            NOM='DEFAUT'
            IL=6
            ENDIF
      RETURN
      END
      
      FUNCTION X01(AK)
            F(T)=AK-T*TANH(T)
            EPS=5.E-6
            ITOUR=0
            XI=0.
            XS=XI
            PAS=AMAX1(AK,SQRT(AK))
   30       XS=XS+PAS
            ITOUR=ITOUR+1
            IF(ITOUR.LE.1000) THEN
            IF(F(XS)*F(XI).GT.0)GOTO 30
            ENDIF
            IITER=0
   10       XM=(XI+XS)*0.5
            IITER=IITER+1
            IF(IITER.GT.1000.OR.ITOUR.GT.1000) THEN
                  WRITE(*,110)ITOUR,IITER
  110             FORMAT(2X,'ERREUR DANS LA RECHERCHE DE LA RACINE',
     1 /2X,'APPROXIMATION =',I5,'   DICHOTOMIE = ',I5)
                  STOP
            ELSE
                  IF(ABS((XS-XI)/XM).GT.EPS) THEN
                        IF(F(XM)*F(XI).LT.0) THEN
                              XS=XM
                        ELSE
                              XI=XM
                        ENDIF
                        GOTO 10
                  ELSE
                        X01=XM
                  ENDIF
            ENDIF
            RETURN
      END
