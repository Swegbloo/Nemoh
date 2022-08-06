!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991             
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr) 2014 
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
MODULE MReadInputFiles

USE CONSTANTS,                  ONLY: II, PI
USE MFileDirectoryList

IMPLICIT NONE
PUBLIC:: Read_NP_GaussQuad,Read_Mechanical_Coefs,Read_FirstOrderLoad,Read_Motion, &
         Read_SourceDistribution,Read_Eps_Zmin
!
TYPE TMech
    REAL,ALLOCATABLE,DIMENSION(:,:) :: MassMat          !Mass-Inertia Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat         !Stifness-Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat_EXT     !Additional Stifness Matrix i.e: mooring
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DampCoefMat_EXT  !Additional damping coefficients
END TYPE

TYPE TLoad1
    COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: excitation
    REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: addedmass
    REAL,ALLOCATABLE,DIMENSION(:,:,:)    :: dampcoef
END TYPE
TYPE TSource
    COMPLEX,ALLOCATABLE,DIMENSION(:,:)   :: ZIGB,ZIGS       !source distribution
END TYPE
!
CONTAINS
        INTEGER FUNCTION Read_NP_GaussQuad(wd)
           CHARACTER(LEN=*) :: wd
           INTEGER       :: NPGQ
           OPEN(10,file=wd//'/input_solver.txt',form='formatted',status='old')
           READ(10,*) NPGQ
           CLOSE(10)
           Read_NP_GaussQuad=NPGQ**2
           RETURN
        END FUNCTION

        FUNCTION Read_Eps_Zmin(wd) result(EPS_ZMIN)
           CHARACTER(LEN=*) :: wd
           REAL       :: EPS_ZMIN
           OPEN(10,file=wd//'/input_solver.txt',form='formatted',status='old')
           READ(10,*)
           READ(10,*) EPS_ZMIN
           CLOSE(10)
           RETURN
        END FUNCTION 


        SUBROUTINE Read_Mechanical_Coefs(wd,Nradiation,MechCoef)
           !input/output
           CHARACTER(LEN=*),            INTENT(IN)::wd
           INTEGER,                     INTENT(IN)::Nradiation
           TYPE(TMech),                 INTENT(OUT)::MechCoef
           !Local
           INTEGER      ::u1,u2,u3,u4,C,I,J
           
           ALLOCATE(MechCoef%MassMat(Nradiation,Nradiation))
           ALLOCATE(MechCoef%StiffMat(Nradiation,Nradiation))
           ALLOCATE(MechCoef%StiffMat_EXT(Nradiation,Nradiation))
           ALLOCATE(MechCoef%DampCoefMat_EXT(Nradiation,Nradiation))
           CALL exist_file(trim(wd)//'/Mechanics/Inertia.dat')
           CALL exist_file(trim(wd)//'/Mechanics/kh.dat')
           CALL exist_file(trim(wd)//'/Mechanics/km.dat')
           CALL exist_file(trim(wd)//'/Mechanics/Badd.dat')
           
           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/Mechanics/Inertia.dat',ACTION='READ')
           OPEN(NEWUNIT=u2,FILE=trim(wd)//'/Mechanics/kh.dat',ACTION='READ')
           OPEN(NEWUNIT=u3,FILE=trim(wd)//'/Mechanics/km.dat',ACTION='READ')
           OPEN(NEWUNIT=u4,FILE=trim(wd)//'/Mechanics/Badd.dat',ACTION='READ')
           
           DO I=1,Nradiation
                READ(u1,*) (MechCoef%MassMat(I,J),J=1,Nradiation)
                READ(u2,*) (MechCoef%StiffMat(I,J),J=1,Nradiation)
                READ(u3,*) (MechCoef%StiffMat_EXT(I,J),J=1,Nradiation)
                READ(u4,*) (MechCoef%DampCoefMat_EXT(I,J),J=1,Nradiation)
           ENDDO
           CLOSE(u1)
           CLOSE(u2)
           CLOSE(u3)
           CLOSE(u4)
        END SUBROUTINE

        SUBROUTINE Read_FirstOrderLoad(wd,Nw,Nbeta,Nintegration,Nradiation,Forces1)
           !Input/Output
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nintegration,Nradiation
           TYPE(TLoad1),                INTENT(OUT):: Forces1
           !Local
           INTEGER                                 :: u1,u2,u3,I,J,K
           REAL,DIMENSION(Nintegration)            :: Amp, Phase
           REAL                                    :: w0

           ALLOCATE(Forces1%excitation(Nw,Nintegration,Nbeta))
           ALLOCATE(Forces1%addedmass(Nw,Nradiation,Nradiation))
           ALLOCATE(Forces1%dampcoef(Nw,Nradiation,Nradiation))
           
           CALL exist_file(trim(wd)//'/results/CA.dat')
           CALL exist_file(trim(wd)//'/results/CM.dat')
           CALL exist_file(trim(wd)//'/results/Fe.dat')

           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/results/CA.dat',ACTION='READ')!Damp. Coef.
           OPEN(NEWUNIT=u2,FILE=trim(wd)//'/results/CM.dat',ACTION='READ')!Added Mass
           READ(u1,*) 
           READ(u2,*) 
           DO I=1,Nw
                READ(u1,*)
                READ(u2,*)
                DO J=1,Nradiation
                   READ(u1,*) ( Forces1%dampcoef(I,J,K),K=1,Nradiation )
                   READ(u2,*) ( Forces1%addedmass(I,J,K),K=1,Nradiation )
                ENDDO
           ENDDO
           CLOSE(u1)
           CLOSE(u2)
           
           OPEN(NEWUNIT=u3,FILE=trim(wd)//'/results/Fe.dat',ACTION='READ')!excitation
           READ(u3,*)
           READ(u3,*)
           DO K=1,Nbeta
                READ(u3,*)
                DO I=1,Nw
                   READ(u3,*) w0,(Amp(J),J=1,Nintegration),(Phase(J),J=1,Nintegration)
                   DO J=1,Nintegration
                      Phase(J)=Phase(J)*PI/180.0
                      Forces1%excitation(I,J,K)=Amp(J)*CEXP(II*Phase(J))
                   ENDDO
                ENDDO
           ENDDO
           CLOSE(u3)
        END SUBROUTINE

        SUBROUTINE Read_Motion(wd,Nw,Nbeta,Nradiation,Motion)
           !Input/Output
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nradiation
           COMPLEX,DIMENSION(Nw,Nradiation,Nbeta),INTENT(OUT):: Motion
           !Local
           INTEGER                                 :: u1,I,J,K
           REAL,DIMENSION(Nradiation)            :: Amp, Phase
           REAL                                    :: w0
           
           CALL exist_file(trim(wd)//'/Motion/RAO.dat')
           OPEN(NEWUNIT=u1,FILE=trim(wd)//'/Motion/RAO.dat',ACTION='READ')!RAO
           READ(u1,*)
           READ(u1,*)
             
           DO K=1,Nbeta
              READ(u1,*)
              DO I=1,Nw
                 READ(u1,*) w0,(Amp(J),J=1,Nradiation),(Phase(J),J=1,Nradiation)
                 DO J=1,Nradiation
                     Phase(J)=Phase(J)*PI/180.0
                     IF (MODULO(J,6)>3.OR.MODULO(J,6)==0) Amp(J)=Amp(J)*PI/180.0
                     Motion(I,J,K)=Amp(J)*CEXP(II*Phase(J))
                 ENDDO
              ENDDO  
           ENDDO
           CLOSE(u1)
        END
        
        SUBROUTINE Read_SourceDistribution(wd,Iw,Nw,Nradiation,Nbeta,Npanels,SourceDistr)
           !INPUT/OUTPUT
           CHARACTER(LEN=*),            INTENT(IN) :: wd
           INTEGER,                     INTENT(IN) :: Nw,Nbeta,Nradiation,Npanels,Iw
           TYPE(TSource),               INTENT(INOUT):: SourceDistr
           !Local
           CHARACTER*5  :: str
           REAL         :: RE,IM
           INTEGER      :: idiffrad,Pbnumber,I,K
           INTEGER      :: u1
           idiffrad=0
           !Radiation problems
           DO K=1,Nradiation
                idiffrad= idiffrad+1
                Pbnumber=(Nbeta+Nradiation)*(Iw-1)+Nbeta+K
                WRITE(str,'(I0.5)') Pbnumber
                CALL exist_file(wd//'/results/sources/sources.'//str//'.dat')
                OPEN(NEWUNIT=u1,FILE=wd//'/results/sources/sources.'//str//'.dat')
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGB(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGS(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                CLOSE(u1)
           ENDDO

           !DIffraction problems
           DO K=1,Nbeta
                idiffrad= idiffrad+1
                Pbnumber=(Nbeta+Nradiation)*(Iw-1)+K
                WRITE(str,'(I0.5)') Pbnumber
                CALL exist_file(wd//'/results/sources/sources.'//str//'.dat')
                OPEN(NEWUNIT=u1,FILE=wd//'/results/sources/sources.'//str//'.dat')
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGB(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                DO I=1,Npanels
                   READ(u1,*) RE,IM
                   SourceDistr%ZIGS(I,idiffrad)=CMPLX(RE,IM)
                ENDDO
                CLOSE(u1)
           ENDDO
        END SUBROUTINE

        SUBROUTINE READ_POTENTIALS_VELOCITIES_BODYWLINE(wd,Nw,Nbeta,NRadiation,NPFlow,      &
                        TotPot,TotVel,RadPot,RadVel,w,k,beta)
          CHARACTER(LEN=*),                         INTENT(IN) ::wd
          INTEGER,                                  INTENT(IN) ::Nw,Nbeta,Nradiation,NPFlow
          COMPLEX,DIMENSION(NPFlow,Nbeta,Nw),       INTENT(OUT)::TotPot
          COMPLEX,DIMENSION(NPFlow,3,Nbeta,Nw),     INTENT(OUT)::TotVel
          COMPLEX,DIMENSION(NPFlow,Nradiation,Nw),  INTENT(OUT)::RadPot
          COMPLEX,DIMENSION(NPFlow,3,Nradiation,Nw),INTENT(OUT)::RadVel
          REAL,DIMENSION(Nbeta),                    INTENT(OUT)::beta
          REAL,DIMENSION(Nw),                       INTENT(OUT)::w,k    !radfreq&wavenumber
          INTEGER                     :: Iw, Ibeta, Irad,Ipanel
          REAL, DIMENSION(NPFlow)     :: REALDATA1,IMDATA1,REALDATA2,IMDATA2,               &
                                         REALDATA3,IMDATA3
          INTEGER                     :: uF1,uF2,uF3,uF4,ILINE
          REAL,DIMENSION(3)           :: wkbeta,wkIrad
          
          OPEN(NEWUNIT=uF1, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotPotFILE,                &
                  STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+2*NPFLOW)
          OPEN(NEWUNIT=uF2, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotVelFILE,                &
                  STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+6*NPFLOW)
          OPEN(NEWUNIT=uF3, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadPotFILE,                &
                  STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+2*NPFLOW)
          OPEN(NEWUNIT=uF4, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadVelFILE,                &
                  STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+6*NPFLOW)

          DO Iw=1,Nw
              DO Ibeta=1,Nbeta
                 ILINE=(Iw-1)*Nbeta+Ibeta
                 READ(uF1,REC=ILINE) wkbeta,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow)
                 TotPot(:,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                 READ(uF2,REC=ILINE) wkbeta,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow),                          &                       
                                (REALDATA2(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA2(Ipanel),Ipanel=1,NPFlow),                          &    
                                (REALDATA3(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA3(Ipanel),Ipanel=1,NPFlow)                             
                 TotVel(:,1,Ibeta,Iw)=CMPLX(REALDATA1,IMDATA1)
                 TotVel(:,2,Ibeta,Iw)=CMPLX(REALDATA2,IMDATA2)
                 TotVel(:,3,Ibeta,Iw)=CMPLX(REALDATA3,IMDATA3)
                 IF (Iw==1) beta(Ibeta)=wkbeta(3)
              ENDDO
                 w(Iw)=wkbeta(1)
                 k(Iw)=wkbeta(2)
              DO Irad=1,NRadiation
                 ILINE=(Iw-1)*Nradiation+Irad
                 READ(uF3,REC=ILINE) wkIrad,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow)
                 RadPot(:,Irad,Iw)=CMPLX(REALDATA1,IMDATA1)
                 READ(uF4,REC=ILINE) wkIrad,                                                &
                                (REALDATA1(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA1(Ipanel),Ipanel=1,NPFlow),                          &                       
                                (REALDATA2(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA2(Ipanel),Ipanel=1,NPFlow),                          &    
                                (REALDATA3(Ipanel),Ipanel=1,NPFlow),                        &
                                (IMDATA3(Ipanel),Ipanel=1,NPFlow)
                 RadVel(:,1,Irad,Iw)=CMPLX(REALDATA1,IMDATA1)
                 RadVel(:,2,Irad,Iw)=CMPLX(REALDATA2,IMDATA2)
                 RadVel(:,3,Irad,Iw)=CMPLX(REALDATA3,IMDATA3)
              ENDDO
          ENDDO
          CLOSE(uF1)
          CLOSE(uF2)
          CLOSE(uF3)
          CLOSE(uF4)

        END SUBROUTINE

        SUBROUTINE READ_GENERALIZED_NORMAL_BODY_dAREA(wd,Npanels,Nintegration,genNormal_dS)
          CHARACTER(LEN=*),                     INTENT(IN)   :: wd
          INTEGER,                              INTENT(IN)   :: Npanels,NIntegration
          REAL, DIMENSION(Nintegration,Npanels),INTENT(INOUT):: genNormal_dS
          !Local
          INTEGER I,J,Ninteg,u
          CALL exist_file(TRIM(wd)//'/Mesh/Integration.dat')

          OPEN(NEWUNIT=u, FILE=TRIM(wd)//'/Mesh/Integration.dat', STATUS='OLD', ACTION='READ')
          READ(u,*) Ninteg
          IF (Ninteg.NE.Nintegration) THEN
            CLOSE(u)
            print*,'Number of rows (Nintegration) in /Mesh/Integration.dat is not correct!'
            STOP
          ENDIF
          DO I = 1, Nintegration
             READ(u, *) (genNormal_dS(I,J), J=1,Npanels)
          END DO
          CLOSE(u)

        END SUBROUTINE

        SUBROUTINE  exist_file(filename)
          CHARACTER(LEN=*),       INTENT(IN) :: filename
          LOGICAL                            ::existfile
          INQUIRE (FILE=filename, EXIST=existfile)       
          IF (.NOT.existfile) THEN
               PRINT*,filename,' data is missing!'
               STOP
          ENDIF
        END SUBROUTINE
        
        
END MODULE
