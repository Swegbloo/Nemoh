MODULE MQSolverPreparation

USE CONSTANTS 
USE MEnvironment,               ONLY: TEnvironment,FunVect_inverseDispersion
USE Elementary_functions,       ONLY: Fun_closest,Fun_MIN,Fun_MAX,              &
                                      CROSS_PRODUCT_COMPLEX
USE MFileDirectoryList 
USE MLogFile
USE linear_interpolation_module
USE MMesh,                      ONLY: Tmesh
USE MFace,                      ONLY: TVFace,TWLine
USE MNemohCal,                  ONLY: TNemCal
USE MReadInputFiles,            ONLY: TMech,TLoad1
IMPLICIT NONE
TYPE TQfreq
     INTEGER                          :: NwQ
     REAL, ALLOCATABLE,DIMENSION(:,:) :: wQ,kQ          !freq for QTF comp
     REAL, ALLOCATABLE,DIMENSION(:,:) :: diffwQ,sumwQ   !diff and sum freq
     INTEGER, ALLOCATABLE,DIMENSION(:)  :: InterpPotSwitch  
END TYPE

TYPE TPotVel
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: TotPot,RadPot ! Total&radiation pot
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: TotVel,RadVel ! Total&radiation vel
END TYPE

CONTAINS



SUBROUTINE PREPARE_POTENTIAL_VELOCITIES(Qfreq,Nw,w,Nbeta,beta, &
                                        NPFlow,Nradiation,datPotVel,datPotVelQ)
        !input/output
        TYPE(TQfreq),                   INTENT(IN)   ::Qfreq
        INTEGER,                        INTENT(IN)   ::Nw,Nbeta
        INTEGER,                        INTENT(IN)   ::NPFlow,Nradiation
        REAL, DIMENSION(Nw),            INTENT(IN)   ::w 
        REAL, DIMENSION(Nbeta),         INTENT(IN)   ::beta 
        TYPE(TPotVel),                  INTENT(INOUT)::datPotVel
        TYPE(TPotVel),                  INTENT(INOUT)::datPotVelQ

        !Local
        INTEGER                    :: Ibeta,Iw1,Iw2,Ipanel,iflag,IwQ
        Type(linear_interp_1d)     :: interpPR,interpVxR,interpVyR,interpVzR
        Type(linear_interp_1d)     :: interpPI,interpVxI,interpVyI,interpVzI
        REAL,DIMENSION(Qfreq%NwQ)  :: TotPotR,TotPotI
        REAL,DIMENSION(3,Qfreq%NwQ):: TotVelR,TotVelI
        
        ALLOCATE(datPotVelQ%TotPot(NPFlow,  Nbeta     ,Qfreq%NwQ))
        ALLOCATE(datPotVelQ%TotVel(NPFlow,3,Nbeta     ,Qfreq%NwQ))
        ALLOCATE(datPotVelQ%RadPot(NPFlow,  Nradiation,Nw))
        ALLOCATE(datPotVelQ%RadVel(NPFlow,3,Nradiation,Nw))

        
        DO Ibeta=1,Nbeta
               IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                   IF(Qfreq%NwQ.NE.Nw) THEN
                        Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                        Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                        datPotVelQ%TotPot(:,   Ibeta,:)=datPotVel%TotPot(:,Ibeta,Iw1:Iw2)
                        datPotVelQ%TotVel(:,:, Ibeta,:)=datPotVel%TotVel(:,:,Ibeta,Iw1:Iw2)
                    ELSE
                        datPotVelQ%TotPot(:,  Ibeta,:)=datPotVel%TotPot(:,Ibeta,:)
                        datPotVelQ%TotVel(:,:,Ibeta,:)=datPotVel%TotVel(:,:,Ibeta,:)
                    ENDIF 
                ELSE
                   DO Ipanel=1,NPFlow
                   !interpolating potential for the wQ rad. frequencies
                   CALL interpPR%initialize(w ,REAL(datPotVel%TotPot(Ipanel,Ibeta,:)),iflag)
                   CALL interpVxR%initialize(w,REAL(datPotVel%TotVel(Ipanel,1,Ibeta,:)),iflag)
                   CALL interpVyR%initialize(w,REAL(datPotVel%TotVel(Ipanel,2,Ibeta,:)),iflag)
                   CALL interpVzR%initialize(w,REAL(datPotVel%TotVel(Ipanel,3,Ibeta,:)),iflag)
                   CALL interpPI%initialize(w ,AIMAG(datPotVel%TotPot(Ipanel,Ibeta,:)),iflag)
                   CALL interpVxI%initialize(w,AIMAG(datPotVel%TotVel(Ipanel,1,Ibeta,:)),iflag)
                   CALL interpVyI%initialize(w,AIMAG(datPotVel%TotVel(Ipanel,2,Ibeta,:)),iflag)
                   CALL interpVzI%initialize(w,AIMAG(datPotVel%TotVel(Ipanel,3,Ibeta,:)),iflag)
                   DO IwQ=1,Qfreq%NwQ
                         CALL interpPR%evaluate(Qfreq%wQ(IwQ,Ibeta), TotPotR(IwQ))
                         CALL interpVxR%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelR(1,IwQ))
                         CALL interpVyR%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelR(2,IwQ))
                         CALL interpVzR%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelR(3,IwQ))
                         CALL interpPI%evaluate(Qfreq%wQ(IwQ,Ibeta), TotPotI(IwQ))
                         CALL interpVxI%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelI(1,IwQ))
                         CALL interpVyI%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelI(2,IwQ))
                         CALL interpVzI%evaluate(Qfreq%wQ(IwQ,Ibeta),TotVelI(3,IwQ))
                   ENDDO
                   datPotVelQ%TotPot(Ipanel,   Ibeta,:)=CMPLX(TotPotR,TotPotI)
                   datPotVelQ%TotVel(Ipanel,:, Ibeta,:)=CMPLX(TotVelR,TotVelI)

                   CALL interpPR%destroy() 
                   CALL interpVxR%destroy()
                   CALL interpVyR%destroy()
                   CALL interpVzR%destroy()
                   CALL interpPI%destroy() 
                   CALL interpVxI%destroy()
                   CALL interpVyI%destroy()
                   CALL interpVzI%destroy()

                   ENDDO
                ENDIF 
        ENDDO
        !keep radiation potential as the calculated potential in NEMOH first order
        !will be taken/interpolated for difference and sum frequency later
        datPotVelQ%RadPot=datPotVel%RadPot
        datPotVelQ%RadVel=datPotVel%RadVel

        !destroy the initial data
        DEALLOCATE(datPotVel%TotPot,datPotVel%TotVel) 
        DEALLOCATE(datPotVel%RadPot,datPotVel%RadVel)
END SUBROUTINE

SUBROUTINE PREPARE_BODY_DISPLACEMENT(Qfreq,Nw,w,Nbeta,Nradiation,NPFlow,Nbodies,      &
                                       Mesh,WLine,InpNEMOHCAL,Motion,BdisplaceQ)
        TYPE(TQfreq),                          INTENT(IN)::Qfreq
        INTEGER,                               INTENT(IN)::Nw,Nbeta,Nradiation
        INTEGER,                               INTENT(IN)::Nbodies,NPFlow
        REAL,DIMENSION(Nw),                    INTENT(IN)::w
        TYPE(TMesh),                           INTENT(IN)::Mesh
        TYPE(TWLine),                          INTENT(IN)::WLine
        TYPE(TNemCal),                         INTENT(IN)::InpNEMOHCAL

        COMPLEX,DIMENSION(Nw,Nbeta,Nradiation),INTENT(IN)::Motion       !RAO
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT):: BdisplaceQ
        !Local
        INTEGER                     :: Ibeta,IdB,Isegline,Iw,Ipanel,iflag
        INTEGER                     :: Iw1,Iw2,IwQ,NPanWL
        REAL,DIMENSION(3)           :: vect_R,XCOG,vect_Rsym
        REAL,DIMENSION(Qfreq%NwQ,3) :: BdispR,BdispI,BdispRsym,BdispIsym
        COMPLEX,DIMENSION(3)        :: vect_Theta     !complex rotation RAO
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: Bdisplace
        Type(linear_interp_1d)      :: interpBdispX_R,interpBdispY_R,interpBdispZ_R
        Type(linear_interp_1d)      :: interpBdispX_I,interpBdispY_I,interpBdispZ_I
        Type(linear_interp_1d)      :: interpBdispXsym_R,interpBdispYsym_R,interpBdispZsym_R
        Type(linear_interp_1d)      :: interpBdispXsym_I,interpBdispYsym_I,interpBdispZsym_I
        ALLOCATE(BdisplaceQ(Qfreq%NwQ,Nbeta,NPFlow,3))
        ALLOCATE(Bdisplace(Nw,Nbeta,NPFlow,3))

        NpanWL=Mesh%Npanels+WLine%NWLineseg
           
        DO Ipanel=1,Mesh%Npanels+WLine%NWLineseg
                IF (Ipanel .LE. Mesh%Npanels) THEN
                  XCOG=InpNEMOHCAL%bodyinput(Mesh%cPanel(Ipanel))%IntCase(4)%Axis(1:3)
                  IdB=6*(Mesh%cPanel(Ipanel)-1)
                  vect_R=Mesh%XM(:,Ipanel)-XCOG
                ELSE
                  Isegline=Ipanel-Mesh%Npanels
                  XCOG=InpNEMOHCAL%bodyinput(Mesh%cPanel(Wline%IndexPanel(Isegline) ))     &
                                                                 %IntCase(4)%Axis(1:3)
                  IdB=6*(Mesh%cPanel(Wline%IndexPanel(Isegline))-1)
                  vect_R=Wline%XM(Isegline,1:3)-XCOG
                ENDIF   

                IF (Mesh%iSym.EQ.1) THEN
                    vect_Rsym=vect_R
                    IF (Ipanel .LE. Mesh%Npanels) THEN
                        vect_Rsym(2)=-Mesh%XM(2,Ipanel)-XCOG(2)
                    ELSE
                        Isegline=Ipanel-Mesh%Npanels
                        vect_Rsym(2)=-Wline%XM(Isegline,2)-XCOG(2)
                    ENDIF   
                ENDIF

                DO Ibeta=1,Nbeta
                   DO Iw=1,Nw
                       vect_Theta=Motion(Iw,Ibeta,IdB+4:IdB+6)
                       Bdisplace(Iw,Ibeta,Ipanel,1:3)=Motion(Iw,Ibeta,IdB+1:IdB+3)              &
                            + CROSS_PRODUCT_COMPLEX(vect_Theta,CMPLX(vect_R,0))
                        IF (Mesh%iSym.EQ.1) THEN
                           Bdisplace(Iw,Ibeta,NpanWL+Ipanel,1:3)=                               &
                                   Motion(Iw,Ibeta,IdB+1:IdB+3)                                 &
                                   + CROSS_PRODUCT_COMPLEX(vect_Theta,CMPLX(vect_Rsym,0))
                        ENDIF
                   ENDDO

                    !Interpolation
                   IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                      IF(Qfreq%NwQ.NE.Nw) THEN
                           Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                           Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                           BdisplaceQ(1:Qfreq%NwQ,Ibeta,Ipanel,1:3)=Bdisplace(Iw1:Iw2,Ibeta,Ipanel,1:3)
                           IF (Mesh%iSym.EQ.1) THEN
                           BdisplaceQ(:,Ibeta,NpanWL+Ipanel,:)=Bdisplace(Iw1:Iw2,Ibeta,NpanWL+Ipanel,:)
                           ENDIF
                       ELSE
                           BdisplaceQ(:,Ibeta,Ipanel,:)=Bdisplace(:,Ibeta,Ipanel,:)
                           IF (Mesh%iSym.EQ.1) THEN
                           BdisplaceQ(:,Ibeta,NpanWL+Ipanel,:)=Bdisplace(:,Ibeta,NpanWL+Ipanel,:)
                           ENDIF
                      ENDIF 
                   
                    ELSE
                         !interpolating displacement for the wQ rad. frequencies
                         CALL interpBdispX_R%initialize(w ,REAL(Bdisplace(:,Ibeta,Ipanel,1)),iflag)
                         CALL interpBdispY_R%initialize(w ,REAL(Bdisplace(:,Ibeta,Ipanel,2)),iflag)
                         CALL interpBdispZ_R%initialize(w ,REAL(Bdisplace(:,Ibeta,Ipanel,3)),iflag)
                         CALL interpBdispX_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,Ipanel,1)),iflag)
                         CALL interpBdispY_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,Ipanel,2)),iflag)
                         CALL interpBdispZ_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,Ipanel,3)),iflag)
      
                         IF (Mesh%iSym.EQ.1) THEN
                            CALL interpBdispXsym_R%initialize(w ,REAL(Bdisplace(:,Ibeta,NpanWL+Ipanel,1)),iflag)
                            CALL interpBdispYsym_R%initialize(w ,REAL(Bdisplace(:,Ibeta,NpanWL+Ipanel,2)),iflag)
                            CALL interpBdispZsym_R%initialize(w ,REAL(Bdisplace(:,Ibeta,NpanWL+Ipanel,3)),iflag)
                            CALL interpBdispXsym_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,NpanWL+Ipanel,1)),iflag)
                            CALL interpBdispYsym_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,NpanWL+Ipanel,2)),iflag)
                            CALL interpBdispZsym_I%initialize(w ,AIMAG(Bdisplace(:,Ibeta,NpanWL+Ipanel,3)),iflag)
                         ENDIF

                         DO IwQ=1,Qfreq%NwQ
                               CALL interpBdispX_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispR(IwQ,1))
                               CALL interpBdispY_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispR(IwQ,2))
                               CALL interpBdispZ_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispR(IwQ,3))
                               CALL interpBdispX_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispI(IwQ,1))
                               CALL interpBdispY_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispI(IwQ,2))
                               CALL interpBdispZ_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispI(IwQ,3))
                               IF (Mesh%iSym.EQ.1) THEN
                               CALL interpBdispXsym_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispRsym(IwQ,1))
                               CALL interpBdispYsym_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispRsym(IwQ,2))
                               CALL interpBdispZsym_R%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispRsym(IwQ,3))
                               CALL interpBdispXsym_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispIsym(IwQ,1))
                               CALL interpBdispYsym_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispIsym(IwQ,2))
                               CALL interpBdispZsym_I%evaluate(Qfreq%wQ(IwQ,Ibeta),BdispIsym(IwQ,3))
                               ENDIF
                         ENDDO
                         BdisplaceQ(:,Ibeta,Ipanel,1)=CMPLX(BdispR(:,1),BdispI(:,1))
                         BdisplaceQ(:,Ibeta,Ipanel,2)=CMPLX(BdispR(:,2),BdispI(:,2))
                         BdisplaceQ(:,Ibeta,Ipanel,3)=CMPLX(BdispR(:,3),BdispI(:,3))
                         IF (Mesh%iSym.EQ.1) THEN
                         BdisplaceQ(:,Ibeta,NpanWL+Ipanel,1)=CMPLX(BdispRsym(:,1),BdispIsym(:,1))
                         BdisplaceQ(:,Ibeta,NPanWL+Ipanel,2)=CMPLX(BdispRsym(:,2),BdispIsym(:,2))
                         BdisplaceQ(:,Ibeta,NPanWL+Ipanel,3)=CMPLX(BdispRsym(:,3),BdispIsym(:,3))
                         ENDIF
                         CALL interpBdispX_R%destroy()
                         CALL interpBdispY_R%destroy()
                         CALL interpBdispZ_R%destroy()
                         CALL interpBdispX_I%destroy()
                         CALL interpBdispY_I%destroy()
                         CALL interpBdispZ_I%destroy()
                         IF (Mesh%iSym.EQ.1) THEN
                         CALL interpBdispXsym_R%destroy()
                         CALL interpBdispYsym_R%destroy()
                         CALL interpBdispZsym_R%destroy()
                         CALL interpBdispXsym_I%destroy()
                         CALL interpBdispYsym_I%destroy()
                         CALL interpBdispZsym_I%destroy()
                         ENDIF 
                   ENDIF
   
                ENDDO
        ENDDO
        DEALLOCATE(Bdisplace)
END SUBROUTINE

SUBROUTINE PREPARE_INERTIA_FORCES(MechCoef,Motion,Nw,Nbeta,Nradiation,Nintegration,&
                w,Qfreq,Forces1,InertiaForceQ)
        !FInertia=[FHydrodynamic1]+[Fhydrostatic1]
        !       =[Fexc1+Ma*ddXidt2+B*dXidt]+[-K*Xi]
        INTEGER,                        INTENT(IN):: Nradiation,Nintegration
        INTEGER,                        INTENT(IN):: Nw,Nbeta
        TYPE(TLoad1),                   INTENT(IN):: Forces1
        TYPE(TMech),                    INTENT(IN):: MechCoef
        REAL,DIMENSION(Nw),             INTENT(IN):: w
        TYPE(TQfreq),                   INTENT(IN):: Qfreq
        COMPLEX,DIMENSION(Nw,Nbeta,Nradiation),INTENT(IN)::Motion       !RAO
        COMPLEX,DIMENSION(Qfreq%NwQ,Nbeta,Nintegration),                          &
                                        INTENT(OUT)::InertiaForceQ
        !local 
        COMPLEX,DIMENSION(Nintegration)          ::HydroDynForce,HydroStatForce
        COMPLEX,DIMENSION(Nw,Nintegration)       ::InertiaForce
        REAL,DIMENSION(Nradiation,Nradiation)    ::StiffMat
        REAL,DIMENSION(Nradiation,Nradiation)    ::AddedMass,DampCoef
        COMPLEX,DIMENSION(Nintegration)          ::ExcitForce
        COMPLEX,DIMENSION(Nintegration)          ::Xi,ddXidt2,dXidt
        INTEGER                                  ::Iw,Ibeta,Iinteg,iflag
        REAL,DIMENSION(Qfreq%NwQ)                ::InForceR,InForceI
        INTEGER                                  ::Iw1,Iw2,IwQ 
        Type(linear_interp_1d)                   :: interpIF_R,interpIF_I
            
        StiffMat=MechCoef%StiffMat+MechCoef%StiffMat_EXT
        DampCoef=MechCoef%DampCoefMat_EXT
        
        DO Ibeta=1,Nbeta
           DO Iw=1,Nw
           ExcitForce   = Forces1%excitation(Iw,:,Ibeta)
           DampCoef     = MechCoef%DampCoefMat_EXT+Forces1%dampcoef(Iw,:,:)
           AddedMass    = Forces1%addedmass(Iw,:,:)
           Xi           = Motion(Iw,Ibeta,:)
           ddXidt2      = -w(Iw)*w(Iw)*Xi
           dXidt        = -II*w(Iw)*Xi
           HydrodynForce= ExcitForce                                           &
                          +MATMUL(CMPLX(AddedMass,0.),ddXidt2)                 &
                          +MATMUL(CMPLX(DampCoef,0.),dXidt)                    
           HydroStatForce=-MATMUL(CMPLX(StiffMat,0.),Xi)
           InertiaForce(Iw,:)=HydroDynForce+HydroStatForce
           ENDDO
            !Interpolation
            IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
               IF(Qfreq%NwQ.NE.Nw) THEN
                    Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                    Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                    InertiaForceQ(:,Ibeta,:)=InertiaForce(Iw1:Iw2,:)
                ELSE
                    InertiaForceQ(:,Ibeta,:)=InertiaForce(:,:)
                ENDIF 
             ELSE
             !interpolating displacement for the wQ rad. frequencies
             DO Iinteg=1,Nintegration
             CALL interpIF_R%initialize(w ,REAL(InertiaForce(:,Iinteg)),iflag)
             CALL interpIF_I%initialize(w ,AIMAG(InertiaForce(:,Iinteg)),iflag)
                DO IwQ=1,Qfreq%NwQ
                       CALL interpIF_R%evaluate(Qfreq%wQ(IwQ,Ibeta),InForceR(IwQ))
                       CALL interpIF_I%evaluate(Qfreq%wQ(IwQ,Ibeta),InForceI(IwQ))
                ENDDO
             InertiaForceQ(:,Ibeta,Iinteg)=CMPLX(InForceR,InForceI)
             CALL interpIF_R%destroy()
             CALL interpIF_I%destroy()
             ENDDO
           ENDIF
        ENDDO
END SUBROUTINE

SUBROUTINE PREPARE_ROTATION_ANGLES(Motion,Nw,Nbeta,Nradiation,&
                Nbodies,w,Qfreq,RotAnglesQ)
        INTEGER,                        INTENT(IN):: Nradiation
        INTEGER,                        INTENT(IN):: Nw,Nbeta,Nbodies
        REAL,DIMENSION(Nw),             INTENT(IN):: w
        TYPE(TQfreq),                   INTENT(IN):: Qfreq
        COMPLEX,DIMENSION(Nw,Nbeta,Nradiation),INTENT(IN)::Motion       !RAO
        COMPLEX,DIMENSION(Qfreq%NwQ,Nbeta,3*Nbodies),                    &
                                        INTENT(OUT)::RotAnglesQ
        !local 
        INTEGER                                  ::Iw,Ibeta,Iinteg,iflag
        INTEGER                                  ::Itheta0,Itheta,Ibody
        REAL,DIMENSION(Qfreq%NwQ)                ::RotAngleR,RotAngleI
        INTEGER                                  ::Iw1,Iw2,IwQ 
        Type(linear_interp_1d)                   ::interpROT_R,interpROT_I
        COMPLEX,DIMENSION(Nw,3)                  ::vect_Theta
            
        DO Ibody=1,Nbodies
           DO Ibeta=1,Nbeta
              Itheta0=6*(Ibody-1)
              DO Iw=1,Nw
              vect_Theta(Iw,:)=Motion(Iw,Ibeta,Itheta0+4:Itheta0+6)
              ENDDO
              Itheta0=3*(Ibody-1)
              !Interpolation
              IF(Qfreq%InterpPotSwitch(Ibeta)==0) THEN
                IF(Qfreq%NwQ.NE.Nw) THEN
                 Iw1=Fun_closest(Nw,w,Qfreq%wQ(1,Ibeta))
                 Iw2=Fun_closest(Nw,w,Qfreq%wQ(Qfreq%NwQ,Ibeta))
                 RotAnglesQ(:,Ibeta,Itheta0+1:Itheta0+3)=vect_Theta(Iw1:Iw2,:)
                ELSE
                 RotAnglesQ(:,Ibeta,Itheta0+1:Itheta0+3)=vect_Theta(:,:)
                ENDIF 
               ELSE
               !interpolating displacement for the wQ rad. frequencies
                DO Itheta=1,3 !theta_x,theta_y,theta_z
                CALL interpROT_R%initialize(w ,REAL(vect_Theta(:,Itheta)),iflag)
                CALL interpROT_I%initialize(w ,AIMAG(vect_Theta(:,Itheta)),iflag)
                 DO IwQ=1,Qfreq%NwQ
                  CALL interpROT_R%evaluate(Qfreq%wQ(IwQ,Ibeta),RotAngleR(IwQ))
                  CALL interpROT_I%evaluate(Qfreq%wQ(IwQ,Ibeta),RotAngleI(IwQ))
                 ENDDO
                RotAnglesQ(:,Ibeta,Itheta0+Itheta)=CMPLX(RotAngleR,RotAngleI)
                CALL interpROT_R%destroy()
                CALL interpROT_I%destroy()
                ENDDO
              ENDIF
           ENDDO
        ENDDO
END SUBROUTINE


SUBROUTINE Discretized_omega_wavenumber_for_QTF               &
                        (Nw,w,kw,NwQ,wQminmax, Nbeta,beta,BdySpeed,    &
                          depth,g,QFreq)
           INTEGER,                         INTENT(IN) ::Nw,NwQ,Nbeta
           REAL,        DIMENSION(2),       INTENT(IN) ::wQminmax
           REAL,        DIMENSION(Nbeta),   INTENT(IN) ::beta
           REAL,                            INTENT(IN) ::BdySpeed,     &
                                                         depth,g
           REAL,        DIMENSION(Nw),      INTENT(IN) ::w,kw
           TYPE(TQfreq),                    INTENT(INOUT)::Qfreq
           
           INTEGER                                    ::Ibeta,Iw
           REAL                                       ::wQmax,wQmin,dwQ
           REAL,DIMENSION(NwQ)                        ::wQ_temp,kQ_temp
           INTEGER                                    ::Ind1,Ind2
           INTEGER                                    ::InterpPotSwitch

           !In case NO forward speed, wQ is not depend on beta     
           !Allocation
           ALLOCATE(Qfreq%wQ(NwQ,Nbeta),Qfreq%kQ(NwQ,Nbeta))
           ALLOCATE(Qfreq%diffwQ(NwQ,Nbeta),Qfreq%sumwQ(NwQ-1,Nbeta))
           ALLOCATE(Qfreq%InterpPotSwitch(Nbeta))

           wQmin=wQminmax(1)
           wQmax=wQminmax(2)

           IF (NwQ.GT.1) THEN
                dwQ=(wQmax-wQmin)/(NwQ-1)
                IF (ABS(w(2)-w(1)-dwQ).LT.0.01) dwQ=w(2)-w(1)
                DO Iw=1,NwQ
                    wQ_temp(Iw)=wQmin+dwQ*(Iw-1)                
                END DO
            ELSE
             WRITE(*,*) &
             "ERROR: QTF can not be run with 1 frequency only!"
             STOP
            END IF
           IF (w(2)-w(1)==dwQ) THEN
                Ind1=Fun_closest(Nw,w,wQmin)
                Ind2=Fun_closest(Nw,w,wQmax)
                kQ_temp=kw(Ind1:Ind2)
                InterpPotSwitch=0
           ELSE
                kQ_temp=FunVect_inverseDispersion(NwQ,wQ_temp,depth,g)!wavenumbers
                InterpPotSwitch=1
           ENDIF
           
           DO Ibeta=1,Nbeta
             Qfreq%wQ(:,Ibeta)=wQ_temp-kQ_temp*BdySpeed*COS(beta(Ibeta))
             IF (BdySpeed>0.) THEN
                Qfreq%kQ(:,Ibeta)=FunVect_inverseDispersion(                 &
                                NwQ,Qfreq%wQ(:,Ibeta),depth,g)!wavenumbers
                Qfreq%InterpPotSwitch(Ibeta)=1
             ELSE
                Qfreq%kQ(:,Ibeta)=kQ_temp
                Qfreq%InterpPotSwitch(Ibeta)=0
                IF (InterpPotSwitch==1) Qfreq%InterpPotSwitch(Ibeta)=1
             ENDIF
             Qfreq%diffwQ(1,Ibeta)=0
             Qfreq%diffwQ(2:NwQ,Ibeta)=Qfreq%wQ(2:NwQ,Ibeta)-Qfreq%wQ(1,Ibeta)
             Qfreq%sumwQ(1:NwQ-1,Ibeta)=Qfreq%wQ(1:NwQ-1,Ibeta)+Qfreq%wQ(1,Ibeta)
             Qfreq%NwQ=NwQ
             IF (Fun_MIN(NwQ-1,Qfreq%diffwQ(2:NwQ,Ibeta))<w(1)) THEN
                    print*,'INPUT ERROR: Min. diff. rad freq < w(1) in QTFpreproc data!'
                    STOP
             ENDIF
             IF (Fun_MAX(NwQ,Qfreq%sumwQ(1:NwQ-1,Ibeta))>w(Nw)) THEN
                    print*,'INPUT ERROR: Max. sum rad freq > w(Nw) in QTFpreproc data!'
                    STOP
             ENDIF
           ENDDO
END SUBROUTINE

SUBROUTINE CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT(Mesh,Nintegration,WLine, &
                                                InpNEMOHCAL,genNormalWLine_dGamma)
        USE MNemohCal,          ONLY:TNemCal
        IMPLICIT NONE
        TYPE(TMesh),                                  INTENT(IN):: Mesh
        TYPE(TWLINE),                                 INTENT(IN):: WLine
        TYPE(TNemCal),                                INTENT(IN):: InpNEMOHCAL
        INTEGER,                                      INTENT(IN):: Nintegration
        REAL, DIMENSION(Nintegration,Wline%NWLineseg*2**Mesh%iSym), &
                                              INTENT(INOUT) ::genNormalWLine_dGamma 
        !local
        INTEGER IdBody, IDMode,c,indsum
        REAL, DIMENSION(Wline%NWLineseg*2**Mesh%iSym) :: NDGamma

        indsum=1 
        DO IdBody=1,InpNEMOHCAL%Nbodies 
           DO IdMode=1,InpNEMOHCAL%bodyinput(IdBody)%NIntegration
           CALL ComputeNDGamma_WLINE(Mesh,WLINE,IdBody,                                 &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%ICase,                    &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Direction(1:3),           &
                InpNEMOHCAL%bodyinput(IdBody)%IntCase(IdMode)%Axis(1:3),NDGamma)
                DO c=1,Wline%NWLineseg*2**Mesh%iSym
                genNormalWLine_dGamma(indsum,c)=NDGamma(c)
                END DO
                indsum=indsum+1
            END DO
        END DO
END SUBROUTINE

!-- SUBROUTINE ComputeNDGamma waterline
    SUBROUTINE ComputeNDGamma_WLINE(Mesh,WLINE,c,iCase,Direction,Axis,NDGamma)  
    USE MMesh
    USE MFace
    IMPLICIT NONE
    TYPE(TMesh)  :: Mesh
    TYPE(TWline) :: Wline
    INTEGER :: c,iCase
    REAL,DIMENSION(3) :: Direction,Axis
    REAL,DIMENSION(*) :: NDGamma
    REAL,DIMENSION(3) :: VEL
    INTEGER :: i,ipanel
    SELECT CASE (iCase)
    CASE (1)
        DO i=1,Wline%NWLineseg
            ipanel=Wline%IndexPanel(i)
            IF (Mesh%cPanel(ipanel).EQ.c) THEN
                VEL(1)=Direction(1)
                VEL(2)=Direction(2)
                VEL(3)=0.
                NDGamma(i)=(Mesh%N(1,ipanel)*VEL(1)+Mesh%N(2,ipanel)*VEL(2)+Mesh%N(3,ipanel)*VEL(3))&
                                *Wline%SegLength(i)
            ELSE
                NDGamma(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(ipanel).EQ.c) THEN
                    VEL(1)=Direction(1)
                    VEL(2)=Direction(2)
                    VEL(3)=0
                    NDGamma(i+WLINE%NWLineseg)=(Mesh%N(1,ipanel)*VEL(1)-Mesh%N(2,ipanel)*VEL(2) &
                                                        +Mesh%N(3,ipanel)*VEL(3))*Wline%Seglength(i)
                 ELSE
                    NDGamma(i+WLINE%NWLineseg)=0.
                 END IF          
            END IF
        END DO
    CASE (2)
        DO i=1,Wline%NWLineseg
            ipanel=Wline%IndexPanel(i)
            IF (Mesh%cPanel(ipanel).EQ.c) THEN
                VEL(1)=Direction(2)*(Mesh%XM(3,ipanel)-Axis(3))-Direction(3)*(Mesh%XM(2,ipanel)-Axis(2))
                VEL(2)=Direction(3)*(Mesh%XM(1,ipanel)-Axis(1))-Direction(1)*(Mesh%XM(3,ipanel)-Axis(3))
                VEL(3)=Direction(1)*(Mesh%XM(2,ipanel)-Axis(2))-Direction(2)*(Mesh%XM(1,ipanel)-Axis(1))  
                NDGamma(i)=(Mesh%N(1,ipanel)*VEL(1)+Mesh%N(2,ipanel)*VEL(2)+Mesh%N(3,ipanel)*VEL(3))&
                                        *Wline%SegLength(i)
            ELSE
                NDGamma(i)=0.
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(ipanel).EQ.c) THEN
                    VEL(1)=Direction(2)*(Mesh%XM(3,ipanel)-Axis(3))-Direction(3)*(-Mesh%XM(2,ipanel)-Axis(2))
                    VEL(2)=Direction(3)*(Mesh%XM(1,ipanel)-Axis(1))-Direction(1)*(Mesh%XM(3,ipanel)-Axis(3))
                    VEL(3)=Direction(1)*(-Mesh%XM(2,ipanel)-Axis(2))-Direction(2)*(Mesh%XM(1,ipanel)-Axis(1))  
                    NDGamma(i+Wline%NWLineseg)=(Mesh%N(1,ipanel)*VEL(1)-Mesh%N(2,ipanel)*VEL(2)             &
                                                                +Mesh%N(3,ipanel)*VEL(3))*Wline%SegLength(i)
                ELSE
                    NDGamma(i+Wline%NWLineseg)=0.
                 END IF 
            END IF
        END DO
    CASE (3)
        WRITE(*,*) 'Error: force case 3 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE


SUBROUTINE WRITE_QTFSOLVERLOGFILE(wd,Nbeta,beta,Qfreq)
       CHARACTER(LEN=*),             INTENT(IN)::wd
       TYPE(TQfreq),                 INTENT(IN)::Qfreq
       INTEGER,                      INTENT(IN)::Nbeta
       REAL,DIMENSION(Nbeta),        INTENT(IN)::beta
       CHARACTER(LEN=1000)                     ::LogTextToBeWritten
       INTEGER Ibeta
            
       WRITE(*, *) ' '
       WRITE(LogTextToBeWritten,*) '----Solver (QTF Module)---'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
       DO Ibeta=1,Nbeta
       WRITE(LogTextToBeWritten,'(A,F7.3,A,I4,A,3(F7.3,A))') ' Beta= ',beta(Ibeta)*180/PI,      &
               ' deg, NFreq= ',Qfreq%NwQ, ', omega = (', Qfreq%wQ(1,Ibeta),':'                  &
               ,Qfreq%wQ(2,Ibeta)-Qfreq%wQ(1,Ibeta),':',Qfreq%wQ(Qfreq%NwQ,Ibeta),') rad/s'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
       ENDDO
       WRITE(LogTextToBeWritten,*) '--------------------------'
       CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
END SUBROUTINE

END MODULE
