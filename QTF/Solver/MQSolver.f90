Module MQSolver
USE MFace,               ONLY:TVFace,TWLine
USE MMesh
USE MQSolverPreparation, ONLY:TPotVel,TQfreq
USE CONSTANTS          , ONLY:II
USE Elementary_functions,ONLY:CROSS_PRODUCT_COMPLEX,Fun_closest,CIH
USE MEnvironment,        ONLY:Fun_Dispersion,TEnvironment
USE linear_interpolation_module
IMPLICIT NONE
CONTAINS
  SUBROUTINE COMPUTATION_QTF_QUADRATIC(Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                       Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,&
                                       rho,g,datPotVelQ,BdisplaceQ,genNormal_dS, &
                                       genNormalWLine_dGamma,wQ,beta,            &
                                       InertiaForceQ,RotAnglesQ,QTFDuok)
          !Input/output
          INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
          INTEGER,                              INTENT(IN) :: Nintegration,NPFlow
          INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nbodies
          REAL,                                 INTENT(IN) :: rho,g
          TYPE(TMesh),                          INTENT(IN) :: Mesh
          TYPE(TVFace),                         INTENT(IN) :: VFace
          TYPE(TWLine),                         INTENT(IN) :: WLine
          TYPE(TPotVel),                        INTENT(IN) :: datPotVelQ 
          COMPLEX,DIMENSION(NwQ,Nbeta,NPFlow,3),INTENT(IN) :: BdisplaceQ
          REAL,DIMENSION(Nintegration,Mesh%Npanels*2**Mesh%Isym),                 &
                                                INTENT(IN) :: genNormal_dS
          REAL,DIMENSION(Nintegration,WLine%NWLineSeg*2**Mesh%Isym),              &
                                                INTENT(IN) :: genNormalWLine_dGamma
          REAL,DIMENSION(NwQ,Nbeta),            INTENT(IN) :: wQ
          REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
          COMPLEX,DIMENSION(NwQ,Nbeta,Nintegration),                              &
                                                INTENT(IN) :: InertiaForceQ
          COMPLEX,DIMENSION(NwQ,Nbeta,3*Nbodies),INTENT(IN):: RotAnglesQ        
          COMPLEX,DIMENSION(Nintegration,2),    INTENT(OUT):: QTFDuok
          !Local
          COMPLEX,DIMENSION(NPFlow)           ::TotPot_Iw1,TotPot_Iw2
          COMPLEX,DIMENSION(NPFlow,3)         ::TotVel_Iw1,TotVel_Iw2
          COMPLEX,DIMENSION(NPFlow,3)         ::Bdisplace_Iw1,Bdisplace_Iw2
          COMPLEX,DIMENSION(Nintegration)     ::InertiaForce_Iw1,InertiaForce_Iw2
          COMPLEX,DIMENSION(3*Nbodies)        ::RotAngles_Iw1,RotAngles_Iw2
          COMPLEX                             ::eta_Iw1,eta_Iw2,zeta_Iw1,zeta_Iw2
          INTEGER                             ::Npanels,Ipanel,Iinteg,NpanWLin
          COMPLEX                             ::quad_M,quad_P
          REAL                                ::w1,w2
          INTEGER                             ::Iwline,Ibody,Itheta0
          INTEGER                             ::Isym
          INTEGER,DIMENSION(2)                ::Ipanelinit,Iwlineinit
          
          w1=wQ(Iw1,Ibeta1)
          w2=wQ(Iw2,Ibeta2)

          TotPot_Iw1=datPotVelQ%TotPot(1:NPFlow,Ibeta1,Iw1)
          TotVel_Iw1=datPotVelQ%TotVel(1:NPFlow,1:3,Ibeta1,Iw1)
  
          TotPot_Iw2=datPotVelQ%TotPot(:,Ibeta2,Iw2)
          TotVel_Iw2=datPotVelQ%TotVel(1:NPFlow,1:3,Ibeta1,Iw2)

          Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,1:NPFlow,1:3)
          Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,1:NPFlow,1:3)
                
          InertiaForce_Iw1=InertiaForceQ(Iw1,Ibeta1,:)
          InertiaForce_Iw2=InertiaForceQ(Iw2,Ibeta2,:)

          RotAngles_Iw1=RotAnglesQ(Iw1,Ibeta1,:)
          RotAngles_Iw2=RotAnglesQ(Iw2,Ibeta2,:)

          Npanels=Mesh%Npanels
          NpanWlin=Npanels+WLine%NWLineSeg
          Ipanelinit=0
          DO Iinteg=1,Nintegration
           QTFDuok(Iinteg,1)=0
           QTFDuok(Iinteg,2)=0
           !integration over body panels
           DO Ipanel=1,Npanels
            IF (Mesh%XM(3, Ipanel)<0.) THEN      !dont compute at lid panels
              DO Isym=1,1+Mesh%Isym
                Ipanelinit=Function_SymIpanel(NpanWlin,Npanels,Isym)
                ! term: product(grad_phi_1,grad_phi2)
                quad_M=DOT_PRODUCT_DIFF_BIHARM(                                 &
                  TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw2(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3)
                quad_P=DOT_PRODUCT_SUM_BIHARM(                                  &
                  TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw2(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                           &
                  TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3)

                QTFDuok(Iinteg,1)=QTFDuok(Iinteg,1)+0.5*rho*                    &
                                quad_M*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
                QTFDuok(Iinteg,2)=QTFDuok(Iinteg,2)+0.5*rho*                    &
                                quad_P*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
              !  ! term: product(displacement,dtgradPhi)
              !  quad_M=DOT_PRODUCT_DIFF_BIHARM(                                 &
              !          Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),                  &
              !          Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),                  &
              !          -II*w1*TotVel_Iw1(Ipanelinit(1)+Ipanel,:),              &
              !          -II*w2*TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3)
              !  quad_P=DOT_PRODUCT_SUM_BIHARM(                                  &
              !          Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),                  &
              !          Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),                  &
              !          -II*w1*TotVel_Iw1(Ipanelinit(1)+Ipanel,:),              &
              !          -II*w2*TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3)
              !  QTFDuok(Iinteg,1)=QTFDuok(Iinteg,1)+rho*                        &
              !                  quad_M*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
              !  QTFDuok(Iinteg,2)=QTFDuok(Iinteg,2)+rho*                        &
              !                  quad_P*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
                ! IF (Iinteg==1) THEN
                ! print*,Ipanel,quad_M*genNormal_dS(Iinteg,Ipanel)*rho,              &
                !                    quad_P*genNormal_dS(Iinteg,Ipanel)*rho
                ! ENDIF
              ENDDO
            ENDIF
           ENDDO
           !integration over waterline segments
          ! DO Iwline=1,WLine%NWLineSeg
          !     DO Isym=1,1+Mesh%Isym
          !      Iwlineinit=Function_SymIwline(NpanWlin,Npanels,WLine%NWLineSeg,Isym)
          !      !term: [eta-zeta]^2
          !      zeta_Iw1=Bdisplace_Iw1(Iwlineinit(1)+Iwline,3)    !vertical Waterline displ.
          !      zeta_Iw2=Bdisplace_Iw2(Iwlineinit(1)+Iwline,3)  
          !      eta_Iw1 =II*w1*TotPot_Iw1(Iwlineinit(1)+Iwline)/g !wave elevation
          !      eta_Iw2 =II*w2*TotPot_Iw2(Iwlineinit(1)+Iwline)/g
          !      quad_M=DOT_PRODUCT_DIFF_BIHARM(                                    &
          !           eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,                            &
          !           eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,1)                 
          !      quad_P=DOT_PRODUCT_SUM_BIHARM(                                     &
          !           eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,                            &
          !           eta_Iw1-zeta_Iw1,eta_Iw2-zeta_Iw2,1)    
          !      QTFDuok(Iinteg,1)=QTFDuok(Iinteg,1)-0.5*rho*g*quad_M               &
          !                      *genNormalWLine_dGamma(Iinteg,Iwlineinit(2)+Iwline)
          !      QTFDuok(Iinteg,2)=QTFDuok(Iinteg,2)-0.5*rho*g*quad_P               &
          !                      *genNormalWLine_dGamma(Iinteg,Iwlineinit(2)+Iwline)
          ! !    IF (Iinteg==1) THEN
          ! !    print*,Iwline,-0.5*rho*g*quad_M*genNormalWLine_dGamma(Iinteg,Iwline),&
          ! !            -0.5*rho*g*quad_P*genNormalWLine_dGamma(Iinteg,Iwline)
          ! !    ENDIF
          !    ENDDO
          ! ENDDO
         ENDDO

        !   ! print*,Iw1,Iw2,QTFDuok(1,1)
        ! !term: Matrix product Rotation matrix and  Inertia Force
        !  DO Ibody=1,Nbodies
        !     Itheta0=3*(Ibody-1)
        !     Iinteg=6*(Ibody-1)
        !     !For diff freq
        !     !R(F_I) for translation modes (1,2,3) of FI
        !     QTFDuok(Iinteg+1:Iinteg+3,1)= QTFDuok(Iinteg+1:Iinteg+3,1)+      &
        !              CROSS_PRODUCT_DIFF_BIHARM(                              &
        !              RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
        !              RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
        !              InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
        !              InertiaForce_Iw2(Iinteg+1:Iinteg+3))                  
        !     !R(F_I) for rotation modes (4,5,6) of FI
        !     QTFDuok(Iinteg+4:Iinteg+6,1)=QTFDuok(Iinteg+4:Iinteg+6,1)+       &
        !              CROSS_PRODUCT_DIFF_BIHARM(                              &
        !              RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
        !              RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
        !              InertiaForce_Iw1(Iinteg+4:Iinteg+6),                    &
        !              InertiaForce_Iw2(Iinteg+4:Iinteg+6))                           
        !     !For sum freq
        !     !R(F_I) for translation modes (1,2,3) of FI
        !     QTFDuok(Iinteg+1:Iinteg+3,2)=QTFDuok(Iinteg+1:Iinteg+3,2)+       &
        !              CROSS_PRODUCT_SUM_BIHARM(                               &
        !              RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
        !              RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
        !              InertiaForce_Iw1(Iinteg+1:Iinteg+3),                    &
        !              InertiaForce_Iw2(Iinteg+1:Iinteg+3))                   
        !     !R(F_I) for rotation modes (4,5,6) of FI
        !     QTFDuok(Iinteg+4:Iinteg+6,2)=QTFDuok(Iinteg+4:Iinteg+6,2)+       &
        !              CROSS_PRODUCT_SUM_BIHARM(                               &
        !              RotAngles_Iw1(Itheta0+1:Itheta0+3),                     &
        !              RotAngles_Iw2(Itheta0+1:Itheta0+3),                     &
        !              InertiaForce_Iw1(Iinteg+4:Iinteg+6),                    &
        !              InertiaForce_Iw2(Iinteg+4:Iinteg+6))   
        !    ! print*, QTFDuok(1,1),QTFDuok(1,2)
        !  ENDDO
        QTFDuok(:,:)=QTFDuok(:,:)/2 !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
  END SUBROUTINE
 
  SUBROUTINE COMPUTATION_QTF_POTENTIAL_BODYFORCE                              &
                                   (Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,       &
                                    Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,&
                                    env,datPotVelQ,BdisplaceQ,genNormal_dS,   &
                                    Nw,w,Qfreq,beta,RotAnglesQ,               &
                                    QTFHasbo)
          !Input/output
          INTEGER,                              INTENT(IN) :: Iw1,Iw2,Ibeta1,Ibeta2
          INTEGER,                              INTENT(IN) :: Nintegration,NPFlow
          INTEGER,                              INTENT(IN) :: NwQ,Nbeta,Nbodies,Nw
          TYPE(TMesh),                          INTENT(IN) :: Mesh
          TYPE(TVFace),                         INTENT(IN) :: VFace
          TYPE(TWLine),                         INTENT(IN) :: WLine
          TYPE(TPotVel),                        INTENT(IN) :: datPotVelQ 
          TYPE(TQfreq),                         INTENT(IN) :: Qfreq
          TYPE(TEnvironment),                   INTENT(IN) :: Env
          COMPLEX,DIMENSION(NwQ,Nbeta,NPFlow,3),INTENT(IN) :: BdisplaceQ
          COMPLEX,DIMENSION(NwQ,Nbeta,3*Nbodies),INTENT(IN):: RotAnglesQ        
          REAL,DIMENSION(Nintegration,Mesh%Npanels*2**Mesh%Isym),                 &
                                                INTENT(IN) :: genNormal_dS
          REAL,DIMENSION(Nw),                   INTENT(IN) :: w
          REAL,DIMENSION(Nbeta),                INTENT(IN) :: beta
          COMPLEX,DIMENSION(Nintegration,2),    INTENT(OUT):: QTFHasbo
          !Local
          COMPLEX,DIMENSION(NPFlow)           ::TotPot_Iw1,TotPot_Iw2
          COMPLEX,DIMENSION(NPFlow,3)         ::TotVel_Iw1,TotVel_Iw2
          COMPLEX,DIMENSION(NPFlow,3)         ::Bdisplace_Iw1,Bdisplace_Iw2
          COMPLEX,DIMENSION(NPFlow,3)         ::dtBdisplace_Iw1,dtBdisplace_Iw2
          COMPLEX,DIMENSION(3*Nbodies)        ::RotAngles_Iw1,RotAngles_Iw2
          COMPLEX                             ::eta_Iw1,eta_Iw2,zeta_Iw1,zeta_Iw2
          INTEGER                             ::Npanels,Ipanel,Iinteg,NpanWLin
          REAL                                ::w1,w2,delw,sumw,k1,k2
          REAL                                ::rho,depth
          INTEGER                             ::Iwline,Ibody,Itheta0
          COMPLEX,DIMENSION(2)                ::Radpot 
          COMPLEX,DIMENSION(3,2)              ::Radvel
          COMPLEX,Dimension(2)                ::Pressure_I
          INTEGER                             ::InterpSwitch
          REAL,DIMENSION(3)                   ::XM,Normal_Vect
          COMPLEX,DIMENSION(2)                ::PHI_I
          COMPLEX,DIMENSION(3,2)              ::GradPHI_I
          COMPLEX,DIMENSION(2)                ::dnPhi_I
          COMPLEX                             ::QB1_M,QB1_P
          COMPLEX                             ::QB2_M,QB2_P
          COMPLEX,DIMENSION(3)                ::RIw1_n0,RIw2_n0
          COMPLEX,DIMENSION(3)                ::dtDisp_GradPhi_I_Iw1, &
                                                dtDisp_GradPhi_I_Iw2
          COMPLEX,DIMENSION(3)                ::dGAMMA_Vect,QB2V_M,QB2V_P
          INTEGER                             ::Isym
          INTEGER,DIMENSION(2)                ::Ipanelinit,Iwlineinit

          rho=Env%rho
          depth=Env%depth
          w1=Qfreq%wQ(Iw1,Ibeta1)
          w2=Qfreq%wQ(Iw2,Ibeta2)
          delw=w1-w2
          sumw=w1+w2
          k1=Qfreq%kQ(Iw1,Ibeta1)
          k2=Qfreq%kQ(Iw2,Ibeta2)

          TotPot_Iw1=datPotVelQ%TotPot(1:NPFlow,Ibeta1,Iw1)
          TotVel_Iw1=datPotVelQ%TotVel(1:NPFlow,1:3,Ibeta1,Iw1)
  
          TotPot_Iw2=datPotVelQ%TotPot(:,Ibeta2,Iw2)
          TotVel_Iw2=datPotVelQ%TotVel(1:NPFlow,1:3,Ibeta1,Iw2)

          Bdisplace_Iw1=BdisplaceQ(Iw1,Ibeta1,1:NPFlow,1:3)
          Bdisplace_Iw2=BdisplaceQ(Iw2,Ibeta2,1:NPFlow,1:3)

          dtBdisplace_Iw1=-II*w1*Bdisplace_Iw1
          dtBdisplace_Iw2=-II*w2*Bdisplace_Iw2

          RotAngles_Iw1=RotAnglesQ(Iw1,Ibeta1,:)
          RotAngles_Iw2=RotAnglesQ(Iw2,Ibeta2,:)
          
          InterpSwitch=Qfreq%InterpPotSwitch(Ibeta1)                      &
                        +Qfreq%InterpPotSwitch(Ibeta2)
          
          Npanels=Mesh%Npanels
          NpanWlin=Npanels+WLine%NWLineSeg
          DO Iinteg=1,Nintegration
           QTFHasbo(Iinteg,1)=0
           QTFHasbo(Iinteg,2)=0
           !integration over body panels
           DO Ipanel=1,Npanels
            IF (Mesh%XM(3, Ipanel)<0.) THEN  !dont compute at lid panels
            XM=Mesh%XM(1:3, Ipanel)
            Normal_Vect=Mesh%N(1:3,Ipanel)
            DO Isym=1,1+Mesh%Isym
              IF (Isym==2) THEN
                 XM(2)=-XM(2)
                 Normal_Vect(2)=-Normal_Vect(2)
              ENDIF
              Ipanelinit=Function_SymIpanel(NpanWlin,Npanels,Isym)

              !interpolating radpot and radvel at delw and sumw        
              Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
                    (Nw,w,datPotVelQ%RadPot(Ipanelinit(1)+Ipanel,Iinteg,:),&
                        InterpSwitch,delw,sumw)
              RadVel(:,1:2)=INTERP_RADIATION_VELOCITY                      &
                  (Nw,w,datPotVelQ%RadVel(Ipanelinit(1)+Ipanel,:,Iinteg,:),&
                        InterpSwitch,delw,sumw)
              !------------------------------------------------------
              !Compute second order Froude-krylov force
              Phi_I=CALC_SECONDORDER_INCOMING_POTENTIAL                     &
                   (Env,w1,w2,k1,k2,beta(Ibeta1),beta(Ibeta2),XM)
              Pressure_I(1)=II*delw*Phi_I(1) !-dtPhi_M
              Pressure_I(2)=II*sumw*Phi_I(2) !-dtPhi_P

              QTFHasbo(Iinteg,1)=QTFHasbo(Iinteg,1)                         &
                 -rho*Pressure_I(1)*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
              QTFHasbo(Iinteg,2)=QTFHasbo(Iinteg,2)                         &
                 -rho*Pressure_I(2)*genNormal_dS(Iinteg,Ipanelinit(2)+Ipanel)
              !------------------------------------------------------
              !Compute second order diffraction force: body force contrib
              !-----------------------------------------------------
              !Compute iw*rho*Integ(dnPhiI*Psi)dS
              GradPhi_I=CALC_SECONDORDER_INCOMING_VELOCITY                  &
                       (k1,k2,beta(Ibeta1),beta(Ibeta2),depth,XM(3),PHI_I) 
              dnPhi_I(1)=DOT_PRODUCT_COMPLEX(                               &
                      GradPhi_I(:,1),CMPLX(Normal_Vect,0),3) !DnPhi^M
              dnPhi_I(2)=DOT_PRODUCT_COMPLEX(                               &
                      GradPhi_I(:,2),CMPLX(Normal_Vect,0),3) !DnPhi^P 
              !
              QTFHasbo(Iinteg,1)=QTFHasbo(Iinteg,1)                         &
                    +II*delw*rho*dnPhi_I(1)*Radpot(1)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2)=QTFHasbo(Iinteg,2)                         &
                    +II*sumw*rho*dnPhi_I(2)*Radpot(2)*Mesh%A(Ipanel)
              !-------------------------------------------------------
              ! term: product(dtBdisplace-grad_Phi,R(n0))*Psi
              ITheta0= INT((Iinteg-1)/6)*3
              RIw1_n0=CROSS_PRODUCT_COMPLEX(                                &
                         RotAngles_Iw1(Itheta0+1:Itheta0+3),                &
                         CMPLX(Normal_Vect,0.))
              RIw2_n0=CROSS_PRODUCT_COMPLEX(                                &
                        RotAngles_Iw2(Itheta0+1:Itheta0+3),                 &
                        CMPLX(Normal_Vect,0.))
              dtDisp_GradPhi_I_Iw1=dtBdisplace_Iw1(Ipanelinit(1)+Ipanel,:)  &
                                  -TotVel_Iw1(Ipanel,:)
              dtDisp_GradPhi_I_Iw2=dtBdisplace_Iw2(Ipanelinit(1)+Ipanel,:)  &
                                  -TotVel_Iw2(Ipanel,:)

              QB1_M=DOT_PRODUCT_DIFF_BIHARM(                             &
                   dtDisp_GradPhi_I_Iw1,dtDisp_GradPhi_I_Iw2,            &
                   RIw1_n0,RIw2_n0,3)
              QB1_P=DOT_PRODUCT_SUM_BIHARM(                              &
                   dtDisp_GradPhi_I_Iw1,dtDisp_GradPhi_I_Iw2,            &
                   RIw1_n0,RIw2_n0,3)
!
              QTFHasbo(Iinteg,1)=QTFHasbo(Iinteg,1)                      &
                    -II*delw*rho*QB1_M*Radpot(1)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2)=QTFHasbo(Iinteg,2)                      &
                    -II*sumw*rho*QB1_P*Radpot(2)*Mesh%A(Ipanel)
              !-------------------------------------------------------
              ! term: product(Grad_PHI,psi).product(Bdisplace,n0) 
              QB2_M=DOT_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(           &
                     TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3),RadVel(:,1),3) &
                    *DOT_PRODUCT_COMPLEX(                                 &
                    Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),                &
                    CMPLX(Normal_vect),3)                                 &
                  +DOT_PRODUCT_COMPLEX(                                   &
                       TotVel_Iw1(Ipanelinit(1)+Ipanel,:),RadVel(:,1),3)  &
                    *DOT_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(          &
                    Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),3),             &
                    CMPLX(Normal_Vect,0.),3)
              QB2_P=DOT_PRODUCT_COMPLEX(                                  &
                     TotVel_Iw2(Ipanelinit(1)+Ipanel,:),RadVel(:,2),3)    &
                    *DOT_PRODUCT_COMPLEX(                                 &
                    Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),                &
                    CMPLX(Normal_vect),3)                                 &
                  +DOT_PRODUCT_COMPLEX(                                   &
                       TotVel_Iw1(Ipanelinit(1)+Ipanel,:),RadVel(:,2),3)  &
                    *DOT_PRODUCT_COMPLEX(                                 &
                    Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),                &
                    CMPLX(Normal_Vect,0.),3) 

              ! term: product(Grad_Psi,Bdisplace).product(gradPhi,n0)
              QB2_M=QB2_M-DOT_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(     &
                     TotVel_Iw2(Ipanelinit(1)+Ipanel,:),3),               &
                     CMPLX(Normal_vect),3)                                &
                     *DOT_PRODUCT_COMPLEX(RadVel(:,1),                    &
                     Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),3)             &
                   -DOT_PRODUCT_COMPLEX(                                  &
                       TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                &
                       CMPLX(Normal_vect),3)                              &
                    *DOT_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(          &
                    Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),3),             &
                    RadVel(:,1),3)
              QB2_P=QB2_P-DOT_PRODUCT_COMPLEX(                            &
                     TotVel_Iw2(Ipanelinit(1)+Ipanel,:),                  &
                     CMPLX(Normal_vect),3)                                &
                     *DOT_PRODUCT_COMPLEX(RadVel(:,2),                    &
                     Bdisplace_Iw1(Ipanelinit(1)+Ipanel,:),3)             &
                   -DOT_PRODUCT_COMPLEX(                                  &
                       TotVel_Iw1(Ipanelinit(1)+Ipanel,:),                &
                       CMPLX(Normal_vect),3)                              &
                    *DOT_PRODUCT_COMPLEX(                                 &
                    Bdisplace_Iw2(Ipanelinit(1)+Ipanel,:),                &
                    RadVel(:,2),3)

              QTFHasbo(Iinteg,1)=QTFHasbo(Iinteg,1)                      &
                    +II*delw*rho*(QB2_M/2)*Mesh%A(Ipanel)
              QTFHasbo(Iinteg,2)=QTFHasbo(Iinteg,2)                      &
                    +II*sumw*rho*(QB2_P/2)*Mesh%A(Ipanel)

              !------------------------------------------------------
              !IF (Iinteg==1) THEN
              !  print*,Ipanel,-rho*Pressure_I(1)*genNormal_dS(Iinteg,Ipanel)
              !ENDIF
              !IF (Iinteg==1 .AND. Isym==2) THEN
              !    print*,Ipanel,dnPhi_I(1)*Radpot(1)*Mesh%A(Ipanel)
              !ENDIF
              !IF (Iinteg==1) THEN
              !    print*,Ipanel,QB2_M/2,QB2_P/2
              !ENDIF
              !  IF (Iinteg==1 .AND. Isym==2) THEN
              !      print*,Ipanel,QB1_M
              !      print*,Ipanel,QB1_P
              !  ENDIF
            ENDDO
            ENDIF
           ENDDO

          ! DO Iwline=1,Wline%NWLineSeg
          !   Normal_Vect=Mesh%N(1:3,Wline%IndexPanel(Iwline))
          !  DO Isym=1,1+Mesh%Isym
          !   Iwlineinit=Function_SymIwline(NpanWlin,Npanels,WLine%NWLineSeg,Isym)
          !   
          !   IF (Isym==2) THEN
          !        Normal_Vect(2)=-Normal_Vect(2)
          !   ENDIF

          !   dGAMMA_Vect(1)=CMPLX(-Normal_Vect(2),0.)
          !   dGAMMA_Vect(2)=CMPLX( Normal_Vect(1),0.)
          !   dGAMMA_Vect(3)=CMPLX( 0.,0.)
          !   dGAMMA_Vect=dGAMMA_Vect*Wline%SegLength(Iwline)

          !   !interpolating radpot and radvel at delw and sumw        
          !   Radpot(1:2)=INTERP_RADIATION_POTENTIAL                       &
          !          (Nw,w,datPotVelQ%RadPot(Iwlineinit(1)+Iwline,Iinteg,:),&
          !              InterpSwitch,delw,sumw)
          !   !-----------------------------------------------------
          !     QB2V_M=CROSS_PRODUCT_COMPLEX(Radpot(1)*                 &
          !            Bdisplace_Iw1(Iwlineinit(1)+Iwline,:)            &
          !           ,COMPLEX_CONJUGATE_VECT(                          &
          !            TotVel_Iw2(Iwlineinit(1)+Iwline,:),3))           &
          !           +CROSS_PRODUCT_COMPLEX(Radpot(1)*                 &
          !            COMPLEX_CONJUGATE_VECT(                          &
          !            Bdisplace_Iw2(Iwlineinit(1)+Iwline,:),3)         &
          !           ,TotVel_Iw1(Iwlineinit(1)+Iwline,:))
          !     QB2V_P=CROSS_PRODUCT_COMPLEX(Radpot(2)*                 &
          !            Bdisplace_Iw1(Iwlineinit(1)+Iwline,:)            &
          !           ,TotVel_Iw2(Iwlineinit(1)+Iwline,:))              &
          !           +CROSS_PRODUCT_COMPLEX(Radpot(2)*                 &
          !            Bdisplace_Iw2(Iwlineinit(1)+Iwline,:)            &
          !           ,TotVel_Iw1(Iwlineinit(1)+Iwline,:))

          !    QTFHasbo(Iinteg,1)=QTFHasbo(Iinteg,1)                    &
          !           -0.5*II*delw*rho                                  &
          !           *DOT_PRODUCT_COMPLEX(QB2V_M,dGAMMA_Vect,3)  
          !    QTFHasbo(Iinteg,2)=QTFHasbo(Iinteg,2)                    &
          !           -0.5*II*sumw*rho                                  &
          !           *DOT_PRODUCT_COMPLEX(QB2V_P,dGAMMA_Vect,3)
          !    !----------------------------------------------------
          !    !IF (Iinteg==1 .AND. Isym==2) THEN
          !    !   print*,Iwline,-0.5*II*delw*rho*DOT_PRODUCT_COMPLEX(QB2V_M,dGAMMA_Vect,3)
          !    !ENDIF
          !  ENDDO
          ! ENDDO
          ENDDO
          !TRANSFER FUNCTION = BICHROMATIC FORCES / 2?
          QTFHasbo(:,:)=QTFHasbo(:,:)/2 
  END SUBROUTINE

  FUNCTION CALC_SECONDORDER_INCOMING_POTENTIAL                         &
                  (Env,w1,w2,k1,k2,beta1,beta2,XMp) result(Phi_I)
          !input/output
          REAL,                INTENT(IN)::w1,w2,k1,k2,beta1,beta2
          REAL,DIMENSION(3),   INTENT(IN)::XMp !XM at a panel
          TYPE(TEnvironment),  INTENT(IN)::Env
          COMPLEX               ::QFI_M,QFI_P   !free surface forcing
          REAL                  ::delw,sumw,OM_M,OM_P,OM_1,OM_2
          REAL,DIMENSION(2)     ::delk_vect,sumk_vect
          REAL,DIMENSION(2)     ::k1_vect,k2_vect
          REAL,DIMENSION(2)     ::XM
          REAL                  ::abs_delk,abs_sumk
          REAL                  ::g,D          
          COMPLEX,DIMENSION(2)  ::Phi_I
          g=Env%g
          D=Env%Depth
          XM(1)=XMp(1)-Env%Xeff
          XM(2)=XMp(2)-Env%Yeff
          XM(3)=XMp(3)
          delw= w1-w2
          sumw= w1+w2
          k1_vect(1)=k1*cos(beta1)
          k1_vect(2)=k1*sin(beta1)
          k2_vect(1)=k2*cos(beta2)
          k2_vect(2)=k2*sin(beta2)
          
          delk_vect(1:2)=k1_vect-k2_vect
          sumk_vect(1:2)=k1_vect+k2_vect
          
          OM_1=FUN_DISPERSION(k1,D,g)
          OM_2=FUN_DISPERSION(k2,D,g)

          QFI_M=II*g*g*exp(II*DOT_PRODUCT(delk_vect,XM(1:2)))*          &
             ( delw/w1/w2*                                              &
             (DOT_PRODUCT(k1_vect,k2_vect)+OM_1**2/g*OM_2**2/g)         &
             +0.5*((k1**2-(OM_1**2/g)**2)/w1-(k2**2-(OM_2**2/g)**2)/w2))
          QFI_P=II*g*g*exp(II*DOT_PRODUCT(sumk_vect,XM(1:2)))*          &
             ( sumw/w1/w2*                                              &
             (DOT_PRODUCT(k1_vect,k2_vect)-OM_1**2/g*OM_2**2/g)         &
             +0.5*((k1**2-(OM_1**2/g)**2)/w1+(k2**2-(OM_2**2/g)**2)/w2))

          abs_delk=SQRT(k1**2+k2**2-2*k1*k2*cos(beta1-beta2))
          abs_sumk=SQRT(k1**2+k2**2+2*k1*k2*cos(beta1-beta2))
          OM_M=FUN_DISPERSION(abs_delk,D,g)
          OM_P=FUN_DISPERSION(abs_sumk,D,g)
          IF (delw.GT.0.) THEN
          Phi_I(1)=QFI_M*CIH(abs_delk,XM(3),D)/(-delw**2+OM_M**2)
          ELSE
          Phi_I(1)=CMPLX(0.,0.)
          ENDIF
          Phi_I(2)=QFI_P*CIH(abs_sumk,XM(3),D)/(-sumw**2+OM_P**2)
  END FUNCTION

  FUNCTION CALC_SECONDORDER_INCOMING_VELOCITY                           &
                           (k1,k2,beta1,beta2,D,ZM,PHI_I) result(GradPhi)
          !input/output
          REAL,                 INTENT(IN)      ::k1,k2,beta1,beta2
          REAL,                 INTENT(IN)      ::D,ZM
          COMPLEX,DIMENSION(2), INTENT(IN)      ::PHI_I

          REAL                                  ::abs_delk,abs_sumk
          COMPLEX,DIMENSION(3,2)                ::GradPhi
          !dxPHI 
          GradPhi(1,1)=II*(k1*cos(beta1)-k2*cos(beta2))*PHI_I(1)
          GradPhi(1,2)=II*(k1*cos(beta1)+k2*cos(beta2))*PHI_I(2)
          !dyPHI
          GradPhi(2,1)=II*(k1*sin(beta1)-k2*sin(beta2))*PHI_I(1)
          GradPhi(2,2)=II*(k1*sin(beta1)+k2*sin(beta2))*PHI_I(2)
          GradPhi(2,1)=II*sin(beta1)*PHI_I(1)
          GradPhi(2,2)=II*sin(beta1)*PHI_I(2)
          !dzPHI
          abs_delk=SQRT(k1**2+k2**2-2*k1*k2*cos(beta1-beta2))
          abs_sumk=SQRT(k1**2+k2**2+2*k1*k2*cos(beta1-beta2))
          GradPhi(3,1)=abs_delk*tanh(abs_delk*(D+ZM))*PHI_I(1)
          GradPhi(3,2)=abs_sumk*tanh(abs_sumk*(D+ZM))*PHI_I(2)
  END FUNCTION
  
  FUNCTION INTERP_RADIATION_POTENTIAL                                     &
                (Nw,w,RadPot,InterpSwitch,delw,sumw) RESULT(RadPotInterp)

          INTEGER,                   INTENT(IN)::Nw,InterpSwitch
          REAL, DIMENSION(Nw),       INTENT(IN)::w
          REAL,                      INTENT(IN)::delw,sumw
          COMPLEX,DIMENSION(Nw),     INTENT(IN)::RadPot
          COMPLEX,DIMENSION(2)                 ::RadPotInterp
          REAL                                 ::RadPotR,RadPotI
          Type(linear_interp_1d)     :: interpR,interpI
          INTEGER                    :: iflag,Iprint
          IF (InterpSwitch==0) THEN
           IF (delw.GT.0.) THEN
                RadPotInterp(1)=RadPot(Fun_closest(Nw,w,delw))
           ELSE
                RadPotInterp(1)=CMPLX(0.,0.)
           ENDIF
           IF (sumw.LE.w(Nw)) THEN
                RadPotInterp(2)=RadPot(Fun_closest(Nw,w,sumw))
           ELSE
                RadPotInterp(2)=CMPLX(0.,0.)
           ENDIF
          ELSE
              CALL interpR%initialize(w ,REAL(RadPot(:)),iflag)
              CALL interpI%initialize(w ,AIMAG(RadPot(:)),iflag)
           IF (delw.GT.0.) THEN
              CALL interpR%evaluate(delw, RadPotR)
              CALL interpI%evaluate(delw, RadPotI)
              RadPotInterp(1)=CMPLX(RadPotR,RadPotI) !for delta omega
           ELSE
              RadPotInterp(1)=CMPLX(0.,0.)
           ENDIF
           IF (sumw.LE.w(Nw)) THEN
              CALL interpR%evaluate(sumw, RadPotR)
              CALL interpI%evaluate(sumw, RadPotI)
              RadPotInterp(2)=CMPLX(RadPotR,RadPotI) !for sum omega
              ELSE
              RadPotInterp(2)=CMPLX(0.,0.)
           ENDIF
           CALL interpR%destroy()
           CALL interpI%destroy()
          ENDIF
  END FUNCTION

  FUNCTION INTERP_RADIATION_VELOCITY                                    &
                (Nw,w,RadVel,InterpSwitch,delw,sumw) RESULT(RadVelInterp)
          INTEGER,                   INTENT(IN)::Nw,InterpSwitch
          REAL, DIMENSION(Nw),       INTENT(IN)::w
          REAL,                      INTENT(IN)::delw,sumw
          COMPLEX,DIMENSION(3,Nw),   INTENT(IN)::RadVel
          COMPLEX,DIMENSION(3,2)               ::RadVelInterp
          REAL                                 ::RadVelR,RadVelI
          Type(linear_interp_1d)     :: interpR,interpI
          INTEGER                    :: I,iflag

          IF (InterpSwitch==0) THEN
             IF (delw.GT.0.) THEN
                RadVelInterp(:,1)=RadVel(:,Fun_closest(Nw,w,delw))
             ELSE
                RadVelInterp(:,1)=CMPLX(0.,0.)
             ENDIF

             IF (sumw.LE.w(Nw)) THEN
                RadVelInterp(:,2)=RadVel(:,Fun_closest(Nw,w,sumw))
             ELSE
                RadVelInterp(:,2)=CMPLX(0.,0.)
             ENDIF
          ELSE
           DO I=1,3      !for Vx,Vy,Vz
             CALL interpR%initialize(w ,REAL(RadVel(I,:)),iflag)
             CALL interpI%initialize(w ,AIMAG(RadVel(I,:)),iflag)
             IF (delw.GT.0.) THEN
              CALL interpR%evaluate(delw, RadVelR)
              CALL interpI%evaluate(delw, RadVelI)
              RadVelInterp(I,1)=CMPLX(RadVelR,RadVelI) !for delta omega
             ELSE
              RadVelInterp(I,1)=CMPLX(0.,0.) !for delta omega
             END IF
             IF (sumw.LE.w(Nw)) THEN
              CALL interpR%evaluate(sumw, RadVelR)
              CALL interpI%evaluate(sumw, RadVelI)
              RadVelInterp(I,2)=CMPLX(RadVelR,RadVelI) !for sum omega
             ELSE
                RadVelInterp(I,2)=CMPLX(0.,0.)
             ENDIF
             CALL interpR%destroy()
             CALL interpI%destroy()
           ENDDO
          ENDIF
  END FUNCTION

  FUNCTION COMPLEX_CONJUGATE_VECT(var,Nvect) RESULT(Conj)
           INTEGER,                  INTENT(IN)::NVect 
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var
           COMPLEX,DIMENSION(Nvect) :: Conj
           Conj=CMPLX(REAL(var),-AIMAG(var))
  END FUNCTION

  FUNCTION DOT_PRODUCT_DIFF_BIHARM(var11,var12,var21,var22,Nvect) RESULT(prod)
           !Dot product for the difference frequency term of biharmonic function
           INTEGER,                  INTENT(IN):: NVect 
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var21,var22
           INTEGER                             :: I
           COMPLEX                             :: prod
           prod=0.5*(DOT_PRODUCT_COMPLEX(      &
                    COMPLEX_CONJUGATE_VECT(var22,Nvect),var11,Nvect)  &
                   +DOT_PRODUCT_COMPLEX(        &
                   COMPLEX_CONJUGATE_VECT(var12,Nvect),var21,Nvect))
  END FUNCTION
  
  FUNCTION DOT_PRODUCT_SUM_BIHARM(var11,var12,var21,var22,Nvect) RESULT(prod)
           !Dot product for the sum frequency term of biharmonic function
           INTEGER,                  INTENT(IN):: NVect 
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(Nvect), INTENT(IN):: var21,var22
           COMPLEX                             :: prod
           prod=0.5*(DOT_PRODUCT_COMPLEX(var22,var11,Nvect)            &
                        +DOT_PRODUCT_COMPLEX(var12,var21,Nvect))
  END FUNCTION

  FUNCTION CROSS_PRODUCT_DIFF_BIHARM(var11,var12,var21,var22) RESULT(prod)
           !CROSS product for the difference frequency term of biharmonic function
           COMPLEX,DIMENSION(3), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(3), INTENT(IN):: var21,var22
           COMPLEX,DIMENSION(3)            :: prod
           prod=0.5*(CROSS_PRODUCT_COMPLEX(var11,COMPLEX_CONJUGATE_VECT(var22,3)) &
                     +CROSS_PRODUCT_COMPLEX(COMPLEX_CONJUGATE_VECT(var12,3),var21))
  END FUNCTION
  
  FUNCTION CROSS_PRODUCT_SUM_BIHARM(var11,var12,var21,var22) RESULT(prod)
           !CROSS product for the sum frequency term of biharmonic function
           COMPLEX,DIMENSION(3), INTENT(IN):: var11,var12
           COMPLEX,DIMENSION(3), INTENT(IN):: var21,var22
           COMPLEX,DIMENSION(3)            :: prod
           prod=0.5*(CROSS_PRODUCT_COMPLEX(var11,var22)                            &
                        +CROSS_PRODUCT_COMPLEX(var12,var21))
  END FUNCTION

  FUNCTION Function_SymIpanel(NpanWlin,Npanels,Isym) RESULT(Ipanelinit)
        INTEGER,            INTENT(IN)  ::NpanWlin,Npanels,Isym
        INTEGER,DIMENSION(2)            ::Ipanelinit
           IF (Isym==1) Ipanelinit(1:2)=0
           IF (Isym==2) THEN 
                   Ipanelinit(1)=NpanWLin
                   Ipanelinit(2)=Npanels
           ENDIF
  END FUNCTION

  FUNCTION Function_SymIwline(NpanWlin,Npanels,Nwline,Isym) RESULT(Iwlineinit)
        INTEGER,            INTENT(IN)  ::NpanWlin,Npanels,Nwline,Isym
        INTEGER,DIMENSION(2)            ::Iwlineinit
           IF (Isym==1) THEN
                   Iwlineinit(1)=Npanels
                   Iwlineinit(2)=0
           ELSEIF (Isym==2) THEN 
                   Iwlineinit(1)=Npanels+NpanWlin
                   Iwlineinit(2)=Nwline
           ENDIF
  END FUNCTION

  FUNCTION DOT_PRODUCT_COMPLEX(var1,var2,Nvect) RESULT(prod)
       INTEGER, INTENT(IN)                  :: Nvect
       COMPLEX,DIMENSION(Nvect), INTENT(IN) :: var1,var2  
       COMPLEX                              ::prod
       !NOTE: DOT_PRODUCT(A,B)=A*.B for A,B complex variables
       !we want to have DOT_PRODUCT(A,B)=A.B
       prod=DOT_PRODUCT(COMPLEX_CONJUGATE_VECT(var1,Nvect),var2)
  END FUNCTION

   
END MODULE

