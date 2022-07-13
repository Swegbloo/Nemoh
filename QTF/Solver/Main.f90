!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991             
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!     THESE DE CHEN XIAO-BO(1988) 
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)  Version 2014 
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!   SOLVER
!--------------------------------------------------------------------------------------
PROGRAM MAIN
!
USE MIdentification
USE MNemohCal,          ONLY:TNemCal,READ_TNEMOHCAL
USE MMesh
USE MFace,              ONLY:TVFace, Prepare_FaceMesh,TWLine,Prepare_Waterline
USE MReadInputFiles,    ONLY:Read_NP_GaussQuad,Read_Mechanical_Coefs,TMech,  &
                             Read_FirstOrderLoad,TLoad1,Read_Motion,TSource, &
                             READ_POTENTIALS_VELOCITIES_BODYWLINE,           &
                             READ_GENERALIZED_NORMAL_BODY_dAREA
USE MEnvironment,       ONLY: TEnvironment,FunVect_inverseDispersion
USE MLogFile
USE MQSolverPreparation !CONTAINS:TQfreq,TpotVel,PREPARE_POTENTIAL_VELOCITIES
                        !PREPARE_BODY_DISPLACEMENT,DISCRETIZED_OMEGA_WAVENUMBER_FOR_QTF
                        !PREPARE_INERTIA_FORCE
                        !CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT
                        !WRITE_QTFSOLVERLOGFILE
USE MQSolver            !CONTAINS:COMPUTATION_QTF_QUADRATIC
USE MQSolverOutputFiles !CONTAINS: OutFileDM,OutFileDP,...
                        !INITIALIZE_OUTPUT_FILE,WRITE_QTF_DATA

IMPLICIT NONE
!
! ------Declaration variables
        TYPE(TID)                       :: ID
        TYPE(TMesh)                     :: Mesh
        TYPE(TNemCal)                   :: inpNEMOHCAL
        TYPE(TEnvironment)              :: Env
        INTEGER                         :: Nw,Nbeta           ! Number of Freq,direction
        REAL, ALLOCATABLE,DIMENSION(:)  :: w,kw,beta          ! vector of freq [rad/s],
                                                              ! wave numbers [rad/m],
                                                              ! direction angle [rad]
        TYPE(TWLine)  :: WLine          ! Waterline 
        TYPE(TVFace)  :: VFace          ! Face of a panel  
        TYPE(TMech)   :: MechCoef       ! Mechanical Coef 
        TYPE(TLoad1)  :: Forces1        ! First order forces
        INTEGER       :: NP_GQ          ! Number of point for Gauss Quad. Integration
        INTEGER       :: Nintegration   ! Number of excitation force integration 
        INTEGER       :: Nradiation     ! Number of radiation problem
        INTEGER       :: Nbodies        ! Number of bodies
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: Motion       ! Complex RAO
        TYPE(TPotVel)                           :: datPotVel    ! data Pot & vel
        TYPE(TPotVel)                           :: datPotVelQ    ! data Pot & vel for QTF
        INTEGER                                 :: I,Ibeta1,Ibeta2,Iinteg,Ipanel
        INTEGER                                 :: Iw1,Iw2,IwQ 
        INTEGER                                 :: NPFlow        ! Number of flow points
        REAL,ALLOCATABLE,DIMENSION(:,:)         :: genNormal_dS  ! generalized Normal on panel
                                                           ! time the area of the panel
        REAL,ALLOCATABLE,DIMENSION(:,:)         :: genNormalWLine_dGamma! generalized Normal 
                                                   !on wLine segment time the segm. length
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:)  :: BDisplaceQ!Body displacement
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: InertiaForceQ
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)    :: RotAnglesQ
        !freq related variables For QTF computation
        TYPE(TQfreq)                            :: Qfreq
        REAL                                    :: winputQ(3), BForwardSpeed,delwiter
        INTEGER                                 :: NwQ    ! Number of wave freq
        !---------
        COMPLEX, ALLOCATABLE,DIMENSION(:,:)   :: QTF_DUOK ,QTF_HASBO   
        COMPLEX, ALLOCATABLE,DIMENSION(:,:)   :: QTF_HASFS,QTF_HASFS_ASYMP
        !

!
!   --- Initialize and read input datas -----------------------------------------------
!
        CALL ReadTID(ID)
        CALL ReadTMesh(Mesh,TRIM(ID%ID)//'/mesh/')  
        CALL READ_TNEMOHCAL(TRIM(ID%ID),InpNEMOHCAL)
        !        
        Env          =InpNEMOHCAL%Env
        Nw           =InpNEMOHCAL%waveinput%NFreq
        Nbeta        =InpNEMOHCAL%waveinput%NBeta
        Nbodies      =InpNEMOHCAL%Nbodies
        Nintegration =InpNEMOHCAL%Nintegtot
        Nradiation   =InpNEMOHCAL%Nradtot
        winputQ      =InpNEMOHCAL%qtfinput%omega
        NwQ          =winputQ(1)
        BForwardSpeed=InpNEMOHCAL%qtfinput%body_forward_speed
        NP_GQ        =Read_NP_GaussQuad(TRIM(ID%ID)) 
        !
        CALL Prepare_FaceMesh(Mesh,NP_GQ,VFace)
        CALL Prepare_Waterline(VFace,Mesh%xy_diameter,Mesh%Npanels,WLine)
        !
        NPFlow   =(Mesh%Npanels+WLine%NWLineSeg)*2**Mesh%Isym !Number of flow point
        !
        !Dynamic Memory allocation
        ALLOCATE(Motion(Nw,Nradiation,Nbeta))
        ALLOCATE(datPotVel%TotPot(NPFlow,Nbeta,Nw))
        ALLOCATE(datPotVel%TotVel(NPFlow,3,Nbeta,Nw))
        ALLOCATE(datPotVel%RadPot(NPFlow,Nradiation,Nw))
        ALLOCATE(datPotVel%RadVel(NPFlow,3,Nradiation,Nw))
        
        ALLOCATE(w(Nw),kw(Nw),beta(Nbeta))
        ALLOCATE(genNormal_dS(Nintegration,Mesh%Npanels*2**Mesh%Isym))
        ALLOCATE(genNormalWLine_dGamma(Nintegration,Wline%NWLineseg*2**Mesh%Isym))
        ALLOCATE(InertiaForceQ(NwQ,Nbeta,Nintegration))
        ALLOCATE(RotAnglesQ(NwQ,Nbeta,3*Nbodies))
        ALLOCATE(QTF_DUOK(Nintegration,2))!2 is for QTF- and QTF+
        ALLOCATE(QTF_HASBO(Nintegration,2))
        !
        !
        CALL Read_Mechanical_Coefs(TRIM(ID%ID),Nradiation,MechCoef)
        CALL Read_FirstOrderLoad(TRIM(ID%ID),Nw,Nbeta,Nintegration,Nradiation,Forces1)
        CALL Read_Motion(TRIM(ID%ID),Nw,Nbeta,Nradiation,Motion)!Complex RAO
        CALL READ_POTENTIALS_VELOCITIES_BODYWLINE(TRIM(ID%ID),Nw,Nbeta,NRadiation,       &
                NPFlow,datPotVel%TotPot,datPotVel%TotVel,                                &
                datPotVel%RadPot,datPotVel%RadVel,w,kw,beta)

        CALL READ_GENERALIZED_NORMAL_BODY_dAREA(TRIM(ID%ID),Mesh%Npanels*2**Mesh%Isym,   &
                                                          Nintegration,genNormal_dS)
        CALL CALC_GENERALIZED_NORMAL_WATERLINE_dSEGMENT(Mesh,Nintegration,WLine,         &
                                                InpNEMOHCAL, genNormalWLine_dGamma)
        CALL Discretized_omega_wavenumber_for_QTF(Nw,w,kw,NwQ,winputQ(2:3),Nbeta,beta,   &
                                                  BForwardSpeed,Env%depth,Env%g,Qfreq) 
        
        CALL WRITE_QTFSOLVERLOGFILE(TRIM(ID%ID),Nbeta,beta,Qfreq)

        CALL PREPARE_POTENTIAL_VELOCITIES(Qfreq,Nw,w,Nbeta,beta,NPFlow,                  &
                Nradiation,datPotVel,datPotVelQ)
        
        CALL PREPARE_BODY_DISPLACEMENT(Qfreq,Nw,w,Nbeta,Nradiation,NPFlow,Nbodies,       &
                                        Mesh,WLine,InpNEMOHCAL,Motion,BdisplaceQ)

        CALL PREPARE_INERTIA_FORCES(MechCoef,Motion,Nw,Nbeta,Nradiation,Nintegration,    &
                                        w,Qfreq,Forces1,InertiaForceQ)
        CALL PREPARE_ROTATION_ANGLES(Motion,Nw,Nbeta,Nradiation,Nbodies,                 &
                                        w,Qfreq,RotAnglesQ)
        !DO Ipanel=1,Mesh%Npanels+WLine%NWLineSeg
        !print*,Ipanel,datPotVelQ%TotPot(Ipanel,1,1)
        !ENDDO
        ! COMPUTE QTF
        CALL  INITIALIZE_OUTPUT_FILES(TRIM(ID%ID))
        DO Ibeta1=1,Nbeta
           DO Ibeta2=1,Nbeta
                WRITE(*,'(A,F7.3,A,F7.3,A)'),'beta1=', beta(Ibeta1)*180/PI,&
                        ', beta2=', beta(Ibeta2)*180/PI, ' [deg]'
                DO IwQ=0,NwQ-1
                    DO Iw1=IwQ+1,NwQ
                        Iw2=Iw1-IwQ
                        IF (Iw1==1 .AND. Iw2==1) THEN
                                delwiter=Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2)
                                WRITE(*,'(A,F7.3,A,F7.3,A)')'w1-w2=',delwiter, ', w1+w2=', &
                                     Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1), ' [rad/s]'
                        ENDIF
                        
                        IF (Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2).GT.1.01*delwiter) THEN
                          delwiter=Qfreq%wQ(Iw1,Ibeta1)-Qfreq%wQ(Iw2,Ibeta2)
                          IF (Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1).LE.w(Nw)) THEN 
                             WRITE(*,'(A,F7.3,A,F7.3,A)'),'w1-w2=',delwiter, ', w1+w2=', &
                                     Qfreq%wQ(Iw2,Ibeta2)+Qfreq%wQ(Iw1,Ibeta1), ' [rad/s]'
                          ELSE
                            WRITE(*,'(A,F7.3,A)'),'w1-w2=',delwiter, ', w1+w2= --NA-- [rad/s]'
                          ENDIF
                        ENDIF

                        CALL COMPUTATION_QTF_QUADRATIC(Iw1,Iw2,Ibeta1,Ibeta2,Nintegration,&
                                Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,Env%rho,Env%g,  &
                                datPotVelQ,BdisplaceQ,genNormal_dS,genNormalWLine_dGamma, &
                                Qfreq%wQ,beta,InertiaForceQ,RotAnglesQ,QTF_DUOK(:,:))

                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileDM,OutFileDP,Nintegration, &
                                Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                beta(Ibeta1),beta(Ibeta2), QTF_DUOK(:,:))

                        CALL COMPUTATION_QTF_POTENTIAL_BODYFORCE(Iw1,Iw2,Ibeta1,Ibeta2,   &
                                Nintegration,Mesh,VFace,WLine,NwQ,Nbeta,NPFlow,Nbodies,   &
                                Env, datPotVelQ,BdisplaceQ,genNormal_dS,                  &
                                Nw, w,Qfreq,beta,RotAnglesQ,                              &
                                QTF_HASBO(:,:))      
                        
                        CALL WRITE_QTF_DATA(TRIM(ID%ID),OutFileHBM,OutFileHBP,Nintegration,&
                                Qfreq%wQ(Iw1,Ibeta1),Qfreq%wQ(Iw2,Ibeta2),                &
                                beta(Ibeta1),beta(Ibeta2), QTF_HASBO(:,:))
                    ENDDO
                ENDDO
           ENDDO
        ENDDO
       
       ! print*,RadVel(900:1000,1,1,1)

! ----- Finalize ---------------------------------------------------------------------------
!       DEALOCATING variables
        DEALLOCATE(w,kw,beta)
        DO I=1,InpNEMOHCAL%Nbodies
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%RadCase)
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%IntCase)
        ENDDO
        DEALLOCATE(inpNEMOHCAL%bodyinput)
        DEALLOCATE(Motion)
        DEALLOCATE(VFace%X,VFace%XM,VFace%N,VFace%A,VFace%tDis)
        DEALLOCATE(VFace%dXdXG_WGQ_per_A,VFace%XM_GQ)
        DEALLOCATE(datPotVelQ%TotPot,datPotVelQ%TotVel) 
        DEALLOCATE(datPotVelQ%RadPot,datPotVelQ%RadVel)
        DEALLOCATE(Qfreq%wQ,Qfreq%kQ)
        DEALLOCATE(Qfreq%diffwQ,Qfreq%sumwQ)
        DEALLOCATE(Qfreq%InterpPotSwitch)
        DEALLOCATE(genNormal_dS,genNormalWLine_dGamma)
        DEALLOCATE(BDisplaceQ)
        DEALLOCATE(InertiaForceQ)
        DEALLOCATE(QTF_DUOK,QTF_HASBO)
       ! DEALLOCATE(QTF_HASFS,QTF_HASFS_ASYMP)

END PROGRAM
