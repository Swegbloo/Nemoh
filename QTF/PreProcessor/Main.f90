!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991             
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)  Version 2014 
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------

PROGRAM Main
!
USE MIdentification
USE MNemohCal,          ONLY:TNemCal,READ_TNEMOHCAL,Discretized_Omega_and_Beta
USE MMesh
USE MFace,              ONLY:TVFace, Prepare_FaceMesh,TWLine,Prepare_Waterline
USE MReadInputFiles,    ONLY:Read_NP_GaussQuad,Read_Mechanical_Coefs,TMech,  &
                             Read_FirstOrderLoad,TLoad1,Read_Motion,TSource, &
                             Read_SourceDistribution
USE M_INITIALIZE_GREEN, ONLY: TGREEN, INITIALIZE_GREEN
USE MQpreprocessor
USE MLogFile               


IMPLICIT NONE
!
!Declaration variables
!
        TYPE(TID)                       :: ID
        TYPE(TMesh)                     :: Mesh
        TYPE(TNemCal)                   :: inpNEMOHCAL
        INTEGER                         :: Nw,Nbeta           ! Number of Freq,direction
        REAL, ALLOCATABLE,DIMENSION(:)  :: w,beta             ! vector of freq [rad/s],
                                                              ! direction angle [rad]
        TYPE(TWLine)  :: WLine          ! Waterline 
        TYPE(TVFace)  :: VFace          ! Face of a panel       
        TYPE(TMech)   :: MechCoef       ! MechCoef 
        TYPE(TLoad1)  :: Forces1        ! First order forces
        INTEGER       :: NP_GQ          ! Number of point for Gauss Quad. Integration
        INTEGER       :: Nintegration   ! Number of excitation force integration 
        INTEGER       :: Nradiation     ! Number of radiation problem
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: Motion
        TYPE(TSource) :: SOURCEDISTR                            !First order NEMOH solution
        TYPE(TGREEN)                         :: IGreen          ! Initial Green variables
        INTEGER       :: I,J,uFile
        CHARACTER(LEN=1000)                :: LogTextToBeWritten
        REAL                               :: tcpu_start
!
!   --- Initialize and read input datas -----------------------------------------------
!
        CALL ReadTID(ID)
        CALL ReadTMesh(Mesh,TRIM(ID%ID)//'/mesh/')  
        CALL READ_TNEMOHCAL(TRIM(ID%ID),InpNEMOHCAL)
        CALL Read_Mechanical_Coefs(TRIM(ID%ID),InpNEMOHCAL%Nbodies,MechCoef)
!        
        Nw          =InpNEMOHCAL%waveinput%NFreq
        Nbeta       =InpNEMOHCAL%waveinput%NBeta
        Nintegration=InpNEMOHCAL%Nintegtot
        Nradiation  =InpNEMOHCAL%Nradtot
        CALL Read_FirstOrderLoad(TRIM(ID%ID),Nw,Nbeta,Nintegration,Nradiation,Forces1)
        ALLOCATE(Motion(Nw,Nradiation,Nbeta))
        CALL Read_Motion(TRIM(ID%ID),Nw,Nbeta,Nradiation,Motion)!RAO
!        
        ALLOCATE(w(Nw),beta(Nbeta))
        CALL Discretized_Omega_and_Beta(1,InpNEMOHCAL%waveinput,Nw,Nbeta,w,beta)
!        
        NP_GQ=Read_NP_GaussQuad(TRIM(ID%ID)) 
!
        CALL Prepare_FaceMesh(Mesh,NP_GQ,VFace)
        CALL Prepare_Waterline(VFace,Mesh%xy_diameter,Mesh%Npanels,WLine)
!
        CALL INITIALIZE_GREEN(VFace,Mesh,InpNEMOHCAL%Env%depth, &
                              WLine%XM,WLine%NWlineseg,IGreen)

!
        CALL WRITE_QTFLOGFILE(TRIM(ID%ID),beta,Nbeta,w,Nw,NP_GQ,                        &
                InpNEMOHCAL%Nbodies,InpNEMOHCAL%Env%depth) 
        CALL START_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE,IdAppend)
        WRITE(LogTextToBeWritten,*) '-------'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

        ALLOCATE(SOURCEDISTR%ZIGB(Mesh%Npanels,Nradiation+Nbeta))
        ALLOCATE(SOURCEDISTR%ZIGS(Mesh%Npanels,Nradiation+Nbeta))
        
        CALL make_directory(TRIM(ID%ID)//'/'//PreprocDir)

! ------Computing potentials and velocities--------------------------------------------------        
        DO I=1,Nw
            CALL Read_SourceDistribution(TRIM(ID%ID),I,Nw,Nradiation,Nbeta,                 &
                                         Mesh%Npanels,SourceDistr)
            CALL COMPUTE_POTENTIALS_AND_VELOCITIES(TRIM(ID%ID),                             &
                                         I,w(I),beta,Nbeta,Nradiation,InpNEMOHCAL%Env,Mesh, &
                                         VFace,WLine,IGreen,SourceDistr,Motion(I,:,:))
        ENDDO
        WRITE(LogTextToBeWritten,*) '-------'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        CALL END_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE)
        WRITE(LogTextToBeWritten,*) '---- DONE ---'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
! ----- Finalize ---------------------------------------------------------------------------
!       DEALOCATING variables
        DEALLOCATE(w,beta)
        DO I=1,InpNEMOHCAL%Nbodies
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%RadCase)
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%IntCase)
        ENDDO
        DEALLOCATE(inpNEMOHCAL%bodyinput)

        DEALLOCATE(VFace%X,VFace%XM,VFace%N,VFace%A,VFace%tDis)
        DEALLOCATE(VFace%dXdXG_WGQ_per_A,VFace%XM_GQ)
        DEALLOCATE(WLine%XM,WLine%SegLength,WLine%IndexPanel)
        DEALLOCATE(MechCoef%MassMat,MechCoef%StiffMat)
        DEALLOCATE(MechCoef%StiffMat_EXT,MechCoef%DampCoefMat_EXT)
        DEALLOCATE(Forces1%addedmass,Forces1%dampcoef)
        DEALLOCATE(Forces1%excitation)
        DEALLOCATE(Motion)
        DEALLOCATE(SOURCEDISTR%ZIGB,SOURCEDISTR%ZIGS)
        DEALLOCATE(IGreen%FSP1,IGreen%FSM1,IGreen%VSP1,IGREEN%VSM1)
        DEALLOCATE(IGreen%FSP1_INF,IGreen%FSM1_INF,IGreen%VSP1_INF,IGREEN%VSM1_INF)
        DEALLOCATE(IGreen%XR,IGreen%XZ)
        DEALLOCATE(IGreen%APD1X,IGREEN%APD2X,IGREEN%APD1Z,IGREEN%APD2Z)
END PROGRAM Main
