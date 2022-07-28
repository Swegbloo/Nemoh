!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991             
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr) 2014 
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
Module MQpreprocessor
USE CONSTANTS
USE MReadInputFiles,            ONLY: TSource
USE MMesh,                      ONLY: TMesh
USE MFace,                      ONLY: TVFace, TWLine
USE MEnvironment,               ONLY: TEnvironment,Fun_inverseDispersion, &
                                      COMPUTE_INC_POTENTIAL_VELOCITY
USE Elementary_functions,       ONLY: X0 !invers of disp. relation
! Green functions
USE M_INITIALIZE_GREEN,         ONLY: TGREEN
USE MInfluenceMatrix,           ONLY: CONSTRUCT_INFLUENCE_MATRIX
!
USE MFileDirectoryList          
USE MLogFile 
!
IMPLICIT NONE


CONTAINS
  SUBROUTINE COMPUTE_POTENTIALS_AND_VELOCITIES(wd,Iw,omega,Vbeta,Nbeta,Nradiation,&
                                              Env,Mesh,VFace,WLine,IGreen,   &
                                              SourceDistr,MotionIw)
        !INPUT/OUTPUT
        CHARACTER(LEN=*),               INTENT(IN)::wd
        INTEGER,                        INTENT(IN)::Iw,Nbeta,Nradiation
        REAL,                           INTENT(IN)::omega       !rad freq w(Iw)
        REAL,DIMENSION(Nbeta),          INTENT(IN)::Vbeta       !Angle vector
        TYPE(TEnvironment),             INTENT(IN)::Env
        TYPE(TMesh),                    INTENT(IN)::Mesh
        TYPE(TVFace),                   INTENT(IN)::VFace
        TYPE(TWLine),                   INTENT(IN)::WLine
        TYPE(TGREEN),                   INTENT(INOUT)::IGreen
        TYPE(TSource),                  INTENT(IN)::SourceDistr 
        COMPLEX,DIMENSION(Nradiation,Nbeta),INTENT(IN):: MotionIw 
        !LOCAL
        INTEGER                                ::NPFLOW,uFile,ILINE 
        REAL                                   :: wavenumber
        COMPLEX, DIMENSION(:,:,:)  , ALLOCATABLE :: S       ! Inf. coef. Integ. of Green func.
        COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE :: GradS   ! Inf. coef. Integ. of Gradient Green func.
        INTEGER                                :: Ibeta,Ipanel,Irad
        COMPLEX,DIMENSION(Mesh%Npanels,Nbeta)  :: ZPGB,ZPGS ! Perturbation sources, B for Body
                                                            ! S for symmetric part
        COMPLEX,ALLOCATABLE,DIMENSION(:)       ::Potential  ! Total Potential=Phi_InC+Phi_Perturb
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)     ::Velocity   ! Total velocity
        COMPLEX,ALLOCATABLE,DIMENSION(:)       ::RadPotential ! Rad pot mod I
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)     ::RadVelocity  ! Rad velocity mod I
        CHARACTER(LEN=1000)                    :: LogTextToBeWritten

        NPFLOW=Mesh%NPanels+WLine%NWlineseg
        ALLOCATE(S(NPFLOW,Mesh%NPanels,2**Mesh%Isym))
        ALLOCATE(GradS(NPFLOW,Mesh%NPanels,3,2**Mesh%Isym))
        ALLOCATE(Potential(NPFLOW*2**Mesh%Isym))      
        ALLOCATE(Velocity(NPFLOW*2**Mesh%Isym,3)) 
        wavenumber=Fun_inverseDispersion(omega,Env%depth,Env%g)  

        CALL CONSTRUCT_INFLUENCE_MATRIX(omega,wavenumber,Env,IGreen,                  &
                        Mesh,VFace,WLine%XM,WLine%NWlineseg,S,GradS)
        
        DO Ibeta=1,Nbeta
           !CONSTRUCT PERTURBATION(Diff+Rad) SINGULAR DISTRIBUTION FOR EACH WAVE DIRECTION
           !assign the diffraction singular source distribution for all panels
           ZPGB(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGB(1:Mesh%Npanels,Nradiation+Ibeta)   
           ZPGS(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGB(1:Mesh%Npanels,Nradiation+Ibeta) 
              !sum the diffraction + radiation singular source distribution for all panels
              DO Irad=1,Nradiation  
              ZPGB(1:Mesh%Npanels,Ibeta)=ZPGB(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGB(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ZPGS(1:Mesh%Npanels,Ibeta)=ZPGS(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGS(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ENDDO
              !--------------------------------------------------------------------------------   
              ! Compute incoming Potential and velocity
              CALL COMPUTE_INC_POTENTIAL_VELOCITY(wavenumber,omega,Vbeta(Ibeta),                &
                                  VFace%XM,Mesh%Npanels,WLine%XM,WLine%NWlineseg,               &
                                  Env,Mesh%ISym,Potential,Velocity)
              !--------------------------------------------------------------------------------
              ! Compute total potential and velocity
              IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                Potential(1:NPFLOW)=Potential(1:NPFLOW)+                                        &
                        MATMUL(S(1:NPFLOW,1:Mesh%Npanels,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vx
                Velocity(1:NPFLOW,1)=Velocity(1:NPFLOW,1)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,1,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vy
                Velocity(1:NPFLOW,2)=Velocity(1:NPFLOW,2)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,2,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vz
                Velocity(1:NPFLOW,3)=Velocity(1:NPFLOW,3)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,3,1),ZPGB(1:Mesh%Npanels,Ibeta))

              ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                 Potential(1:NPFLOW)=Potential(1:NPFLOW)+                                       &
                   MATMUL(S(:, :, 1) + S(:, :, 2), ZPGB(:,Ibeta))/2                             &
                   +MATMUL(S(:, :, 1) - S(:, :, 2), ZPGS(:,Ibeta))/2
                 Potential(NPFLOW+1:2*NPFLOW) =Potential(NPFLOW+1:2*NPFLOW)+                    &
                   MATMUL(S(:, :, 1) - S(:, :, 2), ZPGB(:,Ibeta))/2                             &
                   +MATMUL(S(:, :, 1) + S(:, :, 2), ZPGS(:,Ibeta))/2
                !Vx 
                Velocity(1:NPFLOW,1)=Velocity(1:NPFLOW,1)+                                      &
                   MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,1) =Velocity(NPFLOW+1:2*NPFLOW,1)+                  &
                   MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), ZPGS(:,Ibeta))/2
                !Vy 
                Velocity(1:NPFLOW,2)=Velocity(1:NPFLOW,2)+                                      &
                   MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,2) =Velocity(NPFLOW+1:2*NPFLOW,2)-                 &
                   MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), ZPGB(:,Ibeta))/2               &
                   -MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), ZPGS(:,Ibeta))/2
                !Vz 
                Velocity(1:NPFLOW,3)=Velocity(1:NPFLOW,3)+                                      &
                   MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,3) =Velocity(NPFLOW+1:2*NPFLOW,3)+                  &
                   MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), ZPGS(:,Ibeta))/2
               ! DO Ipanel=1,NPFLOW
               !         print*,Ipanel,Potential(Ipanel)
               ! ENDDO
               ! STOP
              ENDIF
              
              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotPotFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+2*NPFLOW*2**Mesh%Isym)
              ILINE=(Iw-1)*Nbeta+Ibeta
              WRITE(uFile,REC=ILINE) omega,wavenumber,Vbeta(Ibeta),                             &
                             (REAL(Potential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym),            &
                             (AIMAG(Potential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)

              
              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotVelFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+3*2*NPFLOW*2**Mesh%Isym)
              WRITE(uFile,REC=ILINE) omega,wavenumber,Vbeta(Ibeta),                             &
                               ( REAL(Velocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               ( REAL(Velocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               ( REAL(Velocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)
        ENDDO
        DEALLOCATE(Potential,Velocity) 
        
        ! Compute radiation potential and velocity
        ALLOCATE(RadPotential(NPFLOW*2**Mesh%Isym))      
        ALLOCATE(RadVelocity(NPFLOW*2**Mesh%Isym,3)) 
        DO IRad=1,Nradiation
               IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                RadPotential(1:NPFLOW)=                                                         &
                  MATMUL(S(1:NPFLOW,1:Mesh%Npanels,1),SourceDistr%ZIGB(:,Irad))
                !Vx
                RadVelocity(1:NPFLOW,1)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,1,1),SourceDistr%ZIGB(:,Irad))
                !Vy
                RadVelocity(1:NPFLOW,2)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,2,1),SourceDistr%ZIGB(:,Irad))
                !Vz
                RadVelocity(1:NPFLOW,3)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,3,1),SourceDistr%ZIGB(:,Irad))
              ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                 RadPotential(1:NPFLOW)=                                                        &
                   MATMUL(S(:, :, 1) + S(:, :, 2), SourceDistr%ZIGB(:,Irad))/2                  &
                   +MATMUL(S(:, :, 1) - S(:, :, 2),SourceDistr%ZIGS(:,Irad))/2
                 RadPotential(NPFLOW+1:2*NPFLOW) =                                              &
                   MATMUL(S(:, :, 1) - S(:, :, 2), SourceDistr%ZIGB(:,Irad))/2                  &
                   +MATMUL(S(:, :, 1) + S(:, :, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vx 
                RadVelocity(1:NPFLOW,1)=                                                        &
                   MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,1) =                                              &
                   MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vy 
                RadVelocity(1:NPFLOW,2)=                                                        &
                   MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,2) =                                              &
                   -MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), SourceDistr%ZIGB(:,Irad))/2   &
                   -MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vz 
                RadVelocity(1:NPFLOW,3)=                                                        &
                   MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,3) =                                              &
                   MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2),SourceDistr%ZIGS(:,Irad))/2
              ENDIF
              
              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadPotFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+2*NPFLOW*2**Mesh%Isym)
              ILINE=(Iw-1)*Nradiation+Irad
              WRITE(uFile,REC=ILINE) omega,wavenumber,Irad,                                     &
                             (REAL(RadPotential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym),         &
                             (AIMAG(RadPotential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)

              
              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadVelFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3+3*2*NPFLOW*2**Mesh%Isym)
              WRITE(uFile,REC=ILINE) omega,wavenumber,Irad,                                     &
                               ( REAL(RadVelocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               ( REAL(RadVelocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               ( REAL(RadVelocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)

        ENDDO
        IF ((Env%depth == INFINITE_DEPTH).OR.(omega**2*Env%depth/Env%g.GE.20)) THEN
        WRITE(LogTextToBeWritten,'(A,F7.3,A)') 'Omega= ', omega,' [rad/s], Green Fun: Infinite-Depth. DONE!'
        ELSE
        WRITE(LogTextToBeWritten,'(A,F7.3,A)') 'Omega= ', omega,' [rad/s], Green Fun: Finite-Depth. DONE!'
        ENDIF
        CALL WRITE_LOGFILE(trim(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

        DEALLOCATE(RadPotential,RadVelocity) 
        DEALLOCATE(S,GradS)        
  END SUBROUTINE

  SUBROUTINE WRITE_QTFLOGFILE(wd,beta,Nbeta,w,Nw,NP_GQ,Nbodies,depth)
        CHARACTER(LEN=*),             INTENT(IN)::wd
        INTEGER,                      INTENT(IN)::Nbeta,Nw,NP_GQ,Nbodies
        REAL, DIMENSION(Nbeta),       INTENT(IN)::beta
        REAL, DIMENSION(Nw),          INTENT(IN)::w
        REAL,                         INTENT(IN)::depth
        CHARACTER(LEN=1000)                     ::LogTextToBeWritten


        WRITE(*, *) ' '
        WRITE(LogTextToBeWritten,*) '----Pre-Processing (QTF Module)---'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I4,A,3(F7.3,A))') ' NFreq= ', Nw, ', omega = (', w(1),':',w(2)-w(1),':',w(Nw),') rad/s'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        IF (Nbeta>1) THEN
        WRITE(LogTextToBeWritten,'(A,I4,A,3(F7.3,A))') ' Nbeta= ', Nbeta, ', beta  = (',beta(1)*180/PI, &
                ':',(beta(2)-beta(1))*180/PI,':',beta(Nbeta)*180/PI,') deg'
        ELSE
         WRITE(LogTextToBeWritten,'(A,I4,A,F7.3,A)') ' Nbeta= ', Nbeta, ', beta  = ',beta(1)*180/PI,' deg'
        ENDIF
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I3,A)') ' NBodies=', Nbodies, ', NDOF=6'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        IF (depth==INFINITE_DEPTH) THEN
        WRITE(LogTextToBeWritten,*) 'Waterdepth= Infinite'
        ELSE
        WRITE(LogTextToBeWritten,'(A,F7.3,A)') ' Waterdepth= ', depth,' m'
        ENDIF
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
         WRITE(LogTextToBeWritten,'(A,I3)') ' NP Gauss Quadrature Integ.: ', NP_GQ
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

  END SUBROUTINE

 
END MODULE
