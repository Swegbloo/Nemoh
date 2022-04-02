!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!
!   Contributors list:
!   - G. Delhommeau
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------

PROGRAM Main

  USE Constants
  USE MMesh,                ONLY: TMesh,           ReadTMesh
  USE MEnvironment,         ONLY: TEnvironment,    ReadTEnvironment
  USE MBodyConditions,      ONLY: TBodyConditions, ReadTBodyConditions
  USE M_Solver,             ONLY: TSolver,         ReadTSolver, ID_GMRES
  USE MLogFile              !ID 
  ! Preprocessing and initialization
  USE M_INITIALIZE_GREEN,   ONLY: TGREEN, INITIALIZE_GREEN
  USE Elementary_functions, ONLY: X0

  ! Resolution
  USE SOLVE_BEM_DIRECT,     ONLY: SOLVE_POTENTIAL_DIRECT
  ! Post processing and output
  USE OUTPUT,               ONLY: WRITE_DATA_ON_MESH
  USE FORCES,               ONLY: COMPUTE_AND_WRITE_FORCES
  USE KOCHIN,               ONLY: COMPUTE_AND_WRITE_KOCHIN
  USE FREESURFACE,          ONLY: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION

  IMPLICIT NONE

  CHARACTER(LEN=1000)   :: wd             ! Working directory path (max length: 1000 characters, increase if necessary)
  TYPE(TMesh)           :: Mesh           ! Mesh of the floating body
  TYPE(TBodyConditions) :: BodyConditions ! Physical conditions on the floating body
  TYPE(TEnvironment)    :: Env            ! Physical conditions of the environment
  TYPE(TSolver)         :: SolverOpt      ! Solver Option, specified by user in input_solver.txt 
  
  INTEGER                            :: i_problem          ! Index of the current problem
  REAL                               :: omega, wavenumber  ! Wave frequency and wavenumber
  COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB, ZIGS         ! Computed source distribution
  COMPLEX, DIMENSION(:), ALLOCATABLE :: Potential          ! Computed potential

  TYPE(TGREEN)                       :: IGreen              ! Initial Green variables
  
  REAL                               :: tcpu_start
  CHARACTER(LEN=1000)                :: LogTextToBeWritten

  ! Initialization ---------------------------------------------------------------------

  WRITE(*,*) ' '
  WRITE(*,'(A,$)') '  -> Initialisation '

  ! Get working directory from command line argument
  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, wd)
  ELSE
    wd = "."
  END IF

  CALL ReadTMesh(Mesh, TRIM(wd)//'/mesh/')
  ALLOCATE(ZIGB(Mesh%NPanels), ZIGS(Mesh%NPanels))
  ALLOCATE(Potential(Mesh%NPanels*2**Mesh%Isym))

  CALL ReadTBodyConditions            &
  ( BodyConditions,                   &
    Mesh%Npanels*2**Mesh%Isym,        &
    TRIM(wd)//'/Normalvelocities.dat' &
    )

  CALL ReadTEnvironment(Env, file=TRIM(wd)//'/Nemoh.cal')

  call INITIALIZE_GREEN(Mesh,Env%depth,IGreen)

  WRITE(*, *) '. Done !'
  WRITE(*, *) ' '

  ! Solve BVPs and calculate forces ----------------------------------------------------

  WRITE(*, *) ' -> Solve BVPs and calculate forces '
  CALL ReadTSolver(SolverOpt,TRIM(wd))
  WRITE(LogTextToBeWritten,*) 'Linear Solver: ', SolverOpt%SNAME
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
  CALL START_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt',IdAppend)
  WRITE(*, *) ' '

  DO i_problem = 1, BodyConditions%Nproblems
    WRITE(*,'(A,I5,A,I5,A,A,$)') ' Problem ',i_problem,' / ',BodyConditions%Nproblems,' ',CHAR(13)

    omega = BodyConditions%omega(i_problem) ! Wave frequency
    ! Compute wave number
    IF ((Env%depth == INFINITE_DEPTH) .OR. (omega**2*Env%depth/Env%g >= 20)) THEN
      wavenumber = omega**2/Env%g
    ELSE
      wavenumber = X0(omega**2*Env%depth/Env%g)/Env%depth
      ! X0(y) returns the solution of y = x * tanh(x)
    END IF
    !===============
    ! BEM Resolution
    !===============
      CALL SOLVE_POTENTIAL_DIRECT                                              &
      !==========================
      ( Mesh, Env, omega, wavenumber,IGreen,                                   &
        BodyConditions%NormalVelocity(1:Mesh%Npanels*2**Mesh%Isym, i_problem), &
        ZIGB, ZIGS,                                                            &
        Potential(:),SolverOpt,trim(wd))

    !===========================
    ! Post processing and output
    !===========================

    CALL COMPUTE_AND_WRITE_FORCES            &
    !============================
    ( TRIM(wd)//'/mesh/Integration.dat',     &
      Mesh, Env%rho, omega, Potential,       &
      Bodyconditions%Switch_type(i_problem), &
      TRIM(wd)//'/results/Forces.dat'        &
      )
    IF (BodyConditions%Switch_Potential(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_DATA_ON_MESH                                      &
      !=======================
      ( Mesh,                                                      &
        Env%rho*II*omega*Potential(:),                             &
        TRIM(wd)//'/results/pressure.'//string(i_problem)//'.dat'  &
        )
    END IF

    IF (BodyConditions%Switch_Kochin(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_KOCHIN                             &
      !============================
      ( TRIM(wd)//'/mesh/Kochin.dat',                           &
        Mesh, Env, wavenumber, ZIGB, ZIGS,                      &
        TRIM(wd)//'/results/Kochin.'//string(i_problem)//'.dat' &
        )
    END IF
    IF (BodyConditions%Switch_FreeSurface(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION                  &
      !============================================
      ( TRIM(wd)//'/mesh/Freesurface.dat', IGreen,                   &
        Mesh, Env, omega, wavenumber, ZIGB, ZIGS,                    &
        TRIM(wd)//'/results/freesurface.'//string(i_problem)//'.dat' &
        )
    END IF
  END DO
  CALL END_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt')
  WRITE(*,*) '. Done !'
  ! Finalize ---------------------------------------------------------------------------

  DEALLOCATE(ZIGB, ZIGS, Potential)

CONTAINS

  FUNCTION string (i) result (s)
    ! For example 5 -> "00005"
    INTEGER :: i
    CHARACTER(LEN=5) :: s
    WRITE(s, '(I0.5)') i
  END FUNCTION

      
END PROGRAM Main
    

!<<<<<<< HEAD
!!--------------------------------------------------------------------------------------
!!
!!   NEMOH V1.0 - BVP solver - January 2014
!!
!!--------------------------------------------------------------------------------------
!!
!!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!!
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!   You may obtain a copy of the License at
!!
!!       http://www.apache.org/licenses/LICENSE-2.0
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License. 
!!
!!   Contributors list:
!!   - G. Delhommeau
!!   - P. Guével
!!   - J.C. Daubisse
!!   - J. Singh
!!   - A. Babarit  
!!
!!--------------------------------------------------------------------------------------
!!
!    PROGRAM Main
!!   
!    USE MIdentification
!    USE MMesh
!    USE MBodyConditions
!#ifndef GNUFORT
!    USE iflport
!#endif
!    USE SOLVE_BEM
!    USE OUTPUT
!    USE INITIALIZATION
!!    
!    IMPLICIT NONE
!!   ID
!    TYPE(TID)               :: ID
!!   Array sizes
!    INTEGER                 :: NFA,NSYMY    ! Number of panels and symmetry about xOz plane
!!   Body conditions
!    TYPE(TBodyConditions)   :: BodyConditions
!!   Kochin function
!    INTEGER                 :: Ntheta    
!    REAL,DIMENSION(:),ALLOCATABLE :: Theta
!!   Meshes
!    TYPE(TMesh) :: Mesh 
!    TYPE(TMesh) :: MeshFS
!!   Nemoh
!    REAL                    :: T
!    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL,PRESSURE
!    COMPLEX,DIMENSION(:),ALLOCATABLE :: HKochin
!!   Results
!    REAL :: XEFF,YEFF
!    INTEGER :: Nintegration
!    REAL,DIMENSION(:,:),ALLOCATABLE :: NDS
!    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: Force
!    REAL,DIMENSION(:),ALLOCATABLE :: line
!!   Locals
!    REAL                    :: PI
!    INTEGER                 :: c,d,M,l,i,j,k
!    REAL                    :: Discard
!    COMPLEX,PARAMETER       :: II=CMPLX(0.,1.)
!!   CPU TIME
!    REAL                    :: tcpu_start, tcpu_finish, comptime,relcomptime
!    INTEGER,DIMENSION(8)    :: DATETIMEVAL
!    CHARACTER(20)           :: SOLVER_NAME, GreenFun_name
!!
!!   --- Initialisation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!   
!    PI=4.*ATAN(1.)
!    WRITE(*,*) ' '
!    WRITE(*,'(A,$)') '  -> Initialisation ' 
!!   Read case ID
!    CALL ReadTID(ID,'ID.dat') 
!!   Read Mesh
!    CALL ReadTMesh(Mesh,ID)  
!!   Read Body Conditions
!    CALL ReadTBodyConditions(BodyConditions,Mesh%Npanels*2**Mesh%Isym,ID%ID(1:ID%lID)//'/Normalvelocities.dat') 
!!   Initialise Nemoh
!    CALL INITIALIZE(ID,NFA,NSYMY,XEFF,YEFF,Mesh)
!    ALLOCATE(NVEL(NFA*2**NSYMY),PRESSURE(NFA*2**NSYMY))
!!    WRITE(*,'(A,$)') '.'
!!   Initialise Force matrix
!    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Integration.dat')
!    READ(10,*) Nintegration
!    ALLOCATE(NDS(Nintegration,NFA*2**NSYMY))
!    DO i=1,Nintegration
!        READ(10,*) (NDS(i,j),j=1,NFA*2**NSYMY)
!    END DO
!    CLOSE(10)
!!   Initialise Kochin function calculation
!    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Kochin.dat')
!    READ(10,*) Ntheta 
!    ALLOCATE(Theta(Ntheta))
!    IF (Ntheta.GT.0) THEN        
!        DO j=1,Ntheta
!            READ(10,*) Theta(j)
!        END DO        
!    END IF
!    ALLOCATE(HKochin(NTheta))
!    CLOSE(10)
!!   Initialise free surface calculation points
!    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Freesurface.dat')
!    READ(10,*) MeshFS%Npoints,MeshFS%Npanels 
!    IF (MeshFS%Npoints.GT.0) THEN
!        CALL CreateTMesh(MeshFS,MeshFS%Npoints,MeshFS%Npanels,1)
!        DO j=1,MeshFS%Npoints
!            READ(10,*) MeshFS%X(1,j),MeshFS%X(2,j)
!        END DO
!        DO j=1,MeshFS%Npanels
!            READ(10,*) MeshFS%P(1,j),MeshFS%P(2,j),MeshFS%P(3,j),MeshFS%P(4,j)
!        END DO
!    END IF
!    CLOSE(10)
!!   Initialise results table
!    ALLOCATE(Force(Nintegration,Bodyconditions%Nproblems))
!    Force(:,:)=0.
!   ! WRITE(*,*) '. Done !'
!    WRITE(*,*) ' '
!!
!!   --- Solve BVPs and calculate forces -------------------------------------------------------------------------------------------------------------------------------------------------
!!
!    CALL CPU_TIME(tcpu_start)
!    call date_and_time(VALUES=DATETIMEVAL)
!
!  !  OPEN(100,FILE=ID%ID(1:ID%lID)//'/computation_time.txt' ) 
!    WRITE(100,*) 'Date: ',DATETIMEVAL(3),'-',DATETIMEVAL(2),'-',DATETIMEVAL(1)
!    WRITE(100,*) 'Time: ',DATETIMEVAL(5),':',DATETIMEVAL(6),':',DATETIMEVAL(7)
!    WRITE(*,*) ' -> Solve BVPs and calculate forces ' 
!    WRITE(*,*) ' '
!    DO j=1,BodyConditions%Nproblems 
!!         WRITE(*,'(A,I5,A,I5,A,A)',ADVANCE='NO') ' Problem ',j,' / ',BodyConditions%Nproblems,' . . Processing',CHAR(13)
!        WRITE(*,'(A,I5,A,I5,A,A,$)') ' Problem ',j,' / ',BodyConditions%Nproblems,' ',CHAR(13)
!        DO c=1,Mesh%Npanels*2**Mesh%Isym
!            NVEL(c)=BodyConditions%NormalVelocity(c,j)
!        END DO
!!       Solve BVP
!        CALL SOLVE_BVP(j,ID,2.*PI/BodyConditions%Omega(j),NVEL,PRESSURE,BodyConditions%Switch_Kochin(j),NTheta,Theta,HKochin,BodyConditions%Switch_Freesurface(j),MeshFS,BodyConditions%Switch_Potential(j))
!!       Calculate force coefficients
!        DO i=1,Nintegration
!            DO c=1,Mesh%Npanels*2**Mesh%Isym
!                Force(i,j)=Force(i,j)-PRESSURE(c)*NDS(i,c)
!            END DO
!        END DO
!!~         WRITE(*,'(A,I5,A,I5,A,A)',ADVANCE="NO") ' Problem ',j,' / ',BodyConditions%Nproblems,' . . Done !',CHAR(13)
!!~         WRITE(*,*) '. Done !'      
!    END DO    
!    CALL CPU_TIME(tcpu_finish)
!    comptime=tcpu_finish-tcpu_start
!!    relcomptime=comptime/Mesh%Npanels/BodyConditions%Nproblems
!    WRITE(100,*) 'Computation time', comptime, ' [s]'
!    CLOSE(100)
!    WRITE(*,*) 'Computation time=',comptime,' [s]'
!!    WRITE(*,*) 'Relative Comp. time=',relcomptime,'[ s/Npanel/Nproblem]'
!    WRITE(*,*) '. Done !'      
!!    WRITE(*,*) ' ' 
!!    CLOSE(10)
!!    CLOSE(12)
!!
!!   --- Save results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!
!    WRITE(*,*) ' -> Save results ' 
!    WRITE(*,*) ' '
!    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/Forces.dat')
!    ALLOCATE(line(BodyConditions%Nproblems*2))
!    DO c=1,Nintegration
!        DO j=1,BodyConditions%Nproblems
!            IF (Bodyconditions%Switch_type(j).NE.-1.) THEN
!                line(2*j-1)=ABS(Force(c,j))
!                line(2*j)=ATAN2(IMAG(Force(c,j)),REAL(Force(c,j)))
!            ELSE
!                line(2*j-1)=IMAG(Force(c,j))/BodyConditions%Omega(j)
!                line(2*j)=-REAL(Force(c,j))
!            END IF
!        END DO
!        WRITE(10,*) (line(j),j=1,2*BodyConditions%Nproblems)      
!    END DO
!    CLOSE(10)    
!!
!!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!
!   
!    IF (MeshFS%Npoints.GT.0) CALL DeleteTMesh(MeshFS)
!    CALL DEALLOCATE_DATA
!    DEALLOCATE(line)
!    DEALLOCATE(NVEL,PRESSURE)    
!    DEALLOCATE(Force)
!    DEALLOCATE(Theta,HKochin)
!
!    END PROGRAM Main
!!----------------------------------------------------------------

