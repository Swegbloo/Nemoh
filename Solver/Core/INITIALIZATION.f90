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
!   - J. Singh  
!
!--------------------------------------------------------------------------------------
MODULE INITIALIZATION

IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------
    SUBROUTINE INITIALIZE(ID,NF,NSYM,XF,YF,Mesh)
!
    USE MIDENTIFICATION
    USE COM_VAR
    USE PREPARE_MESH
    USE MMesh
!
    IMPLICIT NONE
!   ID
    TYPE(TID) :: ID
!   Geometry
    INTEGER :: NF,NSYM,c,i,M
    REAL :: XF,YF
    TYPE(TMesh) :: Mesh
    CHARACTER(20) :: SOLVER_NAME    
!
!   Read input file and geometry
    OPEN(10,file=ID%ID(1:ID%lID)//'/Nemoh.cal',form='formatted',status='old')
    READ(10,*)
    READ(10,*) RHO
    READ(10,*) G
    READ(10,*) DEPTH
    READ(10,*) XF,YF
    READ(10,*)
    READ(10,*)
     DO c=1,Mesh%Nbodies
        DO i=1,3
            READ(10,*)
        END DO
        READ(10,*) M ! Ndof  
        DO i=1,M
            READ(10,*)
        END DO
        READ(10,*) M !Nforce
        DO i=1,M
            READ(10,*)
        END DO
        READ(10,*) M !Nadditional line
        DO i=1,M
            READ(10,*)
        END DO
    END DO
    DO c=1,9
    READ(10,*)
    END DO
    READ(10,*) Switch_Sources
    CLOSE(10)
    OPEN(10,file=ID%ID(1:ID%lID)//'/input_solver.txt',form='formatted',status='old')
    READ(10,*) Indiq_solver
    CLOSE(10)   
    XEFF=XF
    YEFF=YF
    NFA=Mesh%Npanels
    NP=Mesh%Npoints
    NSYMY=Mesh%Isym
    IF (NSYMY.NE.1) NSYMY=0
    IF (NSYMY.EQ.1) YEFF=0
    NF=NFA
    NSYM=NSYMY
!   Initialise Nemoh
    CALL ALLOCATE_DATA
    w_previous=-1.
    CALL PRE_PROC_MESH(Mesh)
!
    IF (Indiq_solver.EQ.0) THEN
            SOLVER_NAME='GAUSS ELIMINATION'
    ELSE IF (Indiq_solver.EQ.1) THEN
            SOLVER_NAME='LU DECOMPOSITION'
    ELSE 
            SOLVER_NAME='NOT AVAILABLE'
            WRITE(*,*) 'SOLVER 2 is not available'
            STOP
    END IF
   
    WRITE(*,*) 'LINEAR SOLVER: ', SOLVER_NAME
    OPEN(100,FILE=ID%ID(1:ID%lID)//'/computation_time.txt') 
    WRITE(100,*) 'LINEAR SOLVER=', SOLVER_NAME
   ! CLOSE(100)
    END SUBROUTINE INITIALIZE  

END MODULE INITIALIZATION
