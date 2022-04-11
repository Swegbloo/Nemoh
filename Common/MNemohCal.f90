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
!   - R. Kurnia
!
!--------------------------------------------------------------------------------------

        MODULE MNemohCal
         
          USE MEnvironment, only: TEnvironment
          IMPLICIT NONE
      
          PUBLIC      :: READ_TNEMOHCAL

          INTEGER, PARAMETER :: IdRadFreq=1     !unit rad/s
          INTEGER, PARAMETER :: IdFreqHz =2     !unit Hz
          INTEGER, PARAMETER :: IdPeriod =3     !unit s

          TYPE Twaveinput
              INTEGER  :: FreqType
              INTEGER  :: NFreq,NBeta
              REAL     :: Freq1,Freq2,Beta1,Beta2
          END TYPE Twaveinput
           
          TYPE TIRF
              INTEGER  :: Switch
              REAL     :: time_step,duration
          END TYPE TIRF
          
          TYPE TKochin
              INTEGER  :: Switch,Ntheta
              REAL     :: min_theta,max_theta   !degree
          END TYPE TKochin
          
          TYPE TFreesurface
              INTEGER     :: Switch
              INTEGER     :: NX,NY      !Number of points
              REAL        :: Lx,LY      !Length of domain
          END TYPE TFreesurface

          TYPE TOptOutput
              INTEGER            :: Switch_POTENTIAL
              INTEGER            :: Switch_SourceDistr
              TYPE(TIRF)         :: IRF
              TYPE(TKochin)      :: Kochin
              TYPE(TFreesurface) :: Freesurface
          END TYPE TOptOutput

          TYPE TBCase
              INTEGER :: ICase
              REAL,DIMENSION(3) :: Direction,Axis 
          END TYPE TBCase

          TYPE Tbodyinput
               CHARACTER(LEN=100)       :: meshfile
               INTEGER                  :: Npoints,NPanels
               INTEGER                  :: NRadiation
               INTEGER                  :: NIntegration
               TYPE(TBCase),ALLOCATABLE :: RadCase(:)
               TYPE(TBCase),ALLOCATABLE :: IntCase(:)
          END TYPE Tbodyinput

          TYPE TNemCal
              TYPE(Tenvironment)           :: Env
              INTEGER                      :: Nbodies,Nradtot,Nintegtot
              TYPE(Tbodyinput),ALLOCATABLE :: bodyinput(:)
              TYPE(Twaveinput)             :: waveinput
              TYPE(TOptOutput)             :: OptOUTPUT
          END TYPE TNemCal

        CONTAINS

         SUBROUTINE READ_TNEMOHCAL(wd,InpNEMOHCAL)
           
          CHARACTER(LEN=*),     INTENT(IN)      :: wd
          TYPE(TNemCal),        INTENT(OUT)     :: InpNEMOHCAL
          
          !Local var
          INTEGER ufile,I,J,K
          LOGICAL :: exist_dir


          OPEN(NEWUNIT=ufile ,FILE=TRIM(wd)//'/Nemoh.cal')
           READ(ufile,*) !--- Environment ----------------------------!
           READ(ufile,*) InpNEMOHCAL%Env%RHO
           READ(ufile,*) InpNEMOHCAL%Env%G
           READ(ufile,*) InpNEMOHCAL%Env%Depth
           READ(ufile,*) InpNEMOHCAL%Env%Xeff, InpNEMOHCAL%Env%Yeff
           READ(ufile,*) !--- Description of floating bodies-----------!
           READ(ufile,*) InpNEMOHCAL%Nbodies
           !! ALLOCATE body input variable
           ALLOCATE(InpNEMOHCAL%bodyinput(InpNEMOHCAL%Nbodies))
           !!
           InpNEMOHCAL%Nradtot=0
           InpNEMOHCAL%Nintegtot=0
           
           DO I=1,InpNEMOHCAL%Nbodies
             READ(ufile,*) !----Body(I)--------------------------------!
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%meshfile
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%Npoints,           &
                                  InpNEMOHCAL%bodyinput(I)%Npanels
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%NRadiation
             InpNEMOHCAL%Nradtot=InpNEMOHCAL%Nradtot+                  &
                                 InpNEMOHCAL%bodyinput(I)%NRadiation
             ALLOCATE(InpNEMOHCAL%bodyinput(I)%RadCase(                &
                              InpNEMOHCAL%bodyinput(I)%NRadiation))
             DO J=1,InpNEMOHCAL%bodyinput(I)%NRadiation
              READ(ufile,*)InpNEMOHCAL%bodyinput(I)%RadCase(J)%ICase,  &
              (InpNEMOHCAL%bodyinput(I)%RadCase(J)%Direction(K),K=1,3),&
              (InpNEMOHCAL%bodyinput(I)%RadCase(J)%Axis(K),K=1,3)
             END DO
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%NIntegration
             InpNEMOHCAL%Nintegtot=InpNEMOHCAL%Nintegtot+              &
                                 InpNEMOHCAL%bodyinput(I)%NIntegration
             ALLOCATE(InpNEMOHCAL%bodyinput(I)%IntCase(                &
                              InpNEMOHCAL%bodyinput(I)%NIntegration))
             DO J=1,InpNEMOHCAL%bodyinput(I)%NIntegration
              READ(ufile,*)InpNEMOHCAL%bodyinput(I)%IntCase(J)%ICase,  &
              (InpNEMOHCAL%bodyinput(I)%IntCase(J)%Direction(K),K=1,3),&
              (InpNEMOHCAL%bodyinput(I)%IntCase(J)%Axis(K),K=1,3)
             END DO
             READ(ufile,*) !For additional info, ie. for generalize mode, Not implemented yet 
           END DO
           
             READ(ufile,*) !--- Load cases to be solved ---------------!
             READ(ufile,*) InpNEMOHCAL%waveinput%FreqType,             &
                           InpNEMOHCAL%waveinput%NFreq,                &
                           InpNEMOHCAL%waveinput%Freq1,                &
                           InpNEMOHCAL%waveinput%Freq2
             READ(ufile,*) InpNEMOHCAL%waveinput%NBeta,                &
                           InpNEMOHCAL%waveinput%Beta1,                &
                           InpNEMOHCAL%waveinput%Beta2
             READ(ufile,*) !--- Post Processing -----------------------!
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%IRF%Switch,           &
                           InpNEMOHCAL%OptOUTPUT%IRF%time_step,        &
                           InpNEMOHCAL%OptOUTPUT%IRF%duration
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Switch_POTENTIAL
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Kochin%Ntheta,        &
                           InpNEMOHCAL%OptOUTPUT%Kochin%min_theta,     &
                           InpNEMOHCAL%OptOUTPUT%Kochin%max_theta
             IF (InpNEMOHCAL%OptOUTPUT%Kochin%Ntheta.GT.0 ) THEN
                           InpNEMOHCAL%OptOUTPUT%Kochin%Switch=1
             ELSE
                           InpNEMOHCAL%OptOUTPUT%Kochin%Switch=0
             ENDIF

             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Freesurface%Nx,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Ny,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Lx,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Ly
             IF (InpNEMOHCAL%OptOUTPUT%Freesurface%Nx.GT.0 ) THEN
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Switch=1
             ELSE
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Switch=0
             ENDIF
             READ(ufile,*)! ---QTF----
             READ(ufile,*)InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr
          CLOSE(ufile)
        
          IF(InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
            INQUIRE (DIRECTORY=TRIM(wd)//'/results/sources',           &
                       EXIST=exist_dir) 
            IF (.NOT.exist_dir) CALL SYSTEM('mkdir '//TRIM(wd)//       &
                                                    '/results/sources')
          END IF
         END SUBROUTINE

        END MODULE 


           
