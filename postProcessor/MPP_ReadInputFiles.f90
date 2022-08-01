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
MODULE MPP_ReadInputFiles

IMPLICIT NONE

TYPE TMech
    REAL,ALLOCATABLE,DIMENSION(:,:) :: MassMat          !Mass-Inertia Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat         !Stifness-Matrix
    REAL,ALLOCATABLE,DIMENSION(:,:) :: StiffMat_EXT     !Additional Stifness Matrix i.e: mooring
    REAL,ALLOCATABLE,DIMENSION(:,:) :: DampCoefMat_EXT  !Additional damping coefficients
END TYPE TMech
CONTAINS
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
