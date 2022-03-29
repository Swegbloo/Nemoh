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

MODULE MLogFile
        
  IMPLICIT NONE
  
  PUBLIC       ::  WRITE_LOGFILE,START_RECORD_TIME
  !  CPU TIME
  INTEGER,parameter, public  :: IdAppend=1      !Append writing to a file
  INTEGER,parameter, public  :: IdStartLog=0    !Starting a new file
  INTEGER,parameter, public  :: IdprintTerm=1   !print to the terminal
  INTEGER,parameter, public  :: IdNoprintTerm=0 !Not print to the terminal
CONTAINS

  SUBROUTINE WRITE_LOGFILE(logfile,textToBeWritten,FlagW,FlagWTerm)

        CHARACTER(LEN=*) logfile,textToBeWritten
        INTEGER u,FlagW,FlagWTerm
        IF (FlagW==IdStartLog) THEN
                OPEN(NEWUNIT=u, FILE=logfile, ACTION='WRITE')
        ELSE
                OPEN(NEWUNIT=u, FILE=logfile, ACTION='WRITE',POSITION='APPEND')
        ENDIF
        WRITE(u,*) textToBeWritten
        IF (FlagWTerm == IdprintTerm) THEN
        WRITE(*,*) textToBeWritten
        END IF
        CLOSE(u)
  END SUBROUTINE

  SUBROUTINE START_RECORD_TIME(tcpu_start,logfile,FlagW)
        CHARACTER(LEN=*) logfile
        INTEGER,DIMENSION(8)    :: DATETIMEVAL
        CHARACTER(LEN=1000)     :: textToBeWritten
        INTEGER                 :: FlagW
        REAL                    :: tcpu_start

        CALL CPU_TIME(tcpu_start)
        CALL DATE_AND_TIME(VALUES=DATETIMEVAL)
  
        CALL WRITE_LOGFILE(logfile,' STARTING TIME:',FlagW,IdprintTerm)
        WRITE(textToBeWritten,*)  ' Date: ',DATETIMEVAL(3),'-',DATETIMEVAL(2),'-',DATETIMEVAL(1)
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
        WRITE(textToBeWritten,*) ' Time: ',DATETIMEVAL(5),':',DATETIMEVAL(6),':',DATETIMEVAL(7)
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
  END SUBROUTINE        
  
  SUBROUTINE END_RECORD_TIME(tcpu_start,logfile)
        CHARACTER(LEN=*) logfile
        INTEGER,DIMENSION(8)    :: DATETIMEVAL
        CHARACTER(LEN=1000)     :: textToBeWritten
        INTEGER                 :: FlagW
        REAL                    :: tcpu_start,tcpu_finish,comptime
        
        CALL CPU_TIME(tcpu_finish)
        comptime=tcpu_finish-tcpu_start
        WRITE(textToBeWritten,*) 'Computation time', comptime, ' [s]'
        CALL WRITE_LOGFILE(logfile,TRIM(textToBeWritten),IDAppend,IdprintTerm)
  END SUBROUTINE      
END MODULE
