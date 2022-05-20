!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - Postprocessing module
!   MAIN PROGRAM
!   Description:  producing QTF output as the WAMIT QTF output format
!   Contributors list:
!    - Ruddy Kurnia (ECN)
!--------------------------------------------------------------------------------------
PROGRAM Main
!       Used modules
        USE MIdentification
        USE MNemohCal,          ONLY: TNemCal,READ_TNEMOHCAL
        USE MQPostproc,         ONLY: READWRITE_QTFDATA
!
        IMPLICIT NONE
!
!       Variables declaration
        TYPE(TID)       :: ID
        TYPE(TNemCal)   :: inpNEMOHCAL
!
!       Read input files
        CALL ReadTID(ID)
        CALL READ_TNEMOHCAL(TRIM(ID%ID),inpNEMOHCAL)
        CALL READWRITE_QTFDATA(inpNEMOHCAL,TRIM(ID%ID))

END PROGRAM Main
