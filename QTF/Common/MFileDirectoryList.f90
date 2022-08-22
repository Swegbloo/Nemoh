!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!   Contributors list:
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!    List of input/output file
!--------------------------------------------------------------------------------------
MODULE MFileDirectoryList  

CHARACTER(LEN=*), PARAMETER      ::PreprocDir='QTFPreprocOut'
CHARACTER(LEN=*), PARAMETER      ::LogFILE   ='logfileQTF.txt'
CHARACTER(LEN=*), PARAMETER      ::TotPotFILE='TotalPotentialBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::TotVelFILE='TotalVelocityBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::RadPotFILE='RadiationPotentialBodyWLine.bin'
CHARACTER(LEN=*), PARAMETER      ::RadVelFILE='RadiationVelocityBodyWLine.bin'


CHARACTER(LEN=*), PARAMETER      ::  OutQTFDir ='/results/QTF/'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDM ='QTFM_DUOK.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDP ='QTFP_DUOK.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBM='QTFM_HASBO.dat'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBP='QTFP_HASBO.dat'

CHARACTER(LEN=*), PARAMETER      ::  OutFileDM_term ='QTFM_DUOK_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileDP_term ='QTFP_DUOK_term_'

CHARACTER(LEN=*), PARAMETER      ::  OutFileHBM_term='QTFM_HASBO_term_'
CHARACTER(LEN=*), PARAMETER      ::  OutFileHBP_term='QTFP_HASBO_term_'



!Output variables
  !DUOK      : Quadratic QTF
  !HASBO     : Potential QTF-Haskind on Body term
  !HASFS     : Potential QTF-Haskind on Free-Surface in Finite domain
  !HASFS_ASYM: Potential QTF-Haskind on Free-Surface in Asymptotic domain
CONTAINS
   SUBROUTINE  make_directory(dirname)
          CHARACTER(LEN=*),       INTENT(IN) ::dirname
          LOGICAL                            ::existdir
          INQUIRE (DIRECTORY=dirname, EXIST=existdir)       
          IF (.NOT.existdir) CALL SYSTEM('mkdir '//dirname)
   END SUBROUTINE


END MODULE
