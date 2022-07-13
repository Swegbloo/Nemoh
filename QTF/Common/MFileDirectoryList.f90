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

END MODULE
