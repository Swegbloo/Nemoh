!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!    - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!    - Fabien Robaux (EDF/INNOSEA)
!    - G.DELHOMMEAU, THESE DE CHEN XIAO-BO(1988)    
!--------------------------------------------------------------------------------------
  
  PROGRAM Solver2

  INTEGER :: Nbodies, NP, c,M, contrib,Loutduok,Louthasbo,Louthasfs,NW1,NW2,NW,NDOF,nhaskind,dw,npasr
  REAL :: trash,RCEREXT, wmin, wmax   ! wmin, wmax added by RK
  CHARACTER*20 ID
  CHARACTER*255 :: EXEDIR
  CHARACTER*255:: CMD
  INTEGER       :: lengthArg
  logical :: dir_e ! to check if directories exist
  Nbodies = 1 ! Only one body for now
  
  OPEN(UNIT=7,FILE="ID.dat")
  READ(7,*) lID
  READ(7,*) ID
  CLOSE(7)
  
  ! !   check if dir exists
  INQUIRE( DIRECTORY=ID(1:lID)//'/results/QTF', EXIST=dir_e )  !File is changed to Directory by RK
  IF ( .NOT.dir_e ) THEN
	CALL SYSTEM('mkdir '//ID(1:lID)//'/results/QTF')
  END IF
  
  ! !   get exe directory to call exe
  CALL GET_COMMAND_ARGUMENT(0,EXEDIR)
  lengthArg=LEN_TRIM(EXEDIR)-10
  EXEDIR = EXEDIR(1:lengthArg)
  IF (LEN_TRIM(EXEDIR) == 0) THEN
    EXEDIR = './../../../../../bin/' ! For specifying the directory manually        
  END IF
  !WRITE(*,*) EXEDIR
  ! ! 	read number of frequencies
  OPEN(10,FILE=ID(1:lID)//'/Nemoh.cal')
    DO c=1,7
	READ(10,*)
    END DO

    DO c=1,Nbodies
	DO i=1,3
	    READ(10,*)
	END DO
	READ(10,*) M      !Ndof
	DO i=1,M
	    READ(10,*)
	END DO
	READ(10,*) M      !N gen forces
	DO i=1,M
	    READ(10,*)
	END DO
	READ(10,*) M      !N additional lines
	DO i=1,M
	    READ(10,*)    
	END DO
    END DO
    READ(10,*)    
    READ(10,*) NP,wmin,wmax        ! N of wave freq, wmin and wmax
    DO I=1,8
      READ(10,*)
    END DO
    READ(10,*) LQTFP
! ! 	read what contrib to compute
    READ(10,*) contrib
    READ(10,*)
    READ(10,*) Loutduok
    READ(10,*) Louthasbo
    READ(10,*) Louthasfs
    if ( Louthasfs == 1) then
      READ(10,*) NW1,NW2,NDOF
    else
      NW1=0
      NW2=0
      NDOF=0
    endif
  CLOSE(10)
  
  if (NW2 > NW1) then
    NW=NW1
    NW1=NW2
    NW2=NW
  endif
  if (NW1==NW2) then
    NW2=NW1-1
  endif
  
  if (NW1 + NW2 > NP .AND. LQTFP ==1 ) then
    dw=NW1-NW2
    NW1=NP/2
    NW2=max(1,NW1-dw)
  endif
  
  if (Loutduok+Louthasbo+Louthasfs > 0) then
    inquire( DIRECTORY=ID(1:lID)//'/results/QTF/Appendix', exist=dir_e ) ! file replaced with Directory
    if ( .NOT. dir_e ) then
	  CALL SYSTEM('mkdir '//ID(1:lID)//'/results/QTF/Appendix')
    endif
  endif
  
  
  Write(*,*)'Number of frequencies = '
  Write(*,*) NP
  
  
  ICER0=10
  if (contrib>2) then
  OPEN(2,FILE=ID(1:lID)//'/mesh/SF_L12.dat')
  !OPEN(2,FILE=ID(1:lID)//'/mesh/SF_L12_2.dat')
  do i=1,7
    read(2,*)
  enddo
  READ(2,*) trash,trash,trash,trash,trash,trash,npasr
  !READ(2,*) trash,trash,trash,trash,npasr
  CLOSE(2)
  NRCER=npasr-ICER0-10
  else
  NRCER=102
  endif
  
  OPEN(1,FILE='avs.txt')
  write(1,*)NP
  DO 777 I=1,NP
  write(1,*)I
777     CONTINUE
  write(1,*) LQTFP
  write(1,*) contrib
  write(1,*) Loutduok
  write(1,*) Louthasbo
  write(1,*) Louthasfs,NW1,NW2,NDOF
  write(1,*) NRCER
  DO ICER=0,NRCER-1
    write(1,*) ICER+ICER0
  END DO
  CLOSE(1)
  write(*,*) 
  write(*,*) 'duok   ##################################'
  
  CALL SYSTEM(TRIM(EXEDIR)//'/duokap <avs.txt')
  
  IF(contrib>1) THEN
    write(*,*) 
    write(*,*) 'hasbo   ##################################'
    CALL SYSTEM(TRIM(EXEDIR)//'/hasbo <avs.txt')
  ENDIF
  
  IF(contrib>2) THEN
    write(*,*) 
    write(*,*) 'hasfs   ##################################'
    CALL SYSTEM(TRIM(EXEDIR)//'/hasfs <avs.txt')
    write(*,*) 'hasfscalc    #############################'
    CALL SYSTEM(TRIM(EXEDIR)//'/hasfscalc <avs.txt ')
    write(*,*) 
    write(*,*) 'hasfsasym    #############################'
    CALL SYSTEM(TRIM(EXEDIR)//'/asymp <avs.txt')
  ENDIF

! ! Clean old results !
   inquire( file=ID(1:lID)//'/QTF/WRH.dat', exist=dir_e )
   if ( dir_e ) then
        CALL SYSTEM('rm '//ID(1:lID)//'/QTF/WRH.dat')
   endif
   inquire( file=ID(1:lID)//'/avs.txt', exist=dir_e )
   if ( dir_e ) then
       CALL SYSTEM('rm '//ID(1:lID)//'/avs.txt')
   endif
!   CALL SYSTEM('rm '//ID(1:lID)//'/QTF/*.RES 2> /dev/null')
   CALL SYSTEM('rm '//ID(1:lID)//'/QTF/*.wat 2> /dev/null') 

  END
