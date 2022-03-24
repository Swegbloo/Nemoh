!-------------------------------------------------------------------------------------------------
!
!  NEMOH v1.0 - Hydrostatic Calculation - July 2021
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
!   - R Kurnia  
!
!--------------------------------------------------------------------------------------

 PROGRAM Hydrostatic

    USE MEnvironment
    USE MIdentification
    USE MMesh

#ifndef GNUFORT
    USE iflport
#endif

    IMPLICIT NONE
    TYPE(TID) :: ID,DSCRPT              ! Calculation identification data
    TYPE(TMesh) :: Mesh                 ! Mesh data
    TYPE(TEnvironment) :: Environment   ! Environment data  

!   Maillage proprement dit
        INTEGER,PARAMETER :: NFMX=20000 ! Nombre de facettes max
        INTEGER,PARAMETER :: NPMX=20000 ! Nombre de points max
        INTEGER :: Nmailmx              ! Nombre de facettes du maillage std max
!   Maillage du corps
        INTEGER :: NF,NP
        INTEGER,DIMENSION(4,NFMX) :: Facette
	REAL,DIMENSION(NPMX) :: X,Y,Z
!   Partie immergee du maillage
        INTEGER :: NFm,NPm
        INTEGER,DIMENSION(4,NFMX) :: Facettem
	REAL,DIMENSION(NPMX) :: Xm,Ym,Zm
!   Calcul hydrostatique
	REAL DEPLACEMENT,XF,YF,ZF,SF
	REAL,DIMENSION(6,6) :: KH
	REAL :: xG,yG,zG
	REAL :: RHO,G
!   Calcul coque
	REAL,DIMENSION(3,3) :: Icoque
	REAL,DIMENSION(3) :: Gcoque,CDG
        
        INTEGER         :: i,j,Nsym
!
!   --- Initialize and read input datas ----------------------------------------------------------------------------------------
!
    CALL ReadTID(ID)
    CALL ReadTMesh(Mesh,ID) 
    CALL ReadTEnvironment(Environment,TRIM(ID%ID)//'/Nemoh.cal')
    RHO=Environment%RHO
    G=Environment%G
    NP=Mesh%Npoints
    NF=Mesh%Npanels
    Nsym=Mesh%Isym
    OPEN(10,FILE='Mesh.cal')
    READ(10,*) DSCRPT%ID
    DSCRPT%lID=LNBLNK(DSCRPT%ID)
    READ(10,*)
    READ(10,*) 
    READ(10,*) xG,yG,zG
    READ(10,*) 
    CLOSE(10)
    DO j=1,NP
        X(j)=Mesh%X(1,j)-xG
        Y(j)=Mesh%X(2,j)-yG
        Z(j)=Mesh%X(3,j)
    END DO
    DO j=1,NF
        DO i=1,4
        FACETTE(i,j)=Mesh%P(i,j);
        END DO
    END DO
    CALL HYDRO(X,Y,Z,NP,FACETTE,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xm,Ym,Zm,NPm,FACETTEm,NFm,RHO,G)
    DO j=1,NP
                X(j)=X(j)+xG
                Y(j)=Y(j)+yG
    END DO
    IF (Nsym.EQ.1) THEN
        DEPLACEMENT=2.0*DEPLACEMENT
        YF=0.
        SF=2.0*SF
        KH(3,3)=2.*KH(3,3)
        KH(3,4)=0.
        KH(4,3)=0.
        KH(3,5)=2.*KH(3,5)
        KH(5,3)=KH(3,5)
        KH(4,4)=2.*KH(4,4)
        KH(4,5)=0.
        KH(5,4)=0.
        KH(5,5)=2.*KH(5,5)
    END IF
        KH(4,4)=KH(4,4)+deplacement*RHO*G*(ZF-ZG)
        KH(5,5)=KH(5,5)+deplacement*RHO*G*(ZF-ZG)
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/KH.dat')
        DO i=1,6
                WRITE(10,'(6(1X,E14.7))') (KH(i,j),j=1,6)
        END DO
        CLOSE(10)
        write(*,*) ' -> Calculate hull mass and inertia '
        WRITE(*,*) ' '
        CDG(1)=xG
        CDG(2)=yG
        CDG(3)=zG      
        CALL coque(X,Y,Z,NP,facette,NF,Deplacement,Icoque,Gcoque,CDG,Nsym,rho)
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/GC_hull.dat')
        WRITE(10,'(3(1X,E14.7))') Gcoque(1),Gcoque(2),Gcoque(3)
        CLOSE(10)
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Inertia_hull.dat')
        DO i=1,3
                WRITE(10,'(3(1X,E14.7))') (Icoque(i,j),j=1,3)
        END DO
        CLOSE(10)
        WRITE(*,'(A,I3)') '   - Coordinates of buoyancy centre '
        WRITE(*,'(A,F7.3,A)') '     XB = ',XF+xG,'m'
        WRITE(*,'(A,F7.3,A)') '     YB = ',YF+yG,'m'
        WRITE(*,'(A,F7.3,A)') '     ZB = ',ZF,'m'
        WRITE(*,'(A,E14.7,A)') '    - Displacement  = ',DEPLACEMENT,' m^3'
        WRITE(*,'(A,E14.7,A)') '    - Waterplane area  = ',SF, ' m^2'
        WRITE(*,'(A,E14.7,A)') '    - Mass =',DEPLACEMENT*RHO, ' Kg'
        WRITE(*,*) ' '

        IF ((ABS(XF).GT.1E-02).OR.(ABS(YF).GT.1E-02)) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' !!! WARNING !!! '
            WRITE(*,*) ' '
            WRITE(*,'(A,I3)') ' Buoyancy center and gravity center are not vertically aligned. '
            WRITE(*,*) ' This is not an equilibrium position.'
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' XF = ',XF+xG,' XG = ',xG
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' YF = ',YF+yG,' YG = ',yG
        END IF
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Hydrostatics.dat')
        WRITE(10,'(A,F7.3,A,F7.3)') ' XF = ',XF+xG,' - XG = ',xG
        WRITE(10,'(A,F7.3,A,F7.3)') ' YF = ',YF+yG,' - YG = ',yG
        WRITE(10,'(A,F7.3,A,F7.3)') ' ZF = ',ZF,' - ZG = ',zG
        WRITE(10,'(A,E14.7)') ' Displacement = ',DEPLACEMENT
        WRITE(10,'(A,E14.7)') ' Waterplane area = ',SF
        WRITE(10,'(A,E14.7)') ' Mass =',DEPLACEMENT*RHO
        CLOSE(10)



end program Hydrostatic
