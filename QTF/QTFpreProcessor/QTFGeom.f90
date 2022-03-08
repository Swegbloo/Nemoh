!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!   - Gerard Delhommeau
!--------------------------------------------------------------------------------------

MODULE MQTFGeometry

!   writes FA.RES
! construit le contour = la ligne de flottaison du corps

CONTAINS

    SUBROUTINE QTFGEOM(ID,Mesh) 
    
      USE MIdentification
      USE MMesh
      USE QTFCOM_VAR                      ! Global variable declaration
      
      IMPLICIT NONE  
      
      TYPE(TID) :: ID                     ! Calculation identification data
      TYPE(TMesh) :: Mesh                 ! Mesh data
      
      CALL ProcessMesh(Mesh)              ! updating the global variables
   
    END SUBROUTINE QTFGEOM
    
    
    SUBROUTINE ProcessMesh(Mesh)
    
      USE MMesh
      USE QTFCOM_VAR
      
      IMPLICIT NONE  
      TYPE(TMesh) :: Mesh 
      
      INTEGER :: i,j,L,I1,NP1
      REAL :: GL,BL,DNUL,DST
      INTEGER :: KK0(5)
      
      IMX = Mesh%Npanels
      NP = Mesh%Npoints
      NSYMY = Mesh%Isym
      
      I1=IMX+1
      ! LECTURES DES POINTS ET CONNECTIVITES
      DO i=1,NP             ! read nodes
	X(i)=Mesh%X(1,i)
	Y(i)=Mesh%X(2,i)
	Z(i)=Mesh%X(3,i)
      END DO
      !
      DO i=1,IMX	    ! read connectivity node ID in each panel
	M1(I)=Mesh%P(1,i)
	M2(I)=Mesh%P(2,i)
	M3(I)=Mesh%P(3,i)
	M4(I)=Mesh%P(4,i)	
      END DO
      !
      ! CALCUL des quantites A PARTIR DU MAILLAGE (VALIDE LE 17/11/14)
      DO 10 I=1,IMX       !IMX=Npanels 
      CALL AIR(X(M1(I)),Y(M1(I)),Z(M1(I)),X(M2(I)),Y(M2(I)),Z(M2(I)),X(M3(I)),Y(M3(I)),Z(M3(I)),X(M4(I)),Y(M4(I)),Z(M4(I)),AIRE(I),TDIS(I),XM(I),YM(I),ZM(I),P(I),Q(I),R(I))
      ! AIRE surface of panels, XM YM ZM coord of centre panels,Aire is area of the panel, TDIS is max projected Tangen/Binormal to the nodes  in a panel, P Q R are unit normal vector in a panel
      10 CONTINUE          ! this actually the closed loop so I will continue up to IMX
      
      ! CALCUL DE LONGUEURS CARACTERISTIQUES DU CORPS DECRIT
      ! Calculation of characteristic length of the body
      GL=0.
      BL=0.
      NP1=NP-1        ! NP is number of points
      DO 15 I=1,NP1
      I1=I+1
        DO 16 J=I1,NP
        GL=AMAX1(ABS(X(J)-X(I)),GL)
        IF(NSYMY.EQ.0)BL=AMAX1(ABS(Y(J)-Y(I)),BL)
        IF(NSYMY.EQ.1)BL=AMAX1(ABS(Y(J)),ABS(Y(I)),BL)
        16 CONTINUE
      15 CONTINUE
!
      IF(NSYMY.EQ.1)BL=2.*BL     !BL is maximum distance in Y axis between nodes                                                
      GL=AMAX1(GL,BL)            !GL is maximum distance in Y or X axis between nodes                                                 
      ZER1=-GL*1.E-5
!      write(*,"(A,E14.6,A,E14.6,A,E14.6)") 'BL=', BL,' GL=', GL, ' ZER1=',ZER1
!     CALCUL DES GRANDEURS A LA FLOTTAISON                 Calculation of floatation quantities
!     XM,YM : COORDONNEES DU POINT DE CONTROLE DU SEGMENT  XM,YM: Segment checkpoint coordinates
!     ALF : LONGUEUR DU SEGMENT                            Segment length
!     IND : NUMERO DE LA FACETTE DU CORPS ADJACENTE AU SEGMENT Number of facet of the body adjacent to the segment
      IXX=IMX                       !IXX starts as Npanels then will increase
      DO 20 I=1,IMX                 !IMX =Npanels
	IF(ZM(I).LT.ZER1)THEN       !ZER1 ~ 0 so this restrict for wetted panels or ZM< ZER1
	  KK0(1)=M1(I)              !ID of first Node in a panel                                                   
	  KK0(2)=M2(I)
	  KK0(3)=M3(I)                                                               
	  KK0(4)=M4(I)                                                               
	  KK0(5)=KK0(1)
	  DO 25 L=1,4
	    DNUL=SQRT((X(KK0(L+1))-X(KK0(L)))**2+(Y(KK0(L+1))-Y(KK0(L)))**2+(Z(KK0(L+1))-Z(KK0(L)))**2)   !distance of node L+1 and node L                                             
	    IF(DNUL.LT.1.E-5)GOTO 25
	    IF(Z(KK0(L+1))+Z(KK0(L)).GE.ZER1)THEN           ! only for larger than ZERR1  (above mean water level)
	      DST=0.01*DNUL
	      IXX=IXX+1
	      XM(IXX)=(X(KK0(L))+X(KK0(L+1)))*0.5+P(I)*DST  ! updated centroid 
	      YM(IXX)=(Y(KK0(L))+Y(KK0(L+1)))*0.5+Q(I)*DST
	      ALF(IXX-IMX)=DNUL                             ! distance
	      IND(IXX-IMX)=I
	      ZM(IXX)=ZER1
	    ENDIF
	  25 CONTINUE
	ENDIF
      20 CONTINUE
      
      NFFL=IXX-IMX                                         ! Number of waterline/controlled panel 
   !  write(*,"(I5)") NFFL     
    END SUBROUTINE ProcessMesh
    
    SUBROUTINE AIR(X1,Y1,Z1,X2,Y2,Z2,X3,Y3, Z3,X4,Y4,Z4,A,S,XMO,YMO,ZMO,P,Q,R)                                                                    !
!     CALCUL DE LA GEOMETRIE DE LA FACETTE (M1,M2,M3,M4)               !
!     variables are declared implicitly, defined as real   (RK) 
      EPPM=1.E-6
      T1X=X3-X1  !DX,DY,DZ of Nodes 1 and 3
      T1Y=Y3-Y1
      T1Z=Z3-Z1
      T2X=X4-X2  !DX,DY,DZ of Nodes 2 and 4
      T2Y=Y4-Y2
      T2Z=Z4-Z2      
      XNX=T1Y*T2Z-T2Y*T1Z   ! Nx, x component of normal vector, obtained from cross product of 2 lines
      XNY=T1Z*T2X-T2Z*T1X   ! Ny
      XNZ=T1X*T2Y-T2X*T1Y   ! Nz
      XNQ=SQRT(XNX*XNX+XNY*XNY+XNZ*XNZ) ! |N|
      IF(XNQ.LT.EPPM)THEN
      A=0.
      RETURN
      ENDIF
      P=XNX/XNQ              ! unit normal vector (P,Q,R)
      Q=XNY/XNQ
      R=XNZ/XNQ
      XAVER=0.25*(X1+X2+X3+X4)   ! X average
      YAVER=0.25*(Y1+Y2+Y3+Y4)   ! Y average
      ZAVER=0.25*(Z1+Z2+Z3+Z4)   ! Z average
      T1=SQRT(T1X**2+T1Y**2+T1Z**2) ! length
      T1UNX=T1X/T1              ! Tx, x component of unit tangent vector
      T1UNY=T1Y/T1              ! Ty
      T1UNZ=T1Z/T1              ! Tz
      T2UNX=Q*T1UNZ-R*T1UNY     ! x component, Binormal vector B=T x N
      T2UNY=R*T1UNX-P*T1UNZ     ! y comp               
      T2UNZ=P*T1UNY-Q*T1UNX     ! z comp
      XI1=T1UNX*(X1-XAVER)+T1UNY*(Y1-YAVER)+T1UNZ*(Z1-ZAVER)  ! dot product of tangen and the diff of the node to the averaged node
      XI2=T1UNX*(X4-XAVER)+T1UNY*(Y4-YAVER)+T1UNZ*(Z4-ZAVER)  !   
      XI3=T1UNX*(X3-XAVER)+T1UNY*(Y3-YAVER)+T1UNZ*(Z3-ZAVER)
      XI4=T1UNX*(X2-XAVER)+T1UNY*(Y2-YAVER)+T1UNZ*(Z2-ZAVER)
      ETA1=T2UNX*(X1-XAVER)+T2UNY*(Y1-YAVER)+T2UNZ*(Z1-ZAVER) ! dot product of binormal vector and the diff ..
      ETA2=T2UNX*(X4-XAVER)+T2UNY*(Y4-YAVER)+T2UNZ*(Z4-ZAVER)
      ETA4=T2UNX*(X2-XAVER)+T2UNY*(Y2-YAVER)+T2UNZ*(Z2-ZAVER)
      XIO=(XI4*(ETA1-ETA2)+XI2*(ETA4-ETA1))/(3.*(ETA2-ETA4))
      ETAO=-ETA1/3.
      ETAN2=ETA2-ETAO                                          !??? ETAN 1, ETAN3 is not defined? then that are zero
      ETAN4=ETA4-ETAO
      XIN1=XI1-XIO
      XIN3=XI3-XIO
      XMO=XAVER+T1UNX*XIO+T2UNX*ETAO                ! XM0, YM0, ZM0  are coordinate of centroid in a panel
      YMO=YAVER+T1UNY*XIO+T2UNY*ETAO
      ZMO=ZAVER+T1UNZ*XIO+T2UNZ*ETAO
      A=ABS(0.5*(XIN3-XIN1)*(ETAN4-ETAN2))          ! area
      TT1=SQRT(XIN1*XIN1+ETAN1*ETAN1)               ! length of projected tangent in node 1 vector 
      TT2=SQRT(XIN2*XIN2+ETAN2*ETAN2)               ! length of projected binormal in node 2        
      TT3=SQRT(XIN3*XIN3+ETAN3*ETAN3)               ! length of projected tangent
      TT4=SQRT(XIN4*XIN4+ETAN4*ETAN4)               ! length of projected binormal                          
      S=MAX(TT1,TT2,TT3,TT4)                        ! maximum length of projected tangent/binormal to the node
      
      RETURN
      
    END SUBROUTINE AIR
    
END MODULE 

