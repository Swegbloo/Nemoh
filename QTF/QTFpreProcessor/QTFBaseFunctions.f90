!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - November 2014
!   Contributors list:
!   - Gerard Delhommeau
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
!--------------------------------------------------------------------------------------
!  routines initially coming from PROGRAMMES AQUADYN : (15-04-85) P.GUEVEL,J-C.DAUBISSE,G.DELHOMMEAU  !        


MODULE MQTFBaseFunctions 

CONTAINS

  SUBROUTINE VAV(KKK,H,ZER)
  
    USE QTFCOM_VAR
    
    IMPLICIT NONE

    INTEGER:: ISP,IFP
    INTEGER:: KKK,I,J,MK,NJJI,JJ,L,MH,MY,MZ,MJJ
    REAL:: H,ZER,DH,XOI,YOI,ZOI,XGI,YGI,ZGI,DLOG                    
    REAL:: RR(5),DRX(5),DRY(5),DRZ(5)                    
    REAL:: XN(5),YN(5),ZN(5),AIJS(4),VXS(4),VYS(4),VZS(4)
    REAL:: A3J,A6J,A9J,ALDEN,ANL,ANLX,ANLY,ANLZ,ANTX,ANTY,ANTZ
    REAL:: ARG,ASRO,AT,ATX,ATY,ATZ,DAT,DDK,DEN,DENL,DENT,DK,DLOGG
    REAL:: ANT,DNL,DNT,DNTX,DNTY,DNTZ,DR,DS,GY,GYX,GYZ,GZ,PJ,QJ,RJ,RO,SGN,W
    REAL:: GYY,XOJ,YOJ,ZOJ
    
    NJJI=2*(NSYMY+1)
    DH=2*H
    MK=(-1)**(KKK+1)
    DO 10 I=1,IXX
    XOI=XM(I)                                                             
    YOI=YM(I)
    ZOI=MIN(ZM(I),2*ZER)
    DO 24 J=1,IMX
    DO 25 JJ=1,NJJI
    MJJ=(-1)**(JJ+1)
    MY=(-1)**(JJ/3+2)
    MZ=(-1)**(JJ/2+2)
    MH=(1-(-1)**(JJ/2+2))/2
    XOJ=XM(J)                                                               
    YOJ=YM(J)*MY
    ZOJ=ZM(J)*MZ-DH*MH                                               
    A3J=P(J)                                                               
    A6J=Q(J)*MY
    A9J=R(J)*MZ                                                            
    RO=SQRT((XOI-XOJ)**2+(YOI-YOJ)**2+(ZOI-ZOJ)**2)
    IF(RO.GT.7.*TDIS(J))THEN
    AIJS(JJ)=AIRE(J)/RO
    ASRO=AIJS(JJ)/RO**2
    VXS(JJ)=-(XOI-XOJ)*ASRO
    VYS(JJ)=-(YOI-YOJ)*ASRO                                              
    VZS(JJ)=-(ZOI-ZOJ)*ASRO
    ELSE
    AIJS(JJ)=0.
    VXS(JJ)=0.                                                                
    VYS(JJ)=0.                                                                
    VZS(JJ)=0.
    KK(1)=M1(J)
    KK(2)=M2(J)                                                               
    KK(3)=M3(J)                                                               
    KK(4)=M4(J)                                                               
    KK(5)=KK(1)                                                               
    DO 211 L=1,4                                                              
    XN(L)=X(KK(L))                                                            
    YN(L)=Y(KK(L))*MY                                                         
    ZN(L)=Z(KK(L))*MZ-DH*MH                                              
    211 CONTINUE                                                                  
    XN(5)=XN(1)                                                               
    YN(5)=YN(1)                                                               
    ZN(5)=ZN(1)                                                               
    DO 212 L=1,4
    RR(L)=SQRT((XOI-XN(L))**2+(YOI-YN(L))**2+(ZOI-ZN(L))**2)
    DRX(L)=(XOI-XN(L))/RR(L)
    DRY(L)=(YOI-YN(L))/RR(L)
    DRZ(L)=(ZOI-ZN(L))/RR(L)
    212 CONTINUE    
    RR(5)=RR(1)
    DRX(5)=DRX(1)                                                           
    DRY(5)=DRY(1)                                                           
    DRZ(5)=DRZ(1)                                                           
    GZ=(XOI-XOJ)*A3J+(YOI-YOJ)*A6J+(ZOI-ZOJ)*A9J
    DO 29 L=1,4                                                               
    DK=SQRT((XN(L+1)-XN(L))**2+(YN(L+1)-YN(L))**2+(ZN(L+1)-ZN(L))**2)
    IF(DK.GE.1.E-3*TDIS(J))THEN                                             
    PJ=(XN(L+1)-XN(L))/DK                                                     
    QJ=(YN(L+1)-YN(L))/DK                                                     
    RJ=(ZN(L+1)-ZN(L))/DK                                                     
    GYX=A6J*RJ-A9J*QJ                                                     
    GYY=A9J*PJ-A3J*RJ                                                     
    GYZ=A3J*QJ-A6J*PJ                                                     
    GY=(XOI-XN(L))*GYX+(YOI-YN(L))*GYY+(ZOI-ZN(L))*GYZ
    SGN=SIGN(1.,GZ)                                                           
    DDK=2.*DK                                                                 
    ANT=GY*DDK                                                                
    DNT=(RR(L+1)+RR(L))**2-DK*DK+2.*ABS(GZ)*(RR(L+1)+RR(L))          
    ARG=ANT/DNT                                                               
    ANL=RR(L+1)+RR(L)+DK                                                      
    DNL=RR(L+1)+RR(L)-DK                                                      
    DEN=ANL/DNL                                                               
    ALDEN=ALOG(DEN)                                                           
    IF(ABS(GZ).GE.1.E-4*TDIS(J))THEN                                       
    AT=ATAN(ARG)
    ELSE
    AT=0.                                                                     
    ENDIF                                                                  
    AIJS(JJ)=AIJS(JJ)+GY*ALDEN-2.*ABS(GZ)*AT                                  
    DAT=2.*AT*SGN
    ANTX=GYX*DDK
    ANTY=GYY*DDK                                                              
    ANTZ=GYZ*DDK                                                              
    ANLX=DRX(L+1)+DRX(L)                                                      
    ANLY=DRY(L+1)+DRY(L)                                                      
    ANLZ=DRZ(L+1)+DRZ(L)                                                      
    DR=2.*(RR(L+1)+RR(L)+ABS(GZ))                                             
    DS=2.*(RR(L+1)+RR(L))*SGN                                                 
    DNTX=DR*ANLX+A3J*DS                                                     
    DNTY=DR*ANLY+A6J*DS                                                     
    DNTZ=DR*ANLZ+A9J*DS                                                     
    DENL=ANL*DNL                                                              
    DENT=ANT*ANT+DNT*DNT                                                      
    ATX=(ANTX*DNT-DNTX*ANT)/DENT                                              
    ATY=(ANTY*DNT-DNTY*ANT)/DENT                                              
    ATZ=(ANTZ*DNT-DNTZ*ANT)/DENT                                              
    DLOG=(DNL-ANL)/DENL                                                       
    VXS(JJ)=VXS(JJ)+GYX*ALDEN+GY*ANLX*DLOG-2.*ABS(GZ)*ATX-DAT*A3J
    VYS(JJ)=VYS(JJ)+GYY*ALDEN+GY*ANLY*DLOG-2.*ABS(GZ)*ATY-DAT*A6J             
    VZS(JJ)=VZS(JJ)+GYZ*ALDEN+GY*ANLZ*DLOG-2.*ABS(GZ)*ATZ-DAT*A9J            
    ENDIF
    29 CONTINUE
    IF(I.EQ.J.AND.JJ.EQ.1)THEN
    !VXS(1)=VXS(1)-DPI*A3J
    !VYS(1)=VYS(1)-DPI*A6J
    !if(zm(i).GT.ZER)then
    !VZS(1)=VZS(1)+DPI*A9J
    !ELSE
    !VZS(1)=VZS(1)-DPI*A9J
    !ENDIF
		  ELSE
    AIJS(JJ)=AIJS(JJ)*MJJ                                                     
    VXS(JJ)=VXS(JJ)*MJJ
    VYS(JJ)=VYS(JJ)*MJJ
    VZS(JJ)=VZS(JJ)*MJJ
    ENDIF
    ENDIF
    25 CONTINUE                                                                  
    IF(NSYMY.EQ.1)THEN                                                       
    W=AIJS(1)-MK*(AIJS(2)+AIJS(3))+AIJS(4)                                   
    FSP(J)=-W/QPI                                                            
    W=AIJS(1)-MK*(AIJS(2)-AIJS(3))-AIJS(4)                                  
    FSM(J)=-W/QPI
    W=VXS(1)-MK*(VXS(2)+VXS(3))+VXS(4)
    VSXP(J)=-W/QPI
    W=VYS(1)-MK*(VYS(2)+VYS(3))+VYS(4)
    VSYP(J)=-W/QPI                                                            
    W=VZS(1)-MK*(VZS(2)+VZS(3))+VZS(4)                              
    VSZP(J)=-W/QPI
    W=VXS(1)-MK*(VXS(2)-VXS(3))-VXS(4)                               
    VSXM(J)=-W/QPI                                                            
    W=VYS(1)-MK*(VYS(2)-VYS(3))-VYS(4)                             
    VSYM(J)=-W/QPI
    W=VZS(1)-MK*(VZS(2)-VZS(3))-VZS(4)                             
    VSZM(J)=-W/QPI                                                            
    ELSE
    W=AIJS(1)-MK*AIJS(2)
    FSP(J)=-W/QPI                                                            
    FSM(J)=FSP(J)                                                           
    W=VXS(1)-MK*VXS(2)
    VSXP(J)=-W/QPI
    VSXM(J)=VSXP(J)                                                           
    W=VYS(1)-MK*VYS(2)                                                     
    VSYP(J)=-W/QPI                                                            
    VSYM(J)=VSYP(J)                                                           
    W=VZS(1)-MK*VZS(2)                                                     
    VSZP(J)=-W/QPI                                                            
    VSZM(J)=VSZP(J)                                                           
    ENDIF
    24 CONTINUE
    WRITE(8)(FSP(J),J=1,IMX),(FSM(J),J=1,IMX)
    WRITE(8)(VSXP(J),J=1,IMX),(VSYP(J),J=1,IMX),(VSZP(J),J=1,IMX)
    WRITE(8)(VSXM(J),J=1,IMX),(VSYM(J),J=1,IMX),(VSZM(J),J=1,IMX)
    10 CONTINUE
    REWIND 8
    RETURN                                                                    
  
  END SUBROUTINE VAV
  
  
  SUBROUTINE VNS(T,BETA,VA,ZER)
  
    USE QTFCOM_VAR
    
    IMPLICIT NONE
    
    COMPLEX ZIJ(5),CEX(5),GZ(5)
    COMPLEX ZI,C1,C2,ASD,BSD,CSD,ZA0,ZB,ZVS

    REAL :: DPI2,EPS,WH,T,WR,BETA,VA,ZER,TR,GM,BJ,DIJ,AK0,AKK,AKP,AKP2,AKZMAX,ZMIII,BX,YMJJJ,QJJJ,RRR,AKR,ZGAJ,ZZZ,AKZ 
    REAL :: DD, PSURR, XL1, XL2, XL3,XL4,XL5,ZL1,ZL2,ZL3,ZL4,ZL5,F1,F2,F3,F4,F5,PD1Z,PD2Z,EPZ,AKP4,SQ,CSK,SIK,VZ1,VZ2 
    REAL :: PD1X,PD2X,DSK,VR1,VR2,CVX,CVY,TETA,CT,ST,AKAIR,ADPI2,ADPI,AKDPI2,AKDPI
    
    INTEGER:: I,J,NPINTE,IT,NJJ,NP1,I1,JJ,IJUMP,LLL,N,KJ,KI,L

    ZI=(0.,1.)
    DPI2=2.*PI**2
    EPS=0.0001
    WH=DPI/T
    WR=WH-VA*WH**2*COS(BETA)/G
    IF(ABS(WR).LT.1.E-4)THEN
    WRITE(LE,*)'ABS(WR)  = ',ABS(WR),' < 1.E-4'
    !      STOP
    ENDIF
    TR=DPI/WR                                                                 
    WRITE(LE,1978)T,WH,VA,TR
    1978 FORMAT(/1X,'PERIODE = ',F12.6,' SECONDES    PULSATION = ',F12.6,' RAD/S'/1X,'VITESSE D''AVANCE = ',F12.6,' M/S  PERIODE DE RENCONTRE = ', F12.6,' S')
    GM=0.
    NJJ=NSYMY+1
    NP1=NP-1
    DO 7100 I=1,NP1
    I1=I+1
    DO 7200 JJ=1,NJJ
    BJ=(-1.)**(JJ+1)
    DO 7200 J=I1,NP
    DIJ=SQRT((X(J)-X(I))**2+(Y(I)-Y(J)*BJ)**2)
    GM=AMAX1(DIJ,GM)
    7200 CONTINUE
    7100 CONTINUE
    AK0=WR**2/G                                                         
    AKK=AK0*GM
    CALL CINT(AKK,N)
    AKP=AK0/DPI
    AKP2=AK0/DPI2
    AKZMAX=-1.E10
!     PRINT *,' VNS '
    IJUMP=0
    DO 9010 I=1,IXX
    READ(8)(FSP(J),J=1,IMX),(FSM(J),J=1,IMX)
    READ(8)(VSXP(J),J=1,IMX),(VSYP(J),J=1,IMX),(VSZP(J),J=1,IMX)
    READ(8)(VSXM(J),J=1,IMX),(VSYM(J),J=1,IMX),(VSZM(J),J=1,IMX)
    ZMIII=MIN(ZM(I),2*ZER)
    DO 22 JJ=1,NJJ
    BX=(-1)**(JJ+1)
    DO 21 J=1,IMX
    FS1(J,JJ)=0.
    FS2(J,JJ)=0.
    VSX1(J,JJ)=0.
    VSX2(J,JJ)=0.
    VSY1(J,JJ)=0.
    VSY2(J,JJ)=0.
    VSZ1(J,JJ)=0.
    VSZ2(J,JJ)=0.
    IF(ZM(I).LT.ZER.OR.ZM(J).LT.ZER)THEN              
    DO 211 LLL=1,NG
    YMJJJ=YGA(LLL,J)*BX
    QJJJ=Q(J)*BX
    RRR=SQRT((XM(I)-XGA(LLL,J))**2+(YM(I)-YMJJJ)**2)
    AKR=AK0*RRR
    ZGAJ=MIN(ZGA(LLL,J),ZER)
    ZZZ=ZMIII+ZGAJ
    AKZ=AK0*ZZZ
    AKZMAX=MAX(AKZMAX,AKZ)
    IF(AKZ.GT.-1.5E-10)THEN
    IF(IJUMP.NE.1)THEN
    WRITE(LE,*)'AKZ > -1.5 E-10'
    IJUMP=1
    ENDIF
    AKZ=-1.6E-10
    ENDIF
    DD=SQRT(RRR**2+ZZZ**2)
    IF(DD.GT.EPS)THEN
    PSURR=PI/(AK0*DD)**3
    ELSE
    PSURR=0.
    ENDIF
    IF(AKZ.LT.-1.5E-10)THEN
    IF(AKZ.GT.-251.)THEN
    IF(AKR.LT.100)THEN
    KJ=10*(ALOG10(-AKZ)+10.)
    KJ=MAX(KJ,3)
    KJ=MIN(KJ,122)
    IF(AKR.LT.1.)THEN
    KI=10*(ALOG10(AKR+1.E-10)+8)+1
    ELSE
    KI=6*AKR+75
    ENDIF
    KI=MAX(KI,3)
    KI=MIN(KI,674)
    XL1=PL5(XR(KI+2),XR(KI-1),XR(KI  ),XR(KI+1),XR(KI-2),AKR)
    XL2=PL5(XR(KI-2),XR(KI  ),XR(KI+1),XR(KI+2),XR(KI-1),AKR)
    XL3=PL5(XR(KI-1),XR(KI+1),XR(KI+2),XR(KI-2),XR(KI  ),AKR)
    XL4=PL5(XR(KI  ),XR(KI+2),XR(KI-2),XR(KI-1),XR(KI+1),AKR)
    XL5=PL5(XR(KI+1),XR(KI-2),XR(KI-1),XR(KI  ),XR(KI+2),AKR)
    ZL1=PL5(XZ(KJ+2),XZ(KJ-1),XZ(KJ  ),XZ(KJ+1),XZ(KJ-2),AKZ)
    ZL2=PL5(XZ(KJ-2),XZ(KJ  ),XZ(KJ+1),XZ(KJ+2),XZ(KJ-1),AKZ)
    ZL3=PL5(XZ(KJ-1),XZ(KJ+1),XZ(KJ+2),XZ(KJ-2),XZ(KJ  ),AKZ)
    ZL4=PL5(XZ(KJ  ),XZ(KJ+2),XZ(KJ-2),XZ(KJ-1),XZ(KJ+1),AKZ)
    ZL5=PL5(XZ(KJ+1),XZ(KJ-2),XZ(KJ-1),XZ(KJ  ),XZ(KJ+2),AKZ)
    F1=XL1*APD1Z(KI-2,KJ-2)+XL2*APD1Z(KI-1,KJ-2)+XL3*APD1Z(KI  ,KJ-2)+ XL4*APD1Z(KI+1,KJ-2)+XL5*APD1Z(KI+2,KJ-2)
    F2=XL1*APD1Z(KI-2,KJ-1)+XL2*APD1Z(KI-1,KJ-1)+XL3*APD1Z(KI  ,KJ-1)+ XL4*APD1Z(KI+1,KJ-1)+XL5*APD1Z(KI+2,KJ-1)
    F3=XL1*APD1Z(KI-2,KJ  )+XL2*APD1Z(KI-1,KJ  )+XL3*APD1Z(KI  ,KJ  )+ XL4*APD1Z(KI+1,KJ  )+XL5*APD1Z(KI+2,KJ  )
    F4=XL1*APD1Z(KI-2,KJ+1)+XL2*APD1Z(KI-1,KJ+1)+XL3*APD1Z(KI  ,KJ+1)+ XL4*APD1Z(KI+1,KJ+1)+XL5*APD1Z(KI+2,KJ+1)
    F5=XL1*APD1Z(KI-2,KJ+2)+XL2*APD1Z(KI-1,KJ+2)+XL3*APD1Z(KI  ,KJ+2)+ XL4*APD1Z(KI+1,KJ+2)+XL5*APD1Z(KI+2,KJ+2)
    PD1Z=ZL1*F1+ZL2*F2+ZL3*F3+ZL4*F4+ZL5*F5
    F1=XL1*APD2Z(KI-2,KJ-2)+XL2*APD2Z(KI-1,KJ-2)+XL3*APD2Z(KI  ,KJ-2)+ XL4*APD2Z(KI+1,KJ-2)+XL5*APD2Z(KI+2,KJ-2)
    F2=XL1*APD2Z(KI-2,KJ-1)+XL2*APD2Z(KI-1,KJ-1)+XL3*APD2Z(KI  ,KJ-1)+ XL4*APD2Z(KI+1,KJ-1)+XL5*APD2Z(KI+2,KJ-1)
    F3=XL1*APD2Z(KI-2,KJ  )+XL2*APD2Z(KI-1,KJ  )+XL3*APD2Z(KI  ,KJ  )+ XL4*APD2Z(KI+1,KJ  )+XL5*APD2Z(KI+2,KJ  )
    F4=XL1*APD2Z(KI-2,KJ+1)+XL2*APD2Z(KI-1,KJ+1)+XL3*APD2Z(KI  ,KJ+1)+ XL4*APD2Z(KI+1,KJ+1)+XL5*APD2Z(KI+2,KJ+1)
    F5=XL1*APD2Z(KI-2,KJ+2)+XL2*APD2Z(KI-1,KJ+2)+XL3*APD2Z(KI  ,KJ+2)+ XL4*APD2Z(KI+1,KJ+2)+XL5*APD2Z(KI+2,KJ+2)
    PD2Z=ZL1*F1+ZL2*F2+ZL3*F3+ZL4*F4+ZL5*F5
    ELSE
    EPZ=EXP(AKZ)
    AKP4=AKR-PI4
    SQ=SQRT(DPI/AKR)
    CSK=COS(AKP4)
    SIK=SIN(AKP4)
    PD1Z=PSURR*AKZ-PI*EPZ*SQ*SIK
    PD2Z=EPZ*SQ*CSK
    ENDIF
    VZ1=PD1Z-PSURR*AKZ
    VZ2=PD2Z
    ELSE
    PD1Z=PSURR*AKZ
    PD2Z=0.
    VZ1=0.
    VZ2=0.
    ENDIF
    ENDIF
    FS1(J,JJ)=FS1(J,JJ)+XJAC(LLL,J)*PD1Z
    FS2(J,JJ)=FS2(J,JJ)+XJAC(LLL,J)*PD2Z
    IF(RRR.GT.EPS)THEN
    IF(AKZ.LE.-1.5E-10)THEN
    IF(AKZ.GT.-251.)THEN
    IF(AKR.LT.100)THEN
    F1=XL1*APD1X(KI-2,KJ-2)+XL2*APD1X(KI-1,KJ-2)+XL3*APD1X(KI  ,KJ-2)+ XL4*APD1X(KI+1,KJ-2)+XL5*APD1X(KI+2,KJ-2)
    F2=XL1*APD1X(KI-2,KJ-1)+XL2*APD1X(KI-1,KJ-1)+XL3*APD1X(KI  ,KJ-1)+ XL4*APD1X(KI+1,KJ-1)+XL5*APD1X(KI+2,KJ-1)
    F3=XL1*APD1X(KI-2,KJ  )+XL2*APD1X(KI-1,KJ  )+XL3*APD1X(KI  ,KJ  )+ XL4*APD1X(KI+1,KJ  )+XL5*APD1X(KI+2,KJ  )
    F4=XL1*APD1X(KI-2,KJ+1)+XL2*APD1X(KI-1,KJ+1)+XL3*APD1X(KI  ,KJ+1)+ XL4*APD1X(KI+1,KJ+1)+XL5*APD1X(KI+2,KJ+1)
    F5=XL1*APD1X(KI-2,KJ+2)+XL2*APD1X(KI-1,KJ+2)+XL3*APD1X(KI  ,KJ+2)+ XL4*APD1X(KI+1,KJ+2)+XL5*APD1X(KI+2,KJ+2)
    PD1X=ZL1*F1+ZL2*F2+ZL3*F3+ZL4*F4+ZL5*F5
    F1=XL1*APD2X(KI-2,KJ-2)+XL2*APD2X(KI-1,KJ-2)+XL3*APD2X(KI  ,KJ-2)+ XL4*APD2X(KI+1,KJ-2)+XL5*APD2X(KI+2,KJ-2)
    F2=XL1*APD2X(KI-2,KJ-1)+XL2*APD2X(KI-1,KJ-1)+XL3*APD2X(KI  ,KJ-1)+ XL4*APD2X(KI+1,KJ-1)+XL5*APD2X(KI+2,KJ-1)
    F3=XL1*APD2X(KI-2,KJ  )+XL2*APD2X(KI-1,KJ  )+XL3*APD2X(KI  ,KJ  )+ XL4*APD2X(KI+1,KJ  )+XL5*APD2X(KI+2,KJ  )
    F4=XL1*APD2X(KI-2,KJ+1)+XL2*APD2X(KI-1,KJ+1)+XL3*APD2X(KI  ,KJ+1)+ XL4*APD2X(KI+1,KJ+1)+XL5*APD2X(KI+2,KJ+1)
    F5=XL1*APD2X(KI-2,KJ+2)+XL2*APD2X(KI-1,KJ+2)+XL3*APD2X(KI  ,KJ+2)+ XL4*APD2X(KI+1,KJ+2)+XL5*APD2X(KI+2,KJ+2)
    PD2X=ZL1*F1+ZL2*F2+ZL3*F3+ZL4*F4+ZL5*F5
    ELSE
    DSK=0.5/AKR
    PD1X=-PSURR*AKR+PI*EPZ*SQ*(CSK-DSK*SIK)
    PD2X=EPZ*SQ*(SIK+DSK*CSK)
    ENDIF
    VR1=-PD1X-PSURR*AKR
    VR2=-PD2X
    ELSE
    PD1X=-PSURR*AKR
    PD2X=0.
    VR1=0.
    VR2=0.
    ENDIF
    ENDIF
    CVX=(XM(I)-XGA(LLL,J))/RRR
    CVY=(YM(I)-YMJJJ)/RRR
    VSX1(J,JJ)=VSX1(J,JJ)+XJAC(LLL,J)*VR1*CVX
    VSX2(J,JJ)=VSX2(J,JJ)+XJAC(LLL,J)*VR2*CVX
    VSY1(J,JJ)=VSY1(J,JJ)+XJAC(LLL,J)*VR1*CVY
    VSY2(J,JJ)=VSY2(J,JJ)+XJAC(LLL,J)*VR2*CVY
    VSZ1(J,JJ)=VSZ1(J,JJ)+XJAC(LLL,J)*VZ1
    VSZ2(J,JJ)=VSZ2(J,JJ)+XJAC(LLL,J)*VZ2
    ELSE
    PD1X=0.
    PD2X=0.
    VSZ1(J,JJ)=VSZ1(J,JJ)+XJAC(LLL,J)*VZ1
    VSZ2(J,JJ)=VSZ2(J,JJ)+XJAC(LLL,J)*VZ2
    ENDIF
    211 CONTINUE
    ELSE
    FS1(J,JJ)=0.
    FS2(J,JJ)=0.
    VSX1(J,JJ)=0.
    VSX2(J,JJ)=0.
    VSY1(J,JJ)=0.
    VSY2(J,JJ)=0.
    VSZ1(J,JJ)=0.
    VSZ2(J,JJ)=0.
    KK(1)=M1(J)
    KK(2)=M2(J)
    KK(3)=M3(J)
    KK(4)=M4(J)
    KK(5)=KK(1)
    DO 30 IT=1,N
    TETA=QQ(IT)
    CT=COS(TETA)
    ST=SIN(TETA)
    DO 20 L=1,4
    ZIJ(L)=AK0*(ZMIII+Z(KK(L))+ZI*((XM(I)-X(KK(L)))*CT+ (YM(I)-Y(KK(L))*BX)*ST))
    IF(REAL(ZIJ(L)).GT.-25.)THEN
    CEX(L)=CEXP(ZIJ(L))
    ELSE
    CEX(L)=(0.,0.)
    ENDIF
    GZ(L)=ZJ(ZIJ(L),CEX(L))
    20 CONTINUE
    ZIJ(5)=ZIJ(1)
    CEX(5)=CEX(1)
    GZ(5)=GZ(1)
    ZA0=(0.,0.)
    ZB=(0.,0.)
    ZVS=(0.,0.)
    DO 37 L=1,4
    C1=ZIJ(L+1)-ZIJ(L)
    IF(ABS(AIMAG(C1)).LT.EPS.AND.ABS(REAL(C1)).LT.EPS)THEN
    ASD=(GZ(L+1)+GZ(L))*0.5
    BSD=(CEX(L+1)+CEX(L))*0.5
    CSD=ASD-(0.5/ZIJ(L+1)+0.5/ZIJ(L))
    ELSE
    CSD=(GZ(L+1)-GZ(L))/C1
    ASD=(GZ(L+1)-GZ(L)+CLOG(ZIJ(L+1)/ZIJ(L)))/C1
    BSD=(CEX(L+1)-CEX(L))/C1
    ENDIF
    C2=(Q(J)*BX-ZI*R(J)*ST)*(X(KK(L+1))-X(KK(L)))- (P(J)-ZI*R(J)*CT)*(Y(KK(L+1))-Y(KK(L)))*BX
    ZA0=ZA0+ASD*C2
    ZB=ZB+BSD*C2
    ZVS=ZVS+CSD*C2
    37 CONTINUE
    FS1(J,JJ)=FS1(J,JJ)+CQ(IT)*REAL(ZA0)
    FS2(J,JJ)=FS2(J,JJ)+CQ(IT)*REAL(ZB)
    VSX1(J,JJ)=VSX1(J,JJ)-CQ(IT)*CT*AIMAG(ZVS)
    VSX2(J,JJ)=VSX2(J,JJ)-CQ(IT)*CT*AIMAG(ZB)
    VSY1(J,JJ)=VSY1(J,JJ)-CQ(IT)*ST*AIMAG(ZVS)
    VSY2(J,JJ)=VSY2(J,JJ)-CQ(IT)*ST*AIMAG(ZB)
    VSZ1(J,JJ)=VSZ1(J,JJ)+CQ(IT)*REAL(ZVS)
    VSZ2(J,JJ)=VSZ2(J,JJ)+CQ(IT)*REAL(ZB)
    30 CONTINUE                 
    AKAIR=AK0*AIRE(J)                 
    FS1(J,JJ)=FS1(J,JJ)*BX/AKAIR
    FS2(J,JJ)=FS2(J,JJ)*BX/AKAIR
    VSX1(J,JJ)=VSX1(J,JJ)*BX/AKAIR
    VSX2(J,JJ)=VSX2(J,JJ)*BX/AKAIR
    VSY1(J,JJ)=VSY1(J,JJ)*BX/AKAIR
    VSY2(J,JJ)=VSY2(J,JJ)*BX/AKAIR
    VSZ1(J,JJ)=VSZ1(J,JJ)*BX/AKAIR
    VSZ2(J,JJ)=VSZ2(J,JJ)*BX/AKAIR
    ENDIF
    21 CONTINUE
    22 CONTINUE
    IF(NSYMY.EQ.1)THEN
    DO 481 J=1,IMX
    AKAIR=AK0*AIRE(J)
    ADPI2=AKAIR/DPI2
    ADPI=AKAIR/DPI
    SM1(J)=FSM(J)-(FS1(J,1)-FS1(J,2))*ADPI2
    SP1(J)=FSP(J)-(FS1(J,1)+FS1(J,2))*ADPI2
    SM2(J)=-(FS2(J,1)-FS2(J,2))*ADPI                                   
    SP2(J)=-(FS2(J,1)+FS2(J,2))*ADPI
    AKDPI2=ADPI2*AK0
    AKDPI=ADPI*AK0
    VSXP1(J)=VSXP(J)-(VSX1(J,1)+VSX1(J,2))*AKDPI2                             
    VSXM1(J)=VSXM(J)-(VSX1(J,1)-VSX1(J,2))*AKDPI2                            
    VSYP1(J)=VSYP(J)-(VSY1(J,1)+VSY1(J,2))*AKDPI2                           
    VSYM1(J)=VSYM(J)-(VSY1(J,1)-VSY1(J,2))*AKDPI2                            
    VSZP1(J)=VSZP(J)-(VSZ1(J,1)+VSZ1(J,2))*AKDPI2                            
    VSZM1(J)=VSZM(J)-(VSZ1(J,1)-VSZ1(J,2))*AKDPI2
    VSXP2(J)=-(VSX2(J,1)+VSX2(J,2))*AKDPI                                    
    VSXM2(J)=-(VSX2(J,1)-VSX2(J,2))*AKDPI                                    
    VSYP2(J)=-(VSY2(J,1)+VSY2(J,2))*AKDPI                                    
    VSYM2(J)=-(VSY2(J,1)-VSY2(J,2))*AKDPI                                    
    VSZP2(J)=-(VSZ2(J,1)+VSZ2(J,2))*AKDPI
    VSZM2(J)=-(VSZ2(J,1)-VSZ2(J,2))*AKDPI                                   
    481 CONTINUE
    ELSE
    DO 491 J=1,IMX
    AKAIR=AK0*AIRE(J)
    ADPI2=AKAIR/DPI2
    ADPI=AKAIR/DPI
    SP1(J)=FSP(J)-FS1(J,1)*ADPI2                                              
    SM1(J)=SP1(J)
    SP2(J)=-FS2(J,1)*ADPI                                                    
    SM2(J)=SP2(J)
    AKDPI2=ADPI2*AK0
    AKDPI=ADPI*AK0
    VSXP1(J)=VSXP(J)-VSX1(J,1)*AKDPI2
    VSXM1(J)=VSXP1(J)                                                        
    VSYP1(J)=VSYP(J)-VSY1(J,1)*AKDPI2                                        
    VSYM1(J)=VSYP1(J)                                                        
    VSZP1(J)=VSZP(J)-VSZ1(J,1)*AKDPI2
    VSZM1(J)=VSZP1(J)                                                        
    VSXP2(J)=-VSX2(J,1)*AKDPI
    VSXM2(J)=VSXP2(J)                                                        
    VSYP2(J)=-VSY2(J,1)*AKDPI                                             
    VSYM2(J)=VSYP2(J)                                                      
    VSZP2(J)=-VSZ2(J,1)*AKDPI                                             
    VSZM2(J)=VSZP2(J)
    491 CONTINUE
    ENDIF
    WRITE(9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX), (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
    WRITE(9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX), (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
    WRITE(9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX), (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)
    WRITE(9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX), (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
    9010 CONTINUE
    !       PRINT *,'AKZMAX = ',AKZMAX
    RETURN
  END SUBROUTINE VNS
  
  FUNCTION PL5(U1,U2,U3,U4,U5,XU)
  
    PL5 =((XU-U1)*(XU-U2)*(XU-U3)*(XU-U4))/ ((U5-U1)*(U5-U2)*(U5-U3)*(U5-U4))
    RETURN
    
  END FUNCTION
  
  SUBROUTINE VNSF(T,BETA,VA,H,ZER)
  
    USE QTFCOM_VAR
    
    IMPLICIT NONE
    
    REAL :: WH,T,AKH,H,AMH,AMHL,WR,VA,BETA,AM0,TR,GM,BJ,DIJ,AKK,EPS,A,ADPI,ADPI2,COE1,COE2,COE3,COE4,ZMIII,ZER,BX,COF1,COF2,COF3,COF4
    INTEGER :: NJJ,NP1,I,I1,JJ,N,NEXP,J,NEXP1,IJUMP,L,KJ1,KI,KJ2
    REAL :: QJJJ, YMJJJ,RRR,AKR,ZGAJ,ZZZ1,AKZ1,DD1,RR1,PSR1,PSURR1,XL1,XL2,XL3,XL4,XL5,ZL11,ZL21,ZL31,ZL41,ZL51
    REAL :: F1,F2,F3,F4,F5,PD1Z1,PD2Z1,EPZ1,AKP4,SQ,CSK,SIK,VZ1,VZ2 
    REAL :: PD1X,PD2X,DSK,VR1,VR2,CVX,CVY,TETA,CT,ST,AKAIR,AKDPI2,AKDPI
    REAL :: PSK,SCK,VZ11,VZ21,PSR2,PD1Z2,VZ12,VZ22,PD1Z3,PSR3,PD2Z3,VZ13,VZ23,PSR4,PD1Z4,PD2Z4,VZ14,VZ24,ZZZ2,AKZ2,DD2,RR2,PSURR2
    REAL :: PD2Z2,ZL12,ZL22,ZL32,ZL42,ZL52,EPZ2,ZZZ3,AKZ3,DD3,RR3,PSURR3,ZL13,ZL23,ZL33,ZL43,ZL53,ZL531,EPZ3,ZZZ4,AKZ4,DD4,RR4,PSURR4
    REAL::  ZL14,ZL24,ZL34,ZL44,ZL54,QTQQ,STSS,PD1X1,PD2X1,SCDS,SSDS,VR21,PD1X2,PD2X2,VR22,PD1X3,PD2X32,VR23,PD1X4,PD2X4,VR24,ZL1,ZL3,PD2X3,C1V3,C2V3
    REAL :: EPZ4,XPG,YPG,ACT,AQT,ZPG1,ZPG2,ZPG3,ZPG4,RO1,FTS1,ASRO1,VXS1,VYS1,VZS1,RO2,FTS2,ASRO2,VXS2,VYS2,VZS2
    REAL :: RO3,FTS3,ASRO3,VXS3,VYS3,VZS3
    REAL :: RO4,FTS4,ASRO4,VXS4,VYS4,VZS4
    REAL :: OM,DXL,DYL,BKL,FFS1,VX1,VY1,FFS2,VX2,VY2,FFS3,VX3,VY3,VZ3,FFS4,VX4,VY4,VZ4
    INTEGER :: KJ3,KJ4,KE,IT,KL
    REAL :: ZP(4)
    COMPLEX :: ZIJ(4,5),CEX(4,5),GZ(4,5),CL(4,5)
    COMPLEX :: S1,ZV1,S2,ZV2,ZV3,ZV4,Z1,Z0,CS,CD,DS,DD,CL1,CL0,G1,G0
    COMPLEX :: CEX1,CEX0,AUX,ZAM,ZI,ZIRS,ZIRC
    REAL :: XF(5),YF(5),ZF(5)
    
!     PL5(U1,U2,U3,U4,U5,XU)=((XU-U1)*(XU-U2)*(XU-U3)*(XU-U4))/ ((U5-U1)*(U5-U2)*(U5-U3)*(U5-U4))
    
    ZI=(0.,1.)                  ! for taking only imaginary part or making the real part zero
    WH=DPI/T                    ! freq
    AKH=WH**2*H/G               ! kH with deep water dispersion                               
    AMH=X0(AKH)                 !                                             
    AMHL=AMH/H
    AKH=AMH*TANH(AMH)                                                         
    WR=WH-VA*AMH/H*COS(BETA)                                             
    AKH=WR**2*H/G
    AMH=X0(AKH)                 !kH
    AKH=AMH*TANH(AMH)           !w^2H/g=kH tanh(kH)   general dispersion relation                                             
    AM0=AMH/H                   !k          wave number   
    PRINT *,' VNSF '
    IF(ABS(WR).LT.1.E-4)THEN
    WRITE(LE,*)'ABS(WR)  = ',ABS(WR),' < 1.E-4'
    STOP     
    ENDIF    
    TR=DPI/WR                                                                 
    WRITE(LE,1978)T,WH,VA,TR                                            
    1978 FORMAT(/1X,'PERIODE = ',F12.6,' SECONDES    PULSATION = ',F12.6,' RAD/S'/1X,'VITESSE D''AVANCE = ',F12.6,' M/S  PERIODE DE RENCONTRE = ', F12.6,' S')
    WRITE(LE,3001)AKH,AMH                                                     
    3001 FORMAT(1X,'K0 = ',1PE13.6,'  M0 = ',1PE13.6,' ADIMENSIONNALISES PAR RAPPORT A H')                                                           
    IF(AMH-AKH-1.E-03)7106,7106,7101                                          
    7106 WRITE(LE,7102)                                                            
    7102 FORMAT(/5X,'PROFONDEUR QUASI-INFINIE'/5X,'LE PROGRAMME EN PROFONDEUR INFINIE SERAIT PLUS ADAPTE')                  
    GOTO 7104
    7101 IF(AKH-0.1)7103,7103,7104                                                 
    7103 WRITE(LE,7105)                                                            
    7105 FORMAT(/5X,'PROFONDEUR TROP FAIBLE POUR LA LONGUEUR D''ONDE')             
    7104 CONTINUE                                                                  
    GM=0.
    NJJ=NSYMY+1                                 ! NJJ=1 for not symmetric, 2 for symmetric  
    NP1=NP-1                                    ! number of panel
    DO 7100 I=1,NP1
    I1=I+1
        DO 7200 JJ=1,NJJ
        BJ=(-1.)**(JJ+1)
        DO 7200 J=I1,NP
        DIJ=SQRT((X(J)-X(I))**2+(Y(I)-Y(J)*BJ)**2)      ! horizontal distance          
        GM=AMAX1(DIJ,GM)
        7200 CONTINUE
    7100 CONTINUE
    AKK=AM0*GM                                  ! kL  
    CALL CINT(AKK,N)
    CALL LISC(AKH,AMH,NEXP)
    EPS=0.0001
    A=(AMH+AKH)**2/(H*(AMH**2-AKH**2+AKH))
    NEXP1=NEXP+1                                                              
    AMBDA(NEXP1)=0.                                                           
    AR(NEXP1)=2. 
    ADPI2=-A/(8.*PI**2)
    ADPI=-A/(8*PI)
    COE1=ADPI2/AM0
    COE2=ADPI/AM0
    COE3=ADPI2
    COE4=ADPI                                                           
    IJUMP=1
    NJJ=NSYMY+1 
    DO 9010 I=1,IXX     ! IXX number of panel?
        READ(8)(FSP(J),J=1,IMX),(FSM(J),J=1,IMX)
        READ(8)(VSXP(J),J=1,IMX),(VSYP(J),J=1,IMX),(VSZP(J),J=1,IMX)
        READ(8)(VSXM(J),J=1,IMX),(VSYM(J),J=1,IMX),(VSZM(J),J=1,IMX)
        ZMIII=MIN(ZM(I),ZER)
        !       print *,'vnsf zer zm ',i,zer,zmIII
        DO 7122 JJ=1,NJJ
                BX=(-1)**(JJ+1)
                DO 7021 J=1,IMX         !index of panel
                        COF1=COE3*AIRE(J)
                        COF2=COE4*AIRE(J)
                        COF3=AM0*COF1
                        COF4=AM0*COF2                                             
                        FS1(J,JJ)=0.
                        FS2(J,JJ)=0.
                        VSX1(J,JJ)=0.
                        VSX2(J,JJ)=0.
                        VSY1(J,JJ)=0.
                        VSY2(J,JJ)=0.
                        VSZ1(J,JJ)=0.
                        VSZ2(J,JJ)=0.
                        IF(ZM(I).LT.ZER.OR.ZM(J).LT.ZER)THEN    ! ZER=-draft/1000 so near still water level, ZM has negative values 
                                                                ! then this is calculation for panels from ZER to the draft
                               DO 211 L=1,NG                    !number of weighting this is global variable with defined as 1 inthe QTFINIT 
                                QJJJ=BX*Q(J)
                                YMJJJ=BX*YGA(L,J)
                                RRR=SQRT((XM(I)-XGA(L,J))**2+(YM(I)-YMJJJ)**2)
                                AKR=AM0*RRR
                                ZGAJ=MIN(ZGA(L,J),ZER)
                                ZZZ1=ZM(I)+ZGAJ
                                AKZ1=AM0*ZZZ1
                                DD1=SQRT(RRR**2+ZZZ1**2)
                                IF(DD1.GT.EPS)THEN
                                        RR1=AM0*DD1
                                        PSR1=PI/RR1
                                        PSURR1=PI/RR1**3
                                ELSE
                                        PSR1=0.
                                        PSURR1=0.
                                ENDIF
                                IF(AKZ1.GT.-1.5E-10)THEN
                                        IF(IJUMP.NE.1)THEN
                                        WRITE(LE,*)'AKZ1 > -1.5 E-10'
                                        IJUMP=1
                                        ENDIF
                                ELSE                            
                                AKZ1=MIN(AKZ1,-1.6E-10) 
                                        IF(AKZ1.GT.-251.)THEN
                                                IF(AKR.LT.100)THEN
                                                        KJ1=10*(ALOG10(-AKZ1)+10)
                                                        KJ1=MAX(KJ1,3)
                                                        KJ1=MIN(KJ1,122)
                                                        IF(AKR.LT.1.)THEN
                                                                KI=10*(ALOG10(AKR+1.E-10)+8)+1
                                                        ELSE
                                                                KI=6*AKR+75
                                                        ENDIF
                                                        KI=MAX(KI,3)
                                                        KI=MIN(KI,674)
                                                        XL1=PL5(XR(KI+2),XR(KI-1),XR(KI  ),XR(KI+1),XR(KI-2),AKR)
                                                        XL2=PL5(XR(KI-2),XR(KI  ),XR(KI+1),XR(KI+2),XR(KI-1),AKR)
                                                        XL3=PL5(XR(KI-1),XR(KI+1),XR(KI+2),XR(KI-2),XR(KI  ),AKR)
                                                        XL4=PL5(XR(KI  ),XR(KI+2),XR(KI-2),XR(KI-1),XR(KI+1),AKR)
                                                        XL5=PL5(XR(KI+1),XR(KI-2),XR(KI-1),XR(KI  ),XR(KI+2),AKR)
                                                        ZL11=PL5(XZ(KJ1+2),XZ(KJ1-1),XZ(KJ1  ),XZ(KJ1+1),XZ(KJ1-2),AKZ1)
                                                        ZL21=PL5(XZ(KJ1-2),XZ(KJ1  ),XZ(KJ1+1),XZ(KJ1+2),XZ(KJ1-1),AKZ1)
                                                        ZL31=PL5(XZ(KJ1-1),XZ(KJ1+1),XZ(KJ1+2),XZ(KJ1-2),XZ(KJ1  ),AKZ1)
                                                        ZL41=PL5(XZ(KJ1  ),XZ(KJ1+2),XZ(KJ1-2),XZ(KJ1-1),XZ(KJ1+1),AKZ1)
                                                        ZL51=PL5(XZ(KJ1+1),XZ(KJ1-2),XZ(KJ1-1),XZ(KJ1  ),XZ(KJ1+2),AKZ1)
                                                        F1=XL1*APD1Z(KI-2,KJ1-2)+XL2*APD1Z(KI-1,KJ1-2)+ XL3*APD1Z(KI  ,KJ1-2)+XL4*APD1Z(KI+1,KJ1-2)+XL5*APD1Z(KI+2,KJ1-2)
                                                        F2=XL1*APD1Z(KI-2,KJ1-1)+XL2*APD1Z(KI-1,KJ1-1)+ XL3*APD1Z(KI  ,KJ1-1)+XL4*APD1Z(KI+1,KJ1-1)+XL5*APD1Z(KI+2,KJ1-1)
                                                        F3=XL1*APD1Z(KI-2,KJ1  )+XL2*APD1Z(KI-1,KJ1  )+ XL3*APD1Z(KI  ,KJ1  )+XL4*APD1Z(KI+1,KJ1  )+XL5*APD1Z(KI+2,KJ1  )
                                                        F4=XL1*APD1Z(KI-2,KJ1+1)+XL2*APD1Z(KI-1,KJ1+1)+ XL3*APD1Z(KI  ,KJ1+1)+XL4*APD1Z(KI+1,KJ1+1)+XL5*APD1Z(KI+2,KJ1+1)
                                                        F5=XL1*APD1Z(KI-2,KJ1+2)+XL2*APD1Z(KI-1,KJ1+2)+ XL3*APD1Z(KI  ,KJ1+2)+XL4*APD1Z(KI+1,KJ1+2)+XL5*APD1Z(KI+2,KJ1+2)
                                                        PD1Z1=ZL11*F1+ZL21*F2+ZL31*F3+ZL41*F4+ZL51*F5
                                                        F1=XL1*APD2Z(KI-2,KJ1-2)+XL2*APD2Z(KI-1,KJ1-2)+ XL3*APD2Z(KI  ,KJ1-2)+XL4*APD2Z(KI+1,KJ1-2)+XL5*APD2Z(KI+2,KJ1-2)
                                                        F2=XL1*APD2Z(KI-2,KJ1-1)+XL2*APD2Z(KI-1,KJ1-1)+ XL3*APD2Z(KI  ,KJ1-1)+XL4*APD2Z(KI+1,KJ1-1)+XL5*APD2Z(KI+2,KJ1-1)
                                                        F3=XL1*APD2Z(KI-2,KJ1  )+XL2*APD2Z(KI-1,KJ1  )+ XL3*APD2Z(KI  ,KJ1  )+XL4*APD2Z(KI+1,KJ1  )+XL5*APD2Z(KI+2,KJ1  )
                                                        F4=XL1*APD2Z(KI-2,KJ1+1)+XL2*APD2Z(KI-1,KJ1+1)+ XL3*APD2Z(KI  ,KJ1+1)+XL4*APD2Z(KI+1,KJ1+1)+XL5*APD2Z(KI+2,KJ1+1)
                                                        F5=XL1*APD2Z(KI-2,KJ1+2)+XL2*APD2Z(KI-1,KJ1+2)+ XL3*APD2Z(KI  ,KJ1+2)+XL4*APD2Z(KI+1,KJ1+2)+XL5*APD2Z(KI+2,KJ1+2)
                                                        PD2Z1=ZL11*F1+ZL21*F2+ZL31*F3+ZL41*F4+ZL51*F5
                                                        ! write(*,'(E15.4,E15.4,E15.4,E15.4,E15.4)') F1,F2,F3,F4,F5
                                                        ! write(*,'(E15.4,E15.4,E15.4,E15.4,E15.4)') ZL11,ZL21,ZL31,ZL41,ZL51
                                                        ! write(*,'(E15.4,E15.4,E15.4,E15.4,E15.4)') XZ(KJ1+2),XZ(KJ1-1),XZ(KJ1+1),XZ(KJ1-2),AKZ1
                                                        ! write(*,'(E15.4,E15.4,E15.4,E15.4,E15.4)') APD1Z(KI-2,KJ1-2),APD1Z(KI-1,KJ1-1),APD1Z(KI,KJ1),APD1Z(KI+1,KJ1+1),APD1Z(KI+2,KJ1+2)
                                                        ! write(*,'(E15.4,E15.4,E15.4,E15.4,E15.4)') APD2Z(KI-2,KJ1-2),APD2Z(KI-1,KJ1-1),APD2Z(KI,KJ1),APD2Z(KI+1,KJ1+1),APD2Z(KI+2,KJ1+2)
                                                        ! write(*,'(E15.4,E15.4)') PD1Z1,PD2Z1
                                                        ! write(*,*) '-----------------------------------'
                                                ELSE
                                                        EPZ1=EXP(AKZ1)
                                                        AKP4=AKR-PI4                                        
                                                        SQ=SQRT(DPI/AKR)
                                                        CSK=COS(AKP4)
                                                        SIK=SIN(AKP4)
                                                        PSK=PI*SQ*SIK
                                                        SCK=SQ*CSK  
                                                        PD1Z1=PSURR1*AKZ1-PSK*EPZ1
                                                        PD2Z1=EPZ1*SCK
                                                ENDIF     
                                                VZ11=PD1Z1-PSURR1*AKZ1
                                                VZ21=PD2Z1
                                        ELSE
                                                PD1Z1=PSURR1*AKZ1
                                                PD2Z1=0.
                                                VZ11=0.
                                                VZ21=0.
                                        ENDIF    
                                ENDIF
                                PSR2=0.
                                PD1Z2=0.
                                PD2Z2=0.
                                VZ12=0.
                                VZ22=0.
                                PSR3=0.
                                PD1Z3=0.
                                PD2Z3=0.
                                VZ13=0.
                                VZ23=0.
                                PSR4=0.
                                PD1Z4=0.
                                PD2Z4=0.
                                VZ14=0.
                                VZ24=0.
                                ZZZ2=ZGAJ-ZMIII-2*H
                                AKZ2=AM0*ZZZ2
                                DD2=SQRT(RRR**2+ZZZ2**2)
                                IF(DD2.GT.EPS)THEN
                                        RR2=AM0*DD2
                                        PSR2=PI/RR2
                                        PSURR2=PI/RR2**3
                                ELSE
                                        PSR2=0.
                                        PSURR2=0.
                                ENDIF
                                IF(AKZ2.GT.-1.5E-10)THEN
                                        IF(IJUMP.NE.1)THEN
                                                WRITE(LE,*)'AKZ2 > -1.5 E-10'
                                                IJUMP=1
                                        ENDIF
                                ELSE
                                        AKZ2=MIN(AKZ2,-1.6E-10)
                                        IF(AKZ2.GT.-251.)THEN
                                                IF(AKR.LT.100)THEN
                                                        KJ2=10*(ALOG10(-AKZ2)+10)
                                                        KJ2=MAX(KJ2,3)
                                                        KJ2=MIN(KJ2,122)
                                                        ZL12=PL5(XZ(KJ2+2),XZ(KJ2-1),XZ(KJ2  ),XZ(KJ2+1),XZ(KJ2-2),AKZ2)
                                                        ZL22=PL5(XZ(KJ2-2),XZ(KJ2  ),XZ(KJ2+1),XZ(KJ2+2),XZ(KJ2-1),AKZ2)
                                                        ZL32=PL5(XZ(KJ2-1),XZ(KJ2+1),XZ(KJ2+2),XZ(KJ2-2),XZ(KJ2  ),AKZ2)
                                                        ZL42=PL5(XZ(KJ2  ),XZ(KJ2+2),XZ(KJ2-2),XZ(KJ2-1),XZ(KJ2+1),AKZ2)
                                                        ZL52=PL5(XZ(KJ2+1),XZ(KJ2-2),XZ(KJ2-1),XZ(KJ2  ),XZ(KJ2+2),AKZ2)
                                                        F1=XL1*APD1Z(KI-2,KJ2-2)+XL2*APD1Z(KI-1,KJ2-2)+ XL3*APD1Z(KI  ,KJ2-2)+XL4*APD1Z(KI+1,KJ2-2)+XL5*APD1Z(KI+2,KJ2-2)
                                                        F2=XL1*APD1Z(KI-2,KJ2-1)+XL2*APD1Z(KI-1,KJ2-1)+ XL3*APD1Z(KI  ,KJ2-1)+XL4*APD1Z(KI+1,KJ2-1)+XL5*APD1Z(KI+2,KJ2-1)
                                                        F3=XL1*APD1Z(KI-2,KJ2  )+XL2*APD1Z(KI-1,KJ2  )+ XL3*APD1Z(KI  ,KJ2  )+XL4*APD1Z(KI+1,KJ2  )+XL5*APD1Z(KI+2,KJ2  )
                                                        F4=XL1*APD1Z(KI-2,KJ2+1)+XL2*APD1Z(KI-1,KJ2+1)+ XL3*APD1Z(KI  ,KJ2+1)+XL4*APD1Z(KI+1,KJ2+1)+XL5*APD1Z(KI+2,KJ2+1)
                                                        F5=XL1*APD1Z(KI-2,KJ2+2)+XL2*APD1Z(KI-1,KJ2+2)+ XL3*APD1Z(KI  ,KJ2+2)+XL4*APD1Z(KI+1,KJ2+2)+XL5*APD1Z(KI+2,KJ2+2)
                                                        PD1Z2=ZL12*F1+ZL22*F2+ZL32*F3+ZL42*F4+ZL52*F5
                                                        F1=XL1*APD2Z(KI-2,KJ2-2)+XL2*APD2Z(KI-1,KJ2-2)+ XL3*APD2Z(KI  ,KJ2-2)+XL4*APD2Z(KI+1,KJ2-2)+XL5*APD2Z(KI+2,KJ2-2)
                                                        F2=XL1*APD2Z(KI-2,KJ2-1)+XL2*APD2Z(KI-1,KJ2-1)+ XL3*APD2Z(KI  ,KJ2-1)+XL4*APD2Z(KI+1,KJ2-1)+XL5*APD2Z(KI+2,KJ2-1)
                                                        F3=XL1*APD2Z(KI-2,KJ2  )+XL2*APD2Z(KI-1,KJ2  )+ XL3*APD2Z(KI  ,KJ2  )+XL4*APD2Z(KI+1,KJ2  )+XL5*APD2Z(KI+2,KJ2  )
                                                        F4=XL1*APD2Z(KI-2,KJ2+1)+XL2*APD2Z(KI-1,KJ2+1)+ XL3*APD2Z(KI  ,KJ2+1)+XL4*APD2Z(KI+1,KJ2+1)+XL5*APD2Z(KI+2,KJ2+1)
                                                        F5=XL1*APD2Z(KI-2,KJ2+2)+XL2*APD2Z(KI-1,KJ2+2)+ XL3*APD2Z(KI  ,KJ2+2)+XL4*APD2Z(KI+1,KJ2+2)+XL5*APD2Z(KI+2,KJ2+2)
                                                        PD2Z2=ZL12*F1+ZL22*F2+ZL32*F3+ZL42*F4+ZL52*F5
                                                ELSE
                                                        EPZ2=EXP(AKZ2)
                                                        PD1Z2=PSURR2*AKZ2-PSK*EPZ2
                                                        PD2Z2=EPZ2*SCK
                                                ENDIF
                                                VZ12=PD1Z2-PSURR2*AKZ2
                                                VZ22=PD2Z2
                                        ELSE
                                                PD1Z2=PSURR2*AKZ2
                                                PD2Z2=0.
                                                VZ12=0.
                                                VZ22=0.
                                        ENDIF
                                ENDIF
                                ZZZ3=ZMIII-ZGAJ-2*H
                                AKZ3=AM0*ZZZ3
                                DD3=SQRT(RRR**2+ZZZ3**2)
                                IF(DD3.GT.EPS)THEN
                                        RR3=AM0*DD3
                                        PSR3=PI/RR3
                                        PSURR3=PI/RR3**3
                                ELSE
                                        PSR3=0.
                                        PSURR3=0.
                                ENDIF
                                IF(AKZ3.GT.-1.5E-10)THEN
                                        IF(IJUMP.NE.1)THEN
                                                WRITE(LE,*)'AKZ3 > -1.5 E-10'
                                                IJUMP=1
                                        ENDIF
                                ELSE
                                        AKZ3=MIN(AKZ3,-1.6E-10)
                                        IF(AKZ3.GT.-251.)THEN
                                                IF(AKR.LT.100)THEN
                                                        KJ3=10*(ALOG10(-AKZ3)+10)
                                                        KJ3=MAX(KJ3,3)
                                                        KJ3=MIN(KJ3,122)
                                                        ZL13=PL5(XZ(KJ3+2),XZ(KJ3-1),XZ(KJ3  ),XZ(KJ3+1),XZ(KJ3-2),AKZ3)
                                                        ZL23=PL5(XZ(KJ3-2),XZ(KJ3  ),XZ(KJ3+1),XZ(KJ3+2),XZ(KJ3-1),AKZ3)
                                                        ZL33=PL5(XZ(KJ3-1),XZ(KJ3+1),XZ(KJ3+2),XZ(KJ3-2),XZ(KJ3  ),AKZ3)
                                                        ZL43=PL5(XZ(KJ3  ),XZ(KJ3+2),XZ(KJ3-2),XZ(KJ3-1),XZ(KJ3+1),AKZ3)
                                                        ZL53=PL5(XZ(KJ3+1),XZ(KJ3-2),XZ(KJ3-1),XZ(KJ3  ),XZ(KJ3+2),AKZ3)
                                                        F1=XL1*APD1Z(KI-2,KJ3-2)+XL2*APD1Z(KI-1,KJ3-2)+ XL3*APD1Z(KI  ,KJ3-2)+XL4*APD1Z(KI+1,KJ3-2)+XL5*APD1Z(KI+2,KJ3-2)
                                                        F2=XL1*APD1Z(KI-2,KJ3-1)+XL2*APD1Z(KI-1,KJ3-1)+ XL3*APD1Z(KI  ,KJ3-1)+XL4*APD1Z(KI+1,KJ3-1)+XL5*APD1Z(KI+2,KJ3-1)
                                                        F3=XL1*APD1Z(KI-2,KJ3  )+XL2*APD1Z(KI-1,KJ3  )+ XL3*APD1Z(KI  ,KJ3  )+XL4*APD1Z(KI+1,KJ3  )+XL5*APD1Z(KI+2,KJ3  )
                                                        F4=XL1*APD1Z(KI-2,KJ3+1)+XL2*APD1Z(KI-1,KJ3+1)+ XL3*APD1Z(KI  ,KJ3+1)+XL4*APD1Z(KI+1,KJ3+1)+XL5*APD1Z(KI+2,KJ3+1)
                                                        F5=XL1*APD1Z(KI-2,KJ3+2)+XL2*APD1Z(KI-1,KJ3+2)+ XL3*APD1Z(KI  ,KJ3+2)+XL4*APD1Z(KI+1,KJ3+2)+XL5*APD1Z(KI+2,KJ3+2)
                                                        PD1Z3=ZL13*F1+ZL23*F2+ZL33*F3+ZL43*F4+ZL53*F5
                                                        F1=XL1*APD2Z(KI-2,KJ3-2)+XL2*APD2Z(KI-1,KJ3-2)+ XL3*APD2Z(KI  ,KJ3-2)+XL4*APD2Z(KI+1,KJ3-2)+XL5*APD2Z(KI+2,KJ3-2)
                                                        F2=XL1*APD2Z(KI-2,KJ3-1)+XL2*APD2Z(KI-1,KJ3-1)+ XL3*APD2Z(KI  ,KJ3-1)+XL4*APD2Z(KI+1,KJ3-1)+XL5*APD2Z(KI+2,KJ3-1)
                                                        F3=XL1*APD2Z(KI-2,KJ3  )+XL2*APD2Z(KI-1,KJ3  )+ XL3*APD2Z(KI  ,KJ3  )+XL4*APD2Z(KI+1,KJ3  )+XL5*APD2Z(KI+2,KJ3  )
                                                        F4=XL1*APD2Z(KI-2,KJ3+1)+XL2*APD2Z(KI-1,KJ3+1)+ XL3*APD2Z(KI  ,KJ3+1)+XL4*APD2Z(KI+1,KJ3+1)+XL5*APD2Z(KI+2,KJ3+1)
                                                        F5=XL1*APD2Z(KI-2,KJ3+2)+XL2*APD2Z(KI-1,KJ3+2)+ XL3*APD2Z(KI  ,KJ3+2)+XL4*APD2Z(KI+1,KJ3+2)+XL5*APD2Z(KI+2,KJ3+2)
                                                        PD2Z3=ZL13*F1+ZL23*F2+ZL33*F3+ZL43*F4+ZL531*F5
                                                ELSE
                                                        EPZ3=EXP(AKZ3)
                                                        PD1Z3=PSURR3*AKZ3-PSK*EPZ3
                                                        PD2Z3=EPZ3*SCK
                                                ENDIF
                                                VZ13=PD1Z3-PSURR3*AKZ3
                                                VZ23=PD2Z3                
                                        ELSE
                                                PD1Z3=PSURR3*AKZ3
                                                PD2Z3=0.
                                                VZ13=0.
                                                VZ23=0.
                                        ENDIF
                                ENDIF
                                ZZZ4=-ZGAJ-ZMIII-4*H
                                AKZ4=AM0*ZZZ4
                                DD4=SQRT(RRR**2+ZZZ4**2)
                                IF(DD4.GT.EPS)THEN
                                        RR4=AM0*DD4
                                        PSR4=PI/RR4
                                        PSURR4=PI/RR4**3
                                ELSE
                                        PSR4=0.
                                        PSURR4=0.
                                ENDIF
                                IF(AKZ4.GT.-1.5E-10)THEN
                                        IF(IJUMP.NE.1)THEN
                                                WRITE(LE,*)'AKZ4 > -1.5 E-10'
                                                IJUMP=1
                                        ENDIF
                                ELSE
                                        AKZ4=MIN(AKZ4,-1.6E-10)
                                        IF(AKZ4.GT.-251.)THEN
                                                 IF(AKR.LT.100)THEN
                                                        KJ4=10*(ALOG10(-AKZ4)+10)
                                                        KJ4=MAX(KJ4,3)
                                                        KJ4=MIN(KJ4,122)
                                                        ZL14=PL5(XZ(KJ4+2),XZ(KJ4-1),XZ(KJ4  ),XZ(KJ4+1),XZ(KJ4-2),AKZ4)
                                                        ZL24=PL5(XZ(KJ4-2),XZ(KJ4  ),XZ(KJ4+1),XZ(KJ4+2),XZ(KJ4-1),AKZ4)
                                                        ZL34=PL5(XZ(KJ4-1),XZ(KJ4+1),XZ(KJ4+2),XZ(KJ4-2),XZ(KJ4  ),AKZ4)
                                                        ZL44=PL5(XZ(KJ4  ),XZ(KJ4+2),XZ(KJ4-2),XZ(KJ4-1),XZ(KJ4+1),AKZ4)
                                                        ZL54=PL5(XZ(KJ4+1),XZ(KJ4-2),XZ(KJ4-1),XZ(KJ4  ),XZ(KJ4+2),AKZ4)
                                                        F1=XL1*APD1Z(KI-2,KJ4-2)+XL2*APD1Z(KI-1,KJ4-2)+ XL3*APD1Z(KI  ,KJ4-2)+XL4*APD1Z(KI+1,KJ4-2)+XL5*APD1Z(KI+2,KJ4-2)
                                                        F2=XL1*APD1Z(KI-2,KJ4-1)+XL2*APD1Z(KI-1,KJ4-1)+ XL3*APD1Z(KI  ,KJ4-1)+XL4*APD1Z(KI+1,KJ4-1)+XL5*APD1Z(KI+2,KJ4-1)
                                                        F3=XL1*APD1Z(KI-2,KJ4  )+XL2*APD1Z(KI-1,KJ4  )+ XL3*APD1Z(KI  ,KJ4  )+XL4*APD1Z(KI+1,KJ4  )+XL5*APD1Z(KI+2,KJ4  )
                                                        F4=XL1*APD1Z(KI-2,KJ4+1)+XL2*APD1Z(KI-1,KJ4+1)+ XL3*APD1Z(KI  ,KJ4+1)+XL4*APD1Z(KI+1,KJ4+1)+XL5*APD1Z(KI+2,KJ4+1)
                                                        F5=XL1*APD1Z(KI-2,KJ4+2)+XL2*APD1Z(KI-1,KJ4+2)+ XL3*APD1Z(KI  ,KJ4+2)+XL4*APD1Z(KI+1,KJ4+2)+XL5*APD1Z(KI+2,KJ4+2)
                                                        PD1Z4=ZL14*F1+ZL24*F2+ZL34*F3+ZL44*F4+ZL54*F5
                                                        F1=XL1*APD2Z(KI-2,KJ4-2)+XL2*APD2Z(KI-1,KJ4-2)+ XL3*APD2Z(KI  ,KJ4-2)+XL4*APD2Z(KI+1,KJ4-2)+XL5*APD2Z(KI+2,KJ4-2)
                                                        F2=XL1*APD2Z(KI-2,KJ4-1)+XL2*APD2Z(KI-1,KJ4-1)+ XL3*APD2Z(KI  ,KJ4-1)+XL4*APD2Z(KI+1,KJ4-1)+XL5*APD2Z(KI+2,KJ4-1)
                                                        F3=XL1*APD2Z(KI-2,KJ4  )+XL2*APD2Z(KI-1,KJ4  )+ XL3*APD2Z(KI  ,KJ4  )+XL4*APD2Z(KI+1,KJ4  )+XL5*APD2Z(KI+2,KJ4  )
                                                        F4=XL1*APD2Z(KI-2,KJ4+1)+XL2*APD2Z(KI-1,KJ4+1)+ XL3*APD2Z(KI  ,KJ4+1)+XL4*APD2Z(KI+1,KJ4+1)+XL5*APD2Z(KI+2,KJ4+1)
                                                        F5=XL1*APD2Z(KI-2,KJ4+2)+XL2*APD2Z(KI-1,KJ4+2)+ XL3*APD2Z(KI  ,KJ4+2)+XL4*APD2Z(KI+1,KJ4+2)+XL5*APD2Z(KI+2,KJ4+2)
                                                        PD2Z4=ZL14*F1+ZL24*F2+ZL34*F3+ZL44*F4+ZL54*F5
                                                ELSE
                                                        EPZ4=EXP(AKZ4)
                                                        PD1Z4=PSURR4*AKZ4-PSK*EPZ4
                                                        PD2Z4=EPZ4*SCK
                                                ENDIF
                                                VZ14=PD1Z4-PSURR4*AKZ4
                                                VZ24=PD2Z4
                                         ELSE
                                                PD1Z4=PSURR4*AKZ4
                                                PD2Z4=0.
                                                VZ14=0.
                                                VZ24=0.
                                         ENDIF
                                ENDIF
                                    QTQQ=PD1Z1+PD1Z2+PD1Z3+PD1Z4
                                    FS1(J,JJ)=FS1(J,JJ)+COF1*(QTQQ-PSR1-PSR2-PSR3-PSR4)*XJAC(L,J)
                                    STSS=PD2Z1+PD2Z2+PD2Z3+PD2Z4
                                    FS2(J,JJ)=FS2(J,JJ)+COF2*STSS*XJAC(L,J)
                                    IF(RRR.GT.EPS)THEN
                                            IF(AKZ1.LE.-1.5E-10)THEN
                                                    IF(AKZ1.GT.-251.)THEN
                                                            IF(AKR.LT.100)THEN
                                                                    F1=XL1*APD1X(KI-2,KJ1-2)+XL2*APD1X(KI-1,KJ1-2)+ XL3*APD1X(KI  ,KJ1-2)+XL4*APD1X(KI+1,KJ1-2)+XL5*APD1X(KI+2,KJ1-2)
                                                                    F2=XL1*APD1X(KI-2,KJ1-1)+XL2*APD1X(KI-1,KJ1-1)+ XL3*APD1X(KI  ,KJ1-1)+XL4*APD1X(KI+1,KJ1-1)+XL5*APD1X(KI+2,KJ1-1)
                                                                    F3=XL1*APD1X(KI-2,KJ1  )+XL2*APD1X(KI-1,KJ1  )+ XL3*APD1X(KI  ,KJ1  )+XL4*APD1X(KI+1,KJ1  )+XL5*APD1X(KI+2,KJ1  )
                                                                    F4=XL1*APD1X(KI-2,KJ1+1)+XL2*APD1X(KI-1,KJ1+1)+ XL3*APD1X(KI  ,KJ1+1)+XL4*APD1X(KI+1,KJ1+1)+XL5*APD1X(KI+2,KJ1+1)
                                                                    F5=XL1*APD1X(KI-2,KJ1+2)+XL2*APD1X(KI-1,KJ1+2)+ XL3*APD1X(KI  ,KJ1+2)+XL4*APD1X(KI+1,KJ1+2)+XL5*APD1X(KI+2,KJ1+2)
                                                                    PD1X1=ZL11*F1+ZL21*F2+ZL31*F3+ZL41*F4+ZL51*F5
                                                                    F1=XL1*APD2X(KI-2,KJ1-2)+XL2*APD2X(KI-1,KJ1-2)+ XL3*APD2X(KI  ,KJ1-2)+XL4*APD2X(KI+1,KJ1-2)+XL5*APD2X(KI+2,KJ1-2)
                                                                    F2=XL1*APD2X(KI-2,KJ1-1)+XL2*APD2X(KI-1,KJ1-1)+ XL3*APD2X(KI  ,KJ1-1)+XL4*APD2X(KI+1,KJ1-1)+XL5*APD2X(KI+2,KJ1-1)
                                                                    F3=XL1*APD2X(KI-2,KJ1  )+XL2*APD2X(KI-1,KJ1  )+ XL3*APD2X(KI  ,KJ1  )+XL4*APD2X(KI+1,KJ1  )+XL5*APD2X(KI+2,KJ1  )
                                                                    F4=XL1*APD2X(KI-2,KJ1+1)+XL2*APD2X(KI-1,KJ1+1)+ XL3*APD2X(KI  ,KJ1+1)+XL4*APD2X(KI+1,KJ1+1)+XL5*APD2X(KI+2,KJ1+1)
                                                                    F5=XL1*APD2X(KI-2,KJ1+2)+XL2*APD2X(KI-1,KJ1+2)+ XL3*APD2X(KI  ,KJ1+2)+XL4*APD2X(KI+1,KJ1+2)+XL5*APD2X(KI+2,KJ1+2)
                                                                    PD2X1=ZL11*F1+ZL21*F2+ZL31*F3+ZL41*F4+ZL51*F5
                                                            ELSE
                                                                    DSK=0.5/AKR
                                                                    SCDS=PI*SQ*(CSK-DSK*SIK)
                                                                    SSDS=SQ*(SIK+DSK*CSK)
                                                                    PD1X1=-PSURR1*AKR+EPZ1*SCDS
                                                                    PD2X1=EPZ1*SSDS
                                                            ENDIF
                                                            VR21=-PD2X1
                                                    ELSE
                                                            PD1X1=-PSURR1*AKR
                                                            PD2X1=0.
                                                            VR21=0.
                                                    ENDIF
                                            ENDIF
                                            PD1X2=0.
                                            PD2X2=0.
                                            VR22=0.
                                            PD1X3=0.
                                            PD2X32=0.
                                            VR23=0.
                                            PD1X4=0.
                                            PD2X4=0.
                                            VR24=0.
                                            IF(AKZ2.LE.-1.5E-10)THEN
                                                    IF(AKZ2.GT.-251.)THEN
                                                            IF(AKR.LT.100)THEN
                                                                    F1=XL1*APD1X(KI-2,KJ2-2)+XL2*APD1X(KI-1,KJ2-2)+ XL3*APD1X(KI  ,KJ2-2)+XL4*APD1X(KI+1,KJ2-2)+XL5*APD1X(KI+2,KJ2-2)
                                                                    F2=XL1*APD1X(KI-2,KJ2-1)+XL2*APD1X(KI-1,KJ2-1)+ XL3*APD1X(KI  ,KJ2-1)+XL4*APD1X(KI+1,KJ2-1)+XL5*APD1X(KI+2,KJ2-1)
                                                                    F3=XL1*APD1X(KI-2,KJ2  )+XL2*APD1X(KI-1,KJ2  )+ XL3*APD1X(KI  ,KJ2  )+XL4*APD1X(KI+1,KJ2  )+XL5*APD1X(KI+2,KJ2  )
                                                                    F4=XL1*APD1X(KI-2,KJ2+1)+XL2*APD1X(KI-1,KJ2+1)+ XL3*APD1X(KI  ,KJ2+1)+XL4*APD1X(KI+1,KJ2+1)+XL5*APD1X(KI+2,KJ2+1)
                                                                    F5=XL1*APD1X(KI-2,KJ2+2)+XL2*APD1X(KI-1,KJ2+2)+ XL3*APD1X(KI  ,KJ2+2)+XL4*APD1X(KI+1,KJ2+2)+XL5*APD1X(KI+2,KJ2+2)
                                                                    PD1X2=ZL12*F1+ZL22*F2+ZL32*F3+ZL42*F4+ZL52*F5
                                                                    F1=XL1*APD2X(KI-2,KJ2-2)+XL2*APD2X(KI-1,KJ2-2)+ XL3*APD2X(KI  ,KJ2-2)+XL4*APD2X(KI+1,KJ2-2)+XL5*APD2X(KI+2,KJ2-2)
                                                                    F2=XL1*APD2X(KI-2,KJ2-1)+XL2*APD2X(KI-1,KJ2-1)+ XL3*APD2X(KI  ,KJ2-1)+XL4*APD2X(KI+1,KJ2-1)+XL5*APD2X(KI+2,KJ2-1)
                                                                    F3=XL1*APD2X(KI-2,KJ2  )+XL2*APD2X(KI-1,KJ2  )+ XL3*APD2X(KI  ,KJ2  )+XL4*APD2X(KI+1,KJ2  )+XL5*APD2X(KI+2,KJ2  )
                                                                    F4=XL1*APD2X(KI-2,KJ2+1)+XL2*APD2X(KI-1,KJ2+1)+ XL3*APD2X(KI  ,KJ2+1)+XL4*APD2X(KI+1,KJ2+1)+XL5*APD2X(KI+2,KJ2+1)
                                                                    F5=XL1*APD2X(KI-2,KJ2+2)+XL2*APD2X(KI-1,KJ2+2)+ XL3*APD2X(KI  ,KJ2+2)+XL4*APD2X(KI+1,KJ2+2)+XL5*APD2X(KI+2,KJ2+2)
                                                                    PD2X2=ZL12*F1+ZL22*F2+ZL32*F3+ZL42*F4+ZL52*F5
                                                            ELSE
                                                                    PD1X2=-PSURR2*AKR+EPZ2*SCDS
                                                                    PD2X2=EPZ2*SSDS
                                                            ENDIF
                                                            VR22=-PD2X2
                                                    ELSE
                                                            PD1X2=-PSURR2*AKR
                                                            PD2X2=0.
                                                            VR22=0.
                                                    ENDIF
                                            ENDIF
                                            IF(AKZ3.LE.-1.5E-10)THEN
                                                    IF(AKZ3.GT.-251.)THEN
                                                            IF(AKR.LT.100)THEN
                                                                    F1=XL1*APD1X(KI-2,KJ3-2)+XL2*APD1X(KI-1,KJ3-2)+ XL3*APD1X(KI  ,KJ3-2)+XL4*APD1X(KI+1,KJ3-2)+XL5*APD1X(KI+2,KJ3-2)
                                                                    F2=XL1*APD1X(KI-2,KJ3-1)+XL2*APD1X(KI-1,KJ3-1)+ XL3*APD1X(KI  ,KJ3-1)+XL4*APD1X(KI+1,KJ3-1)+XL5*APD1X(KI+2,KJ3-1)
                                                                    F3=XL1*APD1X(KI-2,KJ3  )+XL2*APD1X(KI-1,KJ3  )+ XL3*APD1X(KI  ,KJ3  )+XL4*APD1X(KI+1,KJ3  )+XL5*APD1X(KI+2,KJ3  )
                                                                    F4=XL1*APD1X(KI-2,KJ3+1)+XL2*APD1X(KI-1,KJ3+1)+ XL3*APD1X(KI  ,KJ3+1)+XL4*APD1X(KI+1,KJ3+1)+XL5*APD1X(KI+2,KJ3+1)
                                                                    F5=XL1*APD1X(KI-2,KJ3+2)+XL2*APD1X(KI-1,KJ3+2)+ XL3*APD1X(KI  ,KJ3+2)+XL4*APD1X(KI+1,KJ3+2)+XL5*APD1X(KI+2,KJ3+2)
                                                                    PD1X3=ZL1*F1+ZL23*F2+ZL3*F3+ZL43*F4+ZL53*F5
                                                                    F1=XL1*APD2X(KI-2,KJ3-2)+XL2*APD2X(KI-1,KJ3-2)+ XL3*APD2X(KI  ,KJ3-2)+XL4*APD2X(KI+1,KJ3-2)+XL5*APD2X(KI+2,KJ3-2)
                                                                    F2=XL1*APD2X(KI-2,KJ3-1)+XL2*APD2X(KI-1,KJ3-1)+ XL3*APD2X(KI  ,KJ3-1)+XL4*APD2X(KI+1,KJ3-1)+XL5*APD2X(KI+2,KJ3-1)
                                                                    F3=XL1*APD2X(KI-2,KJ3  )+XL2*APD2X(KI-1,KJ3  )+ XL3*APD2X(KI  ,KJ3  )+XL4*APD2X(KI+1,KJ3  )+XL5*APD2X(KI+2,KJ3  )
                                                                    F4=XL1*APD2X(KI-2,KJ3+1)+XL2*APD2X(KI-1,KJ3+1)+ XL3*APD2X(KI  ,KJ3+1)+XL4*APD2X(KI+1,KJ3+1)+XL5*APD2X(KI+2,KJ3+1)
                                                                    F5=XL1*APD2X(KI-2,KJ3+2)+XL2*APD2X(KI-1,KJ3+2)+ XL3*APD2X(KI  ,KJ3+2)+XL4*APD2X(KI+1,KJ3+2)+XL5*APD2X(KI+2,KJ3+2)
                                                                    PD2X3=ZL13*F1+ZL23*F2+ZL33*F3+ZL43*F4+ZL53*F5
                                                            ELSE
                                                                    PD1X3=-PSURR3*AKR+EPZ3*SCDS
                                                                    PD2X3=EPZ3*SSDS
                                                            ENDIF
                                                            VR23=-PD2X3
                                                    ELSE
                                                            PD1X3=-PSURR3*AKR
                                                            PD2X3=0.
                                                            VR23=0.
                                                    ENDIF
                                            ENDIF
                                            IF(AKZ4.LE.-1.5E-10)THEN
                                                    IF(AKZ4.GT.-251.)THEN
                                                            IF(AKR.LT.100)THEN
                                                                    F1=XL1*APD1X(KI-2,KJ4-2)+XL2*APD1X(KI-1,KJ4-2)+ XL3*APD1X(KI  ,KJ4-2)+XL4*APD1X(KI+1,KJ4-2)+XL5*APD1X(KI+2,KJ4-2)
                                                                    F2=XL1*APD1X(KI-2,KJ4-1)+XL2*APD1X(KI-1,KJ4-1)+ XL3*APD1X(KI  ,KJ4-1)+XL4*APD1X(KI+1,KJ4-1)+XL5*APD1X(KI+2,KJ4-1)
                                                                    F3=XL1*APD1X(KI-2,KJ4  )+XL2*APD1X(KI-1,KJ4  )+ XL3*APD1X(KI  ,KJ4  )+XL4*APD1X(KI+1,KJ4  )+XL5*APD1X(KI+2,KJ4  )
                                                                    F4=XL1*APD1X(KI-2,KJ4+1)+XL2*APD1X(KI-1,KJ4+1)+ XL3*APD1X(KI  ,KJ4+1)+XL4*APD1X(KI+1,KJ4+1)+XL5*APD1X(KI+2,KJ4+1)
                                                                    F5=XL1*APD1X(KI-2,KJ4+2)+XL2*APD1X(KI-1,KJ4+2)+ XL3*APD1X(KI  ,KJ4+2)+XL4*APD1X(KI+1,KJ4+2)+XL5*APD1X(KI+2,KJ4+2)
                                                                    PD1X4=ZL14*F1+ZL24*F2+ZL34*F3+ZL44*F4+ZL54*F5
                                                                    F1=XL1*APD2X(KI-2,KJ4-2)+XL2*APD2X(KI-1,KJ4-2)+ XL3*APD2X(KI  ,KJ4-2)+XL4*APD2X(KI+1,KJ4-2)+XL5*APD2X(KI+2,KJ4-2)
                                                                    F2=XL1*APD2X(KI-2,KJ4-1)+XL2*APD2X(KI-1,KJ4-1)+ XL3*APD2X(KI  ,KJ4-1)+XL4*APD2X(KI+1,KJ4-1)+XL5*APD2X(KI+2,KJ4-1)
                                                                    F3=XL1*APD2X(KI-2,KJ4  )+XL2*APD2X(KI-1,KJ4  )+ XL3*APD2X(KI  ,KJ4  )+XL4*APD2X(KI+1,KJ4  )+XL5*APD2X(KI+2,KJ4  )
                                                                    F4=XL1*APD2X(KI-2,KJ4+1)+XL2*APD2X(KI-1,KJ4+1)+ XL3*APD2X(KI  ,KJ4+1)+XL4*APD2X(KI+1,KJ4+1)+XL5*APD2X(KI+2,KJ4+1)
                                                                    F5=XL1*APD2X(KI-2,KJ4+2)+XL2*APD2X(KI-1,KJ4+2)+ XL3*APD2X(KI  ,KJ4+2)+XL4*APD2X(KI+1,KJ4+2)+XL5*APD2X(KI+2,KJ4+2)
                                                                    PD2X4=ZL14*F1+ZL24*F2+ZL34*F3+ZL44*F4+ZL54*F5
                                                            ELSE
                                                                    PD1X4=-PSURR4*AKR+EPZ4*SCDS
                                                                    PD2X4=EPZ4*SSDS
                                                            ENDIF
                                                            VR24=-PD2X4
                                                    ELSE
                                                            PD1X4=-PSURR4*AKR
                                                            PD2X4=0.
                                                            VR24=0.
                                                    ENDIF  
                                            ENDIF
                                            C1V3=-COF3*(PD1X1+PD1X2+PD1X3+PD1X4)
                                            C2V3=COF4*(VR21+VR22+VR23+VR24)
                                            CVX=(XM(I)-XGA(L,J))/RRR
                                            CVY=(YM(I)-YMJJJ)/RRR
                                            VSX1(J,JJ)=VSX1(J,JJ)+C1V3*CVX*XJAC(L,J)
                                            VSX2(J,JJ)=VSX2(J,JJ)+C2V3*CVX*XJAC(L,J)
                                            VSY1(J,JJ)=VSY1(J,JJ)+C1V3*CVY*XJAC(L,J)
                                            VSY2(J,JJ)=VSY2(J,JJ)+C2V3*CVY*XJAC(L,J)
                                    ELSE
                                            VSX1(J,JJ)=0.                                 
                                            VSX2(J,JJ)=0.                                    
                                            VSY1(J,JJ)=0.                                 
                                            VSY2(J,JJ)=0.
                                    ENDIF
                                    VSZ1(J,JJ)=VSZ1(J,JJ)+COF3*(PD1Z1-PD1Z2+PD1Z3-PD1Z4)*XJAC(L,J)
                                    VSZ2(J,JJ)=VSZ2(J,JJ)+COF4*(VZ21-VZ22+VZ23-VZ24)*XJAC(L,J)
                                    XPG=XM(I)-XGA(L,J)
                                    YPG=YM(I)-YMJJJ
                                    ACT=-0.5*AIRE(J)/QPI
                                    DO 7234 KE=1,NEXP1           
                                            AQT=ACT*AR(KE)
                                            ZPG1=ZMIII-2.*H+H*AMBDA(KE)-ZGAJ
                                            ZPG2=-ZMIII-H*AMBDA(KE)-ZGAJ
                                            ZPG3=-ZMIII-4.*H+H*AMBDA(KE)-ZGAJ
                                            ZPG4=ZMIII+2.*H-H*AMBDA(KE)-ZGAJ
                                            RR1=RRR**2+ZPG1**2
                                            RO1=SQRT(RR1)
                                            IF(RO1.GT.EPS)THEN
                                                    FTS1=AQT/RO1                                                     
                                                    ASRO1=FTS1/RR1                                              
                                            ELSE
                                                    FTS1=0.
                                                    ASRO1=0.
                                            ENDIF
                                            VXS1=-XPG*ASRO1
                                            VYS1=-YPG*ASRO1
                                            VZS1=-ZPG1*ASRO1                                                         
                                            RR2=RRR**2+ZPG2**2
                                            RO2=SQRT(RR2)
                                            IF(RO2.GT.EPS)THEN
                                                    FTS2=AQT/RO2                                                     
                                                    ASRO2=FTS2/RR2                                              
                                            ELSE
                                                    FTS2=0.
                                                    ASRO2=0.
                                            ENDIF
                                            VXS2=-XPG*ASRO2
                                            VYS2=-YPG*ASRO2
                                            VZS2=-ZPG2*ASRO2                                                         
                                            RR3=RRR**2+ZPG3**2
                                            RO3=SQRT(RR3)
                                            IF(RO3.GT.EPS)THEN
                                                    FTS3=AQT/RO3                                                     
                                                    ASRO3=FTS3/RR3                                              
                                            ELSE
                                                    FTS3=0.
                                                    ASRO3=0.
                                            ENDIF
                                            VXS3=-XPG*ASRO3
                                            VYS3=-YPG*ASRO3
                                            VZS3=-ZPG3*ASRO3                                                         
                                            RR4=RRR**2+ZPG4**2
                                            RO4=SQRT(RR4)
                                            IF(RO4.GT.EPS)THEN
                                                    FTS4=AQT/RO4                                                     
                                                    ASRO4=FTS4/RR4                                              
                                            ELSE
                                                    FTS4=0.
                                                    ASRO4=0.
                                            ENDIF
                                            VXS4=-XPG*ASRO4
                                            VYS4=-YPG*ASRO4
                                            VZS4=-ZPG4*ASRO4                                                         
                                            FS1(J,JJ)=FS1(J,JJ)+(FTS1+FTS2+FTS3+FTS4)*XJAC(L,J)
                                            VSX1(J,JJ)=VSX1(J,JJ)+(VXS1+VXS2+VXS3+VXS4)*XJAC(L,J)
                                            VSY1(J,JJ)=VSY1(J,JJ)+(VYS1+VYS2+VYS3+VYS4)*XJAC(L,J)
                                            VSZ1(J,JJ)=VSZ1(J,JJ)+(VZS1-VZS2-VZS3+VZS4)*XJAC(L,J)
                                    7234 CONTINUE
                                211 CONTINUE
                        ELSE    !ZM(I)>ZER OR ZM(J)>ZER  ! this is calculation for panels from 0 to ZER (Near still water level)
                                ! for ireg freq removal treatment?
                            FS1(J,JJ)=0.                                        
                            FS2(J,JJ)=0.                                    
                            VSX1(J,JJ)=0.
                            VSX2(J,JJ)=0.                                         
                            VSY1(J,JJ)=0.                                      
                            VSY2(J,JJ)=0.                                         
                            VSZ1(J,JJ)=0.                                       
                            VSZ2(J,JJ)=0.
                            KK(1)=M1(J)
                            KK(2)=M2(J)
                            KK(3)=M3(J)
                            KK(4)=M4(J)
                            KK(5)=KK(1)
                            DO 30 IT=1,N
                                    TETA=QQ(IT)
                                    CT=COS(TETA)
                                    ST=SIN(TETA)
                                    DO 20 L=1,4
                                            OM=(XM(I)-X(KK(L)))*CT+(YM(I)-BX*Y(KK(L)))*ST
                                            ZIJ(1,L)=AM0*(Z(KK(L))+ZMIII+ZI*OM)
                                            ZIJ(2,L)=AM0*(Z(KK(L))-ZMIII-2*H+ZI*OM)
                                            ZIJ(3,L)=AM0*(ZMIII-Z(KK(L))-2*H+ZI*OM)
                                            ZIJ(4,L)=AM0*(-Z(KK(L))-ZMIII-4*H+ZI*OM)
                                            DO 23 KL=1,4
                                                    IF(REAL(ZIJ(KL,L)).GT.-25.)THEN
                                                            CEX(KL,L)=CEXP(ZIJ(KL,L))
                                                    ELSE
                                                            CEX(KL,L)=(0.,0.)
                                                    ENDIF
                                                    GZ(KL,L)=ZJ(ZIJ(KL,L),CEX(KL,L))
                                                    CL(KL,L)=CLOG(-ZIJ(KL,L))
                                            23 CONTINUE
                                    20 CONTINUE
                                    DO 24 KL=1,4
                                            ZIJ(KL,5)=ZIJ(KL,1)
                                            CEX(KL,5)=CEX(KL,1)
                                            GZ(KL,5)=GZ(KL,1)
                                            CL(KL,5)=CL(KL,1)
                                    24 CONTINUE
                                    S1=(0.,0.)
                                    S2=(0.,0.)
                                    ZV1=(0.,0.)
                                    ZV2=(0.,0.)
                                    ZV3=(0.,0.)
                                    ZV4=(0.,0.)
                                    ZIRS=ZI*R(J)*ST
                                    ZIRC=ZI*R(J)*CT
                                    DO 40 L=1,4
                                            DXL=(X(KK(L+1))-X(KK(L)))
                                            DYL=(Y(KK(L+1))-Y(KK(L)))*BX
                                            DO 50 KL=1,4
                                                    BKL=(-1)**(KL+1)
                                                    IF(KL.LT.3)THEN
                                                            AUX=DXL*(Q(J)*BX-ZIRS)-DYL*(P(J)-ZIRC)
                                                    ELSE
                                                            AUX=DXL*(-Q(J)*BX-ZIRS)-DYL*(-P(J)-ZIRC)
                                                    ENDIF
                                                    Z1=ZIJ(KL,L+1)
                                                    Z0=ZIJ(KL,L)
                                                    CL1=CL(KL,L+1)
                                                    CL0=CL(KL,L)
                                                    G1=GZ(KL,L+1)
                                                    G0=GZ(KL,L)
                                                    CEX1=CEX(KL,L+1)
                                                    CEX0=CEX(KL,L)
                                                    ZAM=Z1-Z0
                                                    IF(ABS(AIMAG(ZAM)).LT.EPS.AND.ABS(REAL(ZAM)).LT.EPS)THEN
                                                            S1=S1+AUX*(G1+G0+CL1+CL0)*0.5
                                                            ZV1=ZV1+AUX*BKL*(G1+G0)*0.5
                                                            ZV3=ZV3+AUX*(G1+G0)*0.5
                                                            S2=S2+AUX*(CEX1+CEX0)*0.5
                                                            ZV2=ZV2+AUX*BKL*(CEX1+CEX0)*0.5
                                                            ZV4=ZV4+AUX*(CEX1+CEX0)*0.5
                                                    ELSE
                                                            S1=S1+AUX*(G1-G0+CL1-CL0+Z1*CL1-Z0*CL0-ZAM)/ZAM
                                                            ZV1=ZV1+AUX*BKL*(G1-G0+CL1-CL0)/ZAM
                                                            ZV3=ZV3+AUX*(G1-G0+CL1-CL0)/ZAM
                                                            S2=S2+AUX*(CEX1-CEX0)/ZAM
                                                            ZV2=ZV2+AUX*BKL*(CEX1-CEX0)/ZAM
                                                            ZV4=ZV4+AUX*(CEX1-CEX0)/ZAM
                                                    ENDIF
                                            50 CONTINUE
                                    40 CONTINUE
                                    FS1(J,JJ)=FS1(J,JJ)+CQ(IT)*REAL(S1)*COE1*BX
                                    FS2(J,JJ)=FS2(J,JJ)+CQ(IT)*REAL(S2)*COE2*BX
                                    VSX1(J,JJ)=VSX1(J,JJ)-CQ(IT)*CT*AIMAG(ZV3)*COE3*BX
                                    VSX2(J,JJ)=VSX2(J,JJ)-CQ(IT)*CT*AIMAG(ZV4)*COE4*BX
                                    VSY1(J,JJ)=VSY1(J,JJ)-CQ(IT)*ST*AIMAG(ZV3)*COE3*BX
                                    VSY2(J,JJ)=VSY2(J,JJ)-CQ(IT)*ST*AIMAG(ZV4)*COE4*BX
                                    VSZ1(J,JJ)=VSZ1(J,JJ)+CQ(IT)*REAL(ZV1)*COE3*BX
                                    VSZ2(J,JJ)=VSZ2(J,JJ)+CQ(IT)*REAL(ZV2)*COE4*BX
                            30 CONTINUE
                            ZP(1)=-ZMIII
                            ZP(2)=ZMIII+2*H
                            ZP(3)=ZMIII-2.*H
                            ZP(4)=-ZMIII-4*H
                            DO 48 L=1,5
                                    XF(L)=X(KK(L))
                                    YF(L)=Y(KK(L))
                                    ZF(L)=Z(KK(L))
                            48 CONTINUE
                            DO 51 KE=1,NEXP1
                                    CALL VSD(XF,YF,ZF,JJ,P(J),Q(J),R(J),AIRE(J),TDIS(J), XM(J),YM(J),ZM(J),XM(I),YM(I),ZP(1)-H*AMBDA(KE),FFS1,VX1,VY1,VZ1)
                                    CALL VSD(XF,YF,ZF,JJ,P(J),Q(J),R(J),AIRE(J),TDIS(J), XM(J),YM(J),ZM(J),XM(I),YM(I),ZP(2)-H*AMBDA(KE),FFS2,VX2,VY2,VZ2)
                                    CALL VSD(XF,YF,ZF,JJ,P(J),Q(J),R(J),AIRE(J),TDIS(J), XM(J),YM(J),ZM(J),XM(I),YM(I),ZP(3)+H*AMBDA(KE),FFS3,VX3,VY3,VZ3)
                                    CALL VSD(XF,YF,ZF,JJ,P(J),Q(J),R(J),AIRE(J),TDIS(J), XM(J),YM(J),ZM(J),XM(I),YM(I),ZP(4)+H*AMBDA(KE),FFS4,VX4,VY4,VZ4)
                                    FS1(J,JJ)=FS1(J,JJ)+(FFS1+FFS2+FFS3+FFS4)*AR(KE)
                                    VSX1(J,JJ)=VSX1(J,JJ)+(VX1+VX2+VX3+VX4)*AR(KE)
                                    VSY1(J,JJ)=VSY1(J,JJ)+(VY1+VY2+VY3+VY4)*AR(KE)
                                    VSZ1(J,JJ)=VSZ1(J,JJ)+(-VZ1+VZ2+VZ3-VZ4)*AR(KE)
                            51 CONTINUE
                        ENDIF
                7021 CONTINUE
        7122 CONTINUE
            IF(NSYMY.EQ.1)THEN
                    DO 481 J=1,IMX  
                            SM1(J)=FSM(J)+FS1(J,1)-FS1(J,2)                                           
                            SP1(J)=FSP(J)+FS1(J,1)+FS1(J,2)                                           
                            SM2(J)=FS2(J,1)-FS2(J,2)                                                  
                            SP2(J)=FS2(J,1)+FS2(J,2)                                                  
                            VSXP1(J)=VSXP(J)+VSX1(J,1)+VSX1(J,2)
                            VSXM1(J)=VSXM(J)+VSX1(J,1)-VSX1(J,2)
                            VSYP1(J)=VSYP(J)+VSY1(J,1)+VSY1(J,2)                                      
                            VSYM1(J)=VSYM(J)+VSY1(J,1)-VSY1(J,2)                                      
                            VSZP1(J)=VSZP(J)+VSZ1(J,1)+VSZ1(J,2)                                      
                            VSZM1(J)=VSZM(J)+VSZ1(J,1)-VSZ1(J,2)                                      
                            VSXP2(J)=VSX2(J,1)+VSX2(J,2)                                              
                            VSXM2(J)=VSX2(J,1)-VSX2(J,2)                                              
                            VSYP2(J)=VSY2(J,1)+VSY2(J,2)                                              
                            VSYM2(J)=VSY2(J,1)-VSY2(J,2)                                              
                            VSZP2(J)=VSZ2(J,1)+VSZ2(J,2)                                              
                            VSZM2(J)=VSZ2(J,1)-VSZ2(J,2)          
                           ! write(*,'(E15.4,E15.4)') FSM(J),FSP(J)    
                    481 CONTINUE
            ELSE
                    DO 491 J=1,IMX
                            SP1(J)=FSP(J)+FS1(J,1)                                                    
                            SM1(J)=SP1(J)                                                             
                            SP2( J)=FS2(J,1)                                                          
                            SM2(J)=SP2(J)                                                             
                            VSXP1(J)=VSXP(J)+VSX1(J,1)
                            VSXM1(J)=VSXP1(J)
                            VSYP1(J)=VSYP(J)+VSY1(J,1)                                                
                            VSYM1(J)=VSYP1(J)                                                         
                            VSZP1(J)=VSZP(J)+VSZ1(J,1)                                                
                            VSZM1(J)=VSZP1(J)                                                         
                            VSXP2(J)=VSX2(J,1)                                                        
                            VSXM2(J)=VSXP2(J)                                                         
                            VSYP2(J)=VSY2(J,1)                                                        
                            VSYM2(J)=VSYP2(J)                                                         
                            VSZP2(J)=VSZ2(J,1)                                                        
                            VSZM2(J)=VSZP2(J)                                                         
                    491 CONTINUE
            ENDIF
            !10 CONTINUE                                           !closed by RK no pair/ or no effect???                       
            WRITE(9)(SP1(J),J=1,IMX),(SM1(J),J=1,IMX), (SP2(J),J=1,IMX),(SM2(J),J=1,IMX)
            WRITE(9)(VSXP1(J),J=1,IMX),(VSXM1(J),J=1,IMX), (VSYP1(J),J=1,IMX),(VSYM1(J),J=1,IMX)
            WRITE(9)(VSZP1(J),J=1,IMX),(VSZM1(J),J=1,IMX), (VSXP2(J),J=1,IMX),(VSXM2(J),J=1,IMX)  
            WRITE(9)(VSYP2(J),J=1,IMX),(VSYM2(J),J=1,IMX), (VSZP2(J),J=1,IMX),(VSZM2(J),J=1,IMX)
          !  WRITE(99,'(<IMX>E14.4)')(SP1(J),J=1,IMX)
    9010 CONTINUE
    RETURN
  END SUBROUTINE VNSF
  
  
  SUBROUTINE VSD(XFT,YFT,ZFT,JJ,P0,Q0,R0,A,T,XG,YG,ZG,XP,YP,ZP,FS,VX,VY,VZ)
  
    USE QTFCOM_VAR
    IMPLICIT NONE 
    
    INTEGER::JJ
    REAL::XG,YG,ZG,XP,YP,ZP,FS,VX,VY,VZ,DLOG
    REAL:: A,T,P0,Q0,R0,XFT(5),YFT(5),ZFT(5)
    INTEGER::MJJ,L
    REAL:: RR(5),DRX(5),DRY(5),DRZ(5)
    REAL:: AIJS,ALDEN,ANL,ANLX,ANLY,ANLZ,ANT,ANTY,ANTX,ANTZ,ARG,ASRO
    REAL:: AT,ATX,ATY,ATZ,DAT,DDK,DEN,DENL,DENT,DK,DLG
    REAL::DNL,DNT,DNTX,DNTY,DNTZ,DR,DS,EPS,GY,GYX,GYZ,GYY,GZ,PJ,QJ 
    REAL::QMP,RJ,RO,SGN,VXS,VZS,VYS

    MJJ=(-1)**(JJ+1)
    QMP=MJJ/(2*QPI)
    EPS=1.E-4
    RO=SQRT((XP-XG)**2+(YP-YG*MJJ)**2+(ZP-ZG)**2)
    GZ=(XP-XG)*P0+(YP-YG*MJJ)*Q0*MJJ+(ZP-ZG)*R0
    IF(RO.GT.7*T)THEN
    FS=-A/(RO*2*QPI)
    ASRO=FS/RO**2
    VX=-(XP-XG)*ASRO
    VY=-(YP-YG*MJJ)*ASRO                                              
    VZ=-(ZP-ZG)*ASRO                                              
    ELSE
    DO 212 L=1,4
    RR(L)=SQRT((XP-XFT(L))**2+(YP-YFT(L)*MJJ)**2+(ZP-ZFT(L))**2)
    DRX(L)=(XP-XFT(L))/RR(L)                                                
    DRY(L)=(YP-YFT(L)*MJJ)/RR(L)                                                
    DRZ(L)=(ZP-ZFT(L))/RR(L)                                                
    212 CONTINUE
    RR(5)=RR(1)
    DRX(5)=DRX(1)                                                           
    DRY(5)=DRY(1)                                                           
    DRZ(5)=DRZ(1)                                                           
    AIJS=0.
    VXS=0.
    VYS=0.
    VZS=0.
    DO 29 L=1,4
    DK=SQRT((XFT(L+1)-XFT(L))**2+(YFT(L+1)-YFT(L))**2+(ZFT(L+1)-ZFT(L))**2)
    IF(DK.GE.1.E-3*T)THEN                                             
    PJ=(XFT(L+1)-XFT(L))/DK                                                     
    QJ=(YFT(L+1)-YFT(L))*MJJ/DK                                                     
    RJ=(ZFT(L+1)-ZFT(L))/DK                                                     
    GYX=Q0*MJJ*RJ-R0*QJ                                                     
    GYY=R0*PJ-P0*RJ                                                     
    GYZ=P0*QJ-Q0*MJJ*PJ                                                     
    GY=(XP-XFT(L))*GYX+(YP-YFT(L)*MJJ)*GYY+(ZP-ZFT(L))*GYZ
    SGN=SIGN(1.,GZ)                                                           
    DDK=2.*DK                                                                 
    ANT=GY*DDK                                                                
    DNT=(RR(L+1)+RR(L))**2-DK*DK+2.*ABS(GZ)*(RR(L+1)+RR(L))
    ANL=RR(L+1)+RR(L)+DK                                                      
    DNL=RR(L+1)+RR(L)-DK                                                      
    DEN=ANL/DNL                                                               
    ALDEN=ALOG(DEN)                                                           
    IF(ABS(GZ).GT.1.E-4*T)THEN
    ARG=ANT/DNT
    AT=ATAN(ARG)
    ELSE
    AT=0.                                                                     
    ENDIF                                                                  
    AIJS=AIJS+GY*ALDEN-2.*ABS(GZ)*AT
    DAT=2.*AT*SGN                                                             
    ANTX=GYX*DDK                                                              
    ANTY=GYY*DDK
    ANTZ=GYZ*DDK                                                              
    ANLX=DRX(L+1)+DRX(L)                                                      
    ANLY=DRY(L+1)+DRY(L)                                                      
    ANLZ=DRZ(L+1)+DRZ(L)                                                      
    DR=2.*(RR(L+1)+RR(L)+ABS(GZ))                                             
    DS=2.*(RR(L+1)+RR(L))*SGN                                                 
    DNTX=DR*ANLX+P0*DS                                                     
    DNTY=DR*ANLY+Q0*MJJ*DS                                                     
    DNTZ=DR*ANLZ+R0*DS                                                     
    DENL=ANL*DNL                                                              
    DENT=ANT*ANT+DNT*DNT                                                      
    ATX=(ANTX*DNT-DNTX*ANT)/DENT                                              
    ATY=(ANTY*DNT-DNTY*ANT)/DENT                                              
    ATZ=(ANTZ*DNT-DNTZ*ANT)/DENT                                              
    DLOG=(DNL-ANL)/DENL
    VXS=VXS+GYX*ALDEN+GY*ANLX*DLOG-2.*ABS(GZ)*ATX-DAT*P0
    VYS=VYS+GYY*ALDEN+GY*ANLY*DLOG-2.*ABS(GZ)*ATY-DAT*Q0*MJJ             
    VZS=VZS+GYZ*ALDEN+GY*ANLZ*DLOG-2.*ABS(GZ)*ATZ-DAT*R0   
    ENDIF        
    29 CONTINUE
    FS=-AIJS*QMP
    VX=-VXS*QMP
    VY=-VYS*QMP
    VZ=-VZS*QMP
    ENDIF
    RETURN
  END SUBROUTINE VSD         
    
  SUBROUTINE CINT(AKK,N)
    
    USE QTFCOM_VAR
    IMPLICIT NONE 
    
    REAL :: Q8(8),CQ8(8),Q12(12),CQ12(12),Q16(16),CQ16(16)
    REAL :: Q24(24),CQ24(24),Q32(32),CQ32(32)
    
    REAL :: AKK
    INTEGER :: N,I,J
    
    
    DATA Q8/.4801449,.3983332,.2627662,.09171732,4*0./
    DATA CQ8/.05061427,.1111905,.1568533,.1813419,4*0./
    DATA Q12/.4907803,.4520586,.3849513,.2936589,.1839157,.06261670,6*0./
    DATA CQ12/.2358766E-1,.5346966E-1,.8003916E-1,.1015837,.1167462,.1245735,6*0./
    DATA Q16/.4947004,.4722875,.4328156,.3777022,.3089381,.2290084,.1408017,.04750625,8*0./
    DATA CQ16/.01357622,.03112676,.04757925,.06231448,.07479799,.08457826,.09130170,.09472530,8*0./
    DATA Q24/.4975936,.4873642,.469137,.4432077,.4100009,.3700621,.3240468,.2727107,.2168967,.1575213,.09555943,.032028446,12*0./
    DATA CQ24/.6170615E-2,.1426569E-1,.2213872E-1,.2964929E-1,.366732E-1,.4309508E-1,.4880932E-1,.5372213E-1,.5775283E-1,.6083523E-1,.6291873E-1,.6396909E-1,12*0./
    DATA Q32/.4986319,.4928057,.4823811,.4674530,.4481605,.4246838,.3972418,.3660910,.3315221,.2938578,.2534499,.2106756,.1659343,.1196436,.07223598,.02415383,16*0./
    DATA CQ32/.350930E-2,.8137197E-2,.1269603E-1,.1713693E-1,.2141794E-1,.2549903E-1,.2934204E-1,.3291111E-1,.3617289E-1,.3909694E-1,.4165596E-1,.4382604E-1,.4558693E-1,.4692219E-1,.4781936E-1,.4827004E-1,16*0./
    IF(AKK-.4)801,801,810
    810 IF(AKK-2.5)802,802,811
    811 IF(AKK-4.)803,803,812
    812 IF(AKK-8.)804,804,813
    813 IF(AKK-25.)805,805,806
    801 N=8
    DO 700 I=1,4
    Q8(I)=Q8(I)
    Q8(9-I)=-Q8(I)
    CQ8(I)=CQ8(I)
    CQ8(9-I)=CQ8(I)
    700 CONTINUE
    DO 110 J=1,N
    QQ(J)=Q8(J)*PI
    CQ(J)=CQ8(J)*PI
    110 CONTINUE
    GOTO 815
    802 N=12
    DO 701 I=1,6
    Q12(I)=Q12(I)
    Q12(13-I)=-Q12(I)
    CQ12(I)=CQ12(I)
    CQ12(13-I)=CQ12(I)
    701 CONTINUE
    DO 120 J=1,N
    QQ(J)=Q12(J)*PI
    CQ(J)=CQ12(J)*PI
    120 CONTINUE
    GOTO 815
    803 N=16
    DO 702 I=1,8
    Q16(I)=Q16(I)
    Q16(17-I)=-Q16(I)
    CQ16(I)=CQ16(I)
    CQ16(17-I)=CQ16(I)
    702 CONTINUE
    DO 130 J=1,N
    QQ(J)=Q16(J)*PI
    CQ(J)=CQ16(J)*PI
    130 CONTINUE
    GOTO 815
    804 N=24
    DO 703 I=1,12
    Q24(I)=Q24(I)
    Q24(25-I)=-Q24(I)
    CQ24(I)=CQ24(I)
    CQ24(25-I)=CQ24(I)
    703 CONTINUE
    DO 140 J=1,N
    QQ(J)=Q24(J)*PI
    CQ(J)=CQ24(J)*PI
    140 CONTINUE
    GOTO 815
    805 N=32
    DO 704 I=1,16
    Q32(I)=Q32(I)
    Q32(33-I)=-Q32(I)
    CQ32(I)=CQ32(I)
    CQ32(33-I)=CQ32(I)
    704 CONTINUE
    DO 150 J=1,N
    QQ(J)=Q32(J)*PI
    CQ(J)=CQ32(J)*PI
    150 CONTINUE
    815 CONTINUE
    WRITE(*,3001)N,AKK
    3001 FORMAT( 2X,'INTEGRATION PAR QG',I2/2X,'K0 ADIMENSIONNEL PAR RAPPORT LA PLUS GRANDE DIMENSION = ',1PE16.6)
    GOTO 807
    806 CONTINUE
    N=51
    IF(AKK.GT.40.)THEN
    WRITE(*,822)
    822 FORMAT(10X,'L''INTEGRATION RISQUE D''ETRE PEU PRECISE')
    N=NPIN
    ENDIF
    DO 160 J=1,N
    QQ(J)=-PI/2.+(J-1.)/(N-1.)*PI
    IF(J-1)161,161,162
    161 CQ(J)=PI/(3.*(N-1.))
    GOTO 160
    162 IF(J-N)163,161,161
    163 IF(MOD(J,2))164,165,164
    164 CQ(J)=2./(3.*(N-1.))*PI
    GOTO 160
    165 CQ(J)=4./(3.*(N-1.))*PI
    160 CONTINUE
    WRITE(*,3002)N,AKK
    3002 FORMAT(2X,'SIMPSON ',I3,' POINTS'/2X,'K0 ADIMENSIONNEL PAR RAPPORT A LA PLUS GRANDE DIMENSION = ',1PE16.6)
    807 CONTINUE
    RETURN
  END SUBROUTINE CINT
  
  
  FUNCTION CH(AK,Z,H)
  IF(AK*H.LE.20)THEN
  CH=COSH(AK*(Z+H))/COSH(AK*H)
  ELSE
  CH=EXP(AK*Z)
  ENDIF
  RETURN
  END FUNCTION
  
  FUNCTION SH(AK,Z,H)
  IF(AK*H.LE.20)THEN
  SH=SINH(AK*(Z+H))/COSH(AK*H)
  ELSE
  SH=EXP(AK*Z)
  ENDIF
  RETURN
  END FUNCTION
  
SUBROUTINE LISC(AK0,AM0,NM)
  
    USE QTFCOM_VAR
    IMPLICIT NONE
    
    
    INTEGER::NM,I,J,NJ,NPP
    REAL:: AK0,AM0                                                              
    REAL:: POL(NEXR),A,B,H                               
    REAL:: S(4*(NEXR-1),NEXR+1),X0(4*(NEXR-1)+1),Y0(4*(NEXR-1)+1)            
    REAL:: SC(NEXR),VR(NEXR),VC(NEXR)
    INTEGER::ISTIR,NMAX,NK,ISOR,NPI,NMO
    REAL:: PRECI,ERMAX,ERMOY,XX,YY,TT,DIF,RT
    COMPLEX:: COM(NEXR)
    REAL :: AMBDAR,T
    
    
    WRITE(LE,7100)                                                            
    7100 FORMAT(/5X,'CARACTERISTIQUES DU LISSAGE'/)                                
    IF(AK0.LT.0.05)THEN
    NM=1
    AMBDAR=-EXP(3.0737-14.269*AM0+132.665*AM0**2-630.831*AM0**3+ 1630.47*AM0**4-1794.54*AM0**5)
    AR(NM)=FF(0.,AK0,AM0)
    AMBDA(NM)=AMBDAR
    GOTO 61
    ENDIF
    PRECI=1.E-02
    ISTIR=0                                                                   
    NMAX=4*(NEXR-1)                                                           
    NK=4                                                                      
    A=-0.1
    B=20.                                                                     
    62 CONTINUE                                                                  
    NM=NK
    NJ=4*NM                                                                   
    7000 CONTINUE                                                                  
    NP=NJ+1                                                                   
    H=(B-A)/NJ                                                                
    DO 10 I=1,NP                                                              
    X0(I)=A+(I-1)*H                                                            
    Y0(I)=FF(X0(I),AK0,AM0)                                                     
    10 CONTINUE                                                                  
    ISOR=0                                                                    
    CALL EXPORS(X0,Y0,NJ,NM,AMBDA,NMAX,S,SC,VR,VC,COM,POL,AR)                
    NPI=2                                                                     
    NMO=NPI*NP-NPI+1                                                          
    ERMAX=0.                                                                  
    ERMOY=0.                                                                  
    DO 20 I=1,NMO                                                             
    XX=(I-1)*B/(NMO-1)                                                        
    YY=FF(XX,AK0,AM0)                                                         
    T=0.                                                                      
    DO 30 J=1,NM                                                              
    RT=AMBDA(J)*XX                                                            
    IF(RT.GT.-20.)THEN
    T=T+AR(J)*EXP(RT)                                                         
    ENDIF
    30 CONTINUE                                                                  
    DIF=YY-T                                                                  
    ERMOY=ERMOY+DIF                                                           
    ERMAX=AMAX1(ERMAX,ABS(DIF))                                               
    IF(ABS(DIF).GT.PRECI)ISOR=1                                               
    20 CONTINUE                                                                  
    ERMOY=ERMOY/NMO                                                           
    WRITE(LE,1111)NM,ERMAX,ERMOY                                              
    1111 FORMAT(5X,I2,' EXPONENTIELLES  ECART MAXI = ',E10.3,'  ECART MOYEN= ',E10.3/)                                                              
    IF(ISTIR.EQ.1)GOTO 61                                                     
    IF(ISOR)63,61,63                                                          
    63 CONTINUE                                                                  
    NK=NK+2                                                                   
    IF(NK-(NEXR-1))62,62,65
    65 WRITE(LE,6500)PRECI,NM                                                    
    6500 FORMAT(/5X,'PRECISION = ',E10.3,'  NON ATTEINTE AVEC ',I2,'  EXPONENTIELLES')                                                               
    61 DO 60 J=1,NM
    WRITE(LE,1100)AR(J),AMBDA(J)                                              
    1100 FORMAT(5X,E16.7,'EXP(',E16.7,')')                                         
    IF(AMBDA(J).GT.0.)STOP                                                    
    60 CONTINUE                                                                  
    RETURN                                                                    
  END SUBROUTINE LISC
  
  
  FUNCTION FF(X,AK,AM)                                                      
    
    FG(T)=(T+AK)*EXP(T)/(T*SINH(T)-AK*COSH(T))-COEF/(T-AM)-2                  
    COEF=(AM+AK)**2/(AM**2-AK**2+AK)                                          
    H=X-AM                                                                    
    TOL=AMAX1(0.05,0.1*AM)
    IF(ABS(H).GT.TOL)THEN                                                  
    FF=FG(X)                                                                  
    ELSE                                                                    
    A=AM-TOL                                                                  
    B=AM                                                                      
    C=AM+TOL                                                                  
    D=FG(A)                                                                   
    E=COEF/(AM+AK)*(AM+AK+1)-(COEF/(AM+AK))**2*AM-2                           
    F=FG(C)                                                                   
    FF=(X-B)*(X-C)*D/((A-B)*(A-C))+(X-C)*(X-A)*E/((B-C)*(B-A))+(X-A)*(X-B)*F/((C-A)*(C-B))                                               
    ENDIF
    RETURN                                                                     
  END FUNCTION
  
  
SUBROUTINE EXPORS(XT,YT,NJ,NM,VCOM,NMAX,S,SC,VR,VC,COM,POL,AR2) 

      IMPLICIT NONE
      
      INTEGER::NJ,NM,NMAX
      REAL:: VCOM(31),POL(31),AR2(31),SC(31),VR(31),VC(31)
      REAL:: S(4*(31-1),31+1),XT(4*(31-1)+1),YT(4*(31-1)+1)           
      COMPLEX:: COM(31) 
      INTEGER::I,J,K,NPP,JJ,II,IJ,MN,NEXP
      INTEGER::IS,IER
      REAL::H,EPS                                                       
                                                     
      NPP=NJ+1                                                                   
      H=(XT(NPP)-XT(1))/NJ                                                         
      K=NPP-NM                                                                  
      DO 2 I=1,K                                                                
      DO 1 J=1,NM                                                               
      JJ=NM-J+I                                                     
 1    S(I,J)=YT(JJ)                                                              
      II=NM+I                                                         
 2    S(I,NM+1)=-YT(II)                                              
      EPS=1.E-20
                                                                
      CALL HOUSRS(S,NMAX,K,NM,1,EPS)          
      DO 5 I=1,NM                                                               
      IJ=NM-I+1                                                                 
 5    SC(IJ)=S(I,NM+1)                                                          
      MN=NM+1                                                                   
      SC(MN)=1.                                                                 
      CALL SPRBM(SC,MN,VR,VC,POL,IS,IER)                              
      DO 6 I=1,NM                                                               
      COM(I)=CMPLX(VR(I),VC(I))                                             
      COM(I)=CLOG(COM(I))/H                                                     
      VR(I)=REAL(COM(I))                                                  
      VC(I)=AIMAG(COM(I))                                  
    6 CONTINUE                                                                  
      I=1                                                                       
      J=0                                                                       
  100 IF(VC(I))110,111,110                                                      
  111 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+1                                                                     
      GO TO 101                                                                 
  110 IF(ABS(VR(I)-VR(I+1))-1.E-5)120,120,121                                   
  120 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+2                                                                     
      GO TO 101                                                                 
  121 J=J+1                                                                     
      VCOM(J)=VR(I)                                                             
      I=I+1                                                                     
  101 IF(I-NM)100,100,102                                                       
  102 NEXP=J                                                                    
      J=0                                                                       
      DO 300 I=1,NEXP                                                           
      J=J+1                                                                     
      IF(VCOM(I).GE.0.)GOTO 301                                                 
      IF(VCOM(I)+20.)301,301,302                                                
  301 J=J-1                                                                     
      GO TO 300                                                                 
  302 VCOM(J)=VCOM(I)                                                           
  300 CONTINUE                                                                  
      NEXP=J                                                                    
      NM=NEXP                                                                 
      CALL MCAS(NEXP,VCOM,XT,YT,NPP,AR2,S,NMAX)                                 
      RETURN                                                                    
      END SUBROUTINE
  
  
  SUBROUTINE HOUSRS(A,NMAX,NL,NC,NS,EPS)
    DIMENSION A(NMAX,*)
    NTC=NC+NS                                                                 
    IF(NC.GT.NL)THEN                                                          
    WRITE(*,3010)                                                            
    3010 FORMAT(' NBRE DE COLONNES > NBRES DE LIGNES')                             
    STOP                                                             
    ENDIF
    DO 13 K=1,NC                                                              
    E=0                                                                       
    DO 1101 I=K,NL                                                            
    E=E+A(I,K)**2                                                             
    1101 CONTINUE                                                                  
    E0=SQRT(E)                                                                
    IF(E0.LT.EPS)THEN
    WRITE(*,201)EPS
    201 FORMAT(1X,'NORME INFERIEURE A ',1PE16.6/)
    STOP                                                                    
    ENDIF
    IF(A(K,K).EQ.0)THEN                                                      
    AR=-E0                                                                    
    ELSE                                                                   
    AR=-SIGN(E0,A(K,K))
    ENDIF
    ETA=AR*(AR-A(K,K))                                                        
    KP1=K+1                                                                   
    DO 10 J=KP1,NTC                                                           
    BA=(A(K,K)-AR)*A(K,J)                                                     
    DO 9 I=KP1,NL                                                            
    BA=BA+A(I,K)*A(I,J)   
    9 CONTINUE                                                   
    A(K,J)=A(K,J)+BA/AR
    DO 11 I=KP1,NL                                                            
    A(I,J)=A(I,J)-A(I,K)*BA/ETA                                               
    11 CONTINUE
    10 CONTINUE                                                                  
    A(K,K)=AR                                                                 
    DO 12 I=KP1,NL                                                            
    12   A(I,K)=0                                                                  
    13   CONTINUE                                                                  
    DO 1006 J=1,NS                                                            
    NCJ=NC+J                                                                  
    A(NC,NCJ)=A(NC,NCJ)/A(NC,NC)                                              
    DO 1005 L=2,NC                                                            
    I1=NC+1-L                                                                 
    M=I1+1                                                                    
    DO 1004 I=M,NC                                                            
    1004 A(I1,NCJ)=A(I1,NCJ)-A(I1,I)*A(I,NCJ)                                      
    1005 A(I1,NCJ)=A(I1,NCJ)/A(I1,I1)                                              
    1006 CONTINUE                                                                  
    RETURN
  END SUBROUTINE HOUSRS
  
  
  SUBROUTINE MCAS(NEXP,TEXP,X,Y,NP,AR,A,NMAX)                            
    
    DIMENSION TEXP(*),Y(*),X(*),AR(*),A(NMAX,*)                         
    EPS=1.E-20                                                                
    DO 1 I=1,NEXP                                                             
    DO 1 J=1,NEXP                                                             
    S=0                                                                       
    DO 3 L=1,NP                                                               
    TT=(TEXP(I)+TEXP(J))*X(L)                                                 
    IF(TT+30)3,4,4                                                            
    4 S=S+EXP(TT)                                                               
    3 CONTINUE                                                                  
    A(I,J)=S                                                                  
    1 CONTINUE                                                                  
    DO 5 I=1,NEXP                                                             
    S=0                                                                       
    DO 6 L=1,NP                                                               
    TTT=TEXP(I)*X(L)                                                          
    IF(TTT+30)6,7,7                                                           
    7 S=S+EXP(TTT)*Y(L)                                                         
    6 CONTINUE                                                                  
    A(I,NEXP+1)=S                                                             
    5 CONTINUE                                                                  
    N=NEXP                                                                    
    M=N+1                                                                     
    CALL HOUSRS(A,NMAX,N,N,1,EPS)                                       
    DO 10 I=1,NEXP                                                            
    10 AR(I)=A(I,NEXP+1)                                                         
    RETURN                                                                    
  END SUBROUTINE MCAS 
  
  
SUBROUTINE SPRBM(C,IC,RR,RC,POL,IR,IER)                                   
    
    USE QTFCOM_VAR
    IMPLICIT NONE 
    
    INTEGER::IC,IR,IER
    REAL:: C(31),RR(31),RC(31),POL(31)                               
    INTEGER::I,J,L,N,LIM,IST
    REAL::A,B,H
    REAL::EPS,Q1,Q2,QQQ(4)

    EPS=1.E-6                                                                 
    LIM=100                                                                   
    IR=IC+1                                                                   
    1 IR=IR-1                                                                   
    IF(IR-1)42,42,2                                                           
    2 IF(C(IR))3,1,3                                                            
    3 IER=0                                                                     
    J=IR                                                                      
    L=0                                                                       
    A=C(IR)                                                                   
    DO 8 I=1,IR                                                               
    IF(L)4,4,7                                                                
    4 IF(C(I))6,5,6                                                             
    5 RR(I)=0.                                                                  
    RC(I)=0.                                                                  
    POL(J)=0.                                                                 
    J=J-1                                                                     
    GO TO 8                                                                   
    6 L=1                                                                       
    IST=I                                                                     
    J=0                                                                       
    7 J=J+1                                                                     
    C(I)=C(I)/A                                                               
    POL(J)=C(I)                                                               
    IF(ABS(POL(J))-1.E27)8,42,42                                              
    8 CONTINUE                                                                  
    Q1=0.                                                                     
    Q2=0.                                                                     
    9 IF(J-2)33,10,14                                                           
    10 A=POL(1)                                                                  
    RR(IST)=-A                                                                
    RC(IST)=0.                                                                
    IR=IR-1                                                                   
    Q2=0.                                                                     
    IF(IR-1)13,13,11                                                          
    11 DO 12 I=2,IR                                                              
    Q1=Q2                                                                     
    Q2=POL(I+1)                                                               
    12 POL(I)=A*Q2+Q1                                                            
    13 POL(IR+1)=A+Q2                                                            
    GO TO 34                                                                  
    14 DO 22 L=1,10                                                              
    N=1                                                                       
    15 QQQ(1)=Q1                                                                   
    QQQ(2)=Q2                                                                   
    CALL SPQFB(POL,J,QQQ,LIM,I)                                                 
    IF(I)16,24,23                                                             
    16 IF(Q1)18,17,18                                                            
    17 IF(Q2)18,21,18                                                            
    18 GOTO(19,20,19,21),N                                                       
    19 Q1=-Q1                                                                    
    N=N+1                                                                     
    GO TO 15                                                                  
    20 Q2=-Q2                                                                    
    N=N+1                                                                     
    GO TO 15                                                                  
    21 Q1=1.+Q1                                                                  
    22 Q2=1.-Q2                                                                  
    IER=3                                                                     
    IR=IR-J                                                                   
    GOTO 45                                                                   
    23 IER=1                                                                     
    24 Q1=QQQ(1)                                                                   
    Q2=QQQ(2)                                                                   
    B=0.                                                                      
    A=0.                                                                      
    I=J                                                                       
    25 H=-Q1*B-Q2*A+POL(I)                                                       
    POL(I)=B                                                                  
    B=A                                                                       
    A=H                                                                       
    I=I-1                                                                     
    IF(I-2)26,26,25                                                           
    26 POL(2)=B                                                                  
    POL(1)=A                                                                  
    L=IR-1                                                                    
    IF(J-L)27,27,29                                                           
    27 DO 28 I=J,L                                                               
    28 POL(I-1)=POL(I-1)+POL(I)*Q2+POL(I+1)*Q1                                   
    29 POL(L)=POL(L)+POL(L+1)*Q2+Q1                                              
    POL(IR)=POL(IR)+Q2                                                        
    H=-.5*Q2                                                                  
    A=H*H-Q1                                                                  
    B=SQRT(ABS(A))                                                            
    IF(A)30,30,31                                                             
    30 RR(IST)=H                                                                 
    RC(IST)=B                                                                 
    IST=IST+1                                                                 
    RR(IST)=H                                                                 
    RC(IST)=-B                                                                
    GO TO 32                                                                  
    31 B=H+SIGN(B,H)                                                             
    RR(IST)=Q1/B                                                              
    RC(IST)=0.                                                                
    IST=IST+1                                                                 
    RR(IST)=B                                                                 
    RC(IST)=0.                                                                
    32 IST=IST+1                                                                 
    J=J-2                                                                     
    GO TO 9                                                                   
    33 IR=IR-1                                                                   
    34 A=0.                                                                      
    DO 38 I=1,IR                                                              
    Q1=C(I)                                                                   
    Q2=POL(I+1)                                                               
    POL(I)=Q2                                                                 
    IF(Q1)35,36,35                                                            
    35 Q2=(Q1-Q2)/Q1                                                             
    36 Q2=ABS(Q2)                                                                
    IF(Q2-A)38,38,37                                                          
    37 A=Q2                                                                      
    38 CONTINUE
    I=IR+1                                                                    
    POL(I)=1.                                                                 
    RR(I)=A                                                                   
    RC(I)=0.                                                                  
    IF(IER)39,39,41                                                           
    39 IF(A-EPS)41,41,40
    40 IER=-1                                                                    
    41 GOTO 45                                                                   
    42 IER=2                                                                     
    IR=0
    45 IF(IER-2)46,47,46                                                         
    47 WRITE(LE,48)IER                                                           
    48 FORMAT(/5X,'IER = ',I3,'  ERREUR DANS SPRBM'/)                            
    STOP                                                                      
    46 RETURN                                                                    
  END SUBROUTINE SPRBM
  
  
  SUBROUTINE SPQFB(C,IC,Q,LIM,IER)                                          
    
    DIMENSION C(*),Q(*)                                                       
    IER=0                                                                     
    J=IC+1                                                                    
    1 J=J-1                                                                     
    IF(J-1) 40,40,2                                                           
    2 IF(C(J)) 3,1,3                                                            
    3 A=C(J)                                                                    
    IF(A-1.) 4,6,4                                                            
    4 DO 5 I=1,J                                                                
    C(I)=C(I)/A                                                               
    IF(ABS(C(I))-1.E27)5,40,40                                                
    5 CONTINUE                                                                  
    6 IF(J-3) 41,38,7                                                           
    7 EPS=1.E-14                                                                
    EPS1=1.E-6                                                                
    L=0                                                                       
    LL=0                                                                      
    Q1=Q(1)                                                                   
    Q2=Q(2)                                                                   
    QQ1=0.                                                                    
    QQ2=0.                                                                    
    AA=C(1)                                                                   
    BB=C(2)                                                                   
    CB=ABS(AA)                                                                
    CA=ABS(BB)                                                                
    IF(CB-CA) 8,9,10                                                          
    8 CC=CB+CB                                                                  
    CB=CB/CA                                                                  
    CA=1.                                                                     
    GO TO 11                                                                  
    9 CC=CA+CA                                                                  
    CA=1.                                                                     
    CB=1.                                                                     
    GO TO 11                                                                  
    10 CC=CA+CA                                                                  
    CA=CA/CB                                                                  
    CB=1.                                                                     
    11 CD=CC*.1                                                                  
    12 A=0.                                                                      
    B=A                                                                       
    A1=A                                                                      
    B1=A                                                                      
    I=J                                                                       
    QQQ1=Q1                                                                   
    QQQ2=Q2                                                                   
    DQ1=HH                                                                    
    DQ2=H                                                                     
    13 H=-Q1*B-Q2*A+C(I)                                                         
    IF(ABS(H)-1.E27)14,42,42                                                  
    14 B=A                                                                       
    A=H                                                                       
    I=I-1                                                                     
    IF(I-1) 18,15,16                                                          
    15 H=0.                                                                      
    16 H=-Q1*B1-Q2*A1+H                                                          
    IF(ABS(H)-1.E27)17,42,42                                                  
    17 C1=B1                                                                     
    B1=A1                                                                     
    A1=H                                                                      
    GO TO 13                                                                  
    18 H=CA*ABS(A)+CB*ABS(B)                                                     
    IF(LL) 19,19,39                                                           
    19 L=L+1                                                                     
    IF(ABS(A)-EPS*ABS(C(1))) 20,20,21                                         
    20 IF(ABS(B)-EPS*ABS(C(2))) 39,39,21                                         
    21 IF(H-CC) 22,22,23                                                         
    22 AA=A                                                                      
    BB=B                                                                      
    CC=H                                                                      
    QQ1=Q1                                                                    
    QQ2=Q2                                                                    
    23 IF(L-LIM) 28,28,24                                                        
    24 IF(H-CD) 43,43,25                                                         
    25 IF(Q(1)) 27,26,27                                                         
    26 IF(Q(2)) 27,42,27                                                         
    27 Q(1)=0.                                                                   
    Q(2)=0.                                                                   
    GO TO 7                                                                   
    28 HH=AMAX1(ABS(A1),ABS(B1),ABS(C1))                                         
    IF(HH) 42,42,29                                                           
    29 A1=A1/HH                                                                  
    B1=B1/HH                                                                  
    C1=C1/HH                                                                  
    H=A1*C1-B1*B1                                                             
    IF(H) 30,42,30                                                            
    30 A=A/HH                                                                    
    B=B/HH                                                                    
    HH=(B*A1-A*B1)/H                                                          
    H=(A*C1-B*B1)/H                                                           
    Q1=Q1+HH                                                                  
    Q2=Q2+H                                                                   
    IF(ABS(HH)-EPS*ABS(Q1)) 31,31,33                                          
    31 IF(ABS(H)-EPS*ABS(Q2)) 32,32,33                                           
    32 LL=1                                                                      
    GO TO 12
    33 IF(L-1)12,12,34                                                           
    34 IF(ABS(HH)-EPS1*ABS(Q1)) 35,35,12                                         
    35 IF(ABS(H)-EPS1*ABS(Q2)) 36,36,12                                          
    36 IF(ABS(QQQ1*HH)-ABS(Q1*DQ1)) 37,44,44                                     
    37 IF(ABS(QQQ2*H)-ABS(Q2*DQ2)) 12,44,44
    38 Q(1)=C(1)                                                                 
    Q(2)=C(2)                                                                 
    Q(3)=0.                                                                   
    Q(4)=0.                                                                   
    GOTO 45                                                                   
    39 Q(1)=Q1                                                                   
    Q(2)=Q2                                                                   
    Q(3)=A                                                                    
    Q(4)=B                                                                    
    GOTO 45                                                                   
    40 IER=-1                                                                    
    GOTO 45                                                                   
    41 IER=-2                                                                    
    GOTO 45                                                                   
    42 IER=-3                                                                    
    GO TO 44                                                                  
    43 IER=1                                                                     
    44 Q(1)=QQ1                                                                  
    Q(2)=QQ2                                                                  
    Q(3)=AA                                                                   
    Q(4)=BB                                                                   
    45 RETURN                                                                    
    END SUBROUTINE SPQFB
    
  FUNCTION X0(AK)
    ! fonction de recherche de la racine de F(x)=w**2*h/g-kh*th(kh) (relation de dispersion)
    F(T)=AK-T*TANH(T)                                                         
    EPS=5.E-6                                                                 
    ITOUR=0                                                                   
    XI=0.                                                                     
    XS=XI                                                                     
    PAS=AMAX1(AK,SQRT(AK))                                                    
    30 XS=XS+PAS                                                                 
    ITOUR=ITOUR+1                                                             
    IF(ITOUR.LE.1000)THEN                                                 
	IF(F(XS)*F(XI).GT.0)GOTO 30                                    
    ENDIF
    IITER=0                                                                   
    10 CONTINUE
    XM=(XI+XS)*0.5                                                            
    IITER=IITER+1                                                             
    IF(IITER.GT.1000.OR.ITOUR.GT.1000)THEN
	WRITE(*,110)ITOUR,IITER                                                  
	110 FORMAT(2X,'ERREUR DANS LA RECHERCHE DE LA RACINE',/2X,'APPROXIMATION =',I5,'   DICHOTOMIE = ',I5)                         
	STOP                                                                      
    ELSE
	IF(ABS((XS-XI)/XM).GT.EPS)THEN
	    IF(F(XM)*F(XI).LT.0)THEN                                                 
		XS=XM                                                                     
	    ELSE                                                                      
		XI=XM                                                                     
	    ENDIF                                                        
	    GOTO 10
	ELSE                                         
	    X0=XM
	ENDIF
    ENDIF
    RETURN                                                                    
  END FUNCTION
  
  
  
  SUBROUTINE SOMGO(X,Y,Z,M1,M2,M3,M4,J,XJAC,XGA,YGA,ZGA,NG)
  
    DIMENSION X(*),Y(*),Z(*),M1(*),M2(*),M3(*),M4(*)
    DIMENSION XJAC(16,*),XGA(16,*),YGA(16,*),ZGA(16,*)
    !
    !     PREPARATION DU CALCUL DES COEFFICIENTS D'INFLUENCE:
    !     DETERMINATION DE LA POSITION DES POINTS DE GAUSS
    !---------------------GAUSS 1 4 9 OU 16--------------------------------
    !     ET DU JACOBIEN POUR CHAQUE FACETTE
    !      NG =1,4,9 OU 16
    !
    DIMENSION XX(16,4),YY(16,4),WG(16,4)
    !
    DATA((XX(I,L),I=1,16),L=1,4)/16*0.,.57735027, .57735027,-.57735027,-.57735027,12*0.,.77459667, .77459667, .77459667,3*0.,-.77459667,-.77459667,-.77459667,7*0.,.86113631, .86113631, .86113631, .86113631,.33998104, .33998104, .33998104, .33998104,-.33998104,-.33998104,-.33998104,-.33998104,-.86113631,-.86113631,-.86113631,-.86113631/
    DATA((YY(I,L),I=1,16),L=1,4)/16*0.,-.57735027,.57735027,-.57735027,.57735027,12*0.,-.77459667,0.,.77459667,-.77459667,0.,.77459667,-.77459667,0.,.77459667,7*0.,-.86113363,-.33998104,.33998104,.86113363,-.86113363,-.33998104,.33998104,.86113363,-.86113363,-.33998104,.33998104,.86113363,-.86113363,-.33998104,.33998104,.86113363/
    DATA((WG(I,L),I=1,16),L=1,4)/1,15*0.,.25,.25,.25,.25,12*0.,.07716049,.12345679,.07716049,.12345679,.19753086,.12345679,.07716049,.12345679,.07716049,7*0.,.30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1,.56712963E-1,.10632333,.10632333,.56712963E-1,.56712963E-1,.10632333,.10632333,.56712963E-1,.30250748E-1,.56712963E-1,.56712963E-1,.30250748E-1/
    !
    IF(NG.EQ.1)THEN
    IMETH=1
    ELSE
    IF(NG.EQ.4)THEN
    IMETH=2
    ELSE
    IF(NG.EQ.9) THEN
    IMETH=3
    ELSE
    IF(NG.EQ.16) THEN
    IMETH=4
    ELSE
    WRITE(*,*)'ARGUMENT INCORRECT DANS SOMGO'
    STOP
    ENDIF
    ENDIF
    ENDIF
    ENDIF
    !
    T1X=X(M3(J))-X(M1(J))
    T1Y=Y(M3(J))-Y(M1(J))
    T1Z=Z(M3(J))-Z(M1(J))
    T2X=X(M4(J))-X(M2(J))
    T2Y=Y(M4(J))-Y(M2(J))
    T2Z=Z(M4(J))-Z(M2(J))
    XNX=T2Y*T1Z-T1Y*T2Z
    XNY=T1X*T2Z-T2X*T1Z
    XNZ=T2X*T1Y-T1X*T2Y
    XNQ=SQRT(XNX**2+XNY**2+XNZ**2)
    AT3=XNX/XNQ
    AT6=XNY/XNQ
    AT9=XNZ/XNQ
    TT=SQRT(T1X**2+T1Y**2+T1Z**2)
    AT1=T1X/TT
    AT4=T1Y/TT
    AT7=T1Z/TT
    AT2=AT6*AT7-AT9*AT4
    AT5=AT9*AT1-AT3*AT7
    AT8=AT3*AT4-AT6*AT1
    XM=0.25*(X(M1(J))+X(M2(J))+X(M3(J))+X(M4(J)))
    YM=0.25*(Y(M1(J))+Y(M2(J))+Y(M3(J))+Y(M4(J)))
    ZM=0.25*(Z(M1(J))+Z(M2(J))+Z(M3(J))+Z(M4(J)))
    XL1=AT1*(X(M1(J))-XM)+AT4*(Y(M1(J))-YM)+AT7*(Z(M1(J))-ZM)
    YL1=AT2*(X(M1(J))-XM)+AT5*(Y(M1(J))-YM)+AT8*(Z(M1(J))-ZM)
    XL2=AT1*(X(M2(J))-XM)+AT4*(Y(M2(J))-YM)+AT7*(Z(M2(J))-ZM)
    YL2=AT2*(X(M2(J))-XM)+AT5*(Y(M2(J))-YM)+AT8*(Z(M2(J))-ZM)
    XL3=AT1*(X(M3(J))-XM)+AT4*(Y(M3(J))-YM)+AT7*(Z(M3(J))-ZM)
    YL3=AT2*(X(M3(J))-XM)+AT5*(Y(M3(J))-YM)+AT8*(Z(M3(J))-ZM)
    XL4=AT1*(X(M4(J))-XM)+AT4*(Y(M4(J))-YM)+AT7*(Z(M4(J))-ZM)
    YL4=AT2*(X(M4(J))-XM)+AT5*(Y(M4(J))-YM)+AT8*(Z(M4(J))-ZM)
    !
    !    DETERMINATION DU JACOBIEN EN CHACUN DES POINTS DE GAUSS
    !
    X41=XL4-XL1
    X32=XL3-XL2
    Y34=YL3-YL4
    Y21=YL2-YL1
    Y41=YL4-YL1
    Y32=YL3-YL2
    X34=XL3-XL4
    X21=XL2-XL1
    DO 15 L=1,NG
    A=(1.-YY(L,IMETH))*X41+(1+YY(L,IMETH))*X32
    B=(1.-XX(L,IMETH))*Y21+(1+XX(L,IMETH))*Y34
    C=(1.-XX(L,IMETH))*X21+(1+XX(L,IMETH))*X34
    D=(1.-YY(L,IMETH))*Y41+(1+YY(L,IMETH))*Y32
    XJAC(L,J)=ABS(A*B-C*D)*WG(L,IMETH)*.25
    15 CONTINUE
    !
    !    COORDONNES ABSOLUES DES POINTS DE GAUSS
    !
    EPSN=0.
    DO 20 L=1,NG
    AA=.25*(1-XX(L,IMETH))*(1-YY(L,IMETH))
    BB=.25*(1-XX(L,IMETH))*(1+YY(L,IMETH))
    CC=.25*(1+XX(L,IMETH))*(1+YY(L,IMETH))
    DD=.25*(1+XX(L,IMETH))*(1-YY(L,IMETH))
    XG1=AA*XL1+BB*XL2+CC*XL3+DD*XL4
    YG1=AA*YL1+BB*YL2+CC*YL3+DD*YL4
    XGA(L,J)=XM+AT1*XG1+AT2*YG1+AT3*EPSN
    YGA(L,J)=YM+AT4*XG1+AT5*YG1+AT6*EPSN
    ZGA(L,J)=ZM+AT7*XG1+AT8*YG1+AT9*EPSN
    20  CONTINUE
    !
    RETURN
  END SUBROUTINE SOMGO
  
  FUNCTION ZJ(Z,CEX)
    
    COMPLEX ZJ, Z,Y,CEX
    IF(REAL(Z)+16.)2,2,1
    2 Y=1./Z                                                                    
    ZJ=Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))                      
    RETURN                                                                    
    1 T=AIMAG(Z)                                                                
    IF(ABS(T)-10.)5,5,6                                                       
    6 ZJ=0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
    IF(T)33,44,44                                                             
    5 CONTINUE
    IF(REAL(Z)+0.5)7,7,8
    8 ZJ=-(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01+Z*(-0.8700861E-03+Z*0.2989204E-03))))))                                                         
    IF(T)33,44,44   
    33 ZJ=ZJ-(0.,3.14159265)*CEX                                  
    RETURN                                                                    
    44 ZJ=ZJ+(0.,3.14159265)*CEX
    RETURN
    7 CONTINUE
    IF(T)3,4,4                                                                
    3 ZJ=((((((( (1.000000,1.3935496E-06)*Z+ (15.82958,-20.14222)) *Z+ (-70.52863,-227.9511))*Z+ (-985.4221,-226.6272))*Z + (-1202.318,1580.907))*Z+ (953.2441,1342.447))*Z + (417.3716,-196.6665))*Z+ (-9.881266,-24.24952))/ (((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,-20.14481))*Z + (-55.66969,-248.1167))*Z+ (-1068.640,-434.4707))*Z + (-2082.250,1522.471))*Z+ (383.3455,2730.378))*Z + (1216.791,351.7189))*Z+ (115.3926,-161.2647))*Z + (-3.777369,-4.510900))-(0.,3.14159265)*CEX
    RETURN
    4 ZJ=((((((( (1.000000,-1.3935496E-06)*Z+ (15.82958,20.14222)) *Z+ (-70.52863,227.9511))*Z+ (-985.4221,226.6272))*Z + (-1202.318,-1580.907))*Z+ (953.2441,-1342.447))*Z + (417.3716,196.6665))*Z+ (-9.881266,24.24952))/ (((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,20.14481))*Z + (-55.66969,248.1167))*Z+ (-1068.640,434.4707))*Z + (-2082.250,-1522.471))*Z+ (383.3455,-2730.378))*Z + (1216.791,-351.7189))*Z+ (115.3926,161.2647))*Z + (-3.777369,4.510900))+(0.,3.14159265)*CEX
    RETURN
  END FUNCTION ZJ
  
  
  
!   FUNCTION IL(NOM)
!   CHARACTER*10 NOM
!   IL=INDEX(NOM,' ')-1
!   IF(IL.LT.0)IL=LEN(NOM)
!   IF(IL.EQ.0)THEN
!   NOM='defaut'
!   IL=6
!   ENDIF
!   RETURN
!   END

!added by RK for preparing second term of the greenfunction
  SUBROUTINE CREK(ID)
    !this is is preparation for constructing Second part of Green Function Snm^2 NEMOH paper Eq 35
    !In this subroutine X discretization (AKR) and Z discretization (AKZ) are constructed base on Eq. 46 in NEMOH paper
    USE QTFCOM_VAR
    USE MIdentification
!
    IMPLICIT NONE
    INTEGER:: I,J,JZ,IR
    REAL:: AKZ,AKR
!   ID
    TYPE(TID) :: ID

    JZ=130
    IR=700
    DO 8006 J=1,JZ
      AKZ=min(10**(J/10.-10),251.)    ! Z discretization Eq 46
      XZ(J)=-AKZ     
8006 CONTINUE
    XR(1)=0.
    DO 8007 I=2,IR
      IF(I.LT.81)THEN
        AKR=10**((I-1.)/10.-8)  ! X discretization Eq 46
      ELSE
        AKR=ABS((I-75.)/6.)
      ENDIF
      XR(I)=AKR
8007 CONTINUE
     DO 8009 J=1,JZ
      DO 8008 I=1,IR
        CALL VNSG(XZ(J),XR(I),I,J)                                           
8008 CONTINUE
8009 CONTINUE

    OPEN(UNIT=44,FILE=ID%ID(1:ID%lID)//'/QTF/GRIN.QAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
    WRITE(44)IR,JZ,(XR(I),I=1,IR),(XZ(J),J=1,JZ)
    DO J=1,JZ
     WRITE(44)(APD1X(I,J),I=1,IR),(APD1Z(I,J),I=1,IR),(APD2X(I,J),I=1,IR),(APD2Z(I,J),I=1,IR)
    ENDDO
    CLOSE(UNIT=44)

    RETURN
!
    END SUBROUTINE
!---------------------------------------------------------                                                                    

      SUBROUTINE VNSG(AKZ,AKR,I,J) 
      USE QTFCOM_VAR
      USE ELEMENTARY_FNS
      IMPLICIT NONE         
      INTEGER:: I,J,IT
      REAL:: AKZ,AKR,CT,ST,CQIT,TETA
      REAL:: QQT(NPIN),CQT(NPIN)
      REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
      COMPLEX:: IM,C1,C2,ZIK,GZ,CEX            
      IM=(0.,1.)                                                                
     ! PI=4.*ATAN(1.)
      CALL COFINT(CQT,QQT)         
      FD1JX=0.                                                              
      FD1JZ=0.                                                              
      FD2JX=0.                                                              
      FD2JZ=0.                                                              
      DO 30 IT=1,NPIN                                           
        TETA=QQT(IT)                                                           
        CQIT=CQT(IT)                                                   
        CT=COS(TETA)
        ST=SIN(TETA)
        ZIK=AKZ+IM*AKR*CT    
        IF(REAL(ZIK)+30.)2,2,1
         2 CEX=(0.,0.)
        GOTO 3 
         1 CEX=CEXP(ZIK)                                               
         3 GZ=GG(ZIK,CEX)                      
        C1=CQIT*(GZ-1./ZIK)
        C2=CQIT*CEX
        FD1JX=FD1JX+CT*AIMAG(C1)                                 
        FD1JZ=FD1JZ+REAL(C1)                                     
        FD2JX=FD2JX+CT*AIMAG(C2)                                 
        FD2JZ=FD2JZ+REAL(C2)                                     
   30 CONTINUE                                                                
      APD1X(I,J)=FD1JX                                               
      APD1Z(I,J)=FD1JZ                                               
      APD2X(I,J)=FD2JX                                               
      APD2Z(I,J)=FD2JZ                                              
      RETURN 
                                                                   
      END SUBROUTINE  
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      SUBROUTINE COFINT(CQT,QQT)  
      USE QTFCOM_VAR
    IMPLICIT NONE 
      INTEGER :: J
      REAL:: QQT(NPIN),CQT(NPIN) 

     ! PI=4.*ATAN(1.)
      DO 160 J=1,NPIN   
      QQT(J)=-PI/2.+(J-1.)/(NPIN-1.)*PI
      IF(J-1)161,161,162
  161 CQT(J)=PI/(3.*(NPIN-1.))
      GOTO 160
  162 IF(J-NPIN)163,161,161
  163 IF(MOD(J,2))164,165,164
  164 CQT(J)=2./(3.*(NPIN-1.))*PI
      GOTO 160
  165 CQT(J)=4./(3.*(NPIN-1.))*PI
  160 CONTINUE

      RETURN                                                                    
      END SUBROUTINE   
!------------------------------------------------------------------

  
END MODULE MQTFBaseFunctions
  
