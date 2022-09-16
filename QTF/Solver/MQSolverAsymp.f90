MODULE MQSOLVERASYMP
!all functions used in ASYMP calculation

USE CONSTANTS          , ONLY:II,CZERO,PI
USE MEnvironment,        ONLY:Fun_Dispersion
USE Elementary_functions,ONLY:CIH,fun_BESSJ,   &
                              Fun_KronDelta
IMPLICIT NONE

CONTAINS
  
  FUNCTION Fun_IDF(w1,w2,k1,k2,delk,sumk,g,D,Rf,NRf,Nbessel,             &
                  Ivartheta1l,Ivartheta2l,IR1l,IR2l)  result(IDF)

     REAL,                              INTENT(IN):: k1,k2,w1,w2,delk,sumk,g,D
     INTEGER,                           INTENT(IN):: NRf, Nbessel
     REAL, DIMENSION(NRf),              INTENT(IN):: Rf
     COMPLEX,DIMENSION(2,NBESSEL+1),    INTENT(IN):: IR1l,IR2l
     COMPLEX,DIMENSION(2,NBESSEL+3),    INTENT(IN):: Ivartheta1l,Ivartheta2l


     COMPLEX,DIMENSION(2,2)                       :: IDF
     COMPLEX,DIMENSION(2)                         :: IDF11,IDF12
     COMPLEX,DIMENSION(2)                         :: IDF21,IDF22
     INTEGER                                      :: Ibessel
     COMPLEX,DIMENSION(2)                         :: sumIRIVartheta_11,sumIRIVartheta_12  
     COMPLEX,DIMENSION(2)                         :: sumIRIVartheta_21,sumIRIVartheta_22  
     
     sumIRIVartheta_11(:)=CZERO
     sumIRIVartheta_12(:)=CZERO
     sumIRIVartheta_21(:)=CZERO
     sumIRIVartheta_22(:)=CZERO

     DO Ibessel=0,Nbessel
        sumIRIVartheta_11=sumIRIVartheta_11                                             &
               +0.5*IR1l(:,Ibessel)*(Ivartheta1l(:,Ibessel-1)+Ivartheta1l(:,Ibessel+1))    
        sumIRIVartheta_12=sumIRIVartheta_12                                             &
               +0.5*IR2l(:,Ibessel)*(Ivartheta2l(:,Ibessel-1)+Ivartheta2l(:,Ibessel+1))
        sumIRIVartheta_21=sumIRIVartheta_21                                             &
               +IR1l(:,Ibessel)*Ivartheta1l(:,Ibessel)   
        sumIRIVartheta_22=sumIRIVartheta_22                                             &
               +IR2l(:,Ibessel)*Ivartheta2l(:,Ibessel) 
     ENDDO
     IDF11(1)=(-II*g*8*PI*SQRT(k2*delk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_11(1)

     IDF11(2)=(-II*g*8*PI*SQRT(k2*sumk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_11(2)

     IDF12(1)=(II*g*8*PI*SQRT(k1*delk)/w2)                                              &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_12(1)

     IDF12(2)=(-II*g*8*PI*SQRT(k1*sumk)/w2)                                             &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_12(2)

     IDF21(1)=(-II*g*8*PI*SQRT(k2*delk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_21(1)

     IDF21(2)=(-II*g*8*PI*SQRT(k2*sumk)/w1)                                             &
               *Fun_Fprofile(k2,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_21(2)

     IDF22(1)=(II*g*8*PI*SQRT(k1*delk)/w2)                                              &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(delk,D,0.)*sumIRIVartheta_22(1)

     IDF22(2)=(-II*g*8*PI*SQRT(k1*sumk)/w2)                                             &
               *Fun_Fprofile(k1,D,0.)*Fun_Fprofile(sumk,D,0.)*sumIRIVartheta_22(2)
     IDF(:,1)=IDF11+IDF12
     IDF(:,2)=IDF21+IDF22
     IF (delk==0.) IDF(1,:)=CZERO
  END FUNCTION

  FUNCTION Fun_KAPPA2_DIFFSUM(k1,k2,w1,w2,D,g) RESULT(KAPPA2)
     REAL,        INTENT(IN):: k1,k2,w1,w2,g,D
     COMPLEX,DIMENSION(2)   :: KAPPA2
     REAL                   :: Nu1,Nu2,SECH2_K1D,SECH2_K2D

     Nu1=Fun_Dispersion(k1,D,g)**2/g 
     Nu2=Fun_Dispersion(k2,D,g)**2/g
     IF (D.EQ.0.) THEN !INFINITE DEPTH
     SECH2_K1D=0
     SECH2_K2D=0
     ELSE              !FINITE DEPTH
     SECH2_K1D=1-tanh(k1*D)**2 
     SECH2_K2D=1-tanh(k2*D)**2
     ENDIF 
     KAPPA2(1)=II*(w1-w2)*Nu1*Nu2+II*w1*w2/g*(k1**2*SECH2_K1D/w1-k2**2*SECH2_K2D/w2)
     KAPPA2(2)=II*(w1+w2)*Nu1*Nu2-II*w1*w2/g*(k1**2*SECH2_K1D/w1+k2**2*SECH2_K2D/w2)

  END FUNCTION 

  FUNCTION Fun_Fprofile(k,D,z) result(F)
     REAL,      INTENT(IN):: k,D,z
     REAL                 :: F
     
     IF(D.EQ.0..OR.(k*(D+z).GT.50)) THEN !INFINITE DEPTH case
         F=exp(k*z)
     ELSE
         F=cosh(k*(D+z))/( (k*D*(1-tanh(k*D)**2)+tanh(k*D))*cosh(k*D) )
     ENDIF
  END FUNCTION

  FUNCTION Fun_IR1l(l,k1,k2,delk,sumk,Rf,NRf) result(IR1l)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     INTEGER               :: eps_l 
     COMPLEX ,DIMENSION(2) :: IR1_0Rext
     COMPLEX ,DIMENSION(2) :: IR1l
     eps_l=Fun_epsilon_l(l)
        
     IR1_0Rext=fun_IR1_0Rext(l,k1,k2,delk,sumk,Rf,NRf)
     IR1l(1)=eps_l*II**l*(fun_gamma(l,delk-k2,k1)-IR1_0Rext(1))
     IR1l(2)=eps_l*II**l*(II*fun_gamma(l,sumk+k2,k1)-IR1_0Rext(2))
     
    !print*,l,IR1_0Rext(:)
    !print*,l,delk-k1,k1,fun_gamma(l,delk-k2,k1)
  END FUNCTION

  FUNCTION Fun_IR1_0Rext(l,k1,k2,delk,sumk,Rf,NRf) result(IR1)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     COMPLEX,DIMENSION(2)           :: IR1   
     REAL                           :: dRf
     INTEGER ::Ir

     IR1=CZERO
     DO Ir=1,NRf-1
        dRf=Rf(Ir+1)-Rf(Ir)
        IR1(1)=IR1(1)+0.5*(                                                     &
                exp(II*(delk-k2)*Rf(Ir))*fun_BESSJ(l,k1*Rf(Ir))                 &
               +exp(II*(delk-k2)*Rf(Ir+1))*fun_BESSJ(l,k1*Rf(Ir+1)) )*dRf
        IR1(2)=IR1(2)+0.5*(                                                     &
              exp(II*( (sumk+k2)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k1*Rf(Ir))        &
             +exp(II*( (sumk+k2)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k1*Rf(Ir+1)) )*dRf
     ENDDO
  END FUNCTION

  FUNCTION Fun_IR2l(l,k1,k2,delk,sumk,Rf,NRf) result(IR2l)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     INTEGER               :: eps_l 
     COMPLEX ,DIMENSION(2) :: IR2_0Rext
     COMPLEX ,DIMENSION(2) :: IR2l
     eps_l=Fun_epsilon_l(l)

     IR2_0Rext=fun_IR2_0Rext(l,k1,k2,delk,sumk,Rf,NRf)
     IR2l(1)=eps_l*CONJG(II**l)*(II*fun_gamma(l,delk+k1,k2)-IR2_0Rext(1))
     IR2l(2)=eps_l*II**l*(II*fun_gamma(l,sumk+k1,k2)-IR2_0Rext(2))

    ! print*,l,IR2_0Rext(:)
     !print*,l,delk+k1,k2,fun_gamma(l,delk+k1,k2)
END FUNCTION

  FUNCTION Fun_IR2_0Rext(l,k1,k2,delk,sumk,Rf,NRf) result(IR2)
     INTEGER,             INTENT(IN):: l,NRf
     REAL,                INTENT(IN):: k1,k2,delk,sumk
     REAL,DIMENSION(NRf), INTENT(IN):: Rf
     COMPLEX,DIMENSION(2)           :: IR2   
     REAL                           :: dRf
     INTEGER ::Ir

     IR2=CZERO
     DO Ir=1,NRf-1
        dRf=Rf(Ir+1)-Rf(Ir)
        IR2(1)=IR2(1)+0.5*(                                                     &
             exp(II*( (delk+k1)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k2*Rf(Ir))         &
            +exp(II*( (delk+k1)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k2*Rf(Ir+1)) )*dRf
        IR2(2)=IR2(2)+0.5*(                                                     &
             exp(II*( (sumk+k1)*Rf(Ir)+PI/2   ))*fun_BESSJ(l,k2*Rf(Ir))         &
            +exp(II*( (sumk+k1)*Rf(Ir+1)+PI/2 ))*fun_BESSJ(l,k2*Rf(Ir+1)) )*dRf
     ENDDO
  END FUNCTION

  FUNCTION PREPARE_KOCHIN_COEFFICIENTS(Isym,Npanels,XM,Apanels,D,Nbessel,k,zig) &
                  result(CmSm)
   INTEGER,                             INTENT(IN) :: Isym,Npanels,Nbessel
   REAL,                                INTENT(IN) :: k,D
   REAL,DIMENSION(Npanels),             INTENT(IN) :: Apanels
   REAL,DIMENSION(3,Npanels),           INTENT(IN) :: XM
   COMPLEX,DIMENSION(Npanels*2**Isym)  ,INTENT(IN) :: zig
        
   COMPLEX,DIMENSION(2,Nbessel+1) :: CmSm
   INTEGER                      :: ll
   DO ll=0,Nbessel
   CmSm(:,ll)=Fun_KochinCoefs_l(Isym,Npanels,XM,Apanels,zig,k,D,ll)
   ENDDO
  END FUNCTION

  FUNCTION Fun_IVartheta1l(Nbessel,beta,CmSmPer,CnSnRad_delk,CnSnRad_sumk,ll)&
                                                         result(Ivartheta)
     INTEGER,                             INTENT(IN) :: ll,Nbessel
     REAL,                                INTENT(IN) :: beta
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CmSmPer
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_delk
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_sumk


     INTEGER                       :: mm,nn
     COMPLEX,DIMENSION(2)          :: Ivartheta
     COMPLEX,DIMENSION(2)          :: CnRad,SnRad
     COMPLEX                       :: CmPer,SmPer

     Ivartheta=CZERO
     DO mm=0,Nbessel
        CmPer=CmSmPer(1,mm)
        SmPer=CmSmPer(2,mm)
        DO nn=0,Nbessel
          CnRad(1)=CnSnRad_delk(1,nn)
          SnRad(1)=CnSnRad_delk(2,nn)
          CnRad(2)=CnSnRad_sumk(1,nn)
          SnRad(2)=CnSnRad_sumk(2,nn)
          !for diff freq
          Ivartheta(1)=Ivartheta(1)                                     &
              +CONJG(CmPer)*CnRad(1)*cos(ll*beta)*Fun_Dlmn(1,ll,mm,nn)  &
              +CONJG(CmPer)*SnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,nn,mm)  &
              +CONJG(SmPer)*CnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,mm,nn)  &
              +CONJG(SmPer)*SnRad(1)*cos(ll*beta)*Fun_Dlmn(0,nn,mm,ll) 
          !for sum freq    
          Ivartheta(2)=Ivartheta(2)                                     &
              +CmPer*CnRad(2)*cos(ll*beta)*Fun_Dlmn(1,ll,mm,nn)         &
              +CmPer*SnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,nn,mm)         &
              +SmPer*CnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,mm,nn)         &
              +SmPer*SnRad(2)*cos(ll*beta)*Fun_Dlmn(0,nn,mm,ll)            
        ENDDO
     ENDDO
  END FUNCTION

  FUNCTION Fun_IVartheta2l(Nbessel,beta,CmSmPer,CnSnRad_delk,CnSnRad_sumk,ll)&
                                                         result(Ivartheta)
     INTEGER,                             INTENT(IN) :: ll,Nbessel
     REAL,                                INTENT(IN) :: beta
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CmSmPer
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_delk
     COMPLEX,DIMENSION(2,Nbessel+1),      INTENT(IN) :: CnSnRad_sumk


     INTEGER                       :: mm,nn
     COMPLEX,DIMENSION(2)          :: Ivartheta
     COMPLEX,DIMENSION(2)          :: CnRad,SnRad
     COMPLEX                       :: CmPer,SmPer

     Ivartheta=CZERO
     DO mm=0,Nbessel
        CmPer=CmSmPer(1,mm)
        SmPer=CmSmPer(2,mm)
        DO nn=0,Nbessel
          CnRad(1)=CnSnRad_delk(1,nn)
          SnRad(1)=CnSnRad_delk(2,nn)
          CnRad(2)=CnSnRad_sumk(1,nn)
          SnRad(2)=CnSnRad_sumk(2,nn)

          Ivartheta(1)=Ivartheta(1)                                     &
              +CmPer*CnRad(1)*cos(ll*beta)*Fun_Dlmn(1,ll,mm,nn)         &
              +CmPer*SnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,nn,mm)         &
              +SmPer*CnRad(1)*sin(ll*beta)*Fun_Dlmn(0,ll,mm,nn)         &
              +SmPer*SnRad(1)*cos(ll*beta)*Fun_Dlmn(0,nn,mm,ll)     
          Ivartheta(2)=Ivartheta(2)                                     &
              +CmPer*CnRad(2)*cos(ll*beta)*Fun_Dlmn(1,ll,mm,nn)         &
              +CmPer*SnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,nn,mm)         &
              +SmPer*CnRad(2)*sin(ll*beta)*Fun_Dlmn(0,ll,mm,nn)         &
              +SmPer*SnRad(2)*cos(ll*beta)*Fun_Dlmn(0,nn,mm,ll)            
        ENDDO
     ENDDO
  END FUNCTION


  FUNCTION Fun_KochinCoefs_l(Isym,Npanels,XM,AreaPanel,sig,k,D,ll) result(ClSl)
        INTEGER,                          INTENT(IN) :: Npanels,ll,Isym
        REAL,                             INTENT(IN) :: k,D
        REAL,DIMENSION(Npanels),          INTENT(IN) :: AreaPanel
        REAL,DIMENSION(3,Npanels),        INTENT(IN) :: XM
        COMPLEX,DIMENSION(Npanels*2*Isym),INTENT(IN) :: sig
        COMPLEX,DIMENSION(2)                         :: ClSl
        COMPlEX                                      :: calcClSl
        INTEGER                                      :: Ipanel
        REAL                                         :: r,alpha
        ClSl=CZERO
        DO Ipanel=1,Npanels
          r=sqrt(XM(1,Ipanel)**2+XM(2,Ipanel)**2)
          alpha=atan2(XM(2,Ipanel),XM(1,Ipanel))
          calcClSl=sig(Ipanel)*CIH(k,XM(3,Ipanel),D)                        &
                  *Fun_epsilon_l(ll)*(-II)**ll*fun_BESSJ(ll,k*r)            &
                  *AreaPanel(Ipanel)
          ClSl(1)=ClSl(1)+calcClSl*cos(ll*alpha)
          ClSl(2)=ClSl(2)+calcClSl*sin(ll*alpha)

          IF (Isym.EQ.1) THEN
            alpha=atan2(-XM(2,Ipanel),XM(1,Ipanel))
            calcClSl=sig(Ipanel)*CIH(k,XM(3,Ipanel),D)                       &
                    *Fun_epsilon_l(ll)*(-II)**ll*fun_BESSJ(ll,k*r)           &
                    *AreaPanel(Ipanel)
            ClSl(1)=ClSl(1)+calcClSl*cos(ll*alpha)
            ClSl(2)=ClSl(2)+calcClSl*sin(ll*alpha)
          ENDIF 
        ENDDO
        ClSl=-ClSl/4/PI
  END FUNCTION
  
  FUNCTION Fun_Dlmn(IDEq,l,m,n) result(Dlmn)
        INTEGER,        INTENT(IN):: IDEq,l,m,n
        REAL                      :: Dlmn

        IF (IDEq.EQ.0) THEN
        !Integ sin(l*vartheta)*sin(m*vartheta)*cos(n*vartheta)dvartheta=0
        Dlmn=(Fun_KronDelta(n,abs(m-l))-Fun_KronDelta(n,m+l))
        ELSE
        !Integ sin(l*vartheta)*sin(m*vartheta)*cos(n*vartheta)dvartheta=0
        Dlmn=(Fun_KronDelta(n,abs(m-l))+Fun_KronDelta(n,m+l))
        ENDIF
        Dlmn=Dlmn*PI/Fun_epsilon_l(n)
  END FUNCTION

  FUNCTION Fun_gamma(m,beta,alpha) result(gam)
     INTEGER, INTENT(IN) :: m
     REAL,    INTENT(IN) :: beta,alpha
     COMPLEX             :: gam

     IF (0.LE.ABS(beta) .AND.ABS(beta).LE.alpha) THEN
        gam=exp(II*m*asin(ABS(beta)/alpha))/SQRT(alpha**2-beta**2)
     ELSEIF(0.LT.alpha .AND. alpha.LT.ABS(beta)) THEN
        gam=II*exp(II*m*PI/2)*(alpha/(ABS(beta)+SQRT(beta**2-alpha**2)))**m  &
              /SQRT(beta**2-alpha**2)
     ENDIF
     IF (beta<0)  gam=CONJG(gam)
  END FUNCTION
 
  
  FUNCTION Fun_epsilon_l(l) result(eps)
     INTEGER, INTENT(IN) :: l
     INTEGER             :: eps
     IF (l.EQ. 0) eps=1
     IF (l.NE. 0) eps=2
  END FUNCTION
!!---------------------------------------------------------------------


END MODULE
