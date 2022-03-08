#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#
#   NEMOH2 - second order (QTF)
#   Contributors list:
#   - Jean-Noel Dory, INNOSEA (jeannoel.dory@gmail.com) --> 07/2015
#
#--------------------------------------------------------------------------------------

from common import *
from loadnemoh import *
import struct
import numpy as np
import os
import shutil


def generate_mesh(nemohpath, NPASTH, Lligne, RM=0.10, RCEREXT=300, RLIGNE=10000):
    ''' Cette fonction genere un maillage de la surface libre
    pour une ligne de flottaison circulaire de RM a RCEREXT avec NPASTH noeuds
    dans la direction radiale.
    Ce maillage est utilise par hasfs.f et Solver2.f90.

    Input:
        NPASTH (int)	nombre de points othoradiaux
        RCEREXT(float)  distance jusqu'a laquelle la surface est maillee (m)
        Lligne (bool)	prise en compte de la ligne ?
        RLIGNE (int)	pour visualisation
        RM 	   (float) 	distance minimale (m) (defaut: a 10cm de la ligne de flottaison)
    Output:
        SF_L12.dat		maillage de la surface au format L12 (lu par NEMOH)
        SF_L12.vtk		maillage au format vtk pour visualisation
    '''
    
    print('NPASTH='+str(NPASTH)+'   RCEREXT='+str(RCEREXT))

    # parametre de calibration
    NPASR0=200		# nombre de pas de calibration
    RCEREXT0=310.		# distance max de calibration
    PUISS=2			# variation de la densitee radiale de noeuds (cf l.73)
    NFASL=300000		# attention NFASL defini dans nemoh --> peut changer
    Mbess=1000.		# argument maximum des fonction de bessel


    changename(nemohpath,'SF_L12','.dat')
    changename(nemohpath,'SF_L12','.vtk')

    # lecture des fichiers de parametrage de nemoh
    w,G,H,namemb=read_Nemoh_cal(nemohpath)
    Wmax=max(w)
    #print(str(nemohpath+namemb[0]))
    nemohpathmeshdat=nemohpath+'/mesh/'
    L12=open(nemohpathmeshdat+namemb[0])
    NSYMY=int(L12.readline().strip().split(' ')[-1])
    NJ=NSYMY+1
    if NJ==1:
        NPASTH*=2
    L12.close()
    ##

    # definition de la taille (radiale) maximale des mailles (dr < lambda_min/30)
    # w^2=k.g.th(k.h)
    KMAX=Wmax**2/G			# nombre d'onde max (deep water)
    if KMAX*H < 0.5:
        KMAX=Wmax/np.sqrt(G*H)	# shalow water
    elif KMAX*H < 2:
        KMAX*=2.			# intermediate water --> on maximise k

    RCEREXT=min(RCEREXT,(Mbess-1.)/KMAX)					# correction du rayon limit : si Rcerext*K> Mbess, asymp de peut pas tourner
    MDR=2*np.pi/30./KMAX
    print('max(dr)= '+str(MDR))
    print('Rext= '+str(RCEREXT))
    ##


    # creation de la variation radiale
    RIND=(RCEREXT0-RM)/NPASR0**PUISS	# rayon caracteristique pour que Ri= RIND.i^2
    IMAX=np.int((MDR/RIND-1)/2)		# indice a partir duquel dr > MDR     (dr=(i+1).RIND)
    NPASR=np.int(np.sqrt((RCEREXT-RM)/RIND))# nombre de point radiaux
    if IMAX>NPASR:	# cas ou on a toujours dr < MDR
        RSL=RM+RIND*np.arange(NPASR+1)**2
    else:		# sinon on prolonge par une serie de noeud avec dr = MDR
        RSL=RM+RIND*np.arange(IMAX)**2
        RSL=np.append(RSL,RSL[-1]+MDR*(1+np.arange((RCEREXT-RSL[-1])/MDR)))
    NPASR=len(RSL)-1
    RCEREXT=RSL[-1]
    print('NPASR = '+str(NPASR)+', RCEREXT = '+str(RCEREXT))
    ##

    RSLL=np.append(RSL,RSL[-1]+MDR*(1+np.arange((RLIGNE-RCEREXT)/MDR)))
    if Lligne:
        NPASL=len(RSLL)-1
    else:
        NPASL=0


    NFACSL=0
    NCONT=0
    NOMBPO=0
    XSL=np.zeros((NPASTH+1)*(NPASR+1))
    YSL=np.zeros((NPASTH+1)*(NPASR+1))
    M1SL=np.zeros((NPASR+2)*NPASTH,dtype=int)
    M2SL=np.zeros((NPASR+2)*NPASTH,dtype=int)
    M3SL=np.zeros(NPASTH*NPASR,dtype=int)
    M4SL=np.zeros(NPASTH*NPASR,dtype=int)
    AIRESL=np.zeros((NPASR+2)*NPASTH)
    CNSL=np.zeros((2,(NPASR+2)*NPASTH))
    XMSL=np.zeros(NPASTH*(NPASR+2)+NPASL)
    YMSL=np.zeros(NPASTH*(NPASR+2)+NPASL)

    ## creation du maillage de la surface
    for IR in range(NPASR+1):
        for ITH in range(NPASTH+1):
            # creation des noeuds
            NUMPT=ITH+1+IR*(NPASTH+1)   #removed '.' after 1 to be an integer by RK
            PSL=RSL[IR]*np.exp(1j*ITH*2*np.pi/(NJ*NPASTH))
            # print('NUMPT='+str(NUMPT)+' PSL='+str(PSL))
            XSL[NUMPT-1]=np.real(PSL)
            YSL[NUMPT-1]=np.imag(PSL)
            NOMBPO=NOMBPO+1
            if (IR != NPASR) and (ITH != NPASTH):
            # creation des mailles (lien entre les noeuds)
                NFACSL=NFACSL+1
                M1SL[NFACSL-1]=NUMPT
                M2SL[NFACSL-1]=NUMPT+1
                M3SL[NFACSL-1]=NUMPT+NPASTH+2
                M4SL[NFACSL-1]=NUMPT+NPASTH+1
    for KKK in range(NFACSL):
        X1=XSL[M1SL[KKK]-1]
        Y1=YSL[M1SL[KKK]-1]
        X2=XSL[M2SL[KKK]-1]
        Y2=YSL[M2SL[KKK]-1]
        X3=XSL[M3SL[KKK]-1]
        Y3=YSL[M3SL[KKK]-1]
        X4=XSL[M4SL[KKK]-1]
        Y4=YSL[M4SL[KKK]-1]
        XMSL[KKK]=0.25*(X1+X2+X3+X4)
        YMSL[KKK]=0.25*(Y1+Y2+Y3+Y4)
        H1=np.sqrt((X3-X2)**2+(Y3-Y2)**2)
        G1=np.sqrt((X2-X1)**2+(Y2-Y1)**2)
        H2=np.sqrt((X4-X1)**2+(Y4-Y1)**2)
        G2=np.sqrt((X4-X3)**2+(Y4-Y3)**2)
        HAUTEUR1=0.5*(H1+H2)
        GRANDEU2=0.5*(G1+G2)
        AIRESL[KKK]=HAUTEUR1*GRANDEU2
    #

    ## creation du contour (normale vers l'exterieur)
    for IR in range(NPASR+1):
        for ITH in range(NPASTH+1):
            NUMPT=ITH+1.+IR*(NPASTH+1)
            if (IR == NPASR) and (ITH != NPASTH):
                NCONT=NCONT+1
                M1SL[NCONT+NFACSL-1]=NUMPT
                M2SL[NCONT+NFACSL-1]=NUMPT+1
                X1=XSL[M1SL[NCONT+NFACSL-1]-1]
                Y1=YSL[M1SL[NCONT+NFACSL-1]-1]
                X2=XSL[M2SL[NCONT+NFACSL-1]-1]
                Y2=YSL[M2SL[NCONT+NFACSL-1]-1]
                XMSL[NCONT+NFACSL-1]=0.5*(X1+X2)
                YMSL[NCONT+NFACSL-1]=0.5*(Y1+Y2)
                G1=np.sqrt((X2-X1)**2+(Y2-Y1)**2)
                AIRESL[NCONT+NFACSL-1]=G1
                CNSL[0,NCONT+NFACSL-1]=(Y2-Y1)/G1
                CNSL[1,NCONT+NFACSL-1]=-(X2-X1)/G1
            elif (IR == 0) and (ITH != NPASTH):
                NCONT=NCONT+1
                M1SL[NCONT+NFACSL-1]=NUMPT
                M2SL[NCONT+NFACSL-1]=NUMPT+1
                X1=XSL[M1SL[NCONT+NFACSL-1]-1]
                Y1=YSL[M1SL[NCONT+NFACSL-1]-1]
                X2=XSL[M2SL[NCONT+NFACSL-1]-1]
                Y2=YSL[M2SL[NCONT+NFACSL-1]-1]
                XMSL[NCONT+NFACSL-1]=0.5*(X1+X2)
                YMSL[NCONT+NFACSL-1]=0.5*(Y1+Y2)
                G1=np.sqrt((X2-X1)**2+(Y2-Y1)**2)
                AIRESL[NCONT+NFACSL-1]=G1
                CNSL[0,NCONT+NFACSL-1]=-(Y2-Y1)/G1
                CNSL[1,NCONT+NFACSL-1]=(X2-X1)/G1
    ##
    NTOT=NFACSL+NCONT


    xyz=np.append(XSL.reshape(NOMBPO,1),YSL.reshape(NOMBPO,1),axis=1)
    xyz=np.append(xyz,np.zeros((NOMBPO,1)),axis=1)
    NN1=np.append(M1SL[:NFACSL].reshape(NFACSL,1)-1,M2SL[:NFACSL].reshape(NFACSL,1)-1,axis=1)
    NN2=np.append(M3SL.reshape(NFACSL,1)-1,M4SL.reshape(NFACSL,1)-1,axis=1)
    NN=np.append(NN1,NN2,axis=1)
    NN=np.append(np.ones((NFACSL,1))*4,NN,axis=1).astype(int)
    meshpath = os.path.join(nemohpath, 'Mesh')
    filename = os.path.join(meshpath, 'SF_L12')

    createVTK(xyz,NN,filename)
    # createL12(xyz,NN,filename,NSYMY)

    # creation de la ligne pour visualisation
    if (NTOT+NPASL>NFASL):
        print(str(NTOT+NPASL)+'>'+str(NFASL))
        print('NTOT+NPASL > NFASL')
        print('STOP')
    for ILIGNE in range(NPASL):
        XMSL[NTOT+ILIGNE]=RSLL[ILIGNE]*np.cos(np.pi/4)
        YMSL[NTOT+ILIGNE]=RSLL[ILIGNE]*np.sin(np.pi/4)

    print('NFACSL: '+str(NFACSL)+'   NCONT: '+str(NCONT)+'  NLIGNE: '+str(NPASL))
    NTOT=NTOT+NPASL
    print('NTOT CALCUL:'+str(NTOT))
    if (NTOT>NFASL):
        print('PLUS QUE NFASL FACETTE DANS LA SL')
    ##


    # creation du fichier de meta data
    FStot=open(filename + '.dat', 'w')
    FStot.write('NSYMY    NFACSL    NCONT    NOMBPO    NLIGNE    RCEREXT    NPASR    NTOT\n')
    FStot.write('XLS[]    YSL[]    : COORDONNEES DES NOEUDS DU MAILLAGE\n')
    FStot.write('XMLS[]    YMSL[]    : COORDONNEES DES CENTRES DES MAILLES\n')
    FStot.write('AIRESL[]    : SURFACE (OU LONGUEUR) DES MAILLES\n')
    FStot.write('CNX[]    CNY[]    : NORMALE AU CONTOUR\n')
    FStot.write('M1SL[]    M2SL[]    : INDICES DES POINTS 1 ET 2 DES MAILLES\n')
    FStot.write('M3SL[]    M4SL[]    : INDICES DES POINTS 3 ET 4 DES MAILLES\n')
    FStot.write(str(NSYMY)+' '+str(NFACSL)+' '+str(NCONT)+' '+str(NOMBPO)+' '+str(NPASL)+' '+str(RCEREXT)+' '+str(NPASR)+' '+str(NTOT)+'\n')
    for i in range(NOMBPO):
        FStot.write(str(XSL[i])+' '+str(YSL[i])+'\n')
    for i in range(NFACSL+NCONT+NPASL):
        FStot.write(str(XMSL[i])+' '+str(YMSL[i])+'\n')
    for i in range(NFACSL+NCONT):
        FStot.write(str(AIRESL[i])+'\n')
    for i in range(NCONT):
        j=i+NFACSL
        FStot.write(str(CNSL[0,j])+' '+str(CNSL[1,j])+'\n')
    for i in range(NFACSL+NCONT):
        FStot.write(str(M1SL[i])+' '+str(M2SL[i])+'\n')
    for i in range(NFACSL):
        FStot.write(str(M3SL[i])+' '+str(M4SL[i])+'\n')
    FStot.close()
    ##

    print('MAILLAGE TERMINE')


if __name__ == '__main__':

    NPASTH = 100		# nombre de points othoradiaux
    RCEREXT = 310.		# distance jusqu'a laquelle la surface est maillee
    Lligne = False		# prise en compte de la ligne?
    RLIGNE = 10000.		# pour visualisation
    RM = 10.01			# distance minimale (a 10cm de la ligne de flottaison)
    f = open('workdir.txt', 'r')
    nemohpath=f.read()

    generate_mesh(nemohpath,NPASTH, Lligne, RM, RCEREXT, RLIGNE)
