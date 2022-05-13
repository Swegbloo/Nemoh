#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#
#   NEMOH2 - second order (QTF) - November 2014
#   Contributors list:
#   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
#--------------------------------------------------------------------------------------

import numpy as np
import math


def loadnemohB(foldername):
    filename = foldername + "/results/CA.dat"
    f = open(filename, 'r')

    L = f.readline().strip().split()
    nw = int(L[4])

    w = np.zeros(nw)
    B = np.zeros((6, 6, nw))

    for i in range(nw):
        w[i] = float(f.readline().strip())
        for l in range(6):
            L = f.readline().strip().split()
            B[l, :, i] = [float(j) for j in L]

    return w, B


def loadnemohMa(foldername):
    filename = foldername + "/results/CM.dat"
    f = open(filename, 'r')

    L = f.readline().strip().split()
    nw = int(L[4])

    w = np.zeros(nw)
    Ma = np.zeros((6, 6, nw))

    for i in range(nw):
        w[i] = float(f.readline().strip())
        for l in range(6):
            L = f.readline().strip().split()
            Ma[l, :, i] = [float(j) for j in L]

    return w, Ma


def loadnemohFe(foldername):
    filename = foldername + "/results/Fe.dat"
    wAP = np.genfromtxt(filename, skip_header=3,encoding='latin1')   #added by RK encoding='latin1'

    w = wAP[:, 0]
    Fea = wAP[:, 1:7]
    Feph = wAP[:, 7:] * math.pi / 180.0

    return w, Fea, Feph


def loadnemohQTF(foldername, qtftype, dofs, contrib="tot"):
    resname = foldername + "/results/QTF/"

    if qtftype == "+":
        qtftype = "QTFP"
    else:
        qtftype = "QTFM"

    basenamepink = resname + qtftype + "_APINK"
    basenameduok = resname + qtftype + "_DUOK"
    basenamehasbo = resname + qtftype + "_HASBO"
    basenamehasfs = resname + qtftype + "_HASFS"
    basenameasymp = resname + qtftype + "_ASYMP"

    filenameR = basenameduok + "_DOF_1_PR.dat"
    w, _ = loadnemohQTFfile(filenameR)

    nw = len(w)
    QTF = np.zeros((6, nw, nw)) * 1j

    #TODO: for now, loading only X and Z DOF because QTFSolver misses some
    # output DOF
    for i in dofs:
        i = i - 1

        if (contrib == "duok" or contrib == "tot"):

            filenameR = basenameduok + stringdofRealPart(i)
            filenameI = basenameduok + stringdofImaginaryPart(i)

            w, QTFCR = loadnemohQTFfile(filenameR)
            w, QTFCI = loadnemohQTFfile(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

        if (contrib == "has" or contrib == "tot"):

            filenameR = basenamehasfs + stringdofRealPart(i)
            filenameI = basenamehasfs + stringdofImaginaryPart(i)

            w, QTFCR, Rext = loadnemohQTFfileRext(filenameR)
            w, QTFCI, Rext = loadnemohQTFfileRext(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

            filenameR = basenameasymp + stringdofRealPart(i)
            filenameI = basenameasymp + stringdofImaginaryPart(i)

            w, QTFCR, Rext = loadnemohQTFfileRext(filenameR)
            w, QTFCI, Rext = loadnemohQTFfileRext(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

        if contrib == "pink":

            filenameR = basenamepink + stringdofRealPart(i)
            filenameI = basenamepink + stringdofImaginaryPart(i)

            w, QTFCR = loadnemohQTFfile(filenameR)
            w, QTFCI = loadnemohQTFfile(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

        if (contrib == "hasbo" or contrib == "tot"):

            filenameR = basenamehasbo + stringdofRealPart(i)
            filenameI = basenamehasbo + stringdofImaginaryPart(i)

            w, QTFCR = loadnemohQTFfile(filenameR)
            w, QTFCI = loadnemohQTFfile(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

        if (contrib == "hasfs" or contrib == "tot"):
            # for comparison with Kim&Yue hasfscalc + hasasymp
            filenameR = basenamehasfs + stringdofRealPart(i)
            filenameI = basenamehasfs + stringdofImaginaryPart(i)

            w, QTFCR, Rext = loadnemohQTFfileRext(filenameR)
            w, QTFCI, Rext = loadnemohQTFfileRext(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

            filenameR = basenameasymp + stringdofRealPart(i)
            filenameI = basenameasymp + stringdofImaginaryPart(i)

            w, QTFCR, Rext = loadnemohQTFfileRext(filenameR)
            w, QTFCI, Rext = loadnemohQTFfileRext(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

    return w, QTF


def loadsubcontribQTF(foldername, nw, qtftype, dofs, contribname, isubcontrib):
    resname = foldername + "/results/QTF/Appendix/"

    QTF = np.zeros((6, nw, nw)) * 1j

    subcontribs = {"duok": ['DADV', 'DDPV', 'DRFI', 'DFIX', 'DSTA', 'DFLO'],
                   "hasbo": ['HVPSI', 'HVN', 'HFLO', 'HVRN', 'HIDN', 'HFRK']}

    if qtftype == "+":
        qtftype = "QTFP"
    else:
        qtftype = "QTFM"

    for k in range(isubcontrib):
        subcontribname = subcontribs[contribname][k]
        basename = resname + qtftype + "_" + subcontribname

        for i in dofs:
            i = i - 1
            filenameR = basename + stringdofRealPart(i)
            filenameI = basename + stringdofImaginaryPart(i)
            w, QTFCR = loadnemohQTFfile(filenameR)
            w, QTFCI = loadnemohQTFfile(filenameI)

            QTF[i, :, :] = QTF[i, :, :] + QTFCR + 1j * QTFCI

    return w, QTF


def stringdofRealPart(i):
    return "_DOF_" + str(i + 1) + "_PR.dat"


def stringdofImaginaryPart(i):
    return "_DOF_" + str(i + 1) + "_PI.dat"


def loadnemohQTFfile(filename):
    # if this file does not exist return empty array
    try:
        # print "file : %s" % filename
        M = np.genfromtxt(filename)
    except:
        print "no such file : %s" % filename
        return None
    w = M[1:, 0]
    QTF = M[1:, 1:]

    return w, QTF


def loadnemohQTFfileRext(filename, icer=-1):
    # if this file does not exist return empty array
    try:
        M = np.genfromtxt(filename)
    except:
        print "no such file : %s" % filename
    w = M[0, 1:]
    nw = len(w)
    nr = M.shape[0] / (nw + 1)
    icer = icer % nr

    QTF = M[icer * (nw + 1) + 1:(icer + 1) * (nw + 1), 1:]
    Rext = M[icer * (nw + 1), 0]

    return w, QTF, Rext


def loadnemohRAO(foldername):
    filename = foldername + "/Motion/RAO.dat"

    M = np.genfromtxt(filename, skip_header=3)

    w = M[:, 0]
    RAOamp = M[:, 1:7]
    RAOphase = M[:, 7:]

    return w, RAOamp, RAOphase


def getQTF2D(win, QTF3D, w0):
    #win : the pulsations of the simulation
    #QTF3D: the QTF data at pulsations win
    #w0: the constant pulsation for which the line is required.

    nw = len(win)
    QTF2D = np.zeros([6, nw])
    for i in range(6):
        for j in range(nw):
            QTF2D[i, j] = np.interp(w0, win, QTF3D[i, j, :])

    return QTF2D


def read_Nemoh_cal(Npath):
    ncal = open(Npath + 'Nemoh.cal')
    ncal.readline()
    rho = float(ncal.readline().strip().split()[0])  # masse volumique
    G = float(ncal.readline().strip().split()[0])  # gravite
    H = float(ncal.readline().strip().split()[0])  # profondeur
    Peff = ncal.readline().strip().split()[
           :2]  # point de definition de la houle
    Xeff = float(Peff[0])
    Yeff = float(Peff[1])
    ncal.readline()
    nb = int(ncal.readline().strip().split()[0])  # nombre de corps
    namemb = []
    ndof = []
    nf = []
    nadd = []
    for b in range(nb):
        ncal.readline()
        namemb.append(ncal.readline().strip().split('\t')[
                          0])  # nom du fichier de maillage du corps b
        ncal.readline()
        ndof.append(int(ncal.readline().strip().split()[
                            0]))  # nombre de degre de liberte du corps b
        for i in range(ndof[b]):
            ncal.readline()
        nf.append(int(ncal.readline().strip().split()[
                          0]))  # nombre de composantes des efforts de
        # sortie pour le corps b
        for i in range(nf[b]):
            ncal.readline()
        nadd.append(int(ncal.readline().strip().split()[0]))
        for i in range(nadd[b]):
            ncal.readline()
    ncal.readline()
    wcal = ncal.readline().strip().split()
   
    wtype= int(wcal[0])
    nw   = int(wcal[1])
    Wmin = float(wcal[2])
    Wmax = float(wcal[3])
    w = np.linspace(Wmin, Wmax, nw)    # pulsations
    bcal = ncal.readline().strip().split()
    nbeta = int(bcal[0])
    bmin = float(bcal[1]) * np.pi / 180
    bmax = float(bcal[2]) * np.pi / 180
    beta = np.linspace(bmin, bmax, nbeta)  # directions
    ncal.close()
    return w, G, H, namemb
