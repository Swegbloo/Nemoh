#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------------
#
#   NEMOH2 - second order (QTF) - January 2015
#   Contributors list:
#   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
#   - Jean-Noel Dory, INNOSEA (jeannoel.dory@gmail.com)          07/2015
# --------------------------------------------------------------------------------------


import cmath
import numpy as np
from plot import *
from common import *
from loadnemoh import *


def writeRAO(w, Xaw, Xphw, foldername):
    """
    w: pulsations
    Xaw: amplitude RAOamp
    Xphw: phase RAO
    foldername: nemoh computation folder
    """

    nw = len(w)

    changename(foldername + '/Motion/', 'RAO', '.dat')
    filename = foldername + "/Motion/RAO.dat"
    f = open(filename, 'w')

    f.write(
        '"VARIABLES="pulsation (rad/s)" "|X| (m/m)" "|Y| (m/m)" "|Z| (m/m)" "|phi| (deg)" "|theta| (deg)" "|psi| (deg)" "ang(x) (deg)" "ang(y) (deg)" "ang(z) (deg)" "ang(phi) (deg)" "ang(theta) (deg)" "ang(psi) (deg)"')
    f.write('\n')

    f.write('Zone t="Corps  1"')
    f.write('\n')

    # TODO: write beta value instead of "XXX"
    beta = "XXX"
    f.write('Zone t="beta =   ' + str(beta) + ' deg",I=    ' + str(nw) + ',F=POINT')
    f.write('\n')

    for i in range(nw):
        f.write(str(w[i]))
        f.write('\t')
        for k in range(6):
            if k < 3:
                f.write(str(Xaw[i, k]))
            else:
                f.write(str(Xaw[i, k] * 180.0 / math.pi))
            f.write('\t')
        
        for k in range(6):
            f.write(str(Xphw[i, k] * 180.0 / math.pi))
            f.write('\t')
        f.write('\n')
    f.close()

def writeQTFcontrib(w, QTFamp, QTFphase, contrib, qtftype, dof, foldername):
    nw = len(w)

    fileparts = [qtftype, contrib, 'DOF', str(dof), 'Amp.dat']
    filenameamp = os.path.join(foldername, '_'.join(fileparts))
    fileparts = [qtftype, contrib, 'DOF', str(dof), 'Phase.dat']
    filenamephase = os.path.join(foldername, '_'.join(fileparts))
    fa = open(filenameamp, 'w')
    fp = open(filenamephase, 'w')

    # first line
    fa.write("-1 ")
    fp.write("-1 ")

    for j in range(nw):
        fa.write(str(w[j]) + " ")
        fp.write(str(w[j]) + " ")

    fa.write("\n")
    fp.write("\n")

    for i in range(nw):
        fa.write(str(w[i]) + " ")
        fp.write(str(w[i]) + " ")
        for j in range(nw):
            fa.write(str(QTFamp[i, j]) + " ")
            fp.write(str(QTFphase[i, j]) + " ")
        fa.write("\n")
        fp.write("\n")

    fa.close()
    fp.close()

def convertComplex(QTF):
    # take module
    QTFAmp = np.abs(QTF)
    # take phase
    QTFPhase = np.angle(QTF, deg=True)
    # this is Aquaplus phase convention. Change it to Nemoh convention.
    QTFPhase += 90
    QTFPhase %= 360

    return QTFAmp, QTFPhase


def computeFrequencyDomainMotionRAO(w, modulesFe, phasesFe, M, Ma, Badd, Bw, Kh, Km=np.zeros((6, 6))):
    """Computes motion RAO from frequency domain results
    This function is made for 1 body with classical 6 DOF

    Inputs:
      w: pulsation values
      modulesFe(iw, i)
      phasesFe(iw, i):
        iw: pulsation index
        i: axis index
      M: mass matrix
      Ma(i,j,iw): added mass matrix
        i axis of measured force on given body
        j dof of motion
        iw: pulsation index
      Badd: Additional damping value
      Bw(i,j,w): damping matrix
        i axis of measured force on given body
        j dof of motion
        w the pulsation
      Kh: hydrostatic stiffness matrix
      Km: mooring stiffness matrix

    Outputs:
      Xmodule[iw, iDOF], Xphase[iw, iDOF]: phase and module of the RAO for each value of w
        iw: pulsation index
        iDOF: DOF index

    """

    nw = len(w)
    A = np.zeros((6, 6), dtype=complex)
    Fe = np.zeros((6, 1), dtype=complex)
    Xmodule = np.zeros((nw, 6))
    Xphase = np.zeros((nw, 6))

    for i in range(nw):
        Ma1 = Ma[:, :, i]
        B1w = Bw[:, :, i]
        w1 = float(w[i])
        # derivation in time domain = x -jw in frequency domain because in nemoh X = X0 exp(-jwt)
        A = -w1 * w1 * (M + Ma1) - 1j * w1 * (Badd + B1w) + (Kh + Km)

        Ai = np.linalg.inv(A)  # inverse matrix A

        for j in range(6):
            Fe[j] = modulesFe[i, j] * cmath.exp(1j * phasesFe[i, j])

        X = np.dot(Ai, Fe)

        for j in range(6):
            x = complex(float(X[j].real), float(X[j].imag))
            Xmodule[i, j] = abs(x)

            Xphase[i, j] = cmath.phase(x) % (2 * math.pi)

    return Xmodule, Xphase
