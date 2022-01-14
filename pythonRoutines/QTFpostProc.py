#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
#
#   NEMOH2 - second order (QTF) - November 2014
#   Contributors list:
#   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
# -----------------------------------------------------------------------------

"""
this script performs the post processing of QTFSolver
it should be translated into Fortran to build an executable "QTFpostProc"
following Nemoh's fashion
"""

import os
import numpy as np
import shutil

from common import QTF_RESULTS_FOLDER
from loadnemoh import loadnemohQTFfile
from tools import writeQTFcontrib, convertComplex


def QTFpostProc(nemohfolder):
    QTFpath = os.path.join(nemohfolder, QTF_RESULTS_FOLDER)
    print QTFpath
    # tmp folder
    # QTFtmppath = os.path.join(QTFpath, 'tmp')

    #os.mkdir(QTFtmppath)

    pr = 'PR'
    pi = 'PI'
    results = {pr: [], pi: []}
    subresults = {pr: [], pi: []}

    for root, dirs, files in os.walk(QTFpath, topdown=True):
        # if dirs != 'Appendix':
        # print files
        # print root
        if not dirs:
            # sub-contributions
            for fname in files:
                if os.path.splitext(fname)[0][-2:] == pi:
                    subresults[pi].append(fname)
                else:
                    subresults[pr].append(fname)
                #filetocopy = os.path.join(root, fname)
                #shutil.move(filetocopy, os.path.join(QTFtmppath, fname))
        else:
            # contributions
            for fname in files:

                if os.path.splitext(fname)[0][-3:] == 'Amp':
                    continue
                if os.path.splitext(fname)[0][-5:] == 'Phase':
                    continue

                print(fname)
                if os.path.splitext(fname)[0][-2:] == pi:
                    results[pi].append(fname)
                else:
                    results[pr].append(fname)
               # filetocopy = os.path.join(root, fname)
               # shutil.move(filetocopy, os.path.join(QTFtmppath, fname))

    # write sub-contribs files
    for index, filename in enumerate(results[pi]):
        filetoread_PI = os.path.join(QTFpath, results[pi][index])
        filetoread_PR = os.path.join(QTFpath, results[pr][index])
        print results[pr][index]
        print filetoread_PR
        print os.path.exists(filetoread_PR)

        w, QTFCR = loadnemohQTFfile(filetoread_PR)
        w, QTFCI = loadnemohQTFfile(filetoread_PI)

        # get properties
        dof = filename.split('_')[3]
        qtftype = filename.split('_')[0]
        contrib = filename.split('_')[1]
        QTF = QTFCR + 1j * QTFCI
        QTFAmp, QTFPhase = convertComplex(QTF)

        writeQTFcontrib(w, QTFAmp, QTFPhase, contrib, qtftype, dof, QTFpath)

        if index == 0:
            nw = len(w)
            QTFtot = np.zeros((6, 3, nw, nw)) * 1j

        if qtftype == 'QTFM':
            if contrib == 'HASBO':
                QTFtot[int(dof)-1, 0, :, :]  += QTF  #DUOK+HASBO
            elif contrib == 'APINK':
                QTFtot[int(dof)-1, 1, :, :]  += QTF  #DUOK+PINKSTER
            else:
                QTFtot[int(dof)-1, 0, :, :]  += QTF  #DUOK+HASBO+HASFS+ASYMP
                QTFtot[int(dof)-1, 1, :, :]  += QTF  #DUOK+PINKSTER+HASFS+ASYMP
        else:
             QTFtot[int(dof)-1, 2, :, :] += QTF

    for index, qtftype in enumerate(['QTFM1','QTFM2','QTFP']):
        # take module
        QTFtotAmp = np.abs(QTFtot[:, index, :, :])
        # take phase
        QTFtotPhase = np.angle(QTFtot[:, index, :, :], deg=True)
        # this is Aquaplus phase convention. Change it to Nemoh convention.
        QTFtotPhase += 90
        QTFtotPhase %= 360
        
        for dof in range(6):
            writeQTFcontrib(w, QTFtotAmp[dof, :, :], QTFtotPhase[dof, :, :], 'TOT', qtftype, dof+1, QTFpath)
            
if __name__== '__main__':
    f = open('workdir.txt', 'r')
    workdir=f.read()
    QTFpostProc(workdir)
