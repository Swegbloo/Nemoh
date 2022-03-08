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

    # tmp folder
    QTFtmppath = os.path.join(QTFpath, 'tmp')
    
    
    os.mkdir(QTFtmppath)

    pr = 'PR'
    pi = 'PI'
    results = {pr: [], pi: []}
    subresults = {pr: [], pi: []}

    for root, dirs, files in os.walk(QTFpath, topdown=True): #os.walk: generates the filenames in a directory tree 
     #   print "file : %s / %s / %s " % (root, dirs, files)
     #   print "------------------------------------"        
       
       #copying all output files to the tmp folder
        # if dirs != 'Appendix':
        if not dirs:
            # sub-contributions
            for fname in files:
                # print "fname : %s" %fname
                if os.path.splitext(fname)[0][-2:] == pi:
                    subresults[pi].append(fname)
                else:
                    subresults[pr].append(fname)
                filetocopy = os.path.join(root, fname)
                shutil.move(filetocopy, os.path.join(QTFtmppath, fname))
        else:
            # contributions
            for fname in files:
                if os.path.splitext(fname)[0][-2:] == pi:
                    results[pi].append(fname)
                else:
                    results[pr].append(fname)
                filetocopy = os.path.join(root, fname)
                shutil.move(filetocopy, os.path.join(QTFtmppath, fname))

    # write sub-contribs files
    for index, filename in enumerate(results[pi]):
        filetoread_PI = os.path.join(QTFtmppath, results[pi][index])
        filetoread_PR = os.path.join(QTFtmppath, results[pr][index])
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
            QTFtot = np.zeros((6, 2, nw, nw)) * 1j

        #  print "qtftype=%s " %qtftype
        if qtftype == 'QTFM':
            QTFtot[int(dof)-1, 0, :, :]  += QTF     #sum the complex QTF ie duok+Pink
        else:
            QTFtot[int(dof)-1, 1, :, :] += QTF

    for index, qtftype in enumerate(['QTFM', 'QTFP']):
        # take module
        QTFtotAmp = np.abs(QTFtot[:, index, :, :])  #amplitude of the QTFtot
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


