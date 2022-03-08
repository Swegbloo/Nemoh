#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#   NEMOH2 - second order (QTF) - January 2015
#   Contributors list:
#   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
#--------------------------------------------------------------------------------------


import numpy as np
from tools import computeFrequencyDomainMotionRAO, writeRAO
from loadnemoh import loadnemohB, loadnemohMa, loadnemohFe
from plot import *


def mechanicalSolver(workdir,IdMove,plot = True):
	
	# load hydrodynamic database
       	w, B = loadnemohB(workdir)
	w, Ma = loadnemohMa(workdir)
	w, Fea, Feph = loadnemohFe(workdir)
	Kh = np.genfromtxt(workdir+"Mechanics/Kh.dat")
	Km = np.genfromtxt(workdir+"Mechanics/Km.dat")
	Badd =np.genfromtxt(workdir+"Mechanics/Badd.dat")
	M =np.genfromtxt(workdir+"Mechanics/Inertia.dat")

	# solve motion
	Xmodule, Xphase = computeFrequencyDomainMotionRAO(w, Fea, Feph, M, Ma, Badd, B, Kh, Km)
	
	# plot RAO
	if plot:
		ifig = getnewfigure()
		plotXY(w,Xmodule*IdMove,handle=ifig)
		showplots()

	# write RAO file
	writeRAO(w,Xmodule*IdMove,Xphase*IdMove,workdir)

if __name__ == '__main__':
        IdMove=1 # 0 to make RAO=0
        f = open('workdir.txt', 'r')
        workdir=f.read()
        
        mechanicalSolver(workdir,IdMove)
