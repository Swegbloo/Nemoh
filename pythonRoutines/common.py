#!/usr/bin/env python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
#
#   NEMOH - System parameters - May 2015
#   Contributors list:
#   - Pierre-Yves WUILLAUME, INNOSEA (pierreyves.wuillaume@innosea.fr)
#
#------------------------------------------------------------------------------

import os

import numpy as np

QTF_RESULTS_FOLDER = os.path.join('results', 'QTF')


def changename(path, name, ext):
    if os.path.exists(path + name + ext):
        print('attention fichier existant (' + name + ext + ')')
        name1 = name + '_old'
        changename(path, name1, ext)
        os.rename(path + name + ext, path + name1 + ext)
        print('fichier renomme en ' + name1 + ext)
    return


def createVTK(xyz, NN, filenameVTK):
    """Create a VTK mesh file
    Inputs:
      xyz: points list
      NN: faces list
      Optionnal:
        AddNormals (default = True): enter False if normals vectors should
        not be calculated and written in the VTK file

    """
    changename('', filenameVTK, '.vtk')
    f = open(filenameVTK + '.vtk', 'w')

    Nnodes = xyz.shape[0]
    nFaces = np.size(NN, 0)

    #header
    f.write('# vtk DataFile Version 2.0\n')
    f.write(
        'hydrodynamic vtk mesh for paraview vizu\nASCII\nDATASET POLYDATA\n')
    f.write('POINTS ' + str(Nnodes) + ' float\n')

    for l in range(Nnodes):
        f.write(str(xyz[l, 0]) + ' ' + str(xyz[l, 1]) + ' ' + str(
            xyz[l, 2]) + ' \n')

    nT = np.sum(NN[:, 0] == 3)
    nQ = np.sum(NN[:, 0] == 4)

    f.write('POLYGONS ' + str(nFaces) + ' ' + str(4 * nT + 5 * nQ) + ' \n')
    for l in range(nFaces):
        if NN[l, 0] == 2:
            f.write(str(NN[l, 0]) + ' ' + str(NN[l, 1]) + ' ' + str(
                NN[l, 2]) + ' ' + ' \n')
        elif NN[l, 0] == 3:
            f.write(str(NN[l, 0]) + ' ' + str(NN[l, 1]) + ' ' + str(
                NN[l, 2]) + ' ' + str(NN[l, 3]) + ' \n')
        elif NN[l, 0] == 4:
            f.write(str(NN[l, 0]) + ' ' + str(NN[l, 1]) + ' ' + str(
                NN[l, 2]) + ' ' + str(NN[l, 3]) + ' ' + str(NN[l, 4]) + ' \n')

    f.close()
    print '>> WriteMesh: createVTK :', filenameVTK + '.vtk', 'created.'

    return


def createL12(xyz, NN, filenameL12, sym=0):
    """Create a L12 mesh file
    Inputs:
      xyz: points list
      NN: faces list
      filenameL12 contains the path
      Optionnal:
        sym (default=0): input 1 if the mesh is symetric (if xyz,
        NN contains only half the points and faces)

    """
    changename('', filenameL12, '.dat')
    f = open(filenameL12 + '.dat', 'w')

    Nnodes = xyz.shape[0]
    nFaces = np.size(NN, 0)

    #header
    f.write('2 ' + str(sym) + ' \n')

    for l in range(Nnodes):
        f.write(str(l + 1) + ' ' + str(xyz[l, 0]) + ' ' + str(
            xyz[l, 1]) + ' ' + str(xyz[l, 2]) + ' \n')

    f.write('0 0 0 0\n')
    for l in range(nFaces):
        if NN[l, 0] == 2:
            f.write(str(NN[l, 1] + 1) + ' ' + str(NN[l, 2] + 1) + ' ' + str(
                NN[l, 2] + 1) + ' ' + str(NN[l, 2] + 1) + ' \n')
        elif NN[l, 0] == 3:
            f.write(str(NN[l, 1] + 1) + ' ' + str(NN[l, 2] + 1) + ' ' + str(
                NN[l, 3] + 1) + ' ' + str(NN[l, 3] + 1) + ' \n')
        elif NN[l, 0] == 4:
            f.write(str(NN[l, 1] + 1) + ' ' + str(NN[l, 2] + 1) + ' ' + str(
                NN[l, 3] + 1) + ' ' + str(NN[l, 4] + 1) + ' \n')
    f.write('0 0 0 0')
    f.close()
    print '>> WriteMesh: createL12 :', filenameL12 + '.dat', 'created.'

    return
