#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#--------------------------------------------------------------------------------------
#
#   NEMOH2 - second order (QTF) - November 2014
#   Contributors list:
#   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)
#--------------------------------------------------------------------------------------


import os
import numpy as np
import math
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

VarName = 'VarName'
VarValues = 'Values'
VarColor = 'Color'
XValues = 'XValues'
VarLinestyle = 'Linestyle'
VarMarker = 'Marker'
XValues = 'XValues'

# Other options of plots
MarkerABFe = 'default'
MarkerMultipleABFe = 'default'
FixedLineWidth = 2
FixedLineWidthEigenfrequencies = 1.5
SaveFormat = "png"
# SavePath
WidthFigure = 15
HeightFigure = 8
TitleSize = 20


def plotXY(X,Y,title = 'Title',  xlabel = 'x' , ylabel = 'y', yname = 'y', handle = 1, 
           save = False, show = False, FixedColor = 'b', savepath = './',
           marker='default',xspan=[float('nan'), float('nan')],yspan=[float('nan'), float('nan')], legloc=1, leg_bbox_to_anchor = 0):
  """ 
  Purpose : Plotting a figure
  
  Inputs:
  - X         : Plotted data for the x-axis
  - Y         : Plotted data for the y-axis
  - title     : Title of the figure
  - xlabel    : Label of the x-axis
  - ylabel    : Label of the y-axis
  - yname     : Name of the curve in the legend
  - handle    : ?
  - save      : Saving or not the figure
  - show      : Showing or not the figure
  - savepath  : Path to save the figure
  - marker    : Marker of the curve
  - xspan     : ?
  - yspan     : ?
  - legloc    : Localisation of the legend
                out of graph (left) = -1
                upper right : 1
                upper left : 2
                lower left : 3
                lower right : 4
                right : 5
                center left : 6
                center right : 7
                lower center : 8
                upper center : 9
                center : 10
  - leg_bbox_to_anchor : ? (left, bottom) or (left, bottom, width height) (Default: 0), outside location of the legend
  - Fixedcolor     : Color of the curve 
  """

  idx = np.argsort(X,0)
  X = np.sort(X,0)

  ncurves = 1

  fig = plt.figure(num = handle,figsize=(WidthFigure,HeightFigure))
  fig.my_title = title


  plt.hold(True)
  ynamei=yname

  xmin = xspan[0]
  if (math.isnan(xmin)):
      xmin = float(np.min(X))
      [xmin0, xmax0] = plt.xlim()
      xmin = min(xmin0,xmin)
  
  xmax = xspan[1]
  if (math.isnan(xmax)):
      xmax = float(np.max(X))
      [xmin0, xmax0] = plt.xlim()
      xmax = max(xmax0,xmax)

  Xc = X[np.nonzero((X>=xmin)*(X<=xmax))]


  for i in range(ncurves):
      Ys = Y[idx]
      Yc = Ys[np.nonzero((X>=xmin)*(X<=xmax))]
      
      if ncurves>1:
          ynamei = yname+"_"+str(i+1)
        
      if (ncurves==1) & (marker!='default'):
          plt.plot(Xc,Yc,marker,label=ynamei,color=FixedColor,linewidth = FixedLineWidth)
        
      else:
          plt.plot(Xc,Yc,label=ynamei,color=FixedColor,linewidth = FixedLineWidth)

  ymin = yspan[0]
  if (math.isnan(ymin)):
    if Y.ndim ==1:
      ymin = float(np.min(Y[np.nonzero((X>=xmin)*(X<=xmax))]))
    else:
      ymin = float(np.min(Y[np.nonzero((X>=xmin)*(X<=xmax)),:]))
    [ymin0, ymax0] = plt.ylim()
    ymin = min(ymin0,ymin)
  
  ymax = yspan[1]
  if (math.isnan(ymax)):
    if Y.ndim ==1:
      ymax = float(np.max(Y[np.nonzero((X>=xmin)*(X<=xmax))]))
    else:
      ymax = float(np.max(Y[np.nonzero((X>=xmin)*(X<=xmax)),:]))
    [ymin0, ymax0] = plt.ylim()
    ymax = max(ymax0,ymax)

  yspan = ymax-ymin
  ymin = ymin-yspan*0.05
  ymax = ymax+yspan*0.05     

  plt.xlim(xmin,xmax)
  plt.ylim(ymin,ymax)

  if legloc>0:
      plt.legend(loc = legloc)
  else:
      plt.legend(bbox_to_anchor=(1.1, 1))
  plt.ylabel(ylabel)
  plt.title(title,y=1.04,size = TitleSize)
  plt.xlabel(xlabel)
  plt.grid(True)    

  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  
  return



def plot3D(X,Y,Z, label="label",handle=1,title="Title",xspan=None, yspan=None, grid=False):
  
  fig = plt.figure(handle)
  fig.my_title = title

  ax = fig.gca(projection='3d')
  idx = range(len(X))
  idy = range(len(Y))

  if xspan is not None:
  
      idx = np.nonzero( (X>=xspan[0])*(X<=xspan[1])  )[0]
      X = X[idx]
  
  if yspan is not None:
  
      idy = np.nonzero( (Y>=yspan[0])*(Y<=yspan[1]) )[0]
      Y = Y[idy]
  
  X0, Y0 = np.meshgrid(X, Y)
  Z = Z[idx[0]:(idx[-1]+1),idy[0]:(idy[-1]+1)]

  if(len(X)>10):
      ix = int(len(X)/10)
  else:
      ix = 1
  if(len(Y)>10):
      iy = int(len(X)/10)
  else:
      iy=1
  if grid:
    ax.plot_wireframe(X0,Y0,Z,label=label,rstride =ix,cstride=iy)
  else:
    ax.plot_surface(X0,Y0,Z,cmap=cm.coolwarm, label=label,rstride =ix,cstride=iy)

  plt.title(title)
 
  
  return

def showplots():
  """Shows the figures
  """
  plt.show()
  return


def closeplots(fig = "all"):
  """Closes the figures
  """
  plt.close(fig)
  return

  
def clearplots():
  plt.clf()
  return


def saveplots(path='./'):
  """Saves the figures in the specified directory
  """
  for i in plt.get_fignums():
    fig=plt.figure(i)
    if os.path.isdir(path)==False:
      os.makedirs(path)
    plt.savefig(path+"/"+fig.my_title+'.png', format='png')
    
  return


def getnewfigure():
  allfig = plt.get_fignums()
  
  if(len(allfig)==0):
    return 1
  else:
    i = allfig[-1]+1
  
  return i
  
