#! /usr/bin/env python

# FieldAtPoint.py - try to calculate the electric field (grad of potential) at arbitrary point
# Forked form PlaneField.py - JMF 2016-01-25 

import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt
import sys #for argv

## Input section (define the plane with 3 points, fractional coordinates)
a_point = [0, 0, 0]
b_point = [1, 0, 1]
c_point = [0, 1, 0]

#LOCPOT.CsPbI3_cubic  LOCPOT.CsPbI3_distorted  LOCPOT.MAPI_pseudocubic
#input_file = 'LOCPOT.CsPbI3_distorted'
input_file = sys.argv[1]
print("Input file ",input_file)

#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = pot.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
## Get the gradiens (Field), if required.
## Comment out if not required, due to compuational expense.
grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
#------------------------------------------------------------------


##------------------------------------------------------------------
## Get the equation for the plane
## This is the section for plotting on a user defined plane; 
## uncomment commands if this is the option that you want.
##------------------------------------------------------------------

## Convert the fractional points to grid points on the density surface
a = pot.numbers_2_grid(a_point,NGX,NGY,NGZ)
b = pot.numbers_2_grid(b_point,NGX,NGY,NGZ)
c = pot.numbers_2_grid(c_point,NGX,NGY,NGZ)
plane_coeff = pot.points_2_plane(a,b,c)

## Get the gradients
XY = np.multiply(grad_x,grad_y)
grad_mag = np.multiply(XY,grad_z)

## Create the plane
xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_y)
## Plot the surface
plt.contour(xx,yy,grd,20)
plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

# CsPbI3 - distorted
PB_X=0.469972*NGX
PB_Y=0.530081*NGY
PB_Z=0.468559*NGZ

# CsPbI3 - perfect cubic
#PB_X=0.5*NGX
#PB_Y=0.5*NGY
#PB_Z=0.5*NGZ

# MAPBI3 - pseudo cubic distorted
#PB_X=0.476171*NGX
#PB_Y=0.500031*NGY
#PB_Z=0.475647*NGZ

# Read out massive grad table, in {x,y,z} components
print(grad_x[PB_X][PB_Y][PB_Z],grad_y[PB_X][PB_Y][PB_Z],grad_z[PB_X][PB_Y][PB_Z])
# Norm of electric field at this point
print(np.linalg.norm([ grad_x[PB_X][PB_Y][PB_Z],grad_y[PB_X][PB_Y][PB_Z],grad_z[PB_X][PB_Y][PB_Z] ]))

# OK, let's try this with a Spectral method (FFT)
# JMF - Not currently working; unsure of data formats, need worked example
from scipy import fftpack
V_FFT=fftpack.fftn(grid_pot[:,:,:])
V_deriv=fftpack.diff(grid_pot[:,:,:]) #V_FFT,order=1)

# Standard catch all to drop into ipython at end of script for variable inspection etc.
from IPython import embed; embed() # End on an interactive ipython console to inspect variables etc.

