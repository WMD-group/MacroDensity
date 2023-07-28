#! /usr/bin/env python

'''
WARNING: THIS TOOL IS STILL UNDER DEVELOPMENT. KNOWN BUGS ARE PRESENT.
'''

'''
Inputs:
(define the plane with 3 points)
a_point = [0, 0, 0]
b_point = [1, 0, 1]
c_point = [0, 1, 0]
'''

import math
import numpy as np
import matplotlib.pyplot as plt
from macrodensity.density import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average, volume_average
from macrodensity.beta_tools import create_plotting_mesh,points_2_plane

#------------------------------------------------------------------
# Get the potential
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
cutoff_varience = tolerance
## Get the gradiens (Field), if required.
## Comment out if not required, due to compuational expense.
grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
## Convert the fractional points to grid points on the density surface
a = pot.numbers_2_grid(a_point,NGX,NGY,NGZ)
b = pot.numbers_2_grid(b_point,NGX,NGY,NGZ)
c = pot.numbers_2_grid(c_point,NGX,NGY,NGZ)
plane_coeff = pot.points_2_plane(a,b,c)

## Get the gradients
XY = np.multiply(grad_x,grad_y)
grad_mag = np.multiply(XY,grad_z)

## Create the plane
xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
## Plot the surface
plt.contourf(xx,yy,grd,V)
plt.show()

##------------------------------------------------------------------
## Plotting a planar average (Field/potential) throughout the material
##------------------------------------------------------------------
## FIELDS
planar = pot.planar_average(grad_x,NGX,NGY,NGZ)
## POTENTIAL
planar = pot.planar_average(grid_pot,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
macro  = pot.macroscopic_average(planar,4.80,resolution_z)
plt.plot(planar)
plt.plot(macro)
plt.savefig('Planar.eps')
plt.show()
