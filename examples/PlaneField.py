#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt


## Input section (define the plane with 3 points, fractional coordinates)
a_point = [0, 0, 0]
b_point = [1, 0, 1]
c_point = [0, 1, 0]

input_file = 'LOCPOT.slab'

#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)
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
a = md.numbers_2_grid(a_point,NGX,NGY,NGZ)
b = md.numbers_2_grid(b_point,NGX,NGY,NGZ)
c = md.numbers_2_grid(c_point,NGX,NGY,NGZ)
plane_coeff = md.points_2_plane(a,b,c)

## Get the gradients
XY = np.multiply(grad_x,grad_y)
grad_mag = np.multiply(XY,grad_z)

## Create the plane
xx,yy,grd =  md.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
## Plot the surface
plt.contour(xx,yy,grd,1)
plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

