import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt






#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density('LOCPOT')
vector_a,vector_b,vector_c = pot.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
## Get the gradiens (Field), if required.
## Comment out if not required, due to compuational expense.
#grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
#------------------------------------------------------------------


##------------------------------------------------------------------
## Get the equation for the plane
## This is the section for plotting on a user defined plane; 
## uncomment commands if this is the option that you want.
##------------------------------------------------------------------
## Input section (define the plane with 3 points)
#a_point = [0, 0, 0]
#b_point = [1, 0, 1]
#c_point = [0, 1, 0]

## Convert the fractional points to grid points on the density surface
#a = pot.numbers_2_grid(a_point,NGX,NGY,NGZ)
#b = pot.numbers_2_grid(b_point,NGX,NGY,NGZ)
#c = pot.numbers_2_grid(c_point,NGX,NGY,NGZ)
#plane_coeff = pot.points_2_plane(a,b,c)

## Get the gradients
#XY = np.multiply(grad_x,grad_y)
#grad_mag = np.multiply(XY,grad_z)

## Create the plane
#xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
## Plot the surface
#plt.contourf(xx,yy,grd,V)
#plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Plotting a planar average (Field/potential) throughout the material
##------------------------------------------------------------------
## FIELDS
#planar = pot.planar_average(grad_x,NGX,NGY,NGZ)
## POTENTIAL
#planar = pot.planar_average(grid_pot,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
#macro  = pot.macroscopic_average(planar,4.8,resolution_z)
#plt.plot(planar)
#plt.plot(macro)
#plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
# Getting the average potential in a single cube of arbitrary size
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
#cube = [2,2,2]
## origin defines the bottom left point of the cube the "0,0,0" point
#origin = [1,1,1]
## travelled; do not alter this variable
#travelled = [0,0,0]
## Uncomment the lines below to do the business
#cube_potential, cube_var = pot.cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
#print "Potential            Variance"
#print "--------------------------------"
#print cube_potential,"   ", cube_var
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Plotting the average in a moving cube along a vector
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
#cube = [2,2,2]
## vector is the vector you wish to travel along
#vector = [1,1,0]
## cube defines the origin of the line in units of mesh points (NGX/Y/Z)
#origin = [140,0,140]
## magnitude defines the length of the line, in units of mesh points (NGX/Y/Z)
#magnitude = 280
## IF YOU WANT TO PLOT THE POTENTIAL:
#cubes_potential = pot.cuboid_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
#plt.plot(cubes_potential)
#plt.show()
#fp = open('ElectricPotential.dat','w')
#np.savetxt(fp,cubes_potential)

## IF YOU WANT TO PLOT THE FIELD MAGNITUDE ALSO:
## Get the gradients (of the field, if required)
#XY = np.multiply(grad_x,grad_y)
#grad_mag = np.multiply(XY,grad_z)
#cubes_field = pot.cuboid_average(grad_mag,cube,origin,vector,NGX,NGY,NGZ,magnitude)
#plt.plot(cubes_field)
#plt.show()
#ff = open('ElectricField.dat','w')
#np.savetxt(ff,cubes_field)

##------------------------------------------------------------------
##------------------------------------------------------------------
