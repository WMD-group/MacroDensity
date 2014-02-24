import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt


## Input section (define the plane with 3 points)
#a_point = [0, 0, 0]
#b_point = [1, 0, 1]
#c_point = [0, 1, 0]



# Define the plane
#a = np.array([0,0,0])
#b = np.array([18,0,0])
#c = np.array([0,28,0])




# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density('LOCPOT')
vector_a,vector_b,vector_c = pot.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
print resolution_x,resolution_y,resolution_z
grid_pot = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
## Get the gradiens (Field), if required.
## Comment out if not required, due to compuational expense.
#grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
#------------------------------------------------------------------

# Convert the fractional points to grid points on the density surface
#a = pot.numbers_2_grid(a_point,NGX,NGY,NGZ)
#b = pot.numbers_2_grid(b_point,NGX,NGY,NGZ)
#c = pot.numbers_2_grid(c_point,NGX,NGY,NGZ)

##------------------------------------------------------------------
## Get the equation for the plane
## This is the section for plotting on a user defined plane; 
## uncomment commands if this is the option that you want.
##------------------------------------------------------------------
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
## Plotting a planar average throughout the material
##------------------------------------------------------------------
#planar = pot.planar_average(grad_x,NGX,NGY,NGZ)
planar = pot.planar_average(grid_pot,NGX,NGY,NGZ)
macro  = pot.macroscopic_average(planar,4.8,resolution_z)
plt.plot(planar)
plt.plot(macro)
plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
# Getting the average potential in a single cube of arbitrary size
##------------------------------------------------------------------
#cube = [2,2,2]
#origin = [1,1,1]
#travelled = [0,0,0]
#cube_potential, cube_var = pot.cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
#print "Potential            Variance"
#print "--------------------------------"
#print cube_potential,"   ", cube_var
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Plotting the average in a moving cube along a vector
##------------------------------------------------------------------
#cube = [2,2,2]
#origin = [0,0,0]
#vector = [1,1,0]
#origin = [140,0,140]
#vector = [0,1,0]
#magnitude = 280
#cubes_field = pot.cuboid_average(grad_mag,cube,origin,vector,NGX,NGY,NGZ,magnitude)
#plt.plot(cubes_field)
#plt.show()
#ff = open('ElectricField.dat','w')
#np.savetxt(ff,cubes_field)
#cubes_potential = pot.cuboid_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
#plt.plot(cubes_potential)
#plt.show()
#fp = open('ElectricPotential.dat','w')
#np.savetxt(fp,cubes_potential)
##------------------------------------------------------------------
##------------------------------------------------------------------
