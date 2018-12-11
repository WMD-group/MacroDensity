#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import izip




#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density('LOCPOT.slab')
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
## Input section (define the plane with 3 points)
#a_point = [0, 0, 0]
#b_point = [1, 0, 1]
#c_point = [0, 1, 0]

## Convert the fractional points to grid points on the density surface
#a = md.numbers_2_grid(a_point,NGX,NGY,NGZ)
#b = md.numbers_2_grid(b_point,NGX,NGY,NGZ)
#c = md.numbers_2_grid(c_point,NGX,NGY,NGZ)
#plane_coeff = md.points_2_plane(a,b,c)

## Get the gradients
#XY = np.multiply(grad_x,grad_y)
#grad_mag = np.multiply(XY,grad_z)

## Create the plane
#xx,yy,grd =  md.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
## Plot the surface
#plt.contourf(xx,yy,grd,V)
#plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Plotting a planar average (Field/potential) throughout the material
##------------------------------------------------------------------
## FIELDS
#planar = md.planar_average(grad_x,NGX,NGY,NGZ)
## POTENTIAL
#planar = md.planar_average(grid_pot,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
#macro  = md.macroscopic_average(planar,4.80,resolution_z)
#plt.plot(planar)
#plt.plot(macro)
#plt.savefig('Planar.eps')
#plt.show()
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
# Getting the average potential in a single cube of arbitrary size
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
#cube = [2,2,2]
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
#origin = [0,0,0]
## travelled; do not alter this variable
#travelled = [0,0,0]
## Uncomment the lines below to do the business
#volume_average, cube_var = md.volume_average(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
#print "Potential            Variance"
#print "--------------------------------"
#print volume_average,"   ", cube_var
##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Plotting the average in a moving cube along a vector
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
cube = [2,2,2]
## vector is the vector you wish to travel along
vector = [1,1,0]
## cube defines the origin of the line in units of mesh points (NGX/Y/Z)
origin = [0.5,0,0.5]
## magnitude defines the length of the line, in units of mesh points (NGX/Y/Z)
magnitude = 280
## IF YOU WANT TO PLOT THE POTENTIAL:
cubes_potential = md.travelling_volume_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
abscissa = md.vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)
plt.plot(abscissa, cubes_potential)
plt.xlabel("$z (\AA)$")
plt.ylabel("Potential (eV)")
#plt.show()
#fp = open('ElectricPotential.dat','w')
#np.savetxt(fp,cubes_potential)

## IF YOU WANT TO PLOT THE FIELD MAGNITUDE ALSO:
## Get the gradients (of the field, if required)
grad_mag = md.gradient_magnitude(grad_x,grad_y,grad_z)
cubes_field = md.travelling_volume_average(grad_mag,cube,origin,vector,NGX,NGY,NGZ,magnitude)
abscissa = md.vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)
plt.plot(abscissa, cubes_field)
plt.xlabel("$z (\AA)$")
plt.ylabel("Field $(eV/\AA)$")
#plt.show()
with open('ElectricField.dat','wb')as f:
    writer = csv.writer(f)
    writer.writerows(izip(abscissa, cubes_field))
#ff = open('ElectricField.dat','w')
#np.savetxt(ff,cubes_field)

##------------------------------------------------------------------
##------------------------------------------------------------------

##------------------------------------------------------------------
## Getting the potentials for a group of atoms, in this case the Os
## NOTE THIS REQUIRES ASE to be available https://wiki.fysik.dtu.dk/ase/index.html
##------------------------------------------------------------------
##------------------------------------------------------------------
#import ase                # Only add this if want to read in coordinates
#from ase.io import write  # Only add this if want to read in coordinates
#from ase.io import vasp   # Only add this if want to read in coordinates

#coords = ase.io.vasp.read_vasp('POSCAR')
#scaled_coords = coords.get_scaled_positions()
#ox_coords = []
#i = -1
#for atom in coords:
#    i = i + 1
#    if atom.get_symbol() == "O":
#        ox_coords.append(scaled_coords[i])
#grid_position = np.zeros(shape=(3))
#potentials_list = []
#i = 0
#num_bins = 20
#for coord in ox_coords:
#    i = i + 1
#    grid_position[0] = int(coord[0]*NGX)
#    grid_position[1] = int(coord[1]*NGY)
#    grid_position[2] = int(coord[2]*NGZ)
#    cube = [5,5,3]    # The size of the cube x,y,z in units of grid resolution.
#    origin = [grid_position[0]-2,grid_position[1]-2,grid_position[2]-1]
#    travelled = [0,0,0] # Should be left as it is.
#    volume_average, cube_var = pot.volume_average(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
#    potentials_list.append(volume_average)
#n, bins, patches = plt.hist(potentials_list, num_bins,normed=100, facecolor='#6400E1', alpha=0.5)
#plt.xlabel('Hartree potential (V)',fontsize = 22)
#plt.ylabel('% of O centres',fontsize = 22)
#plt.savefig('Potentials.png',dpi=300)
#plt.show()
