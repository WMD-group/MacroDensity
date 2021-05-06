#! /usr/bin/env python
import macrodensity as md
import matplotlib.pyplot as plt

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

##------------------------------------------------------------------
## Plotting the average in a moving cube along a vector
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
cube = [2, 2, 2]
## vector is the vector you wish to travel along
vector = [1, 1, 0]
## cube defines the origin of the line in units of mesh points (NGX/Y/Z)
origin = [0.5, 0, 0.5]
## magnitude defines the length of the line, in units of mesh points (NGX/Y/Z)
magnitude = 280
## IF YOU WANT TO PLOT THE POTENTIAL:
cubes_potential = md.cuboid_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
abscissa = md.vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)
plt.plot(abscissa, cubes_potential)
plt.xlabel("$z (\AA)$")
plt.ylabel("Potential (eV)")
