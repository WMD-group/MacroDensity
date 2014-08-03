#! /usr/bin/env python
import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt

input_file = 'LOCPOT.MiL'
cube_size = [2,2,2]    # This size is in units of mesh points
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
cube_origin = [0,0,0]
# No need to alter anything after here
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
#------------------------------------------------------------------

##------------------------------------------------------------------
# Getting the average potential in a single cube of arbitrary size
##------------------------------------------------------------------
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)
cube = cube_size
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
origin = cube_origin
## travelled; do not alter this variable
travelled = [0,0,0]
## Uncomment the lines below to do the business
cube_potential, cube_var = pot.cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
print "Potential            Variance"
print "--------------------------------"
print cube_potential,"   ", cube_var
##------------------------------------------------------------------
##------------------------------------------------------------------

