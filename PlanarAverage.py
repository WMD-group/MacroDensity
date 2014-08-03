#! /usr/bin/env python
import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt

input_file = 'LOCPOT.slab'
lattice_vector = 4.75
output_file = 'planar.dat'
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
## POTENTIAL
planar = pot.planar_average(grid_pot,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
macro  = pot.macroscopic_average(planar,lattice_vector,resolution_z)
plt.plot(planar)
plt.plot(macro)
plt.savefig('Planar.eps')
plt.show()
np.savetxt(output_file,planar)
##------------------------------------------------------------------
