#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
import ase.io.cube

input_file = 'cube_002_total_density.cube'
lattice_vector = 4.75
output_file = 'planar.dat'
# No need to alter anything after here
#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
#vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(input_file)
#vector_a,vector_b,vector_c,av,bv,cv = pot.matrix_2_abc(Lattice)
#resolution_x = vector_a/NGX
#resolution_y = vector_b/NGY
#resolution_z = vector_c/NGZ
#grid_pot, electrons = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
potential, atoms = ase.io.cube.read_cube(input_file,read_data=True)
vector_a = np.linalg.norm(atoms.cell[1])
vector_b = np.linalg.norm(atoms.cell[1])
vector_c = np.linalg.norm(atoms.cell[2])
NGX = len(potential)
NGY = len(potential[0])
NGZ = len(potential[0][0])
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
print NGX,NGY,NGZ
#------------------------------------------------------------------
## POTENTIAL
planar = md.planar_average(potential,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
macro  = md.macroscopic_average(planar,lattice_vector,resolution_z)
plt.plot(planar)
plt.plot(macro)
plt.savefig('Planar.eps')
#plt.show()
np.savetxt(output_file,planar)
##------------------------------------------------------------------
