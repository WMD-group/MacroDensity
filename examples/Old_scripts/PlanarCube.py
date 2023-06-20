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
## POTENTIAL
planar = md.planar_average(potential,NGX,NGY,NGZ)
## MACROSCOPIC AVERAGE
macro  = md.macroscopic_average(planar,lattice_vector,resolution_z)
plt.plot(planar)
plt.plot(macro)
plt.savefig('Planar.eps')
np.savetxt(output_file,planar)
