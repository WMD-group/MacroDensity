#! /usr/bin/env python
import sys
sys.path.insert(0, '../')
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

print md.__file__

input_file = 'gulp.out'
output_file = 'planar.dat'
md.read_gulp_potential(input_file)
new_resolution = 3000
lattice_vector = 3.0

# No need to alter anything after here
#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
pot, NGX, NGY, NGZ, Lattice = md.read_gulp_potential(input_file)
vector_a, vector_b, vector_c, av, bv, cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot = md.density_2_grid_gulp(pot, NGX, NGY, NGZ)
#------------------------------------------------------------------
## POTENTIAL PLANAR AVERAGE
planar = md.planar_average(grid_pot, NGX, NGY, NGZ)
np.savetxt(output_file, planar)
## MACROSCOPIC AVERAGE
new_abscissa = np.linspace(0, NGZ - 1, new_resolution)
f = interp1d(range(NGZ), planar, kind='cubic')
interpolated_potential = [f(i) for i in new_abscissa]
#macro  = md.macroscopic_average(planar, lattice_vector, vector_c/new_resolution)
## PLOT
plt.plot(interpolated_potential)
plt.show()
plt.plot(planar)
plt.show()
