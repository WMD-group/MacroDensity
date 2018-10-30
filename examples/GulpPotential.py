#! /usr/bin/env python
import sys
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt

print md.__file__

input_file = 'gulp.out'
md.read_gulp_potential(input_file)

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
grid_pot, electrons = md.density_2_grid(pot, NGX, NGY, NGZ)
#------------------------------------------------------------------
## POTENTIAL
planar = md.planar_average(grid_pot, NGX, NGY, NGZ)
## MACROSCOPIC AVERAGE
plt.plot(planar)
plt.savefig('Planar.eps')
plt.show()
