#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt

"""
Calculates the Spherical Average around a given point.

Inputs:
cube_size = size of the cube in units of FFT mesh points (NGX/Y/Z)
cube_origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
input_file = VASP LOCPOT input filename to be read (DEFAULT = 'LOCPOT')

Outputs:
cube_potential, cube_variance (Terminal)
"""

### INPUT PARAMETERS
input_file = 'LOCPOT.MiL'
cube_size = [2,2,2]
cube_origin = [0,0,0]
### END INPUT PARAMETERS

# No need to alter anything after here
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)
cube = cube_size
origin = cube_origin
volume_average, cube_var = md.volume_average(origin, cube, grid_pot, NGX, NGY, NGZ)
print "Potential            Variance"
print "--------------------------------"
print volume_average,"   ", volume_var
