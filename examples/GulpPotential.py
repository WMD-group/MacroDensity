#! /usr/bin/env python

'''
Planar and macroscopic average with interpolation scheme for GULP outputs

Inputs:
lattice_vector = 3.0
input_file = Name of GULP input file
output_file = Name of output data file
img_file = Name of output image file
new_resolution = Total number of points for the interpolated planar avarage 

Outputs:
planar average, macroscopic average, interpolated planar average
.csv data file containing: planar average, macroscopic average, interpolated planar average
image file plotting .csv data
'''

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from macrodensity.density import matrix_2_abc, planar_average, macroscopic_average,density_2_grid_gulp,read_gulp_potential

## INPUT SECTION
lattice_vector = 3.0
input_file = 'gulp.out'
output_file = 'GulpPotential.csv'
img_file = 'GulpPotential.png'
new_resolution = 3000
## END INPUT SECTION

# GET POTENTIAL
pot, NGX, NGY, NGZ, Lattice = read_gulp_potential(input_file)
vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot = density_2_grid_gulp(pot, NGX, NGY, NGZ)

## POTENTIAL PLANAR AVERAGE
planar = planar_average(grid_pot, NGX, NGY, NGZ)
np.savetxt(output_file, planar)

## MACROSCOPIC AVERAGE
new_abscissa = np.linspace(0, NGZ - 1, new_resolution)
f = interp1d(range(NGZ), planar, kind='cubic')
interpolated_potential = [f(i) for i in new_abscissa]
macro  = macroscopic_average(planar, lattice_vector, vector_c/new_resolution)

## PLOTTING
fig, ax1 = plt.subplots()
ax1.set_ylabel('V/V')
ax1.set_xlabel('Grid Position')
ax1.plot(planar,linestyle = ' ',marker = 'o')
ax1.plot(macro)
ax2 = ax1.twiny()
ax2.plot(interpolated_potential)
plt.savefig(img_file)

## SAVING
df = pd.DataFrame.from_dict({'Planar':planar,'Macroscopic':macro,'Interpolated':interpolated_potential},orient='index')
df = df.transpose()
df.to_csv(output_file)
