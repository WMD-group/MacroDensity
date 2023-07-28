#! /usr/bin/env python

'''
Planar and macroscopic average calculation

Inputs:
input_file = input filename to be read (must be in .cube format)
lattice_vector = Repeating unit over which the potential is averaged to get the macroscopic average (Angstroms)
output file = name of output data file
img_file = name of output image file

Outputs:
planar average, macroscopic average, interpolated planar average
.csv data file containing: planar average, macroscopic average, interpolated planar average
image file plotting .csv data
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from macrodensity.density import read_vasp_density, matrix_2_abc, density_2_grid, planar_average, macroscopic_average

## INPUT SECTION
input_file = 'LOCPOT'
lattice_vector = 5.41
output file = 'PlanarAverage.csv')
img_file = 'PlanarAverage.png')
## END INPUT SECTION

# GETTING POTENTIAL
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

## PLANAR AVERAGE
planar = planar_average(grid_pot,NGX,NGY,NGZ)

## MACROSCOPIC AVERAGE
macro = macroscopic_average(planar,lattice_vector,resolution_z)

## PLOTTING
plt.ylabel('V/V')
plt.xlabel('Grid Position')
plt.plot(planar)
plt.plot(macro)
plt.savefig(img_file)

## SAVING
df = pd.DataFrame.from_dict({'Planar':planar,'Macroscopic':macro},orient='index')
df = df.transpose()
df.to_csv(output_file)
