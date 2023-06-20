#! /usr/bin/env python

'''
Planar and macroscopic average for cube files

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
import math
import numpy as np
import matplotlib.pyplot as plt
import ase.io.cube
import pandas as pd
from macrodensity.density_tools import planar_average, macroscopic_average

## INPUT SECTION
input_file = "cube_001_spin_density.cube"
lattice_vector = 4.75
output file = 'PlanarCube.csv')
img_file = 'PlanarCube.png')
## END INPUT SECTION

# GETTING POTENTIAL
potential, atoms = ase.io.cube.read_cube_data(input_file)
vector_a = np.linalg.norm(atoms.cell[1])
vector_b = np.linalg.norm(atoms.cell[1])
vector_c = np.linalg.norm(atoms.cell[2])
NGX = len(potential)
NGY = len(potential[0])
NGZ = len(potential[0][0])
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ

## PLANAR AVERAGE
planar = planar_average(potential,NGX,NGY,NGZ)

## MACROSCOPIC AVERAGE
macro  = macroscopic_average(planar,lattice_vector,resolution_z)

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
