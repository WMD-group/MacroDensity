#! /usr/bin/env python

'''
Electrostatic potential at atomic sites

Inputs:
potential_file = VASP LOCPOT
coordinate_file = The coordinates file NOTE This must be in vasp 4 format
species = The species whose on-site potential you are interested in (string)
sample_cube = The size of the sampling cube in units of mesh points (NGX/Y/Z)
output file = name of output data file
img_file = name of output image file

Outputs:
.png histogram output
.csv data output
'''
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import ase
import pandas as pd
from ase.io import write
from ase.io import vasp
from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average,volume_average

## INPUT SECTION
potential_file = 'LOCPOT')
coordinate_file = 'POSCAR')
species = "Zn"
sample_cube = [5,5,5]
output file = 'OnSitePotential.csv'
img_file = 'OnSitePotential.png'
## END INPUT SECTION

## GETTING POTENTIALS
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(potential_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
coords = ase.io.vasp.read_vasp(coordinate_file)
scaled_coords = coords.get_scaled_positions()
symbols = coords.get_chemical_symbols()
ox_coords = []

for i, atom in enumerate(coords):
    if symbols[i] == species:
        ox_coords.append(scaled_coords[i])
grid_position = np.zeros(shape=(3))
potentials_list = []
i = 0
num_bins = 20
for coord in ox_coords:
    i = i + 1
    grid_position[0] = coord[0]
    grid_position[1] = coord[1]
    grid_position[2] = coord[2]
    cube = sample_cube
    origin = [grid_position[0]-2,grid_position[1]-2,grid_position[2]-1]
    travelled = [0,0,0]
    cube_potential, cube_var = volume_average(origin,cube,grid_pot,NGX,NGY,NGZ)
    potentials_list.append(cube_potential)

## PLOTTING
n, bins, patches = plt.hist(potentials_list, num_bins, facecolor='#6400E1', alpha=0.5)
plt.xlabel('Hartree potential (V)')
plt.ylabel('% of centres')
plt.savefig(img_file)

## SAVING
df = pd.DataFrame.from_dict({'Potential':potentials_list},orient='index')
df = df.transpose()
df.to_csv(output_file)
