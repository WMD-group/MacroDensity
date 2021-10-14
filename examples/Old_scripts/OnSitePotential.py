#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import izip
import ase
from ase.io import write
from ase.io import vasp


potential_file = 'LOCPOT'
coordinate_file = 'POSCAR'
species = "O"
sample_cube = [5,5,5]


# Nothing below here should require changing

vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density(potential_file)
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)
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
    volume_average, cube_var = md.volume_average(origin, cube, grid_pot, NGX, NGY, NGZ)
    potentials_list.append(volume_average)
n, bins, patches = plt.hist(potentials_list, num_bins,normed=100, facecolor='#6400E1', alpha=0.5)
plt.xlabel('Hartree potential (V)',fontsize = 22)
plt.ylabel('% of centres',fontsize = 22)
plt.savefig('Potentials.png',dpi=300)
plt.show()