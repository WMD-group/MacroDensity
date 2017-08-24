#! /usr/bin/env python
import NewPotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import izip


potential_file = 'LOCPOT' # The file with VASP output for potential
coordinate_file = 'POSCAR' # The coordinates file NOTE NOTE This must be in vasp 4 format 
species = "O"  # The species whose on-site potential you are interested in 
sample_cube = [5,5,5] # The size of the sampling cube in units of mesh points (NGX/Y/Z)

# Nothing below here should require changing
#------------------------------------------------------------------
# Get the potential
# This section should not be altered
#------------------------------------------------------------------
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(potential_file)
vector_a,vector_b,vector_c,av,bv,cv = pot.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ)
## Get the gradiens (Field), if required.
## Comment out if not required, due to compuational expense.
grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
#------------------------------------------------------------------

##------------------------------------------------------------------
## Getting the potentials for a group of atoms, in this case the Os
## NOTE THIS REQUIRES ASE to be available https://wiki.fysik.dtu.dk/ase/index.html
##------------------------------------------------------------------
##------------------------------------------------------------------
import ase                # Only add this if want to read in coordinates
from ase.io import write  # Only add this if want to read in coordinates
from ase.io import vasp   # Only add this if want to read in coordinates

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
    cube = sample_cube    # The size of the cube x,y,z in units of grid resolution.
    origin = [grid_position[0]-2,grid_position[1]-2,grid_position[2]-1]
    travelled = [0,0,0] # Should be left as it is.
    cube_potential, cube_var = pot.cube_potential(origin,travelled,cube,grid_pot,NGX,NGY,NGZ)
    potentials_list.append(cube_potential)
n, bins, patches = plt.hist(potentials_list, num_bins,normed=100, facecolor='#6400E1', alpha=0.5)
plt.xlabel('Hartree potential (V)',fontsize = 22)
plt.ylabel('% of centres',fontsize = 22)
plt.savefig('Potentials.png',dpi=300)
plt.show()
