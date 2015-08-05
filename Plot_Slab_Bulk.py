import PotentialModule as pot
import math
import numpy as np
import matplotlib.pyplot as plt
import sys


### THE INPUT PARAMETERS

input_file_1='CHGCAR-Slab'
input_file_2='CHGCAR-Bulk'
absc_1 = 2
absc_2 = 2
### END INPUT PARAMETERS
#------------------------------------------------------------------
#   READING
# Get the two potentials and change them to a planar average.
# This section should not be altered
#------------------------------------------------------------------
## SLAB
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(input_file_1)
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = pot.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
Volume = pot.get_volume(vec_a,vec_b,vec_c)
grid_pot_slab, electrons_slab = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
# Save the lattce vectors for use later
Vector_A = [vec_a,vec_b,vec_c]
##----------------------------------------------------------------------------------
## CONVERT TO PLANAR DENSITIES
##----------------------------------------------------------------------------------
slab = pot.planar_average_charge(grid_pot_slab,NGX,NGY,NGZ,Vector_A)
slab = pot.add_abscissa(slab,np.linalg.norm(Vector_A[absc_1]))
np.savetxt('Slab.dat',slab)
# Save the length of this for later to re-adjust dz in the Areal density
#----------------------------------------------------------------------------------
## BULK
vasp_pot, NGX, NGY, NGZ, Lattice = pot.read_vasp_density(input_file_2)
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = pot.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
Volume = pot.get_volume(vec_a,vec_b,vec_c)
grid_pot_bulk, electrons_bulk = pot.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
# Save the lattce vectors for use later
Vector_B = [vec_a,vec_b,vec_c]
##----------------------------------------------------------------------------------
## CONVERT TO PLANAR DENSITIES
##----------------------------------------------------------------------------------
bulk = pot.planar_average_charge(grid_pot_bulk,NGX,NGY,NGZ,Vector_B)
# Save the length of this for later to re-adjust dz in the Areal density
#----------------------------------------------------------------------------------
# FINISHED READING
#----------------------------------------------------------------------------------
# Create 2D arrays with the length from the 1D arrays of charge
bulk = pot.add_abscissa(bulk,np.linalg.norm(Vector_B[absc_2]))
#----------------------------------------------------------------------------------
# MOVE THE SLAB TO AVOID PBC
#----------------------------------------------------------------------------------
plt.plot(slab[:,0],slab[:,1])
np.savetxt('Bulk.dat',bulk)
