#! /usr/bin/env python
import macrodensity as md
import matplotlib.pyplot as plt

cube = [2, 2, 2]
vector = [1, 1, 0]
origin = [0.5, 0, 0.5]
magnitude = 280

# No need to alter anything after here
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density('LOCPOT.slab')
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)

## PLOTING
cubes_potential = md.cuboid_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
abscissa = md.vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)
plt.plot(abscissa, cubes_potential)
plt.xlabel("$z (\AA)$")
plt.ylabel("Potential (eV)")
