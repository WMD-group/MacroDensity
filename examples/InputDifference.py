#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt


'''
WARNING: THIS TOOL IS STILL UNDER DEVELOPMENT. KNOWN BUGS ARE PRESENT.
'''


#------------------------------------------------------------------
#   READING
# Get the two potentials and change them to a planar average.
# This section should not be altered
#------------------------------------------------------------------
# SLAB
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density('CHGCAR.Slab')
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = md.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
Volume = md.get_volume(vec_a,vec_b,vec_c)
grid_pot_slab, electrons_slab = md.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
# Save the lattce vectors for use later
Vector_A = [vec_a,vec_b,vec_c]
#----------------------------------------------------------------------------------
# CONVERT TO PLANAR DENSITIES
#----------------------------------------------------------------------------------
planar_slab = md.planar_average(grid_pot_slab,NGX,NGY,NGZ)
# BULK
vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density('CHGCAR.Bulk')
mag_a,mag_b,mag_c,vec_a,vec_b,vec_c = md.matrix_2_abc(Lattice)
resolution_x = mag_a/NGX
resolution_y = mag_b/NGY
resolution_z = mag_c/NGZ
# Save the lattce vectors for use later
Vector_B = [vec_a,vec_b,vec_c]
Volume = md.get_volume(vec_a,vec_b,vec_c)
#----------------------------------------------------------------------------------
# CONVERT TO PLANAR DENSITIES
#----------------------------------------------------------------------------------
grid_pot_bulk, electrons_bulk = md.density_2_grid(vasp_pot,NGX,NGY,NGZ,True,Volume)
planar_bulk = md.planar_average(grid_pot_bulk,NGX,NGY,NGZ)
#----------------------------------------------------------------------------------
# FINISHED READING
#----------------------------------------------------------------------------------
# GET RATIO OF NUMBERS OF ELECTRONS
#----------------------------------------------------------------------------------
elect_ratio = int(electrons_slab/electrons_bulk)
#----------------------------------------------------------------------------------
# SPLINE THE TWO GENERATING A DISTANCE ON ABSCISSA
#----------------------------------------------------------------------------------
slab, bulk = md.matched_spline_generate(planar_slab,planar_bulk,Vector_A[2],Vector_B[1])
#----------------------------------------------------------------------------------
# EXTEND THE BULK POTENTIAL TO MATCH THE ELECTRON NUMBERS
#----------------------------------------------------------------------------------
bulk = md.extend_potential(bulk,elect_ratio,Vector_B[1])
#----------------------------------------------------------------------------------
# MATCH THE RESOLUTIONS OF THE TWO
#----------------------------------------------------------------------------------
slab, bulk = md.match_resolution(slab, bulk)
plt.plot(bulk[:,1])
plt.show()
#----------------------------------------------------------------------------------
# TRANSLATE THE BULK POTENTIAL TO GET OVERLAP
#----------------------------------------------------------------------------------
bulk_trans = md.translate_grid(bulk, 3.13,True,np.dot(Vector_B[1],elect_ratio),0.42)
bulk_trans = md.translate_grid(bulk_trans, 6.57,False,np.dot(Vector_B[1],elect_ratio),0.42)
slab_trans = md.translate_grid(slab, 6.5653,True,Vector_A[2])
#potential_difference = md.diff_potentials(pot_slab_orig,bulk_extd,10,40,tol=0.04)

##------------------------------------------------------------------
## SET THE CHARGE DENSITY TO ZERO OUTSIDE THE BULK
bulk_vacuum = md.bulk_vac(bulk_trans, slab_trans)
plt.plot(bulk_vacuum[:,0],bulk_vacuum[:,1])
plt.plot(slab_trans[:,0],slab_trans[:,1],)
plt.show()
# GET THE DIFFERENCES (within a numerical tolerence)
difference = md.subs_potentials(slab_trans,bulk_vacuum,tol=0.01)
difference = md.spline_generate(difference)
plt.plot(difference[:,0],difference[:,1])

#plt.plot(difference[:,0],difference[:,1])
plt.show()
