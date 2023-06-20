#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import izip

cube = [2,2,2]
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
origin = [0,0,0]
cutoff_varience = 1E-4



vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density('LOCPOT.slab')
vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)

vacuum = []
non_vacuum = []

for i in range(0,NGX,cube[0]):
    print float(i)/NGX
    for j in range(0,NGY,cube[1]):
	for k in range(0,NGZ,cube[2]):
	    origin = [float(i)/NGX,float(j)/NGY,float(k)/NGZ]
            volume_average, cube_var = md.voulme_average(origin, cube, grid_pot, NGX, NGY, NGZ)
	    if cube_var <= cutoff_varience:
		vacuum.append(origin)
	    else:
		non_vacuum.append(origin)
print "Number of vacuum cubes: ", len(vacuum)
print "Number of non-vacuum cubes: ", len(non_vacuum)

print "Percentage of vacuum cubes: ",(float(len(vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.)
print "Percentage of non-vacuum cubes: ",(float(len(non_vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.)
