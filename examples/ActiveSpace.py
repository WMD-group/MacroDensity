#! /usr/bin/env python

'''
Distinguish plateau regions in the electrostatic potential

Inputs:
cube_size = size of the cube in units of FFT mesh points (NGX/Y/Z)
cube_origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
tolerance = threshold below which the electrostatic potential is considered to be plateaued
input_file = VASP LOCPOT input filename to be read

Outputs:
Percentage of vaccum vs non-vacuum cubes
'''

import math
import numpy as np
from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, volume_average

## INPUT SECTION
cube_size = [2,2,2]
cube_origin = [0.5,0.5,0.5]
tolerance = 1E-4
input_file = 'LOCPOT'
## END INPUT SECTION

## GETTING POTENTIAL
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
cutoff_varience = tolerance
grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
travelled = [0,0,0]

## DISTNGUISHING VACCUM FROM NON_VACUUM
vacuum = []
non_vacuum = []
for i in range(0,NGX,cube_size[0]):
    for j in range(0,NGY,cube_size[1]):
        for k in range(0,NGZ,cube_size[2]):
            sub_origin = [float(i)/NGX,float(j)/NGY,float(k)/NGZ]
            cube_pot, cube_var = volume_average(origin=sub_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])
            if cube_var <= cutoff_varience:
                vacuum.append(sub_origin)
            else:
                non_vacuum.append(sub_origin)

## PRINTING
print("Number of vacuum cubes: ", len(vacuum))
print("Number of non-vacuum cubes: ", len(non_vacuum))
print("Percentage of vacuum cubes: ",(float(len(vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
print("Percentage of non-vacuum cubes: ",(float(len(non_vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
