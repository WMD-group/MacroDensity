#! /usr/bin/env python

"""
Calculates the Spherical Average around a given point.

Inputs:
cube_size = size of the cube in units of FFT mesh points (NGX/Y/Z)
cube_origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
input_file = VASP LOCPOT input filename to be read

Outputs:
cube_potential, cube_variance (Terminal)
"""
from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, volume_average

## INPUT SECTION
cube_size = [2,2,2]
cube_origin = [0.5,0.5,0.5]
input_file = 'LOCPOT')
## END INPUT SECTION

## GETTING POTENTIAL
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

cube = cube_size
origin = cube_origin
travelled = [0,0,0]
cube_pot, cube_var = volume_average(origin=cube_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])

## PRINTING
print("Potential            Variance")
print("--------------------------------")
print(cube_pot,"   ", cube_var)
