#! /usr/bin/env python

'''
Electrostatic potential plot spanning a vector acros the unit cell

Inputs:
cube = size of the cube in units of FFT mesh points (NGX/Y/Z)
origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
vector = vector across which the unit cell is traversed (hkl convention)
magnitude = length travelled along the selected vector in units of FFT mesh points (NGX/Y/Z)
input_file = VASP LOCPOT input filename to be read
output file = name of output data file
img_file = name of output image file

Outputs:
averaged electrostatic potential for the set cube size (list)
.csv file containing the above data
.png file presenting the above data
'''
from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, vector_2_abscissa, travelling_volume_average
import matplotlib.pyplot as plt
import pandas as pd

## INPUT SECTION
cube = [1,1,1]
origin = [0.17,0.17,0.17]
vector = [1,1,1]
magnitude = 16
input_file = 'LOCPOT'
output file ='MovingCube.csv'
img_file ='MovingCube.png'
## END INPUT SECTION

## GETTING POTENTIAL
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
cubes_potential = travelling_volume_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
abscissa = vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)

## PLOTTING
plt.plot(abscissa, cubes_potential)
plt.xlabel("$z (\AA)$")
plt.ylabel("Potential (eV)")
plt.savefig(img_file)

##SAVING
df = pd.DataFrame.from_dict({'Potential':cubes_potential},orient='index')
df = df.transpose()
df.to_csv(output_file)
