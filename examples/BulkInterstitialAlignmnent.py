#! /usr/bin/env python

'''
Alignment of the band edges with the interstitial bulk reference potential.

Inputs:
intersices = Positions of the pores/interstices ([[interstice1],[interstice2],...])
outcar = VASP OUTCAR input filename
locpot = VASP LOCPOT input filename
cube_size = a cube defined by LOCPOT FFT mesh points

Output:
Aligned Valence Band, Aligned Conduction Band, Interstitial variances
'''

from macrodensity.density import read_vasp_density, matrix_2_abc, density_2_grid, volume_average
from macrodensity.vasp import get_band_extrema

## INPUT SECTION
interstices = ([0.5,0.5,0.5],[0.25,0.25,0.25])
outcar = 'OUTCAR'
locpot = 'LOCPOT'
cube_size = [2,2,2]
## END INPUT SECTION

## GETTING POTENTIAL
vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(locpot,quiet=True)
vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

## GETTING BAND EDGES
band_extrema = get_band_extrema(outcar)
VB_eigenvalue = band_extrema[0]
CB_eigenvalue = band_extrema[1]

## CALCULATING REFERENCE STATE
interstitial_potentials = []
interstitial_variances = []
for interstice in interstices:
    locpot_extract = volume_average(origin=interstice,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ)
    interstitial_potentials.append(locpot_extract[0])
    interstitial_variances.append(locpot_extract[1])

## CALCULATING ALIGNED BAND ENERGIES
sum_interstitial_potential = 0
for ele in interstitial_potentials:
    sum_interstitial_potential += ele
average_interstitial_potential = sum_interstitial_potential/len(interstitial_potentials)
VB_aligned = round(VB_eigenvalue - average_interstitial_potential,2)
CB_aligned = round(CB_eigenvalue - average_interstitial_potential,2)

## PRINTING
print("Reading band edges from file: "+str(outcar))
print("Reading potential from file: "+str(locpot))
print("Interstital variances: "+str(interstitial_variances))
print("VB_aligned      CB_aligned")
print("--------------------------------")
print(VB_aligned,"         ",CB_aligned)
