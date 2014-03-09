#!/bin/bash
################################################################################
# Copyright Keith Butler(2014) #
# #
# This file NewPotentialModule.py is free software: you can #
# redistribute it and/or modify it under the terms of the GNU General Public #
# License as published by the Free Software Foundation, either version 3 of #
# the License, or (at your option) any later version. #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for #
# more details. #
# You should have received a copy of the GNU General Public License along with #
# this program. If not, see <http://www.gnu.org/licenses/>. #
# #
################################################################################

import numpy
import numpy as np
import math


def macroscopic_average(potential,periodicity,resolution):
    """Getting the macroscopic average of potential
    Args:
         potential : array containig the electrostaticpotential/charge density
	 periodicity : real number; the period over which to average
	 resolution : the grid resolution in the direction of averaging
    Returns:
	 macro_average : array with the macroscopically averaged values"""

    macro_average = np.zeros(shape=(len(potential)))
    period_points = int((periodicity/resolution))
    print periodicity, resolution, period_points*periodicity
# Re-arrange so that period points divides evenly by resolution
    for i in range(len(potential)):
	for j in range(i-int(period_points/2),i+int(period_points/2)):
	    if j < 0:
	    	macro_average[i] = macro_average[i]+potential[j+len(potential)]
	    elif j >= len(potential):
	    	macro_average[i] = macro_average[i]+potential[j-len(potential)]
	    else:
	    	macro_average[i] = macro_average[i]+potential[j]
	macro_average[i] = macro_average[i]/period_points

    print ("Average of the average = ",numpy.average(macro_average))
    return macro_average

def cube_potential(origin,travelled,cube,Grid,nx,ny,nz):
    """Populates the sampling cube with the potential required"""

    potential_cube = np.zeros(shape=(cube[0],cube[1],cube[2]))
    for x in range(0,cube[0]):
        for y in range(0,cube[1]):
    	    for z in range(0,cube[2]):
# Assign the values of coordinates in the original grid
		xv = origin[0]+travelled[0]+x
		yv = origin[1]+travelled[1]+y
		zv = origin[2]+travelled[2]+z
# Minimum image convention
	    	zv = zv - nz*round(zv/nz)
	    	yv = yv - ny*round(yv/ny)
	    	xv = xv - nx*round(xv/nx)
        	potential_cube[x,y,z] = Grid[xv,yv,zv]

    return potential_cube.mean(), np.var(potential_cube)
#------------------------------------------------------------------------------------------

def cuboid_average(Grid,cube,origin,vector,nx,ny,nz,magnitude):
   """Calculates the average in a cube defined by size cube(a,b,c), beginning at origin and 
    reavelling as far as magnitude."""

   plotting_average = np.zeros(shape=(magnitude))
   i = 0
   while i < magnitude:
 	travelled = np.multiply(i, vector) 
    	plotting_average[i], varience = cube_potential(origin,travelled,cube,Grid,nx,ny,nz)
	i = i + 1

   return plotting_average 
#------------------------------------------------------------------------------------------

def planar_average(Grid,nx,ny,nz):
    """Calculate the average in a given plane for the full length of the normal;
    e.g. the full length of z in the xy plane."""
    axis = raw_input("Which axis do you wish to plot along?(x,y,z)LOWER CASE!! ")
    if axis == 'x':
	x_plane = np.zeros(shape=(ny,nz))
	Average = np.zeros(shape=(nx))
        for x_value in range(nx):
	    x_plane[:,:] = Grid[x_value,:,:]
            Average[x_value] = x_plane.mean()
    if axis == 'y':
	Average = np.zeros(shape=(ny))
	y_plane = np.zeros(shape=(nx,nz))
        for y_value in range(ny):
	    y_plane[:,:] = Grid[:,y_value,:]
            Average[y_value] = y_plane.mean()
    if axis == 'z':
	Average = np.zeros(shape=(nz))
	z_plane = np.zeros(shape=(nx,ny))
        for z_value in range(nz):
	    z_plane[:,:] = Grid[:,:,z_value]
            Average[z_value] = z_plane.mean()

    return Average
#------------------------------------------------------------------------------------------

def create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad):
    """Create the mesh of points for a contour plot"""
    xx, yy = np.mgrid[0:NGX,0:NGY]
    x = 0
    grd = np.zeros(shape=(NGX,NGY))
    while x < NGX:
        y = 0
        while y < NGY:
            z = 0
            while z < NGZ:
                z_value = (plane_coeff[3]-plane_coeff[0]*x-plane_coeff[1]*y)/plane_coeff[2]
		while z_value >= NGZ-1:
		    z_value = z_value - NGZ
                grd[x,y] = grad[x,y,z_value]
                z = z + 1
            y = y + 1
        x = x + 1
    return xx,yy,grd
#------------------------------------------------------------------------------------------

def numbers_2_grid(a,NGX,NGY,NGZ):
    """Takes a point (in fractional coordinates) and converts it to a VASP grid
    point based on the NGX/Y/Z values."""
    a_grid = np.zeros(shape=(3))
    a_grid[0] = round(float(a[0])*NGX)
    a_grid[1] = round(float(a[1])*NGY)
    a_grid[2] = round(float(a[2])*NGZ)

    return a_grid
#------------------------------------------------------------------------------------------

def matrix_2_abc(Lattice):
    """Thke the VASP lattice and convert to the a,b,c,alpha,beta,gamma format"""
    a = np.sqrt(Lattice[0,0]**2+Lattice[0,1]**2+Lattice[0,2]**2)
    b = np.sqrt(Lattice[1,0]**2+Lattice[1,1]**2+Lattice[1,2]**2)
    c = np.sqrt(Lattice[2,0]**2+Lattice[2,1]**2+Lattice[2,2]**2)

    return a,b,c
#------------------------------------------------------------------------------------------

def read_vasp_density(FILE):
    """Generic reading of CHGCAR LOCPOT etc files from VASP"""
    f = open(FILE,"r")
    lines = f.readlines()
    f.close()
# Get Header information
    i = 0
    lattice = np.zeros(shape=(3,3))
    for line in lines:
     inp = line.split()
     if inp == []:
      continue
     if len(inp) > 0:
      i = i+1
     if i == 2:
      scale_factor = float(inp[0])
     if i >= 3 and i < 6:
      lattice[i-3,:]=inp[:]
     if i == 6:
      num_species=len(inp)
      species=inp
     if i == 7:
      num_type=inp
      j = 0
      while (j < num_species):
       num_type[j-1] = int(num_type[j-1])
       j = j + 1
      num_atoms=sum(num_type)
     if i == 8:
      coord_type = inp
    
    for i in range(2):
     for j in range(2):
      lattice[i,j] = lattice[i,j]*scale_factor
# Restart reading to get the coordinates...it's just easier this way!
    i=0
    Coordinates = numpy.zeros(shape=(num_atoms,3))
    for line in lines:
     inp = line.split()
     if len(inp) > 0:
      i = i + 1
     if i >= 9 and i <= num_atoms+8 and len(inp) > 0:
      Coordinates[i-9,0] = float(inp[0])
      Coordinates[i-9,1] = float(inp[1])
      Coordinates[i-9,2] = float(inp[2])
# Now get the info about the charge grid
    i = 0
    for line in lines:
     inp = line.split()
     if len(inp) > 0:
      i = i + 1
     if i == num_atoms + 9:
      NGX = int(inp[0])
      NGY = int(inp[1])
      NGZ = int(inp[2])
      k = 0
      Potential = numpy.zeros(shape=(NGX*NGY*NGZ))
# Read in the potential data
     if i > num_atoms + 9 and i < num_atoms + 10 + NGX*NGY*NGZ/5:
      Potential[k]   = inp[0]
      Potential[k+1] = inp[1]
      Potential[k+2] = inp[2]
      Potential[k+3] = inp[3]
      Potential[k+4] = inp[4]
      k = k + 5
      if math.fmod(k,100000) == 0:
       print "Reading potetial, at point", k

    print	"BBBB		OOOO		OOOO		MMMMM	"
    print	"BBBB		OOOO		OOOO		MMMMM	"
    print	"BBBB		OOOO		OOOO		MMMMM	"
    print	"B  B	        OOOO		OOOO		MMMMM	"
    print	"B  B	        O  O		O  O		MMMMM	"
    print	"B  B	        O  O		O  O		MMMMM	"
    print	"B  B	        O  O		O  O		MMMMM	"
    print	"B  B	        O  O		O  O		MMMMM	"
    print	"BBBB	        O  O            O  O		M M M	"
    print	"BBBB	        O  O		O  O		M M M	"
    print	"BBBB	        O  O		O  O		M M M	"
    print	"B  B	        O  O		O  O		M M M	"
    print	"B  B	        O  O		O  O		M M M	"
    print	"B  B	        O  O		O  O		M M M	"
    print	"B  B	        O  O		O  O		M M M	"
    print	"B  B	        OOOO    	OOOO		M M M	"
    print	"BBBB            OOOO	        OOOO		M M M	"
    print	"BBBB            OOOO	        OOOO		M M M	"
    print	"BBBB            OOOO	        OOOO		M M M	"
 

    print ("Average of the potential = ",numpy.average(Potential))
    f.close()
    return Potential, NGX, NGY, NGZ, lattice
#------------------------------------------------------------------------------------------
def density_2_grid(Density,nx,ny,nz):
    """Convert the Potetnial list to a grid for ease of manipulation"""
    l = 0   
    Potential_grid = np.zeros(shape=(nx,ny,nz))
    for k in range(nz):
	for j in range(ny):
	    for i in range(nx):
		Potential_grid[i,j,k] = Density[l]
		l = l + 1

    return Potential_grid
#------------------------------------------------------------------------------------------

def points_2_plane(a,b,c):
    """define a plane based on 3 points
       method as outlined on http://www.had2know.com/academics/equation-plane-through-3-points.html
    """

    coefficients = np.zeros(shape=(4))

    ca = c - a
    ba = b - a
    normal = np.cross(ca,ba)
    d = normal[0]*a[0] + normal[1]*a[1] + normal[2]*a[2]
    for i in 0, 1, 2:
	coefficients[i] = normal[i]
    coefficients[3] = d
    return coefficients 

#------------------------------------------------------------------------------------------
def get_third_coordinate(plane_coeff,NGX,NGY):
    """Given the two sides of a plane, calculate the values of the 'plane' based on the
       equation eariler """

    zz = []
    i = j = 0
    while i <= NGX:
	i = i + 1
	j = 0
	while j <= NGY:
	    j = j + 1
	    rounded = round(((plane_coeff[0]*j+plane_coeff[1]*i) / plane_coeff[2]))
	    standard = ((plane_coeff[0]*j+plane_coeff[1]*i)/plane_coeff[2])
 	    if rounded == standard:   # Is it a whole number?
		zz.append(-(plane_coeff[0]*i+plane_coeff[1]*j)/plane_coeff[2])

    return zz
	 
#------------------------------------------------------------------------------------------
