#! /usr/bin/env python

# Copyright Keith Butler(2014) #
# #
# This file MacroDensity.beta_tools.py is free software: you can #
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
# Note all tools in here are development version only and do not get tested in 
# the build - TREAT WITH CAUTION.
################################################################################

import numpy
import numpy as np
import math
from scipy import interpolate

#------------------------------------------------------------------------------------------
def subs_potentials(A,B,tol):
    """Difference between two sets of data of the same length
    Args:
	A/B: the arrays (2D)
  	tol: the tolerence of the difference
    Returns:
	C: a new aaray (2D)
    """

    C = A
    for i in range(len(A)):
    	C[i,0] = A[i,0]
    	if abs(A[i,1] - B[i,1]) <= tol:
            C[i,1] = 0
    	else:
            C[i,1] = A[i,1] - B[i,1]

    return C

#------------------------------------------------------------------------------------------
def bulk_vac(bulk, slab):
    """ This sets the electron density to zero in the regions of the projected 'bulk' which correspond to the vacuum region of the slab calculation.
    Args:
	bulk : 2D array.
	vacuum : 2D array.
    Returns:
	new_bulk : 2D array with vacuum zeros included.""" 
    new_bulk = np.zeros(shape=(len(slab),2))
    i = -1
    for s_pot in slab:
	i = i + 1
	found = False
	for j in range(len(bulk)):
	    if s_pot[0] <= bulk[j,0] and s_pot[0] > bulk[j-1,0]:
		new_bulk[i,:] = bulk[j,:]
		found = True
     	if found == False:
	    new_bulk[i,0] = s_pot[0]
	    new_bulk[i,1] = 0

    return new_bulk
#------------------------------------------------------------------------------------------

def match_resolution(A,B):
    """Match the resolutions of two data-sets, given their range
    Args:
	A/B: two 2D arrays
    Returns:
	A_new/B_new : two new 2D arrays
    """

    np.append(A,A[0,:])
    np.append(B,B[0,:])
    resolution_a = (max(A[:,0])-min(A[:,0]))/len(A)
    resolution_b = (max(B[:,0])-min(B[:,0]))/len(B)
    new_resolution = min(resolution_a,resolution_b)/3
# Generate the function f for each spline
    f_a = interpolate.interp1d(A[:,0],A[:,1],kind='cubic')
    f_b = interpolate.interp1d(B[:,0],B[:,1],kind='cubic')
# Generate the new abscissa values, at new_resolution
    abscissa_a = np.arange(0,max(A[:,0]),new_resolution)
    abscissa_b = np.arange(0,max(B[:,0]),new_resolution)
# New datasets
    A_new = np.zeros(shape=(len(abscissa_a),2))
    B_new = np.zeros(shape=(len(abscissa_b),2))
    A_new[:,0] = abscissa_a
    B_new[:,0] = abscissa_b
    A_new[:,1] = f_a(abscissa_a)
    B_new[:,1] = f_b(abscissa_b)

    return A_new,B_new

#------------------------------------------------------------------------------------------
def spline_generate(A,new_res_factor):
    """Generate a spline of the data in a 2D array
    Args:
	A: 2D array
	new_res_factor: the factor by which to scale the resolution
    Returns:
	B: A spline of the data
    """
    resolution = (A[len(A)-1,0]-A[0,0])*new_res_factor/len(A)
    array_a = np.arange(min(A[:,0]),max(A[:,0]),resolution)
    f_a = interpolate.interp1d(A[:,0],A[:,1],kind='cubic')
    #ius = interpolate.InterpolatedUnivariateSpline(A[:,0],A[:,1])
    S = f_a(array_a)
    B = np.zeros(shape=(len(A)/new_res_factor,2))
    for i in range(len(B)):
    	B[i,0] = i*resolution+A[0,0]
    	B[i,1] = S[i]

    return B
#------------------------------------------------------------------------------------------

def matched_spline_generate(A,B, V_A, V_B):
    """Create 2D splines of 2 datasets, with an x-axis of units AA
    Args:
	A/B: The two datasets to match up 1-D arrays.
	V_A/B: The vectors of the direction of plotting.
    Returns:
	A/B_new : the new 2D Splined datasets.
    """

# Convert vectors to magnitude
    length_A = np.sqrt(V_A.dot(V_A))
    length_B = np.sqrt(V_B.dot(V_B))
# Determine the new resolution to plot at; twice the highest existing resolution
    res_a = length_A/(len(A))
    res_b = length_B/(len(B))
    new_resolution = (min(res_a,res_b))
# Create an array containing the indices of each potential point 0,1,2,....N
    array_a = np.arange(0,len(A))
    array_b = np.arange(0,len(B))
# Generate the function f for each spline
    f_a = interpolate.interp1d(array_a,A,kind='cubic')
    f_b = interpolate.interp1d(array_b,B,kind='cubic')
# Generate new arrays with the same resolution
    limits_a_new = np.arange(0,len(A))
    limits_b_new = np.arange(0,len(B))
# Make the arrays
    A_new = f_a(limits_a_new)
    B_new = f_b(limits_b_new)
# Convert to 2D arrays with AA in the first column
    TD_A = np.zeros(shape=(len(A_new),2))
    TD_B = np.zeros(shape=(len(B_new),2))
    res_a = length_A/float(len(A_new))
    res_b = length_B/float(len(B_new))
    for i in range(len(A_new)):
	TD_A[i,1] = A[i]
	TD_A[i,0] = i*res_a
    for i in range(len(B_new)):
	TD_B[i,1] = B[i]
	TD_B[i,0] = i*res_b
    return TD_A, TD_B

#------------------------------------------------------------------------------------------
def scissors_shift(potential,delta):
    """Scissors shifts a full potential by delta
    Args:
	potential: a 2D array
	delta: a real number
    Returns:
	shifted_potential: a 2D array
    """
    shifted_potential = potential
    for i in range(len(potential)):
	shifted_potential[i,0] = potential[i,0]
	shifted_potential[i,1] = potential[i,1] - delta

    return shifted_potential

#------------------------------------------------------------------------------------------
def extend_potential(potential,extension,vector):
    """Takes a potential and expands it.
    Args:
   	potential: 2D array, the potential to be expanded.
	extension: integer, the number of times to extend the potential.
 	vector: 1D array, the vector along which you re epanding.
    Returns:
	extended_potential: the 2D array, extended n times"""
    extended_potential = np.zeros(shape=(int(extension*len(potential)),2))
    idx = 0
    diff = np.sqrt(vector.dot(vector))
    increment = diff/len(potential[:,0])
    for i in range(int(extension)):
	for j in range(len(potential)):
    	    extended_potential[idx,0]=potential[j,0]+i*diff
	    extended_potential[idx,1] = potential[j,1]
	    idx = idx + 1

    if int(extension) != extension:            # For non-integer extensions
	i = i + 1
	over_shoot = extension - int(extension) 
	for j in range(int(len(potential)*over_shoot)):
   	    extended_potential[idx,0]=potential[j,0]+i*(max(potential[:,0])-min(potential[:,0]))+increment*i
	    extended_potential[idx,1] = potential[j,1]
            idx = idx + 1

    return extended_potential
#------------------------------------------------------------------------------------------
def sort_potential(potential):
    """I had to write my own method to sort a 2D arry by the first row...I don't know why
    Args:
 	potential: 2-D array.
    Returns:
	sorted_potential the same array sorted by the first column.
    """

    idx = sorted(potential[:,0])
    sorted_potential = potential.copy()
    for i in range(len(idx)):
        for j in range(len(potential)):
            if potential[j,0] == idx[i]:
                sorted_potential[i,0] = idx[i]
                sorted_potential[i,1] = potential[j,1]

    return sorted_potential
#------------------------------------------------------------------------------------------
def diff_potentials(potential_a, potential_b,start,end,tol=0.04):
    """ Gets the difference betweeen two potentials, returns a 1D array
    Args:
   	potential_a/b: 2D arrays
	start/end : the start and finish coordinates 
	tol : the tolerence for the coordinated being the same for subtraction
    Returns:
	new_potential: 2D array
    """
    resolution = potential_a[0,0] - potential_b[0,1]
    new_potential = np.zeros(shape=((start-end)/resolution,2))
    
    for i in range(len(potential_a)):
	if potential_a[i,0] >= start and potential_a[i,0] <= end:
	    for j in range(len(potential_b)):
		if abs(potential_b[j,0] - potential_a[i,0]) <= tol:
	    	    new_potential[i,1] = potential_a[i,1] - potential_b[i,1]
	    	    new_potential[i,0] = potential_a[i,0]	
	    

    return new_potential
#------------------------------------------------------------------------------------------
def translate_grid(potential, translation, periodic=False, vector=[0,0,0],boundary_shift=0.0):
    """Move the potential around, useful for overlapping two plots, or
    shifting a slab to the centre of the plot, this works if periodic=True

    The new index in the potential arrays [:,0] gives a coordinate in grid points.
    Args:
	potential : 2D array containig the electrostaticpotential/charge density
	translation : the distance to move it
	periodic : boolean, perform wrapping of coordinates
	vector: the vector along which you are transorming only required for periodic = True
	boundary_shift : real, number of AA to shift the location of the periodic boundary
    Returns:
	new_potential_trans : 2D array containig the translated electrostaticpotential/charge density 
    """
    new_potential_trans = np.zeros((len(potential),2))
    new_potential_orig = np.zeros((len(potential),2))
    length = np.sqrt(vector.dot(vector))


    for i in range(len(potential)):
    	new_potential_trans[i,0] = potential[i,0] + translation
	new_potential_trans[i,1] = potential[i,1]
        if periodic == True:
	    new_potential_trans[i,0] = new_potential_trans[i,0] - length*int((new_potential_trans[i,0]+boundary_shift)/length)

    if periodic == True:
# Sort the numbers out if you have done periodic wrapping
    	sorted_potential_trans = sort_potential(new_potential_trans)
    else:
	sorted_potential_trans = new_potential_trans

    #print sorted_potential_trans

    return sorted_potential_trans

#------------------------------------------------------------------------------------------

def create_plotting_mesh(NGX,NGY,NGZ,pc,grad):
    """Create the mesh of points for a contour plot
    pc: coefficients of the plane equation.
    """

    if pc[0] == 0 and pc[1] == 0: a = NGX; b = NGY; p = 'zzo'; c = int(pc[3] / pc[2]) - 1
    if pc[0] == 0 and pc[2] == 0: a = NGX; b = NGZ; p = 'zoz'; c = int(pc[3] / pc[1]) - 1
    if pc[1] == 0 and pc[2] == 0: a = NGY; b = NGZ; p = 'ozz'; c = int(pc[3] / pc[0]) - 1
    plane = np.zeros(shape=(a,b))
    for x in range(a):
        for y in range(b):
	    if p == 'zzo':
	        plane[x,y] = grad[x,y,c]

    return plane

#------------------------------------------------------------------------------------------

def read_cube_density(FILE):
    f = open(FILE,"r")
    lines = f.readlines()
    f.close()
    lattice = np.zeros(shape=(3,3))
    for line in lines:
	inp = line.split()
	if inp == []:
	    continue
	if len(inp) == 4:
	    natms = inp[0]
	    
#------------------------------------------------------------------------------------------

def points_2_plane(a,b,c):
    """define a plane based on 3 points
       method as outlined on http://www.had2know.com/academics/equation-plane-through-3-points.html
    """

    coefficients = np.zeros(shape=(4))

    ca = c - a
    ba = b - a
    normal = np.cross(ba,ca)
    #GCD_coeff = GCD_List(normal)
    #normal = [x / GCD_coeff for x in normal]
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
