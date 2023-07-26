###############################################################################
# Copyright Keith Butler(2014)                                                #
#                                                                             #
# This file MacroDensity.beta_tools.py is free software: you can              #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.                         #
# This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                               #
# You should have received a copy of the GNU General Public License along with#
# this program. If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                             #
# Note all tools in here are development version only and do not get tested in#
# the build - TREAT WITH CAUTION.                                             #
###############################################################################

from __future__ import division
import numpy as np
from scipy import interpolate

#------------------------------------------------------------------------------
def subs_potentials(A,B,tol):
    """
    Subtract potentials between two datasets based on a tolerance value.

    Parameters:
        A (numpy.ndarray): The first dataset containing potential values in the format (x, potential).
        B (numpy.ndarray): The second dataset containing potential values in the format (x, potential).
        tol (float): The tolerance value for potential subtraction.

    Returns:
        C (numpy.ndarray): The resulting dataset containing the subtracted potentials in the format (x, potential).

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> B = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> tolerance = 1e-2
        >>> C = subs_potentials(A, B, tolerance)
        >>> print(C)
    """
    C = A
    for i in range(len(A)):
        C[i,0] = A[i,0]
        if abs(A[i,1] - B[i,1]) <= tol:
            C[i,1] = 0
        else:
            C[i,1] = A[i,1] - B[i,1]

    return C

#------------------------------------------------------------------------------
def bulk_vac(bulk, slab):
    """
    Subtract potentials between a bulk dataset and a slab dataset based on their positions.

    Parameters:
        bulk (numpy.ndarray): The dataset containing bulk potential values in the format (x, potential).
        slab (numpy.ndarray): The dataset containing slab potential values in the format (x, potential).

    Returns:
        new_bulk (numpy.ndarray): The resulting bulk dataset with matching positions of the slab dataset.

    Example:
        >>> bulk = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> slab = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> result = bulk_vac(bulk, slab)
        >>> print(result)

    """
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
#------------------------------------------------------------------------------

def match_resolution(A,B):
    """
    Match the resolutions of two datasets by cubic spline interpolation.

    Parameters:
        A (numpy.ndarray): The first dataset containing potential values in the format (x, potential).
        B (numpy.ndarray): The second dataset containing potential values in the format (x, potential).

    Returns:
        A_new (numpy.ndarray): The first dataset with matched resolution and interpolated values.
        B_new (numpy.ndarray): The second dataset with matched resolution and interpolated values.

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> B = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> result_A, result_B = match_resolution(A, B)
        >>> print(result_A)
        >>> print(result_B)
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

    return A_new, B_new

#------------------------------------------------------------------------------
def spline_generate(A,new_res_factor):
    """
    Generate a new dataset with higher resolution using cubic spline interpolation.

    Parameters:
        A (numpy.ndarray): The dataset containing potential values in the format (x, potential).
        new_res_factor (float): The factor by which to increase the resolution.

    Returns:
        B (numpy.ndarray): The new dataset with higher resolution and interpolated values.

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> new_res_factor = 2
        >>> result = spline_generate(A, new_res_factor)
        >>> print(result)
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
#------------------------------------------------------------------------------

def matched_spline_generate(A,B, V_A, V_B):
    """
    Generate matched datasets with the same resolution using cubic spline interpolation.

    Parameters:
        A (numpy.ndarray): The first dataset containing potential values.
        B (numpy.ndarray): The second dataset containing potential values.
        V_A (numpy.ndarray): Vector information for the first dataset.
        V_B (numpy.ndarray): Vector information for the second dataset.

    Returns:
        TD_A (numpy.ndarray): The first dataset with matched resolution and interpolated values.
        TD_B (numpy.ndarray): The second dataset with matched resolution and interpolated values.

    Example:
        >>> A = np.array([1, 2, 3, 4])
        >>> B = np.array([5, 6, 7, 8])
        >>> V_A = np.array([1, 2, 3])
        >>> V_B = np.array([4, 5, 6])
        >>> result_TD_A, result_TD_B = matched_spline_generate(A, B, V_A, V_B)
        >>> print(result_TD_A)
        >>> print(result_TD_B)

    """
    # Convert vectors to magnitude
    length_A = np.sqrt(V_A.dot(V_A))
    length_B = np.sqrt(V_B.dot(V_B))
    # Determine new resolution to plot at; twice highest existing resolution
    res_a = length_A/(len(A))
    res_b = length_B/(len(B))
    new_resolution = (min(res_a,res_b))
    # Create an array containing indices of each potential point 0,1,2,....N
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

    #--------------------------------------------------------------------------
def scissors_shift(potential,delta):
    """
    Shift the potential values by a constant amount.

    Parameters:
        potential (numpy.ndarray): The dataset containing potential values in the format (x, potential).
        delta (float): The constant value to shift the potentials.

    Returns:
        shifted_potential (numpy.ndarray): The resulting dataset with shifted potential values.

    Example:
        >>> potential = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> delta = 0.5
        >>> result = scissors_shift(potential, delta)
        >>> print(result)

    """
    shifted_potential = potential
    for i in range(len(potential)):
        shifted_potential[i,0] = potential[i,0]
        shifted_potential[i,1] = potential[i,1] - delta

    return shifted_potential

#------------------------------------------------------------------------------
def extend_potential(potential,extension,vector):
    """
    Extend a dataset by duplicating potential values along a specified vector direction.

    Parameters:
        potential (numpy.ndarray): The dataset containing potential values in the format (x, potential).
        extension (float): The extension factor specifying how many times to extend the dataset.
        vector (list): The vector specifying the direction along which to extend the dataset.

    Returns:
        extended_potential (numpy.ndarray): The resulting dataset with extended potential values.

    Example:
        >>> potential = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> extension = 2
        >>> vector = np.array([1, 1, 1])
        >>> result = extend_potential(potential, extension, vector)
        >>> print(result)

    """
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
            extended_potential[idx,0] = (potential[j,0] +
                                         i*(max(potential[:,0]) -
                                         min(potential[:,0])) +
                                         increment*i)
            extended_potential[idx,1] = potential[j,1]
            idx = idx + 1

    return extended_potential
#------------------------------------------------------------------------------

def sort_potential(potential):
    """
    Sort the dataset based on the potential values in ascending order.

    Parameters:
        potential (numpy.ndarray): The dataset containing potential values in the format (x, potential).

    Returns:
        sorted_potential (numpy.ndarray): The sorted dataset based on the potential values.

    Example:
        >>> potential = np.array([[3, 4], [1, 2], [2, 3], [0, 1]])
        >>> result = sort_potential(potential)
        >>> print(result)

    """
    idx = sorted(potential[:,0])
    sorted_potential = potential.copy()
    for i in range(len(idx)):
        for j in range(len(potential)):
            if potential[j,0] == idx[i]:
                sorted_potential[i,0] = idx[i]
                sorted_potential[i,1] = potential[j,1]

    return sorted_potential
#------------------------------------------------------------------------------

def diff_potentials(potential_a, potential_b,start,end,tol=0.04):
    """
    Subtract potential values between two datasets within a specified range.

    Parameters:
        potential_a (numpy.ndarray): The first dataset containing potential values in the format (x, potential).
        potential_b (numpy.ndarray): The second dataset containing potential values in the format (x, potential).
        start (float): The starting position for potential subtraction.
        end (float): The ending position for potential subtraction.
        tol (float, optional): The tolerance value for potential comparison. Default is 0.04.

    Returns:
        new_potential (numpy.ndarray): The resulting dataset with subtracted potential values.

    Example:
        >>> potential_a = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> potential_b = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> start = 0.5
        >>> end = 2.5
        >>> tol = 0.02
        >>> result = diff_potentials(potential_a, potential_b, start, end, tol)
        >>> print(result)

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
#------------------------------------------------------------------------------

def translate_grid(potential, translation, periodic=False,
                   vector=[0,0,0], boundary_shift=0.0):
    """
    Translates the grid points of a given potential by a specified translation along the vector direction.

    Parameters:
        potential (numpy.ndarray): Array containing potential data with shape (N, 2), where N is the number of grid points.
        translation (float): The amount of translation to apply to the grid points along the specified vector direction.
        periodic (bool, optional): Whether to apply periodic boundary conditions. Default is False.
        vector (list, optional): The direction vector for translation. Default is [0, 0, 0].
        boundary_shift (float, optional): The amount of shift to consider when applying periodic boundary conditions. Default is 0.0.

    Returns:
        numpy.ndarray: An array containing the translated potential data with shape (N, 2).

    Example:
        # Sample potential data
        >>> potential = np.array([[0.0, 1.0], [0.5, 2.0], [1.0, 3.0]])

        # Translate the grid by 0.2 along the x-direction
        >>> translated_potential = translate_grid(potential, 0.2)

        >>> print(translated_potential)
    """
    new_potential_trans = np.zeros((len(potential),2))
    new_potential_orig = np.zeros((len(potential),2))
    length = np.sqrt(vector.dot(vector))


    for i in range(len(potential)):
        new_potential_trans[i,0] = potential[i,0] + translation
        new_potential_trans[i,1] = potential[i,1]
        if periodic == True:
            new_potential_trans[i,0] = (
                new_potential_trans[i,0] -
                length * int((new_potential_trans[i,0]+boundary_shift)/length))

    if periodic == True:
    # Sort the numbers out if you have done periodic wrapping
        sorted_potential_trans = sort_potential(new_potential_trans)
    else:
        sorted_potential_trans = new_potential_trans

    return sorted_potential_trans
#------------------------------------------------------------------------------

def create_plotting_mesh(NGX,NGY,NGZ,pc,grad):
    """
    Creates a plotting mesh based on the given grid data and plane coefficients.

    Parameters:
        NGX (int): Number of grid points along the x-direction.
        NGY (int): Number of grid points along the y-direction.
        NGZ (int): Number of grid points along the z-direction.
        pc (numpy.ndarray): Array containing plane coefficients with shape (4,).
        grad (numpy.ndarray): Array containing gradient data with shape (NGX, NGY, NGZ).

    Returns:
        numpy.ndarray: A 2D array representing the plotting mesh with shape (a, b), where 'a' and 'b' depend on the plane direction.

    Example:
        # Sample grid data and plane coefficients
        >>> NGX, NGY, NGZ = 10, 10, 10
        >>> pc = np.array([0, 0, 1, 5])
        >>> grad = np.random.rand(NGX, NGY, NGZ)

        # Create the plotting mesh
        >>> plotting_mesh = create_plotting_mesh(NGX, NGY, NGZ, pc, grad)

        >>> print(plotting_mesh)
    """
    if pc[0] == 0 and pc[1] == 0:
        a = NGX; b = NGY; p = 'zzo'; c = int(pc[3] / pc[2]) - 1
    if pc[0] == 0 and pc[2] == 0:
        a = NGX; b = NGZ; p = 'zoz'; c = int(pc[3] / pc[1]) - 1
    if pc[1] == 0 and pc[2] == 0:
        a = NGY; b = NGZ; p = 'ozz'; c = int(pc[3] / pc[0]) - 1
    plane = np.zeros(shape=(a,b))
    for x in range(a):
        for y in range(b):
            if p == 'zzo':
                plane[x,y] = grad[x,y,c]

    return plane
#------------------------------------------------------------------------------

def read_cube_density(FILE):
    """
    Reads a cube density file and extracts relevant information.

    Parameters:
        FILE (str): The path to the cube density file.

    Returns:
        numpy.ndarray: A 3x3 numpy array representing the lattice.

    Example:
        >>> file_path = 'path/to/your/cube_density_file.cube'

        # Read the cube density file and get the lattice
        >>> lattice = read_cube_density(file_path)
        >>> print(lattice)
    """
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
#------------------------------------------------------------------------------

def points_2_plane(a,b,c):
    """
    Calculates the plane coefficients from three points in space.

    Parameters:
        a (numpy.ndarray): First point with shape (3,).
        b (numpy.ndarray): Second point with shape (3,).
        c (numpy.ndarray): Third point with shape (3,).

    Returns:
        numpy.ndarray: An array containing the plane coefficients with shape (4,).

    Example:
        # Sample points in space
        >>> a = np.array([0, 0, 0])
        >>> b = np.array([1, 0, 0])
        >>> c = np.array([0, 1, 0])

        # Calculate plane coefficients
        >>> plane_coefficients = points_2_plane(a, b, c)
        >>> print(plane_coefficients)
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
#------------------------------------------------------------------------------

def get_third_coordinate(plane_coeff,NGX,NGY):
    """
    Computes the third coordinate of the plane for given plane coefficients.

    Parameters:
        plane_coeff (numpy.ndarray): An array containing the plane coefficients with shape (4,).
        NGX (int): Number of grid points along the x-direction.
        NGY (int): Number of grid points along the y-direction.

    Returns:
        list: A list of third coordinates for the plane.

    Example:
        # Sample plane coefficients and grid dimensions
        >>> plane_coeff = np.array([1, 1, 1, 5])
        >>> NGX, NGY = 10, 10

        # Calculate the third coordinate of the plane
        >>> third_coordinates = get_third_coordinate(plane_coeff, NGX, NGY)
        >>> print(third_coordinates)
    """
    zz = []
    i = j = 0
    while i <= NGX:
        i = i + 1
        j = 0
        while j <= NGY:
            j = j + 1
            rounded = round(((plane_coeff[0]*j+plane_coeff[1]*i) /
                             plane_coeff[2]))
            standard = ((plane_coeff[0]*j+plane_coeff[1]*i) /
                        plane_coeff[2])
            if rounded == standard:   # Is it a whole number?
                zz.append(-(plane_coeff[0]*i+plane_coeff[1]*j)/plane_coeff[2])

    return zz
