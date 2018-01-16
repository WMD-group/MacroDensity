###############################################################################
# Copyright Keith Butler(2014)                                                #
#                                                                             #
# This file MacroDensity.density_tools.py is free software: you can           #
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
################################################################################

from __future__ import print_function
from functools import reduce
import math

import numpy
import numpy as np
from scipy import interpolate

#------------------------------------------------------------------------------
def gradient_magnitude(gx, gy, gz):
    """Converts the separate gradient magnitudes to a single magnitude
    Args:
        gx/y/z : fields in x y and z directions 2D array
    Returns:
        grad_mag : gradient of fields at each point"""

    grad_mag = gx
    for i in range(gx.shape[0]):
        for j in range(gy.shape[1]):
            for k in range(gz.shape[2]):
                grad_mag[i,j,k] = np.sqrt(gx[i,j,k]**2 +
                                          gy[i,j,k]**2 +
                                          gz[i,j,k]**2)

    return grad_mag
#------------------------------------------------------------------------------
def vector_2_abscissa(vector, magnitude, dx, dy, dz):
    """Converts a vector with a magnitude given in units of grid density
    (NGX/Y/Z) to AA for plotting
    Args:
        vector : the vector along which the line is being plotted [(3x1) array]
        magnitude : the number of steps that were taken along that vector
            [Integer]
        dx/y/z: the resolution of the density grid in AA-1 [Real]
    Returns:
        abscissa : the values for plotting on the abscissa in AA [1D array]
    """
    vec_mag = np.linalg.norm([vector[0] * dx, vector[1] * dy, vector[2] * dz])
    abscissa = [i * vec_mag for i in range(magnitude)]

    return np.asarray(abscissa)
#------------------------------------------------------------------------------
def number_in_field(gradients, cutoff):
    """Get number of grid elements with a field magnitude greater than cutoff
    Args:
        gradients: the grid of field gradients (Real(ngx,ngy,ngz))
        cutoff: the value above which tocout them (Real)
    Returns:
        number_of_elements: the number satisfying the condition (Integer)
    """
    number_of_elements = 0
    for element in np.nditer(gradients):
        if element >= cutoff:
            number_of_elements += 1

    return number_of_elements
#------------------------------------------------------------------------------
def element_vol(vol, nx, ny, nz):
    """Calculates the volume of each of the elements on the grid.
    Args:
        vol: the cell volume (real)
        x : the number of grid points in each direction (real)
    Returns:
        ele_vol : the volume (real)
    """
    number_of_elements = nx * ny * nz
    ele_vol = vol / number_of_elements

    return ele_vol

#------------------------------------------------------------------------------
def one_2_2d(Array, resolution, vector):
    """Converts the 1d potential array to 2D with angstroms in A[0]
    Args:
        Array: 1D array
        resolution: density of sampling of distance (1/AA)
        vector: The vector of the direction of sampling
    Returns
        New_array: 2D array
    """
    length = np.sqrt(vector.dot(vector))
    New_array = np.zeros(shape=(len(Array) - 1, 2))
    resolution = length / len(Array)
    for i in range(len(Array) - 1):
        New_array[i,0] = i*resolution
        New_array[i,1] = Array[i]

    return New_array
#------------------------------------------------------------------------------

def macroscopic_average(potential, periodicity, resolution):
    """Getting the macroscopic average of potential
    Args:
        potential : array containig the electrostaticpotential/charge density
        periodicity : real number; the period over which to average
        resolution : the grid resolution in the direction of averaging
    Returns:
        macro_average : array with the macroscopically averaged values"""

    macro_average = np.zeros(shape=(len(potential)))
    period_points = int((periodicity/resolution))
    # Re-arrange so that period points divides evenly by resolution
    for i in range(len(potential)):
        for j in range(i - int(period_points / 2),
                       i + int(period_points / 2)):
            if j < 0:
                macro_average[i] = (macro_average[i] +
                                    potential[j + len(potential)])
            elif j >= len(potential):
                macro_average[i] = (macro_average[i] +
                                    potential[j - len(potential)])
            else:
                macro_average[i] = macro_average[i] + potential[j]
        macro_average[i] = macro_average[i] / period_points

    print("Average of the average = ", numpy.average(macro_average))
    return macro_average
#------------------------------------------------------------------------------

def cube_potential(origin, travelled, cube, Grid, nx, ny, nz):
    """Populates the sampling cube with the potential required"""

    # Recalc the origin as grid point coordinates
    n_origin = np.zeros(shape=(3))
    n_origin[0] = int(origin[0]*nx)
    n_origin[1] = int(origin[1]*ny)
    n_origin[2] = int(origin[2]*nz)
    potential_cube = np.zeros(shape=(cube[0],cube[1],cube[2]))
    for x in range(0,cube[0]):
        for y in range(0,cube[1]):
            for z in range(0,cube[2]):
                # Assign the values of coordinates in the original grid
                xv = int(n_origin[0]+travelled[0]+x)
                yv = int(n_origin[1]+travelled[1]+y)
                zv = int(n_origin[2]+travelled[2]+z)
                # Minimum image convention
                zv = int(zv - nz*round(zv/nz))
                yv = int(yv - ny*round(yv/ny))
                xv = int(xv - nx*round(xv/nx))
                potential_cube[x,y,z] = Grid[xv,yv,zv]

    return potential_cube.mean(), np.var(potential_cube)
#------------------------------------------------------------------------------

def cuboid_average(Grid, cube, origin, vector, nx, ny, nz, magnitude):
   """Calculates the average in a cube defined by size cube(a,b,c), beginning
    at origin and travelling as far as magnitude."""

   plotting_average = np.zeros(shape=(magnitude))
   i = 0
   while i < magnitude:
         travelled = np.multiply(i, vector)
         plotting_average[i], varience = cube_potential(origin, travelled,
                                                        cube, Grid,
                                                        nx, ny, nz)
         i = i + 1

   return plotting_average
#------------------------------------------------------------------------------

def planar_average(Grid, nx, ny, nz, axis='z'):
    """Calculate the average in a given plane for the full length of the
    normal; e.g. the full length of z in the xy plane."""
    if axis == 'x':
        x_plane = np.zeros(shape=(ny, nz))
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
#------------------------------------------------------------------------------
def get_volume(a,b,c):
    """Calculate the volume of the cell from lattice vectors
    Args:
        a/b/c: vectors of the lattice edges
    """
    volume = np.dot(a,np.cross(b,c))

    return volume
#------------------------------------------------------------------------------

def numbers_2_grid(a,NGX,NGY,NGZ):
    """Takes a point (in fractional coordinates) and converts it to a VASP grid
    point based on the NGX/Y/Z values."""
    a_grid = np.zeros(shape=(3))
    a_grid[0] = round(float(a[0])*NGX)
    a_grid[1] = round(float(a[1])*NGY)
    a_grid[2] = round(float(a[2])*NGZ)

    return a_grid
#------------------------------------------------------------------------------

def matrix_2_abc(Lattice):
    """The the VASP lattice and convert to the a,b,c,alpha,beta,gamma format"""

    a = np.sqrt(Lattice[0,0]**2+Lattice[0,1]**2+Lattice[0,2]**2)
    b = np.sqrt(Lattice[1,0]**2+Lattice[1,1]**2+Lattice[1,2]**2)
    c = np.sqrt(Lattice[2,0]**2+Lattice[2,1]**2+Lattice[2,2]**2)

    a_vec = Lattice[0,:]
    b_vec = Lattice[1,:]
    c_vec = Lattice[2,:]

    return a,b,c,a_vec,b_vec,c_vec
#------------------------------------------------------------------------------
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

    for i in range(3):
       for j in range(3):
           lattice[i,j] = lattice[i,j] * scale_factor
    # Restart reading to get the coordinates...it's just easier this way!
    i = 0
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
            Potential = numpy.zeros(shape=(NGX * NGY * NGZ))
            # Read in the potential data
            upper_limit =  (int(NGX * NGY * NGZ / 5) +
                            np.mod(NGX * NGY * NGZ, 5))
        if i > (num_atoms + 9) and i < (num_atoms + 10 + upper_limit):
            for m in range(len(inp)):
                Potential[k + m] = inp[m]
            k = k + 5
            if math.fmod(k, 100000) == 0:
                print("Reading potential at point", k)

    print("BBBB       OOOO        OOOO        MMMMM   ")
    print("BBBB       OOOO        OOOO        MMMMM   ")
    print("BBBB       OOOO        OOOO        MMMMM   ")
    print("B  B       OOOO        OOOO        MMMMM   ")
    print("B  B       O  O        O  O        MMMMM   ")
    print("B  B       O  O        O  O        MMMMM   ")
    print("B  B       O  O        O  O        MMMMM   ")
    print("B  B       O  O        O  O        MMMMM   ")
    print("BBBB       O  O        O  O        M M M   ")
    print("BBBB       O  O        O  O        M M M   ")
    print("BBBB       O  O        O  O        M M M   ")
    print("B  B       O  O        O  O        M M M   ")
    print("B  B       O  O        O  O        M M M   ")
    print("B  B       O  O        O  O        M M M   ")
    print("B  B       O  O        O  O        M M M   ")
    print("B  B       OOOO        OOOO        M M M   ")
    print("BBBB       OOOO        OOOO        M M M   ")
    print("BBBB       OOOO        OOOO        M M M   ")
    print("BBBB       OOOO        OOOO        M M M   ")

    print("Average of the potential = ", numpy.average(Potential))
    f.close()
    return Potential, NGX, NGY, NGZ, lattice
#------------------------------------------------------------------------------
def density_2_grid(Density, nx, ny, nz, Charge=False, Volume=1):
    """Convert the Potetnial list to a grid for ease of manipulation
    Args:
        Density: Array of the output from a VAsp calulation charge/potential
        nx,y,z : Number of mesh points in x/y/z
        Charge : Boolean, is it charge or potential (charge needs to be
            normalised by vol)
        Volume : The lattice vectors, only required for normalising charge.
     Returns:
        Potential_grid: the (normalised) quantity on a mesh
        total_electrons : the number of electrons in the system
        """
    l = 0
    Potential_grid = np.zeros(shape=(nx,ny,nz))
    total_electrons = 0
    is_CHGCAR = True
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                Potential_grid[i,j,k] = Density[l] / Volume
                if Charge == True:
                    # Convert the charge density to a number of electrons
                    point_volume = Volume / (nx*ny*nz)
                    Potential_grid[i,j,k] = Potential_grid[i,j,k]*point_volume
                total_electrons = total_electrons + Density[l]
                l = l + 1
    if Charge == True:
        print("Total electrons: ", total_electrons / (nx * ny * nz))
    total_electrons = total_electrons / (nx * ny * nz)
    return Potential_grid, total_electrons
#------------------------------------------------------------------------------
def GCD(a,b):
    """ The Euclidean Algorithm """
    a = abs(a)
    b = abs(b)
    while a:
        a, b = (b % a), a
    return b
#------------------------------------------------------------------------------
def GCD_List(list):
    """ Finds the GCD of numbers in a list.
    Input: List of numbers you want to find the GCD of
            E.g. [8, 24, 12]
    Returns: GCD of all numbers
    """
    return reduce(GCD, list)
#------------------------------------------------------------------------------
