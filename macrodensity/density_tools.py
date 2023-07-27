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
###############################################################################

from __future__ import division, print_function

import math
from functools import reduce
from itertools import chain

import numpy
import numpy as np


#------------------------------------------------------------------------------
def gradient_magnitude(gx, gy, gz):
    """
    Calculate the magnitude of the gradient at each point in a 3D field.

    Parameters:
        gx (numpy.ndarray): Gradient along the x-axis.
        gy (numpy.ndarray): Gradient along the y-axis.
        gz (numpy.ndarray): Gradient along the z-axis.

    Returns:
        numpy.ndarray: 3D array representing the magnitude of the gradient at each point.

    Example:
        >>> gx = np.random.rand(3, 3, 3)
        >>> gy = np.random.rand(3, 3, 3)
        >>> gz = np.random.rand(3, 3, 3)
        >>> grad_magnitude = gradient_magnitude(gx, gy, gz)
        >>> print("Gradient Magnitude:")
        >>> print(grad_magnitude)

    """
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
    """
    Convert a 3D vector to an array of abscissa values.

    Parameters:
        vector (tuple or list): 3D vector represented as (x, y, z).
        magnitude (int or float): Magnitude of the vector.
        dx (int or float): Spacing along the x-axis.
        dy (int or float): Spacing along the y-axis.
        dz (int or float): Spacing along the z-axis.

    Returns:
        numpy.ndarray: 1D array containing abscissa values based on the vector and spacing.

    Example:
        >>> vector = (1, 2, 3)
        >>> magnitude = 5.0
        >>> dx, dy, dz = 0.1, 0.2, 0.3
        >>> abscissa_array = vector_2_abscissa(vector, magnitude, dx, dy, dz)
        >>> print("Abscissa Array:")
        >>> print(abscissa_array)
        
    """
    vec_mag = np.linalg.norm([vector[0] * dx, vector[1] * dy, vector[2] * dz])
    abscissa = [i * vec_mag for i in range(magnitude)]

    return np.asarray(abscissa)
#------------------------------------------------------------------------------

def number_in_field(gradients, cutoff):
    """
    Count the number of elements in a field that have a value greater than or equal to the cutoff.

    Parameters:
        gradients (numpy.ndarray): 3D array representing the field.
        cutoff (int or float): Threshold value for counting elements.

    Returns:
        int: Number of elements in the field satisfying the cutoff condition.

    Example:
        >>> gradients_field = np.random.rand(4, 4, 4)
        >>> cutoff_value = 0.5
        >>> num_elements_above_cutoff = number_in_field(gradients_field, cutoff_value)
        >>> print("Number of Elements Above Cutoff:", num_elements_above_cutoff)
    """
    number_of_elements = 0
    for element in np.nditer(gradients):
        if element >= cutoff:
            number_of_elements += 1

    return number_of_elements
#------------------------------------------------------------------------------

def element_vol(vol, nx, ny, nz):
    """
    Calculate the volume of each element in a 3D grid.

    Parameters:
        vol (int or float): Total volume of the 3D grid.
        nx (int): Number of elements along the x-axis.
        ny (int): Number of elements along the y-axis.
        nz (int): Number of elements along the z-axis.

    Returns:
        float: Volume of each individual element in the grid.

    Example:
        >>> volume = 10.0
        >>> nx, ny, nz = 5, 5, 5
        >>> element_volume = element_vol(volume, nx, ny, nz)
        >>> print("Volume of Each Element:", element_volume)

    """
    number_of_elements = nx * ny * nz
    ele_vol = vol / number_of_elements

    return ele_vol

#------------------------------------------------------------------------------

def one_2_2d(Array, resolution, vector):
    """
    Transform a 1D array to a 2D array with abscissa values based on the given resolution and vector.

    Parameters:
        Array (numpy.ndarray): 1D array to be transformed.
        resolution (int or float): Spacing between abscissa values.
        vector (numpy.ndarray): 3D vector used for the transformation.

    Returns:
        numpy.ndarray: 2D array with abscissa values and the corresponding Array values.

    Example:
        >>> Array = np.random.rand(10)
        >>> resolution = 0.5
        >>> vector = np.array([1, 2, 3])
        >>> transformed_array = one_2_2d(Array, resolution, vector)
        >>> print("Transformed Array:")
        >>> print(transformed_array)
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
    """
    Calculate the macroscopic average of a 1D potential field with periodicity.

    Parameters:
        potential (numpy.ndarray): 1D array representing the potential field.
        periodicity (int or float): Periodicity of the field.
        resolution (int or float): Spacing between potential data points.

    Returns:
        numpy.ndarray: 1D array containing the macroscopic average of the potential field.

    Example:
        >>> potential = np.random.rand(20)
        >>> periodicity = 2.0
        >>> resolution = 0.1
        >>> macro_avg_result = macroscopic_average(potential, periodicity, resolution)
        >>> print("Macroscopic Average Result:")
        >>> print(macro_avg_result)
    """
    macro_average = np.zeros(shape=(len(potential)))
    period_points = int((periodicity/resolution))
    # Period points must be even
    if period_points % 2 != 0:
        period_points = period_points + 1

    length = len(potential)
    for i in range(length):
        start = i - int(period_points / 2)
        end = i + int(period_points / 2)
        if start < 0:
            start = start + length
            macro_average[i] = macro_average[i] + sum(potential[0:end]) + sum(potential[start:length])
            macro_average[i] = macro_average[i] / period_points
        elif end >= length:
            end = end - length
            macro_average[i] = macro_average[i] + sum(potential[start:length]) + sum(potential[0:end])
            macro_average[i] = macro_average[i] / period_points
        else:
            macro_average[i] = macro_average[i] + sum(potential[start:end]) / period_points

    print("Average of the average = ", numpy.average(macro_average))

    return macro_average
#------------------------------------------------------------------------------

def volume_average(origin, cube, grid, nx, ny, nz, travelled=[0, 0, 0]):
    """
    Calculate the volume average and variance of a cube in a 3D grid.

    Parameters:
        origin (tuple or list): Coordinates of the origin point.
        cube (tuple or list): Dimensions of the cube (x, y, z).
        grid (numpy.ndarray): 3D array representing the data grid.
        nx (int): Number of points along the x-axis in the grid.
        ny (int): Number of points along the y-axis in the grid.
        nz (int): Number of points along the z-axis in the grid.
        travelled (tuple or list, optional): Distance travelled from the origin in each direction (x, y, z). Default is [0, 0, 0].

    Returns:
        tuple: A tuple containing the volume average and variance of the cube.

    Example:
        >>> origin = (0.5, 0.5, 0.5)
        >>> cube = (3, 3, 3)
        >>> grid = np.random.rand(10, 10, 10)
        >>> nx, ny, nz = 10, 10, 10
        >>> travelled = [1, 2, 3]
        >>> avg, variance = volume_average(origin, cube, grid, nx, ny, nz, travelled)
        >>> print("Volume Average:", avg)
        >>> print("Variance:", variance)

    """
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
                potential_cube[x,y,z] = grid[int(xv),int(yv),int(zv)]

    return potential_cube.mean(), np.var(potential_cube)
#------------------------------------------------------------------------------

def travelling_volume_average(grid, cube, origin, vector, nx, ny, nz, magnitude):
   """
    Calculate the volume average at multiple positions along a given vector.

    Parameters:
        grid (numpy.ndarray): 3D array representing the data grid.
        cube (tuple or list): Dimensions of the cube (x, y, z).
        origin (tuple or list): Coordinates of the origin point.
        vector (tuple or list): 3D vector representing the direction of travel.
        nx (int): Number of points along the x-axis in the grid.
        ny (int): Number of points along the y-axis in the grid.
        nz (int): Number of points along the z-axis in the grid.
        magnitude (int): Number of positions to travel along the vector.

    Returns:
        numpy.ndarray: 1D array containing the volume averages at each position along the vector.

    Example:
        >>> vector = (0.1, 0.2, 0.3)
        >>> magnitude = 5
        >>> travelling_avg = travelling_volume_average(grid, cube, origin, vector, nx, ny, nz, magnitude)
        >>> print("Travelling Volume Average:")
        >>> print(travelling_avg)
    """
   plotting_average = np.zeros(shape=(magnitude))
   i = 0
   while i < magnitude:
         travelled = np.multiply(i, vector)
         plotting_average[i], varience = volume_average(origin,
                                                        cube, grid,
                                                        nx, ny, nz, travelled)
         i = i + 1

   return plotting_average
#------------------------------------------------------------------------------

def planar_average(grid, nx, ny, nz, axis='z'):
    """
    Calculate the planar average of a 3D grid along a specified axis.

    Parameters:
        grid (numpy.ndarray): 3D array representing the data grid.
        nx (int): Number of points along the x-axis in the grid.
        ny (int): Number of points along the y-axis in the grid.
        nz (int): Number of points along the z-axis in the grid.
        axis (str, optional): Axis along which to calculate the average ('x', 'y', or 'z'). Default is 'z'.

    Returns:
        numpy.ndarray: 1D array containing the planar average along the specified axis.

    Example:
        >>> axis = 'z'
        >>> planar_avg = planar_average(grid, nx, ny, nz, axis)
        >>> print("Planar Average along axis", axis)
        >>> print(planar_avg)
    """
    if axis == 'x':
        x_plane = np.zeros(shape=(ny, nz))
        Average = np.zeros(shape=(nx))
        for x_value in range(nx):
            x_plane[:,:] = grid[x_value,:,:]
            Average[x_value] = x_plane.mean()
    if axis == 'y':
        Average = np.zeros(shape=(ny))
        y_plane = np.zeros(shape=(nx,nz))
        for y_value in range(ny):
            y_plane[:,:] = grid[:,y_value,:]
            Average[y_value] = y_plane.mean()
    if axis == 'z':
        Average = np.zeros(shape=(nz))
        z_plane = np.zeros(shape=(nx,ny))
        for z_value in range(nz):
            z_plane[:,:] = grid[:,:,z_value]
            Average[z_value] = z_plane.mean()

    return Average
#------------------------------------------------------------------------------

def get_volume(a,b,c):
    """
    Calculate the volume of a parallelepiped defined by three vectors a, b, and c.

    Parameters:
        a (numpy.ndarray): 1D array representing vector a.
        b (numpy.ndarray): 1D array representing vector b.
        c (numpy.ndarray): 1D array representing vector c.

    Returns:
        float: Volume of the parallelepiped defined by the three vectors.

    Example:
        >>> a = np.array([1, 0, 0])
        >>> b = np.array([0, 1, 0])
        >>> c = np.array([0, 0, 1])
        >>> volume = get_volume(a, b, c)
        >>> print("Volume of parallelepiped:", volume)
    """
    volume = np.dot(a,np.cross(b,c))

    return volume
#------------------------------------------------------------------------------

def numbers_2_grid(a,NGX,NGY,NGZ):
    """
    Convert fractional coordinates to grid point coordinates.

    Parameters:
        a (tuple or list): Fractional coordinates (x, y, z).
        NGX (int): Number of grid points along the x-axis.
        NGY (int): Number of grid points along the y-axis.
        NGZ (int): Number of grid points along the z-axis.

    Returns:
        numpy.ndarray: 1D array containing the grid point coordinates (x, y, z).

    Example:
        >>> fractional_coords = [0.3, 0.4, 0.5]
        >>> NGX, NGY, NGZ = 10, 10, 10
        >>> grid_coords = numbers_2_grid(fractional_coords, NGX, NGY, NGZ)
        >>> print("Grid Point Coordinates:", grid_coords)
    """
    a_grid = np.zeros(shape=(3))
    a_grid[0] = round(float(a[0])*NGX)
    a_grid[1] = round(float(a[1])*NGY)
    a_grid[2] = round(float(a[2])*NGZ)

    return a_grid
#------------------------------------------------------------------------------

def matrix_2_abc(Lattice):
    """
    Extract lattice parameters and vectors from a 3x3 matrix representing a lattice.

    Parameters:
        Lattice (numpy.ndarray): 3x3 matrix representing the lattice.

    Returns:
        tuple: A tuple containing the lattice parameters a, b, c and lattice vectors a_vec, b_vec, c_vec.

    Example:
        >>> Lattice = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        >>> a, b, c, a_vec, b_vec, c_vec = matrix_2_abc(Lattice)
        >>> print("Lattice parameters:", a, b, c)
        >>> print("Lattice vectors:")
        >>> print(a_vec)
        >>> print(b_vec)
        >>> print(c_vec)
    """
    a = np.sqrt(Lattice[0,0]**2+Lattice[0,1]**2+Lattice[0,2]**2)
    b = np.sqrt(Lattice[1,0]**2+Lattice[1,1]**2+Lattice[1,2]**2)
    c = np.sqrt(Lattice[2,0]**2+Lattice[2,1]**2+Lattice[2,2]**2)

    a_vec = Lattice[0,:]
    b_vec = Lattice[1,:]
    c_vec = Lattice[2,:]

    return a,b,c,a_vec,b_vec,c_vec
#------------------------------------------------------------------------------

def _print_boom(quiet=False):
    if not quiet:
        print("\n")
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

def read_vasp_density(FILE, use_pandas=None, quiet=False):
    """
    Read density data from a VASP CHGCAR-like file.

    Parameters:
        FILE (str): Path to the CHGCAR-like file.
        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files).
                                             If False, use Numpy. If None, automatically use Pandas if available.
                                             Default is None.
        quiet (bool, optional): If True, suppress print statements during reading. Default is False.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - numpy.ndarray: 3x3 array representing the lattice vectors.

    Example:
        >>> FILE = "path/to/your/CHGCAR_file"
        >>> potential_data, NGX, NGY, NGZ, lattice = read_vasp_density(FILE)
        >>> print("Potential Data:")
        >>> print(potential_data)
        >>> print("Number of grid points along x, y, z axes:", NGX, NGY, NGZ)
        >>> print("Lattice Vectors:")
        >>> print(lattice)
    """
    # Get Header information by reading a line at a time

    if use_pandas:
        from pandas import read_table as pandas_read_table
    elif use_pandas is None:
        try:
            from pandas import read_table as pandas_read_table
            use_pandas = True
        except ImportError:
            use_pandas = False

    print("Reading header information...")
    with open(FILE, "r") as f:
        _ = f.readline()
        scale_factor = float(f.readline())

        lattice = np.zeros(shape=(3,3))
        for row in range(3):
            lattice[row] = [float(x) for x in f.readline().split()]
        lattice = lattice * scale_factor

        num_species = len(f.readline().split())
        num_type = [int(x) for x in f.readline().split()]
        num_atoms = sum(num_type)
        coord_type = f.readline().strip()

        coordinates = numpy.zeros(shape=(num_atoms, 3))
        for atom_i in range(num_atoms):
            coordinates[atom_i] = [float(x) for x in f.readline().split()]

        # Skip blank line
        _ = f.readline()

        NGX, NGY, NGZ = [int(x) for x in f.readline().split()]

        if use_pandas:
            print("Reading 3D data using Pandas...")
            skiprows = 10 + num_atoms
            readrows = int(math.ceil(NGX * NGY * NGZ / 5))

            dat = pandas_read_table(FILE, delim_whitespace=True,
                                    skiprows=skiprows, header=None,
                                    nrows=readrows)
            Potential = dat.iloc[:readrows, :5].values.flatten()
            remainder = (NGX * NGY * NGZ) % 5
            if remainder > 0:
                Potential = Potential[:(-5 + remainder)]

        else:
            print("Reading 3D data...")
            Potential = (f.readline().split()
                             for i in range(int(math.ceil(NGX * NGY * NGZ / 5))))
            Potential = numpy.fromiter(chain.from_iterable(Potential), float)

    ##_print_boom(quiet=quiet)
    if not quiet:
        print("Average of the potential = ", numpy.average(Potential))

    return Potential, NGX, NGY, NGZ, lattice
#------------------------------------------------------------------------------

def _read_partial_density(FILE, use_pandas, num_atoms, NGX, NGY, NGZ, spin=0):
    """
    This function is used internally within the read_casp_parchg, reading partial density data from a VASP-PARCHG file.

    Args:
        FILE (str): Path to the CHGCAR-like file.
        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files).
                                             If False, use Numpy. If None, automatically use Pandas if available.
                                             Default is None.
        num_atoms (int): Total number of atoms in the system.
        NGX (int): Number of grid points along the x-axis.
        NGY (int): Number of grid points along the y-axis.
        NGZ (int): Number of grid points along the z-axis.
        spin (int, optional): If 0, read the first spin channel (total density).
                              If 1, read the second spin channel (spin-up or spin-down).
                              Default is 0.

    Returns:
        numpy.ndarray: 1D array representing the partial density data for the specified spin channel.
    """
    print("PANDAS:", use_pandas)
    if use_pandas:
        from pandas import read_table as pandas_read_table
    elif use_pandas is None:
        try:
            from pandas import read_table as pandas_read_table
            use_pandas = True
        except ImportError:
            use_pandas = False


    with open(FILE, "r") as f:
        _ = f.readline()
        scale_factor = float(f.readline())

        lattice = np.zeros(shape=(3,3))
        for row in range(3):
            lattice[row] = [float(x) for x in f.readline().split()]
        lattice = lattice * scale_factor

        num_species = len(f.readline().split())
        num_type = [int(x) for x in f.readline().split()]
        num_atoms = sum(num_type)
        coord_type = f.readline().strip()

        coordinates = numpy.zeros(shape=(num_atoms, 3))
        for atom_i in range(num_atoms):
            coordinates[atom_i] = [float(x) for x in f.readline().split()]

        # Skip blank line
        _ = f.readline()

        NGX, NGY, NGZ = [int(x) for x in f.readline().split()]

        if use_pandas:
            print("Reading 3D data using Pandas...")
            skiprows = 10 + num_atoms + spin * \
                          (math.ceil(NGX * NGY * NGZ / 10) + 2)
            readrows = int(math.ceil(NGX * NGY * NGZ / 10))

            dat = pandas_read_table(FILE, delim_whitespace=True,
                                    skiprows=skiprows, header=None,
                                    nrows=readrows)
            density = dat.iloc[:readrows, :10].values.flatten()
            remainder = (NGX * NGY * NGZ) % 10
            if remainder > 0:
                density = density[:(-10 + remainder)]
        else:
            print("Reading 3D data...")
            density = (f.readline().split()
                             for i in range(int(math.ceil(NGX * NGY * NGZ / 10))))
            density = numpy.fromiter(chain.from_iterable(density), float)

    return density

#------------------------------------------------------------------------------
def read_vasp_parchg(FILE, use_pandas=None, quiet=False, spin=False):
    """
    Read density data or spin-polarized partial density data from a VASP PARCHG-like file.

    Parameters:
        FILE (str): Path to the PARCHG-like file.
        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files).
                                             If False, use Numpy. If None, automatically use Pandas if available.
                                             Default is None.
        quiet (bool, optional): If True, suppress print statements during reading. Default is False.
        spin (bool, optional): If True, read spin-polarized partial densities. Default is False.

    Returns:
        tuple: A tuple containing:
            - list or numpy.ndarray: List containing numpy arrays representing the density data for each spin channel,
                                     or numpy.ndarray for the total density if spin is False.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - numpy.ndarray: 3x3 array representing the lattice vectors.

    Example:
        >>> FILE = "path/to/PARCHG-like-file"
        >>> density, NGX, NGY, NGZ, lattice = read_vasp_parchg(FILE)
        >>> if isinstance(density, list):
            >>> print("Spin-polarized Partial Densities:")
            >>> for i, spin_density in enumerate(density):
                >>> print(f"Spin {i+1} Density Data:", spin_density)
        >>> else:
            >>> print("Total Density Data:", density)
        >>> print("Number of Grid Points (NGX, NGY, NGZ):", NGX, NGY, NGZ)
        >>> print("Lattice Vectors:")
        >>> print(lattice)
    """
    # Get Header information by reading a line at a time

    print("Reading header information...")
    with open(FILE, "r") as f:
        _ = f.readline()
        scale_factor = float(f.readline())

        lattice = np.zeros(shape=(3,3))
        for row in range(3):
            lattice[row] = [float(x) for x in f.readline().split()]
        lattice = lattice * scale_factor

        num_species = len(f.readline().split())
        num_type = [int(x) for x in f.readline().split()]
        num_atoms = sum(num_type)
        coord_type = f.readline().strip()

        coordinates = numpy.zeros(shape=(num_atoms, 3))
        for atom_i in range(num_atoms):
            coordinates[atom_i] = [float(x) for x in f.readline().split()]

        # Skip blank line
        _ = f.readline()

        NGX, NGY, NGZ = [int(x) for x in f.readline().split()]

        if not spin:
            density = _read_partial_density(FILE, use_pandas, num_atoms, NGX, NGY, NGZ)
        else:
            densities = []
            densities.append(_read_partial_density(FILE, use_pandas, num_atoms, NGX, NGY, NGZ
                , spin=0))
            densities.append(_read_partial_density(FILE, use_pandas, num_atoms, NGX, NGY, NGZ
                , spin=1))
            alpha = densities[0] + densities[1]
            beta = densities[0] - densities[1]
            density = [alpha, beta]
    ##_print_boom(quiet=quiet)

    return density, NGX, NGY, NGZ, lattice

def read_vasp_density_classic(FILE):
    """
    Read density data from a classic VASP-style file.

    Parameters:
        FILE (str): Path to the density file.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - numpy.ndarray: 3x3 array representing the lattice vectors.
    Example:
        >>> FILE = "path/to/classic-VASP-density-file"
        >>> potential, NGX, NGY, NGZ, lattice = read_vasp_density_classic(FILE)
        >>> print("Potential Data:", potential)
        >>> print("Number of Grid Points (NGX, NGY, NGZ):", NGX, NGY, NGZ)
        >>> print("Lattice Vectors:")
        >>> print(lattice)
    """
    with open(FILE, "r") as f:
        lines = f.readlines()
    return _read_vasp_density_fromlines(lines)
#------------------------------------------------------------------------------
def _read_vasp_density_fromlines(lines):
    """
    Read density data from a list of lines (classic VASP-style format). This function is used internally within the read_vasp_density_classic function (add hyperlink to read_vasp_density_classic function)

    Parameters:
        lines (list): List of lines from the density file.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - numpy.ndarray: 3x3 array representing the lattice vectors.
    """
    i, j, k = 0, 0, 0
    NGX, NGY, NGZ = 0, 0, 0

    lattice = np.zeros(shape=(3,3))
    upper_limit, num_species, scale_factor = 0, 0, 0
    num_atoms = 1 # First test needs to fail until headers have been read
    Potential, Coordinates = np.zeros(1), np.zeros(1)

    for line in lines:
        inp = line.split()

        if inp == []:
            continue
        else:
            i += 1
        if i > (num_atoms + 9) and i < (num_atoms + 10 + upper_limit):
            for m, val in enumerate(inp):
                Potential[k + m] = val
            k = k + 5
            if math.fmod(k, 100000) == 0:
                print("Reading potential at point", k)
        elif i == 2:
            scale_factor = float(inp[0])
        elif i >= 3 and i < 6:
            lattice[i-3,:]=inp[:]
        elif i == 6:
            num_species = len(inp)
            species = inp
        elif i == 7:
            num_type = inp
            num_atoms = sum(int(x) for x in num_type)
        elif i == 8:
            coord_type = inp
            Coordinates = numpy.zeros(shape=(num_atoms,3))
        elif i >= 9 and i <= num_atoms + 8:
            Coordinates[i-9,0] = float(inp[0])
            Coordinates[i-9,1] = float(inp[1])
            Coordinates[i-9,2] = float(inp[2])
        elif i == num_atoms + 9:
            NGX = int(inp[0])
            NGY = int(inp[1])
            NGZ = int(inp[2])
            Potential = numpy.zeros(shape=(NGX * NGY * NGZ))
            # Read in the potential data
            upper_limit = (int(NGX * NGY * NGZ / 5) +
                           np.mod(NGX * NGY * NGZ, 5))

    ##_print_boom(quiet=quiet)
    print("Average of the potential = ", numpy.average(Potential))

    lattice = lattice * scale_factor

    return Potential, NGX, NGY, NGZ, lattice
#------------------------------------------------------------------------------

def density_2_grid(Density, nx, ny, nz, Charge=False, Volume=1):
    """
    Convert density data to a 3D grid.

    Parameters:
        Density (numpy.ndarray): 1D array representing the density data.
        nx (int): Number of grid points along the x-axis.
        ny (int): Number of grid points along the y-axis.
        nz (int): Number of grid points along the z-axis.
        Charge (bool, optional): If True, convert charge density to the number of electrons.
                                 Default is False.
        Volume (int or float, optional): Volume of the grid cell. Used to convert charge density to electrons.
                                         Default is 1.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: 3D array representing the potential grid.
            - float: Total number of electrons in the grid (if Charge is True).
    
    Example:
        >>> Density = np.random.rand(NGX * NGY * NGZ)  # Replace this with actual density data
        >>> nx, ny, nz = NGX, NGY, NGZ
        >>> Charge = False  # Set to True if Density represents charge density
        >>> Volume = 1.0  # Volume of the grid cell (if Charge is True)
        >>> potential_grid, total_electrons = density_2_grid(Density, nx, ny, nz, Charge, Volume)
        >>> print("Potential Grid:")
        >>> print(potential_grid)
        >>> if Charge:
            >>> print("Total Electrons:", total_electrons)
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

def density_2_grid_gulp(Density, nx, ny, nz):
    """
    Convert density data to a 3D grid in the GULP format.

    Parameters:
        Density (numpy.ndarray): 1D array representing the density data.
        nx (int): Number of grid points along the x-axis.
        ny (int): Number of grid points along the y-axis.
        nz (int): Number of grid points along the z-axis.

    Returns:
        numpy.ndarray: 3D array representing the potential grid in GULP format.

    Example:
        >>> Density = np.random.rand(nx * ny * nz)  # Replace this with actual density data
        >>> potential_grid_gulp = density_2_grid_gulp(Density, nx, ny, nz)
        >>> print("Potential Grid in GULP format:")
        >>> print(potential_grid_gulp)
    """
    l = 0
    Potential_grid = np.zeros(shape=(nx,ny,nz))
    total_electrons = 0
    is_CHGCAR = True
    for k in range(nx):
        for j in range(ny):
            for i in range(nz):
                Potential_grid[k,j,i] = Density[l]
                l = l + 1
    return Potential_grid

#------------------------------------------------------------------------------
def read_gulp_potential(gulpfile='gulp.out'):
    """
    Read electrostatic potential data from a GULP output file.

    Parameters:
        gulpfile (str, optional): Path to the GULP output file (gulp.out). Default is 'gulp.out'.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: 1D array representing the electrostatic potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - numpy.ndarray: 3x3 array representing the Cartesian lattice vectors.
    
    Example:
        >>> gulpfile = 'path/to/gulp.out'
        >>> potential, NGX, NGY, NGZ, lattice = read_gulp_potential(gulpfile)
        >>> print("Electrostatic Potential Data:", potential)
        >>> print("Number of Grid Points (NGX, NGY, NGZ):", NGX, NGY, NGZ)
        >>> print("Lattice Vectors:")
        >>> print(lattice)
    """
    potential = []

    try:
        file_handle=open(gulpfile)
    except IOError:
        print("File not found or path is incorrect")

    lines = file_handle.readlines()
    for n, line in enumerate(lines):
        if line.rfind('Cartesian lattice vectors') > -1:
            lattice = np.zeros(shape=(3, 3))
            for r in range(3):
                lattice[r] = lines[n + 2 + r].split()
            break

    for n, line in enumerate(lines):
        if line.rfind('Electrostatic potential on a grid') > -1:
            NGX = int(lines[n + 3].split()[3])
            NGY = int(lines[n + 3].split()[5])
            NGZ = int(lines[n + 3].split()[7])
            break

    for n, line in enumerate(lines):
        if line.rfind('Electrostatic potential on a grid') > -1:
            for k in reversed(range(9, NGX*NGY*NGZ + 9)):
                potential.append(float(lines[n + k].split()[3]))


    return np.asarray(potential), NGX, NGY, NGZ, lattice


#------------------------------------------------------------------------------

def GCD(a,b):
    """
    Compute the Greatest Common Divisor (GCD) of two integers a and b.

    Parameters:
        a (int): First integer.
        b (int): Second integer.

    Returns:
        int: The Greatest Common Divisor of a and b.
    
    Example:
        >>> a = 36
        >>> b = 48
        >>> gcd = GCD(a, b)
        >>> print("GCD of", a, "and", b, "is:", gcd)
    """
    a = abs(a)
    b = abs(b)
    while a:
        a, b = (b % a), a
    return b
#------------------------------------------------------------------------------

def GCD_List(list):
    """
    Compute the Greatest Common Divisor (GCD) of a list of integers.

    Parameters:
        lst (list): List of integers.

    Returns:
        int: The Greatest Common Divisor of the elements in the list.
    
    Example:
        >>> numbers = [24, 36, 60]
        >>> gcd = GCD_List(numbers)
        >>> print("GCD of", numbers, "is:", gcd)
    """
    return reduce(GCD, list)
#------------------------------------------------------------------------------
def inverse_participation_ratio(density):
    """
    Calculate the inverse participation ratio (IPR) for a given density.

    Parameters:
        density (list or numpy.ndarray): List or 1D array representing the density data.

    Returns:
        float: The inverse participation ratio value.
    
    Example:
        >>> density = np.array([0.2, 0.4, 0.6, 0.8])
        >>> ipr = inverse_participation_ratio(density)
        >>> print("Inverse Participation Ratio (IPR) for the density:", ipr)
    """
    sq = sum(i**2 for i in density)
    fr = sum(i**4 for i in density)
    ifr = 1 / (len(density) * fr)
    isq = 1 / (len(density) * sq)
    return fr / sq**2
