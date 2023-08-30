"""
macrodensity.density contains functions to calculate the electron density of a structure, as well as functions to calculate the electrostatic potential and the electrostatic potential gradients of a structure. 
"""

from __future__ import division, print_function

import numpy as np


def gradient_magnitude(gx: np.ndarray, gy: np.ndarray, gz: np.ndarray) -> np.ndarray:
    """
    Calculate the magnitude of the gradient at each point in a 3D field.

    Parameters:
        gx (np.ndarray): Gradient along the x-axis.

        gy (np.ndarray): Gradient along the y-axis.

        gz (np.ndarray): Gradient along the z-axis.

    Returns:
        np.ndarray: 3D array representing the magnitude of the gradient at each point.

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


def number_in_field(gradients: np.ndarray, cutoff: float) -> int:
    """
    Count the number of elements in a field that have a value greater than or equal to the cutoff.

    Parameters:
        gradients (np.ndarray): 3D array representing the field.

        cutoff (float): Threshold value for counting elements.

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


def element_vol(vol: float, nx: int, ny: int, nz: int) -> float:
    """
    Calculate the volume of each element in a 3D grid.

    Parameters:
        vol (float): Total volume of the 3D grid.

        nx (int): Number of elements along the x-axis.

        ny (int): Number of elements along the y-axis.

        nz (int): Number of elements along the z-axis.

    Returns:
        float: volume of each individual element in the grid.

    Example:
        >>> volume = 10.0
        >>> nx, ny, nz = 5, 5, 5
        >>> element_volume = element_vol(volume, nx, ny, nz)
        >>> print("volume of Each Element:", element_volume)

    """
    number_of_elements = nx * ny * nz
    ele_vol = vol / number_of_elements

    return ele_vol


def one_2_2d(array: np.ndarray, resolution: float, vector: np.ndarray) -> np.ndarray:
    """
    Transform a 1D array to a 2D array with abscissa values based on the given resolution and vector.

    Parameters:
        array (np.ndarray): 1D array to be transformed.

        resolution (float): Spacing between abscissa values.

        vector (np.ndarray): 3D vector used for the transformation.

    Returns:
        np.ndarray: 2D array with abscissa values and the corresponding Array values.

    Example:
        >>> Array = np.random.rand(10)
        >>> resolution = 0.5
        >>> vector = np.array([1, 2, 3])
        >>> transformed_array = one_2_2d(Array, resolution, vector)
        >>> print("Transformed Array:")
        >>> print(transformed_array)
    """
    length = np.sqrt(vector.dot(vector))
    new_array = np.zeros(shape=(len(array) - 1, 2))
    resolution = length / len(array)
    for i in range(len(array) - 1):
        new_array[i,0] = i*resolution
        new_array[i,1] = array[i]
    return new_array


def macroscopic_average(
    potential: np.ndarray, periodicity: float, resolution: float
) -> np.ndarray:
    """
    Calculate the macroscopic average of a 1D potential field with periodicity.

    Parameters:
        potential (np.ndarray): 1D array representing the potential field.

        periodicity (float): Periodicity of the field.

        resolution (float): Spacing between potential data points.

    Returns:
        np.ndarray: 1D array containing the macroscopic average of the potential field.

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

    print("Average of the average = ", np.average(macro_average))

    return macro_average


def volume_average(origin: tuple, cube: tuple, grid: np.ndarray, nx: int, ny: int, nz: int, travelled: list=[0, 0, 0]) -> tuple:
    """
    Calculate the volume average and variance of a cube in a 3D grid.

    Parameters:
        origin (tuple): Coordinates of the origin point.

        cube (tuple): Dimensions of the cube (x, y, z).

        grid (np.ndarray): 3D array representing the data grid.

        nx (int): Number of points along the x-axis in the grid.

        ny (int): Number of points along the y-axis in the grid.

        nz (int): Number of points along the z-axis in the grid.

        travelled (list, optional): Distance travelled from the origin in each direction (x, y, z). Default is [0, 0, 0].

    Returns:
        tuple: A tuple containing the volume average and variance of the cube.

    Example:
        >>> origin = (0.5, 0.5, 0.5)
        >>> cube = (3, 3, 3)
        >>> grid = np.random.rand(10, 10, 10)
        >>> nx, ny, nz = 10, 10, 10
        >>> travelled = [1, 2, 3]
        >>> avg, variance = volume_average(origin, cube, grid, nx, ny, nz, travelled)
        >>> print("volume Average:", avg)
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


def travelling_volume_average(grid: np.ndarray, cube: tuple, origin: tuple, vector: list, nx: int, ny: int, nz: int, magnitude: int) -> np.ndarray:
   """
    Calculate the volume average at multiple positions along a given vector.

    Parameters:
        grid (np.ndarray): 3D array representing the data grid.

        cube (tuple): Dimensions of the cube (x, y, z).

        origin (tuple): Coordinates of the origin point.

        vector (list): 3D vector representing the direction of travel.

        nx (int): Number of points along the x-axis in the grid.

        ny (int): Number of points along the y-axis in the grid.

        nz (int): Number of points along the z-axis in the grid.

        magnitude (int): Number of positions to travel along the vector.

    Returns:
        np.ndarray: 1D array containing the volume averages at each position along the vector.

    Example:
        >>> vector = (0.1, 0.2, 0.3)
        >>> magnitude = 5
        >>> travelling_avg = travelling_volume_average(grid, cube, origin, vector, nx, ny, nz, magnitude)
        >>> print("Travelling volume Average:")
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


def planar_average(grid: np.ndarray, nx: int, ny: int, nz: int, axis: str='z') -> np.ndarray:
    """
    Calculate the planar average of a 3D grid along a specified axis.

    Parameters:
        grid (np.ndarray): 3D array representing the data grid.

        nx (int): Number of points along the x-axis in the grid.

        ny (int): Number of points along the y-axis in the grid.

        nz (int): Number of points along the z-axis in the grid.

        axis (str, optional): Axis along which to calculate the average ('x', 'y', or 'z'). Default is 'z'.

    Returns:
        np.ndarray: 1D array containing the planar average along the specified axis.

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


def get_volume(a: np.ndarray,b: np.ndarray,c: np.ndarray) -> float:
    """
    Calculate the volume of a parallelepiped defined by three vectors a, b, and c.

    Parameters:
        a (np.ndarray): 1D array representing vector a.

        b (np.ndarray): 1D array representing vector b.

        c (np.ndarray): 1D array representing vector c.

    Returns:
        float: volume of the parallelepiped defined by the three vectors.

    Example:
        >>> a = np.array([1, 0, 0])
        >>> b = np.array([0, 1, 0])
        >>> c = np.array([0, 0, 1])
        >>> volume = get_volume(a, b, c)
        >>> print("volume of parallelepiped:", volume)
    """
    volume = np.dot(a,np.cross(b,c))

    return volume


def numbers_2_grid(a: tuple,NGX: int,NGY: int,NGZ: int) -> np.ndarray:
    """
    Convert fractional coordinates to grid point coordinates.

    Parameters:
        a (tuple): Fractional coordinates (x, y, z).

        NGX (int): Number of grid points along the x-axis.

        NGY (int): Number of grid points along the y-axis.

        NGZ (int): Number of grid points along the z-axis.

    Returns:
        np.ndarray: 1D array containing the grid point coordinates (x, y, z).

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


# TODO: Update variables here to be lower case following python convention
def density_2_grid(
    density: np.ndarray, 
    nx: int, 
    ny: int, 
    nz: int, 
    charge: bool=False, 
    volume: float=1, 
    Format: str = 'VASP'
) -> tuple:
    """
    Convert density data to a 3D grid.

    Parameters:
        density (np.ndarray): 1D array representing the density data.
        nx (int): Number of grid points along the x-axis.
        ny (int): Number of grid points along the y-axis.
        nz (int): Number of grid points along the z-axis.
        charge (bool, optional): If True, convert charge density to the number of electrons. Default is False.
        volume (float, optional): volume of the grid cell. Used to convert charge density to electrons. Default is 1.
        Format (str, optional): Format of the density data (e.g., 'VASP', 'GULP'). Default is 'VASP'.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: 3D array representing the potential grid.
            - float: Total number of electrons in the grid (if charge is True).

    Example:
        >>> density = np.random.rand(NGX * NGY * NGZ)  # Replace this with actual density data
        >>> nx, ny, nz = NGX, NGY, NGZ
        >>> charge = False  # Set to True if density represents charge density
        >>> volume = 1.0  # volume of the grid cell (if charge is True)
        >>> potential_grid, total_electrons = density_2_grid(density, nx, ny, nz, charge, volume)
        >>> print("Potential Grid:")
        >>> print(potential_grid)
        >>> if charge:
                print("Total Electrons:", total_electrons)
    """
    l = 0
    Potential_grid = np.zeros(shape=(nx, ny, nz))
    
    if Format.lower() == "gulp":
        for k in range(nx):
            for j in range(ny):
                for i in range(nz):
                    Potential_grid[k,j,i] = density[l]
                    l = l + 1
        return Potential_grid
    
    elif Format.lower() == "vasp":
        total_electrons = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    Potential_grid[i,j,k] = density[l]/volume
                    if charge == True:
                        # Convert the charge density to a number of electrons
                        point_volume = volume / (nx*ny*nz)
                        Potential_grid[i,j,k] = Potential_grid[i,j,k]*point_volume
                    total_electrons = total_electrons + density[l]
                    l = l + 1
                    
        total_electrons = total_electrons / (nx * ny * nz)
        if charge == True:
            print("Total electrons: ", total_electrons)    
        return Potential_grid, total_electrons
    
    else:
        raise ValueError("Invalid Format. Format must be 'VASP' or 'GULP'.")


def inverse_participation_ratio(density: np.ndarray) -> float:
    """
    Calculate the inverse participation ratio (IPR) for a given density.

    Parameters:
        density (np.ndarray): List or 1D array representing the density data.

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


def planar_average_charge(
    grid: np.ndarray,
    nx: int,
    ny: int,nz: int,
    vector: np.ndarray
) -> np.ndarray:

    a, b = 0, 0
    axis = ""
    a_vec, b_vec, c_vec = vector[0], vector[1], vector[2]

    print(a_vec)
    print(b_vec)
    print(c_vec)
    
    if (a_vec == 0).all() and (b_vec == 0).all():
        a, b = nx, ny
        axis = 'z'
        c = int(c_vec[2]) - 1
    elif (a_vec == 0).all() and (c_vec == 0).all():
        a, b = nx, nz
        axis = 'y'
        c = int(b_vec[1]) - 1
    elif (b_vec == 0).all() and (c_vec == 0).all():
        a, b = ny, nz
        axis = 'x'
        c = int(a_vec[0]) - 1
    else:
        raise ValueError("Invalid vector coefficients. Cannot determine plane direction.")

    average_charge = planar_average(grid, a, b, c, axis)

    return average_charge
