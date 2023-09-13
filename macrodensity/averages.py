"""
This module contains functions for calculating different types of averages from the dataset.
"""
from __future__ import division, print_function

import numpy as np

from macrodensity.io import read_gulp_potential, read_vasp_density
from macrodensity.utils import density_2_grid, matrix_2_abc


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
    period_points = int((periodicity / resolution))
    # Period points must be even
    if period_points % 2 != 0:
        period_points = period_points + 1

    length = len(potential)
    for i in range(length):
        start = i - int(period_points / 2)
        end = i + int(period_points / 2)
        if start < 0:
            start = start + length
            macro_average[i] = (
                macro_average[i]
                + sum(potential[0:end])
                + sum(potential[start:length])
            )
            macro_average[i] = macro_average[i] / period_points
        elif end >= length:
            end = end - length
            macro_average[i] = (
                macro_average[i]
                + sum(potential[start:length])
                + sum(potential[0:end])
            )
            macro_average[i] = macro_average[i] / period_points
        else:
            macro_average[i] = (
                macro_average[i] + sum(potential[start:end]) / period_points
            )

    print("Average of the average = ", np.average(macro_average))

    return macro_average


def volume_average(
    origin: tuple,
    cube: tuple,
    grid: np.ndarray,
    nx: int,
    ny: int,
    nz: int,
    travelled: list = [0, 0, 0],
) -> tuple:
    """
    Calculate the volume average and variance of a cube in a 3D grid.

    Parameters:
        origin (tuple): Coordinates of the origin point.

        cube (tuple): Dimensions of the cube (x, y, z).

        grid (np.ndarray): 3D array representing the data grid.

        nx (int): Number of points along the x-axis in the grid.

        ny (int): Number of points along the y-axis in the grid.

        nz (int): Number of points along the z-axis in the grid.

        travelled (list, optional): Distance travelled from the origin in each
            direction (x, y, z). Default is [0, 0, 0].

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
    n_origin[0] = int(origin[0] * nx)
    n_origin[1] = int(origin[1] * ny)
    n_origin[2] = int(origin[2] * nz)
    potential_cube = np.zeros(shape=(cube[0], cube[1], cube[2]))
    for x in range(0, cube[0]):
        for y in range(0, cube[1]):
            for z in range(0, cube[2]):
                # Assign the values of coordinates in the original grid
                xv = int(n_origin[0] + travelled[0] + x)
                yv = int(n_origin[1] + travelled[1] + y)
                zv = int(n_origin[2] + travelled[2] + z)
                # Minimum image convention
                zv = int(zv - nz * round(zv / nz))
                yv = int(yv - ny * round(yv / ny))
                xv = int(xv - nx * round(xv / nx))
                potential_cube[x, y, z] = grid[int(xv), int(yv), int(zv)]

    return potential_cube.mean(), np.var(potential_cube)


def spherical_average(
    cube_size: list,
    cube_origin: list,
    input_file: str = "LOCPOT",
    print_output: bool = True,
) -> (float, float):
    """
    Calculate the volume average of the electronic potential within a spherical region.

    This function calculates the volume average of the electronic potential within a
    spherical region defined by a specific size and origin.
    The size of the spherical region is specified by the cube_size
    parameter, which determines the number of mesh points along each direction (NGX/Y/Z).
    The origin of the sphere is given by the cube_origin parameter,
    specified in fractional coordinates. The function reads the electronic potential data
    from the specified input file (e.g., LOCPOT) and calculates the potential and variance
    within the spherical region.

    Parameters:
        cube_size (:obj:`list`): The size of the spherical region
        in units of mesh points (NGX/Y/Z).

        cube_origin (:obj:`list`): The origin of the spherical region
        in fractional coordinates.

        input_file (:obj:`str`, optional): The filename of the file containing
        the electronic potential (e.g., LOCPOT). Default is 'LOCPOT'.

        print_output (:obj:`bool`, optional): If True, the function prints
        the calculated potential and variance. Default is True.

    Returns:
        :obj:`tuple`: A tuple containing the volume-averaged potential and
        the variance within the spherical region.

    Outputs:
        cube_potential, cube_variance
    """
    ## GETTING POTENTIAL
    if "cube" in input_file:
        cube_pot, NGX, NGY, NGZ, lattice = cube.read_cube_data(input_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot, electrons = density_2_grid(
            cube_pot, NGX, NGY, NGZ, config="CUBE"
        )
    elif (
        "vasp" in input_file
        or "LOCPOT" in input_file
        or "CHGCAR" in input_file
    ):
        vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot, electrons = density_2_grid(
            vasp_pot, NGX, NGY, NGZ, config="VASP"
        )
    elif "gulp" in input_file or ".out" in input_file:
        gulp_pot, NGX, NGY, NGZ, lattice = read_gulp_potential(input_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot, electrons = density_2_grid(
            gulp_pot, NGX, NGY, NGZ, config="GULP"
        )
    else:
        raise ValueError(
            "Invalid input file. File must be in VASP, GULP, or CUBE configuration."
        )
    cube = cube_size
    origin = cube_origin
    travelled = [0, 0, 0]
    cube_pot, cube_var = volume_average(
        origin=cube_origin,
        cube=cube_size,
        grid=grid_pot,
        nx=NGX,
        ny=NGY,
        nz=NGZ,
        travelled=[0, 0, 0],
    )

    ## PRINTING
    if print_output == True:
        print("Potential            Variance")
        print("--------------------------------")
        print(cube_pot, "   ", cube_var)
    return cube_pot, cube_var


def travelling_volume_average(
    grid: np.ndarray,
    cube: tuple,
    origin: tuple,
    vector: list,
    nx: int,
    ny: int,
    nz: int,
    magnitude: int,
) -> np.ndarray:
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
        np.ndarray: 1D array containing the volume averages at each position
            along the vector.

    Example:
        >>> vector = (0.1, 0.2, 0.3)
        >>> magnitude = 5
        >>> travelling_avg = travelling_volume_average(grid, cube, origin, vector,
            nx, ny, nz, magnitude)
        >>> print("Travelling volume Average:")
        >>> print(travelling_avg)
    """
    plotting_average = np.zeros(shape=(magnitude))
    i = 0
    while i < magnitude:
        travelled = np.multiply(i, vector)
        plotting_average[i], varience = volume_average(
            origin, cube, grid, nx, ny, nz, travelled
        )
        i = i + 1

    return plotting_average


def planar_average(
    grid: np.ndarray, nx: int, ny: int, nz: int, axis: str = "z"
) -> np.ndarray:
    """
    Calculate the planar average of a 3D grid along a specified axis.

    Parameters:
        grid (np.ndarray): 3D array representing the data grid.

        nx (int): Number of points along the x-axis in the grid.

        ny (int): Number of points along the y-axis in the grid.

        nz (int): Number of points along the z-axis in the grid.

        axis (str, optional): Axis along which to calculate the average
            ('x', 'y', or 'z'). Default is 'z'.

    Returns:
        np.ndarray: 1D array containing the planar average along the specified axis.

    Example:
        >>> axis = 'z'
        >>> planar_avg = planar_average(grid, nx, ny, nz, axis)
        >>> print("Planar Average along axis", axis)
        >>> print(planar_avg)
    """
    if axis == "x":
        x_plane = np.zeros(shape=(ny, nz))
        average = np.zeros(shape=(nx))
        for x_value in range(nx):
            x_plane[:, :] = grid[x_value, :, :]
            average[x_value] = x_plane.mean()
    if axis == "y":
        average = np.zeros(shape=(ny))
        y_plane = np.zeros(shape=(nx, nz))
        for y_value in range(ny):
            y_plane[:, :] = grid[:, y_value, :]
            average[y_value] = y_plane.mean()
    if axis == "z":
        average = np.zeros(shape=(nz))
        z_plane = np.zeros(shape=(nx, ny))
        for z_value in range(nz):
            z_plane[:, :] = grid[:, :, z_value]
            average[z_value] = z_plane.mean()

    return average
