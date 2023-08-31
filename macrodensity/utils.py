"""General utility functions."""

from functools import reduce

import numpy as np


def matrix_2_abc(lattice: np.ndarray) -> (float, float, float, float, float, float):
    """
    Extract lattice parameters and vectors from a 3x3 matrix representing a lattice.

    Parameters:
        lattice (np.ndarray): 3x3 matrix representing the lattice.

    Returns:
        tuple: A tuple containing the lattice parameters a, b, c and lattice vectors a_vec, b_vec, c_vec.

    Example:
        >>> lattice = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
        >>> a, b, c, a_vec, b_vec, c_vec = matrix_2_abc(lattice)
        >>> print("lattice parameters:", a, b, c)
        >>> print("lattice vectors:")
        >>> print(a_vec)
        >>> print(b_vec)
        >>> print(c_vec)
    """
    a = np.sqrt(lattice[0,0]**2+lattice[0,1]**2+lattice[0,2]**2)
    b = np.sqrt(lattice[1,0]**2+lattice[1,1]**2+lattice[1,2]**2)
    c = np.sqrt(lattice[2,0]**2+lattice[2,1]**2+lattice[2,2]**2)

    a_vec = lattice[0,:]
    b_vec = lattice[1,:]
    c_vec = lattice[2,:]

    return a,b,c,a_vec,b_vec,c_vec


def vector_2_abscissa(vector:list, magnitude: float, dx: float, dy: float, dz: float) -> np.ndarray:
    """
    Convert a 3D vector to an array of abscissa values.

    Parameters:
        vector (list): 3D vector represented as (x, y, z).

        magnitude (float): Magnitude of the vector.

        dx (float): Spacing along the x-axis.

        dy (float): Spacing along the y-axis.

        dz (float): Spacing along the z-axis.

    Returns:
        np.ndarray: 1D array containing abscissa values based on the vector and spacing.

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


def GCD(a: int,b: int) -> int:
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


def GCD_List(list: list) -> int:
    """
    Compute the Greatest Common Divisor (GCD) of a list of integers.

    Parameters:
        list (list): List of integers.

    Returns:
        int: The Greatest Common Divisor of the elements in the list.
    
    Example:
        >>> numbers = [24, 36, 60]
        >>> gcd = GCD_List(numbers)
        >>> print("GCD of", numbers, "is:", gcd)
    """
    return reduce(GCD, list)


def get_volume(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
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
    volume = np.dot(a, np.cross(b, c))

    return volume


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


def numbers_2_grid(a: tuple, NGX: int, NGY: int, NGZ: int) -> np.ndarray:
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


def points_2_plane(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Calculates the plane coefficients from three points in space.

    Parameters:
        a (numpy.ndarray): First point with shape (3,).

        b (numpy.ndarray): Second point with shape (3,).

        c (numpy.ndarray): Third point with shape (3,).

    Returns:
        numpy.ndarray: An array containing the plane coefficients with shape (4,).

    Example:
        >>> # Sample points in space
        >>> a = np.array([0, 0, 0])
        >>> b = np.array([1, 0, 0])
        >>> c = np.array([0, 1, 0])
        >>> # Calculate plane coefficients
        >>> plane_coefficients = points_2_plane(a, b, c)
        >>> print(plane_coefficients)
    """
    coefficients = np.zeros(shape=(4))

    ca = c - a
    ba = b - a
    normal = np.cross(ba,ca)
    d = -np.dot(normal,a)
    A, B, C = normal[0], normal[1], normal[2]
    D = d
    coefficients = np.array([A, B, C, D])
    return coefficients


def subs_potentials(A: np.ndarray, B: np.ndarray, tol: float) -> np.ndarray:
    """
    Subtract potentials between two datasets based on a tolerance value.

    Parameters:
        A (:obj:`numpy.ndarray`): The first dataset containing potential values in the format (x, potential).

        B (:obj:`numpy.ndarray`): The second dataset containing potential values in the format (x, potential).

        tol (:obj:`float`): The tolerance value for potential subtraction.

    Returns:
        C (:obj:`numpy.ndarray`): The resulting dataset containing the subtracted potentials in the format (x, potential).

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