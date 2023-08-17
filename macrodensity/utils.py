"""General utility functions."""

from functools import reduce
import numpy as np

def matrix_2_abc(Lattice: np.ndarray) -> (float, float, float, float, float, float):
    """
    Extract lattice parameters and vectors from a 3x3 matrix representing a lattice.

    Parameters:
        Lattice (np.ndarray): 3x3 matrix representing the lattice.

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