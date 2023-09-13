"""
macrodensity.density contains functions to calculate the electron
density of a structure, as well as functions to calculate the electrostatic potential
and the electrostatic potential gradients of a structure.
"""

from __future__ import division, print_function

import numpy as np

from macrodensity.io import read_vasp_density
from macrodensity.utils import matrix_2_abc


def gradient_magnitude(
    gx: np.ndarray, gy: np.ndarray, gz: np.ndarray
) -> np.ndarray:
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
                grad_mag[i, j, k] = np.sqrt(
                    gx[i, j, k] ** 2 + gy[i, j, k] ** 2 + gz[i, j, k] ** 2
                )

    return grad_mag


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
