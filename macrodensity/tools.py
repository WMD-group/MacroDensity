#! /usr/bin/env python
""" 
macrodensity.tools contains functions to read and manipulate the electronic 
density data from a material.
"""

from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from macrodensity.density import (
    density_2_grid,
    travelling_volume_average,
    volume_average,
)
from macrodensity.io import get_band_extrema, read_vasp_density
from macrodensity.utils import matrix_2_abc, vector_2_abscissa


def bulk_interstitial_alignment(
    interstices: list,
    outcar: str = "OUTCAR",
    locpot: str = "LOCPOT",
    cube_size: list = [2, 2, 2],
    print_output: bool = True,
) -> (float, float, float):
    """
    Calculate the aligned band energies for a bulk material with interstitial sites.

    This function calculates the aligned valence band (VB) and conduction band (CB) energies
    by considering the effect of interstitial sites on the electronic structure of the bulk material.

    Parameters:
        interstices (:obj:`list` of tuples): A list of tuples representing the
        coordinates of the interstitial sites for which the aligned band energies
        will be calculated.

        outcar (:obj:`str`, optional): The filename of the OUTCAR file containing
        electronic band structure information. Default is "OUTCAR".

        locpot (:obj:`str`, optional): The filename of the LOCPOT file containing the
        electronic density information. Default is "LOCPOT".

        cube_size (:obj:`list` of int, optional): The size of the cube (in grid points)
        around each interstitial site used to calculate the local potential.
        Default is [2, 2, 2].

        print_output (:obj:`bool`, optional): Whether to print the intermediate and final results.
        Default is True.

    Returns:
        :obj:`tuple`: A tuple containing the aligned VB energy, aligned CB energy,
        and a list of interstitial variances.
        The variances represent the deviation of the potential from the
        reference state at each interstitial site.

    Output:
        Aligned Valence Band, Aligned Conduction Band, Interstitial variances
    """

    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(locpot, quiet=True)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    grid_pot, electrons = density_2_grid(
        vasp_pot, NGX, NGY, NGZ, config="VASP"
    )

    ## GETTING BAND EDGES
    band_extrema = get_band_extrema(outcar)
    VB_eigenvalue = band_extrema[0]
    CB_eigenvalue = band_extrema[1]

    ## CALCULATING REFERENCE STATE
    interstitial_potentials = []
    interstitial_variances = []
    for interstice in interstices:
        locpot_extract = volume_average(
            origin=interstice,
            cube=cube_size,
            grid=grid_pot,
            nx=NGX,
            ny=NGY,
            nz=NGZ,
        )
        interstitial_potentials.append(locpot_extract[0])
        interstitial_variances.append(locpot_extract[1])

    ## CALCULATING ALIGNED BAND ENERGIES
    sum_interstitial_potential = 0
    for ele in interstitial_potentials:
        sum_interstitial_potential += ele
    average_interstitial_potential = sum_interstitial_potential / len(
        interstitial_potentials
    )
    VB_aligned = round(VB_eigenvalue - average_interstitial_potential, 2)
    CB_aligned = round(CB_eigenvalue - average_interstitial_potential, 2)

    ## PRINTING
    if print_output:
        print("Reading band edges from file: " + str(outcar))
        print("Reading potential from file: " + str(locpot))
        print("Interstital variances: " + str(interstitial_variances))
        print("VB_aligned      CB_aligned")
        print("--------------------------------")
        print(VB_aligned, "         ", CB_aligned)

    return (VB_aligned, CB_aligned, interstitial_variances)


def _find_active_space(
    cube_size: list,
    cube_origin: list,
    tolerance: float = 1e-4,
    input_file="LOCPOT",
    print_output=True,
    show_plot=True,
) -> tuple:
    """
    Plot the active space (vacuum and non-vacuum regions) based on potential variations.

    This function analyzes the potential variations within the specified cubes of the
    given size and determines whether each cube belongs to the vacuum or non-vacuum
    region based on the provided tolerance.
    This function also plots the cube potentials of vacuum and non vacuum cubes.

    Parameters:
        cube_size (list of int): The size of the cubes in units of mesh points (NGX/Y/Z)
        for analysis.

        cube_origin (list of float): The starting point (origin) of the cubes in
        fractional coordinates (range [0, 1]).

        tolerance (float, optional): The cutoff variance value to
        distinguish vacuum from non-vacuum cubes. Default is 1E-4.

        input_file (str, optional): The file with VASP output for potential.
        Default is 'LOCPOT'.

        print_output (bool, optional): Whether to print the analysis results.
        Default is True.

    Returns:
        dict: A dictionary containing the potentials for the vacuum and non-vacuum regions.

    Note:
        The function calculates the potential variation within each cube and
        compares it to the tolerance value. Cubes with potential variations below the
        tolerance are considered vacuum regions, while others are non-vacuum regions.

    Example:
        >>> cube_size = [2, 2, 2]
        >>> cube_origin = [0.0, 0.0, 0.0]
        >>> find_active_space(cube_size, cube_origin, tolerance=1E-5)

    """
    # Get potential
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    grid_pot, electrons = density_2_grid(
        vasp_pot, NGX, NGY, NGZ, config="VASP"
    )
    cutoff_variance = tolerance

    # Distinguish vacuum from non-vacuum
    vacuum = []
    vac_pot = []
    non_vacuum = []
    nvac_pot = []

    for i in range(0, NGX, cube_size[0]):
        for j in range(0, NGY, cube_size[1]):
            for k in range(0, NGZ, cube_size[2]):
                sub_origin = [float(i) / NGX, float(j) / NGY, float(k) / NGZ]
                cube_pot, cube_var = volume_average(
                    origin=sub_origin,
                    cube=cube_size,
                    grid=grid_pot,
                    nx=NGX,
                    ny=NGY,
                    nz=NGZ,
                    travelled=[0, 0, 0],
                )
                if cube_var <= cutoff_variance:
                    vacuum.append(sub_origin)
                    vac_pot.append(cube_pot)
                else:
                    non_vacuum.append(sub_origin)
                    nvac_pot.append(cube_pot)

    if print_output:
        len_vac, len_non_vac = len(vacuum), len(non_vacuum)
        print("Number of vacuum cubes: ", len_vac)
        print("Number of non-vacuum cubes: ", len_non_vac)
        print(
            "Percentage of vacuum cubes: ",
            round((len_vac / (len_vac + len_non_vac) * 100.0), 1),
            "%",
        )
        print(
            "Percentage of non-vacuum cubes: ",
            round((len_non_vac / (len_vac + len_non_vac) * 100.0), 1),
            "%",
        )

    return {
        "Vacuum": vacuum,
        "Vacuum Potential": vac_pot,
        "Non-vacuum": non_vacuum,
        "Non-vacuum Potential": nvac_pot,
    }


def bulk_vac(bulk: np.ndarray, slab: np.ndarray) -> np.ndarray:
    """
    Subtract potentials between a bulk dataset and a slab dataset based on their positions.

    Parameters:
        bulk (:obj:`numpy.ndarray`): The dataset containing bulk potential values in the format (x, potential).

        slab (:obj:`numpy.ndarray`): The dataset containing slab potential values in the format (x, potential).

    Returns:
        new_bulk (:obj:`numpy.ndarray`): The resulting bulk dataset with matching positions of the slab dataset.

    Example:
        >>> bulk = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> slab = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> result = bulk_vac(bulk, slab)
        >>> print(result)

    """
    new_bulk = np.zeros(shape=(len(slab), 2))
    i = -1
    for s_pot in slab:
        i = i + 1
        found = False
        for j in range(len(bulk)):
            if s_pot[0] <= bulk[j, 0] and s_pot[0] > bulk[j - 1, 0]:
                new_bulk[i, :] = bulk[j, :]
                found = True
        if found == False:
            new_bulk[i, 0] = s_pot[0]
            new_bulk[i, 1] = 0

    return new_bulk


def match_resolution(A: np.ndarray, B: np.ndarray) -> tuple:
    """
    Match the resolutions of two datasets by cubic spline interpolation.

    Parameters:
        A (:obj:`numpy.ndarray`): The first dataset containing potential values in the format (x, potential).

        B (:obj:`numpy.ndarray`): The second dataset containing potential values in the format (x, potential).

    Returns:
        A_new (:obj:`numpy.ndarray`): The first dataset with matched resolution and interpolated values.

        B_new (:obj:`numpy.ndarray`): The second dataset with matched resolution and interpolated values.

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> B = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> result_A, result_B = match_resolution(A, B)
        >>> print(result_A)
        >>> print(result_B)
    """
    np.append(A, A[0, :])
    np.append(B, B[0, :])
    resolution_a = (max(A[:, 0]) - min(A[:, 0])) / len(A)
    resolution_b = (max(B[:, 0]) - min(B[:, 0])) / len(B)
    new_resolution = min(resolution_a, resolution_b) / 3
    # Generate the function f for each spline
    f_a = interp1d(A[:, 0], A[:, 1], kind="cubic")
    f_b = interp1d(B[:, 0], B[:, 1], kind="cubic")
    # Generate the new abscissa values, at new_resolution
    abscissa_a = np.arange(0, max(A[:, 0]), new_resolution)
    abscissa_b = np.arange(0, max(B[:, 0]), new_resolution)
    # New datasets
    A_new = np.zeros(shape=(len(abscissa_a), 2))
    B_new = np.zeros(shape=(len(abscissa_b), 2))
    A_new[:, 0] = abscissa_a
    B_new[:, 0] = abscissa_b
    A_new[:, 1] = f_a(abscissa_a)
    B_new[:, 1] = f_b(abscissa_b)

    return A_new, B_new


def spline_generate(A: np.ndarray, new_res_factor: int) -> np.ndarray:
    """
    Generate a new dataset with higher resolution using cubic spline interpolation.

    Parameters:
        A (:obj:`numpy.ndarray`): The dataset containing potential values in the format (x, potential).

        new_res_factor (:obj:`int`): The factor by which to increase the resolution.

    Returns:
        B (:obj:`numpy.ndarray`): The new dataset with higher resolution and interpolated values.

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> new_res_factor = 2
        >>> result = spline_generate(A, new_res_factor)
        >>> print(result)
    """
    resolution = (A[len(A) - 1, 0] - A[0, 0]) * new_res_factor / len(A)
    array_a = np.arange(min(A[:, 0]), max(A[:, 0]), resolution)
    f_a = interp1d(A[:, 0], A[:, 1], kind="cubic")
    # ius = interpolate.InterpolatedUnivariateSpline(A[:,0],A[:,1])
    S = f_a(array_a)
    B = np.zeros(shape=(len(A) // new_res_factor, 2))
    for i in range(len(B)):
        B[i, 0] = i * resolution + A[0, 0]
        B[i, 1] = S[i]

    return B


def matched_spline_generate(
    A: np.ndarray, B: np.ndarray, V_A: np.ndarray, V_B: np.ndarray
) -> tuple:
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
    res_a = length_A / (len(A))
    res_b = length_B / (len(B))
    new_resolution = min(res_a, res_b)
    # Create an array containing indices of each potential point 0,1,2,....N
    array_a = np.arange(0, len(A))
    array_b = np.arange(0, len(B))
    # Generate the function f for each spline
    f_a = interp1d(array_a, A, kind="cubic")
    f_b = interp1d(array_b, B, kind="cubic")
    # Generate new arrays with the same resolution
    limits_a_new = np.arange(0, len(A))
    limits_b_new = np.arange(0, len(B))
    # Make the arrays
    A_new = f_a(limits_a_new)
    B_new = f_b(limits_b_new)
    # Convert to 2D arrays with AA in the first column
    TD_A = np.zeros(shape=(len(A_new), 2))
    TD_B = np.zeros(shape=(len(B_new), 2))
    res_a = length_A / float(len(A_new))
    res_b = length_B / float(len(B_new))
    for i in range(len(A_new)):
        TD_A[i, 1] = A[i]
        TD_A[i, 0] = i * res_a
    for i in range(len(B_new)):
        TD_B[i, 1] = B[i]
        TD_B[i, 0] = i * res_b
    return TD_A, TD_B


def scissors_shift(potential: np.ndarray, delta: float) -> np.ndarray:
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
        shifted_potential[i, 0] = potential[i, 0]
        shifted_potential[i, 1] = potential[i, 1] - delta

    return shifted_potential


def extend_potential(
    potential: np.ndarray, extension: float, vector: list
) -> np.ndarray:
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
    extended_potential = np.zeros(shape=(int(extension * len(potential)), 2))
    idx = 0
    diff = np.sqrt(vector.dot(vector))
    increment = diff / len(potential[:, 0])
    for i in range(int(extension)):
        for j in range(len(potential)):
            extended_potential[idx, 0] = potential[j, 0] + i * diff
            extended_potential[idx, 1] = potential[j, 1]
            idx = idx + 1

    if int(extension) != extension:  # For non-integer extensions
        i = i + 1
        over_shoot = extension - int(extension)
        for j in range(int(len(potential) * over_shoot)):
            extended_potential[idx, 0] = (
                potential[j, 0]
                + i * (max(potential[:, 0]) - min(potential[:, 0]))
                + increment * i
            )
            extended_potential[idx, 1] = potential[j, 1]
            idx = idx + 1

    return extended_potential


def sort_potential(potential: np.ndarray) -> np.ndarray:
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
    idx = sorted(potential[:, 0])
    sorted_potential = potential.copy()
    for i in range(len(idx)):
        for j in range(len(potential)):
            if potential[j, 0] == idx[i]:
                sorted_potential[i, 0] = idx[i]
                sorted_potential[i, 1] = potential[j, 1]

    return sorted_potential


def diff_potentials(
    potential_a: np.ndarray,
    potential_b: np.ndarray,
    start: float,
    end: float,
    tol: float = 0.04,
) -> np.ndarray:
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
    resolution = potential_a[0, 0] - potential_b[0, 1]
    new_potential = np.zeros(shape=((start - end) / resolution, 2))

    for i in range(len(potential_a)):
        if potential_a[i, 0] >= start and potential_a[i, 0] <= end:
            for j in range(len(potential_b)):
                if abs(potential_b[j, 0] - potential_a[i, 0]) <= tol:
                    new_potential[i, 1] = potential_a[i, 1] - potential_b[i, 1]
                    new_potential[i, 0] = potential_a[i, 0]

    return new_potential


def subs_potentials(A: np.ndarray, B: np.ndarray, tol: float) -> np.ndarray:
    """
    Subtract potentials between two datasets based on a tolerance value.

    Parameters:
        A (:obj:`numpy.ndarray`): The first dataset containing potential
        values in the format (x, potential).

        B (:obj:`numpy.ndarray`): The second dataset containing potential
        values in the format (x, potential).

        tol (:obj:`float`): The tolerance value for potential subtraction.

    Returns:
        C (:obj:`numpy.ndarray`): The resulting dataset containing the subtracted
        potentials in the format (x, potential).

    Example:
        >>> A = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
        >>> B = np.array([[0, 1], [1, 3], [2, 2], [3, 4]])
        >>> tolerance = 1e-2
        >>> C = subs_potentials(A, B, tolerance)
        >>> print(C)
    """
    C = A
    for i in range(len(A)):
        C[i, 0] = A[i, 0]
        if abs(A[i, 1] - B[i, 1]) <= tol:
            C[i, 1] = 0
        else:
            C[i, 1] = A[i, 1] - B[i, 1]

    return C


def translate_grid(
    potential: np.ndarray,
    translation: float,
    periodic: bool = False,
    vector: list = [0, 0, 0],
    boundary_shift: float = 0.0,
) -> np.ndarray:
    """
    Translates the grid points of a given potential by a specified translation
    along the vector direction.

    Parameters:
        potential (numpy.ndarray): Array containing potential data with shape (N, 2),
        where N is the number of grid points.

        translation (float): The amount of translation to apply to the grid points
        along the specified vector direction.

        periodic (bool, optional): Whether to apply periodic boundary conditions.
        Default is False.

        vector (list, optional): The direction vector for translation.
        Default is [0, 0, 0].

        boundary_shift (float, optional): The amount of shift to consider when
        applying periodic boundary conditions. Default is 0.0.

    Returns:
        numpy.ndarray: An array containing the translated potential data with shape (N, 2).

    Example:
        >>> # Sample potential data
        >>> potential = np.array([[0.0, 1.0], [0.5, 2.0], [1.0, 3.0]])
        >>> # Translate the grid by 0.2 along the x-direction
        >>> translated_potential = translate_grid(potential, 0.2)
        >>> print(translated_potential)
    """
    new_potential_trans = np.zeros((len(potential), 2))
    length = np.sqrt(vector.dot(vector))

    for i in range(len(potential)):
        new_potential_trans[i, 0] = potential[i, 0] + translation
        new_potential_trans[i, 1] = potential[i, 1]
        if periodic == True:
            new_potential_trans[i, 0] = new_potential_trans[
                i, 0
            ] - length * int(
                (new_potential_trans[i, 0] + boundary_shift) / length
            )

    if periodic == True:
        # Sort the numbers out if you have done periodic wrapping
        sorted_potential_trans = sort_potential(new_potential_trans)
    else:
        sorted_potential_trans = new_potential_trans

    return sorted_potential_trans


def create_plotting_mesh(
    NGX: int, NGY: int, NGZ: int, pc: np.ndarray, grad: np.ndarray
) -> np.ndarray:
    """
    Creates a plotting mesh based on the given grid data and plane coefficients.

    Parameters:
        NGX (int): Number of grid points along the x-direction.

        NGY (int): Number of grid points along the y-direction.

        NGZ (int): Number of grid points along the z-direction.

        pc (numpy.ndarray): Array containing plane coefficients with shape (4,).

        grad (numpy.ndarray): Array containing gradient data with shape (NGX, NGY, NGZ).

    Returns:
        numpy.ndarray: A 2D array representing the plotting mesh with shape (a, b),
        where 'a' and 'b' depend on the plane direction.

    Example:
        >>> # Sample grid data and plane coefficients
        >>> NGX, NGY, NGZ = 10, 10, 10
        >>> pc = np.array([0, 0, 1, 5])
        >>> grad = np.random.rand(NGX, NGY, NGZ)
        >>> # Create the plotting mesh
        >>> plotting_mesh = create_plotting_mesh(NGX, NGY, NGZ, pc, grad)
        >>> print(plotting_mesh)
    """

    a = 0
    b = 0
    c = 0
    p = ""

    if pc[0] == 0 and pc[1] == 0:
        a = NGX
        b = NGY
        p = "zzo"
        c = int(pc[3] / pc[2]) - 1
    elif pc[0] == 0 and pc[2] == 0:
        a = NGX
        b = NGZ
        p = "zoz"
        c = int(pc[3] / pc[1]) - 1
    elif pc[1] == 0 and pc[2] == 0:
        a = NGY
        b = NGZ
        p = "ozz"
        c = int(pc[3] / pc[0]) - 1
    plane = np.zeros(shape=(a, b))
    for x in range(a):
        for y in range(b):
            if p == "zzo":
                plane[x, y] = grad[x, y, c]
            else:
                pass

    return plane
