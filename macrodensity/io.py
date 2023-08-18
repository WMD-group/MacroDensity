"""
Module containing functions to read output files
from VASP, GULP and FHI-AIMS.
"""

import math
from itertools import chain

import numpy as np


def read_gulp_potential(gulpfile: str='gulp.out') -> tuple:
    """
    Read electrostatic potential data from a GULP output file.

    Parameters:
        gulpfile (str, optional): Path to the GULP output file (gulp.out). Default is 'gulp.out'.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: 1D array representing the electrostatic potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - np.ndarray: 3x3 array representing the Cartesian lattice vectors.
    
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


def read_cube_density(FILE: str) -> np.ndarray:
    """
    Reads a cube density file and extracts relevant information.

    Parameters:
        FILE (str): The path to the cube density file.

    Returns:
        numpy.ndarray: A 3x3 numpy array representing the lattice.

    Example:
        >>> file_path = 'path/to/your/cube_density_file.cube'
        >>> # Read the cube density file and get the lattice
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
            

def _read_partial_density(FILE: str, use_pandas: bool, num_atoms: int, NGX: int, NGY: int, NGZ: int, spin: int=0) -> np.ndarray:
    """
    This function is used internally within the read_casp_parchg, reading partial density data from a VASP-PARCHG file.

    Parameters:
        FILE (str): Path to the CHGCAR-like file.

        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files). 
            If False, use np. If None, automatically use Pandas if available. Default is None.

        num_atoms (int): Total number of atoms in the system.

        NGX (int): Number of grid points along the x-axis.

        NGY (int): Number of grid points along the y-axis.

        NGZ (int): Number of grid points along the z-axis.

        spin (int, optional): If 0, read the first spin channel (total density). If 1, read the second spin channel (spin-up or spin-down). Default is 0.

    Returns:
        np.ndarray: 1D array representing the partial density data for the specified spin channel.
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

        coordinates = np.zeros(shape=(num_atoms, 3))
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
            density = np.fromiter(chain.from_iterable(density), float)

    return density


def _read_vasp_density_fromlines(lines: list) -> tuple:
    """
    Read density data from a list of lines (classic VASP-style format). This function is used internally within the read_vasp_density_classic function (add hyperlink to read_vasp_density_classic function)

    Parameters:
        lines (list): List of lines from the density file.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - np.ndarray: 3x3 array representing the lattice vectors.
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
            Coordinates = np.zeros(shape=(num_atoms,3))
        elif i >= 9 and i <= num_atoms + 8:
            Coordinates[i-9,0] = float(inp[0])
            Coordinates[i-9,1] = float(inp[1])
            Coordinates[i-9,2] = float(inp[2])
        elif i == num_atoms + 9:
            NGX = int(inp[0])
            NGY = int(inp[1])
            NGZ = int(inp[2])
            Potential = np.zeros(shape=(NGX * NGY * NGZ))
            # Read in the potential data
            upper_limit = (int(NGX * NGY * NGZ / 5) +
                           np.mod(NGX * NGY * NGZ, 5))

    print("Average of the potential = ", np.average(Potential))

    lattice = lattice * scale_factor

    return Potential, NGX, NGY, NGZ, lattice


def read_vasp_density_classic(FILE: str) -> tuple:
    """
    Read density data from a classic VASP-style file.

    Parameters:
        FILE (str): Path to the density file.

    Returns:
        tuple: A tuple containing:
            - np.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - np.ndarray: 3x3 array representing the lattice vectors.
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


def read_vasp_parchg(FILE: str, use_pandas: bool=None, quiet: bool=False, spin: bool=False) -> tuple:
    """
    Read density data or spin-polarized partial density data from a VASP PARCHG-like file.

    Parameters:
        FILE (str): Path to the PARCHG-like file.

        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files). If False, use np. If None, automatically use Pandas if available. Default is None.

        quiet (bool, optional): If True, suppress print statements during reading. Default is False.

        spin (bool, optional): If True, read spin-polarized partial densities. Default is False.

    Returns:
        tuple: A tuple containing:
            - list or np.ndarray: List containing np arrays representing the density data for each spin channel, or np.ndarray for the total density if spin is False.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - np.ndarray: 3x3 array representing the lattice vectors.

    Example:
        >>> FILE = "path/to/PARCHG-like-file"
        >>> density, NGX, NGY, NGZ, lattice = read_vasp_parchg(FILE)
        >>> if isinstance(density, list):
                print("Spin-polarized Partial Densities:")
                for i, spin_density in enumerate(density):
                    print(f"Spin {i+1} Density Data:", spin_density)
        >>> else:
                print("Total Density Data:", density)
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

        coordinates = np.zeros(shape=(num_atoms, 3))
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

    return density, NGX, NGY, NGZ, lattice


def read_vasp_density(FILE: str, use_pandas: bool=None, quiet: bool=False) -> tuple:
    """
    Read density data from a VASP CHGCAR-like file.

    Parameters:
        FILE (str): Path to the CHGCAR-like file.

        use_pandas (bool or None, optional): If True, use Pandas to read 3D data (recommended for large files). If False, use np. If None, automatically use Pandas if available. Default is None.

        quiet (bool, optional): If True, suppress print statements during reading. Default is False.

    Returns:
        tuple: 
            - np.ndarray: 1D array representing the potential data.
            - int: Number of grid points along the x-axis (NGX).
            - int: Number of grid points along the y-axis (NGY).
            - int: Number of grid points along the z-axis (NGZ).
            - np.ndarray: 3x3 array representing the lattice vectors.

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

        coordinates = np.zeros(shape=(num_atoms, 3))
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
            Potential = np.fromiter(chain.from_iterable(Potential), float)

    if not quiet:
        print("Average of the potential = ", np.average(Potential))

    return Potential, NGX, NGY, NGZ, lattice          
            

def get_band_extrema(input_file: str)->list:
    '''
    Get the valence band maximum and conduction band minimum from VASP OUTCAR.

    This function reads the VASP OUTCAR file and extracts the valence band maximum (VBM) and
    conduction band minimum (CBM). It also checks for partial occupancy and prints a warning
    message if found.

    Args:
        input_file (str): The path to the VASP OUTCAR file.

    Returns:
        list: A list containing the valence band maximum (VBM) and conduction band minimum (CBM).
              list[0] = VBM, list[1] = CBM.

    Example:
        >>> input_file = 'path/to/OUTCAR'
        >>> band_extrema = get_band_extrema(input_file)
        >>> print("Valence Band Maximum (VBM):", band_extrema[0])
        >>> print("Conduction Band Minimum (CBM):", band_extrema[1])
    '''
    lines = open(input_file, 'r').readlines()
    for line in lines:
        if line.rfind('NKPTS') > -1:
            nkpts = int(line.split()[3])
        if line.rfind('ISPIN') > -1:
            ispin = int(line.split()[2])
        if line.rfind('NELECT') > -1:
            nelect = float(line.split()[2])
    if ispin == 1:
        top_band = int(nelect/2)
    else:
        top_band = int(nelect)

    vbm = []
    cbm = []
    for i, line in enumerate(lines):
        if line.rfind('No.') > -1:
            vbm.append(lines[i + top_band].split()[1])
            cbm.append(lines[i + top_band + 1].split()[1])
            if (float(lines[i + top_band].split()[2]) != 1.00 and
                float(lines[i + top_band].split()[2]) != 2.000):
                print('Partial occupancy, be aware!',
                      lines[i + top_band].split()[2])

    return [float(max(vbm)), float(min(cbm))]