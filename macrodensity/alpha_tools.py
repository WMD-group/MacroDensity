#! /usr/bin/env python

def bulk_interstitial_alignment(interstices,outcar="OUTCAR",locpot="LOCPOT",cube_size=[2,2,2],print_output=True):
    """
    Calculate the aligned band energies for a bulk material with interstitial sites.

    This function calculates the aligned valence band (VB) and conduction band (CB) energies
    by considering the effect of interstitial sites on the electronic structure of the bulk material.

    Parameters:
        interstices (list of tuples): A list of tuples representing the coordinates of the interstitial sites
                                     for which the aligned band energies will be calculated.

        outcar (str, optional): The filename of the OUTCAR file containing electronic band structure information.
                                Default is "OUTCAR".

        locpot (str, optional): The filename of the LOCPOT file containing the electronic density information.
                                Default is "LOCPOT".

        cube_size (list of int, optional): The size of the cube (in grid points) around each interstitial site
                                           used to calculate the local potential. Default is [2, 2, 2].

        print_output (bool, optional): Whether to print the intermediate and final results. Default is True.

    Returns:
        tuple: A tuple containing the aligned VB energy, aligned CB energy, and a list of interstitial variances.
               The variances represent the deviation of the potential from the reference state at each interstitial site.

    Example:
        >>> interstices = [(0.25, 0.25, 0.25), (0.5, 0.5, 0.5), (0.75, 0.75, 0.75)]
        VB_aligned, CB_aligned, interstitial_variances = bulk_interstitial_alignment(interstices)
    """
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, volume_average
    from macrodensity.vasp_tools import get_band_extrema
    
    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(locpot,quiet=True)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    ## GETTING BAND EDGES
    band_extrema = get_band_extrema(outcar)
    VB_eigenvalue = band_extrema[0]
    CB_eigenvalue = band_extrema[1]

    ## CALCULATING REFERENCE STATE
    interstitial_potentials = []
    interstitial_variances = []
    for interstice in interstices:
        locpot_extract = volume_average(origin=interstice,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ)
        interstitial_potentials.append(locpot_extract[0])
        interstitial_variances.append(locpot_extract[1])

    ## CALCULATING ALIGNED BAND ENERGIES
    sum_interstitial_potential = 0
    for ele in interstitial_potentials:
        sum_interstitial_potential += ele
    average_interstitial_potential = sum_interstitial_potential/len(interstitial_potentials)
    VB_aligned = round(VB_eigenvalue - average_interstitial_potential,2)
    CB_aligned = round(CB_eigenvalue - average_interstitial_potential,2)

    ## PRINTING
    if print_output == True:
        print("Reading band edges from file: "+str(outcar))
        print("Reading potential from file: "+str(locpot))
        print("Interstital variances: "+str(interstitial_variances))
        print("VB_aligned      CB_aligned")
        print("--------------------------------")
        print(VB_aligned,"         ",CB_aligned)

    return VB_aligned, CB_aligned, interstitial_variances

#------------------------------------------------------------------------------

def plot_active_space(cube_size,cube_origin,tolerance=1E-4,input_file='LOCPOT',print_output=True):
    '''
    Plot the active space (vacuum and non-vacuum regions) based on potential variations.

    This function analyzes the potential variations within the specified cubes of the given size
    and determines whether each cube belongs to the vacuum or non-vacuum region based on the provided tolerance.

    Parameters:
        cube_size (list of int): The size of the cubes in units of mesh points (NGX/Y/Z) for analysis.
        cube_origin (list of float): The starting point (origin) of the cubes in fractional coordinates (range [0, 1]).
        tolerance (float, optional): The cutoff variance value to distinguish vacuum from non-vacuum cubes. Default is 1E-4.
        input_file (str, optional): The file with VASP output for potential. Default is 'LOCPOT'.
        print_output (bool, optional): Whether to print the analysis results. Default is True.

    Returns:
        tuple: A tuple containing the number of cubes identified as vacuum and non-vacuum regions.

    Note:
        The function calculates the potential variation within each cube and compares it to the tolerance value.
        Cubes with potential variations below the tolerance are considered vacuum regions, while others are non-vacuum regions.

    Example:
        >>> cube_size = [2, 2, 2]
        >>> cube_origin = [0.0, 0.0, 0.0]
        >>> plot_active_space(cube_size, cube_origin, tolerance=1E-5)
        Number of vacuum cubes:  20
        Number of non-vacuum cubes:  4
        Percentage of vacuum cubes:  83.33333333333334
        Percentage of non-vacuum cubes:  16.666666666666664
        (20, 4)
    '''
    
    import math
    import numpy as np
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, volume_average

    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    cutoff_varience = tolerance
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    travelled = [0,0,0]

    ## DISTNGUISHING VACCUM FROM NON_VACUUM
    vacuum = []
    non_vacuum = []
    for i in range(0,NGX,cube_size[0]):
        for j in range(0,NGY,cube_size[1]):
            for k in range(0,NGZ,cube_size[2]):
                sub_origin = [float(i)/NGX,float(j)/NGY,float(k)/NGZ]
                cube_pot, cube_var = volume_average(origin=sub_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])
                if cube_var <= cutoff_varience:
                    vacuum.append(sub_origin)
                else:
                    non_vacuum.append(sub_origin)
    if print_output == True:
        print("Number of vacuum cubes: ", len(vacuum))
        print("Number of non-vacuum cubes: ", len(non_vacuum))
        print("Percentage of vacuum cubes: ",(float(len(vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
        print("Percentage of non-vacuum cubes: ",(float(len(non_vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
    return len(vacuum), len(non_vacuum)

#------------------------------------------------------------------------------

def plot_gulp_potential(lattice_vector,input_file='gulp.out',output_file='GulpPotential.csv',img_file='GulpPotential.png',new_resolution = 3000):
    '''
    Plot GULP potential analysis results.

    This function reads GULP potential data from the specified input file and performs a planar average
    as well as a macroscopic average. It then generates a plot comparing the planar and macroscopic averages,
    along with an interpolated curve for better visualization.

    Parameters:
        lattice_vector (float): The repeating unit over which the potential is averaged to get the macroscopic average (Angstroms).
        input_file (str, optional): The filename of the GULP output file containing potential data. Default is 'gulp.out'.
        output_file (str, optional): Name of the output data file to store the planar average. Default is 'GulpPotential.csv'.
        img_file (str, optional): Name of the output image file for the generated plot. Default is 'GulpPotential.png'.
        new_resolution (int, optional): The number of grid points used for the interpolated curve. Default is 3000.

    Returns:
        tuple: A tuple containing the planar average, macroscopic average, and the interpolated potential curve.

    Example:
        >>> lattice_vector = 20.0
        >>> input_file = 'gulp.out'
        >>> output_file = 'GulpPotential.csv'
        >>> img_file = 'GulpPotential.png'
        >>> new_resolution = 5000
        >>> planar_avg, macro_avg, interpolated_potential = plot_gulp_potential(lattice_vector, input_file, output_file, img_file, new_resolution)
        ... (plot generated and data saved to 'GulpPotential.png' and 'GulpPotential.csv')
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy.interpolate import interp1d
    from macrodensity.density_tools import matrix_2_abc, planar_average, macroscopic_average,density_2_grid_gulp,read_gulp_potential

    # GET POTENTIAL
    pot, NGX, NGY, NGZ, Lattice = read_gulp_potential(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot = density_2_grid_gulp(pot, NGX, NGY, NGZ)

    ## POTENTIAL PLANAR AVERAGE
    planar = planar_average(grid_pot, NGX, NGY, NGZ)
    np.savetxt(output_file, planar)

    ## MACROSCOPIC AVERAGE
    new_abscissa = np.linspace(0, NGZ - 1, new_resolution)
    f = interp1d(range(NGZ), planar, kind='cubic')
    interpolated_potential = [f(i) for i in new_abscissa]
    macro  = macroscopic_average(planar, lattice_vector, vector_c/new_resolution)

    ## PLOTTING
    fig, ax1 = plt.subplots()
    ax1.set_ylabel('V/V')
    ax1.set_xlabel('Grid Position')
    ax1.plot(planar,linestyle = ' ',marker = 'o')
    ax1.plot(macro)
    ax2 = ax1.twiny()
    ax2.plot(interpolated_potential)
    plt.savefig(img_file)

    ## SAVING
    df = pd.DataFrame.from_dict({'Planar':planar,'Macroscopic':macro,'Interpolated':interpolated_potential},orient='index')
    df = df.transpose()
    df.to_csv(output_file)
    return planar, macro, interpolated_potential

#------------------------------------------------------------------------------

def plot_on_site_potential(species,sample_cube,potential_file='LOCPOT',coordinate_file='POSCAR',output_file='OnSitePotential.csv',img_file='OnSitePotential.png'):
    '''
    Plot on-site electrostatic potential for a specific species.

    This function reads the electronic potential from the specified VASP output file (LOCPOT) and the atomic coordinates
    from the POSCAR file. It then calculates the on-site electrostatic potential for the specified species and generates
    a histogram to visualize its distribution across the sample cube.

    Parameters:
        species (str): The chemical symbol of the species whose on-site potential is of interest.
        sample_cube (list of int): The size of the sampling cube in units of mesh points (NGX/Y/Z).
        potential_file (str, optional): The filename of the VASP output file containing the electronic potential (LOCPOT).
                                       Default is 'LOCPOT'.
        coordinate_file (str, optional): The filename of the POSCAR file containing atomic coordinates.
                                         Default is 'POSCAR'.
        output_file (str, optional): Name of the output data file to store the on-site potential values. Default is 'OnSitePotential.csv'.
        img_file (str, optional): Name of the output image file for the histogram plot. Default is 'OnSitePotential.png'.

    Returns:
        list: A list containing the on-site electrostatic potential values for the specified species.

    Example:
        >>> species = 'O'
        >>> sample_cube = [2, 2, 2]
        >>> potential_file = 'LOCPOT'
        >>> coordinate_file = 'POSCAR'
        >>> output_file = 'OnSitePotential.csv'
        >>> img_file = 'OnSitePotential.png'
        >>> on_site_potential = plot_on_site_potential(species, sample_cube, potential_file, coordinate_file, output_file, img_file)
        ... (plot generated and on-site potential data saved to 'OnSitePotential.png' and 'OnSitePotential.csv')
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import csv
    import ase
    import pandas as pd
    from ase.io import write
    from ase.io import vasp
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average,volume_average

    ## GETTING POTENTIALS
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(potential_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    coords = ase.io.vasp.read_vasp(coordinate_file)
    scaled_coords = coords.get_scaled_positions()
    symbols = coords.get_chemical_symbols()
    ox_coords = []

    for i, atom in enumerate(coords):
        if symbols[i] == species:
            ox_coords.append(scaled_coords[i])
    grid_position = np.zeros(shape=(3))
    potentials_list = []
    i = 0
    num_bins = 20
    for coord in ox_coords:
        i = i + 1
        grid_position[0] = coord[0]
        grid_position[1] = coord[1]
        grid_position[2] = coord[2]
        cube = sample_cube
        origin = [grid_position[0]-2,grid_position[1]-2,grid_position[2]-1]
        travelled = [0,0,0]
        cube_potential, cube_var = volume_average(origin,cube,grid_pot,NGX,NGY,NGZ)
        potentials_list.append(cube_potential)

    ## PLOTTING
    n, bins, patches = plt.hist(potentials_list, num_bins, facecolor='#6400E1', alpha=0.5)
    plt.xlabel('Hartree potential (V)')
    plt.ylabel('% of centres')
    plt.savefig(img_file)

    ## SAVING
    df = pd.DataFrame.from_dict({'Potential':potentials_list},orient='index')
    df = df.transpose()
    df.to_csv(output_file)
    return potentials_list

#------------------------------------------------------------------------------

def plot_planar_average(lattice_vector,input_file='LOCPOT',output_file='PlanarAverage.csv',img_file='PlanarAverage.png'):
    '''
    Plot planar and macroscopic average of the electronic potential.

    This function reads the electronic potential from the specified VASP output file (LOCPOT) and calculates the planar
    and macroscopic averages of the potential along the lattice vector direction. It then generates two plots: one for
    the planar average and another for the macroscopic average. The function also saves the data to a CSV file.

    Parameters:
        lattice_vector (float): The repeating unit over which the potential is averaged to get the macroscopic average (in Angstroms).
        input_file (str, optional): The filename of the VASP output file containing the electronic potential (LOCPOT).
                                    Default is 'LOCPOT'.
        output_file (str, optional): Name of the output data file to store the planar and macroscopic average data.
                                     Default is 'PlanarAverage.csv'.
        img_file (str, optional): Name of the output image file for the plots. Default is 'PlanarAverage.png'.

    Returns:
        tuple: A tuple containing two arrays:
            - The planar average of the electronic potential.
            - The macroscopic average of the electronic potential.

    Example:
        >>> lattice_vector = 15.0
        >>> input_file = 'LOCPOT'
        >>> output_file = 'PlanarAverage.csv'
        >>> img_file = 'PlanarAverage.png'
        >>> planar_avg, macro_avg = plot_planar_average(lattice_vector, input_file, output_file, img_file)
        ... (plots generated and planar and macroscopic average data saved to 'PlanarAverage.csv')
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, planar_average, macroscopic_average

    # GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    ## PLANAR AVERAGE
    planar = planar_average(grid_pot,NGX,NGY,NGZ)

    ## MACROSCOPIC AVERAGE
    macro = macroscopic_average(planar,lattice_vector,resolution_z)

    ## PLOTTING
    plt.ylabel('V/V')
    plt.xlabel('Grid Position')
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig(img_file)

    ## SAVING
    df = pd.DataFrame.from_dict({'Planar':planar,'Macroscopic':macro},orient='index')
    df = df.transpose()
    df.to_csv(output_file)
    return planar, macro

#------------------------------------------------------------------------------

def plot_planar_cube(input_file,lattice_vector,output_file='PlanarCube.csv',img_file='PlanarCube.png'):
    '''
    Plot planar and macroscopic average of the electronic potential from a cube file.

    This function reads the electronic potential from the specified cube file and calculates the planar and macroscopic
    averages of the potential along the lattice vector direction. It then generates two plots: one for the planar average
    and another for the macroscopic average. The function also saves the data to a CSV file.

    Parameters:
        input_file (str): The filename of the cube file containing the electronic potential.
        lattice_vector (float): The repeating unit over which the potential is averaged to get the macroscopic average (in Angstroms).
        output_file (str, optional): Name of the output data file to store the planar and macroscopic average data.
                                     Default is 'PlanarCube.csv'.
        img_file (str, optional): Name of the output image file for the plots. Default is 'PlanarCube.png'.

    Returns:
        tuple: A tuple containing two arrays:
            - The planar average of the electronic potential.
            - The macroscopic average of the electronic potential.

    Example:
        >>> input_file = 'potential.cube'
        >>> lattice_vector = 15.0
        >>> output_file = 'PlanarCube.csv'
        >>> img_file = 'PlanarCube.png'
        >>> planar_avg, macro_avg = plot_planar_cube(input_file, lattice_vector, output_file, img_file)
        ... (plots generated and planar and macroscopic average data saved to 'PlanarCube.csv')
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    import ase.io.cube
    import pandas as pd
    from macrodensity.density_tools import planar_average, macroscopic_average

    # GETTING POTENTIAL
    potential, atoms = ase.io.cube.read_cube_data(input_file)
    vector_a = np.linalg.norm(atoms.cell[1])
    vector_b = np.linalg.norm(atoms.cell[1])
    vector_c = np.linalg.norm(atoms.cell[2])
    NGX = len(potential)
    NGY = len(potential[0])
    NGZ = len(potential[0][0])
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ

    ## PLANAR AVERAGE
    planar = planar_average(potential,NGX,NGY,NGZ)

    ## MACROSCOPIC AVERAGE
    macro  = macroscopic_average(planar,lattice_vector,resolution_z)

    ## PLOTTING
    plt.ylabel('V/V')
    plt.xlabel('Grid Position')
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig(img_file)

    ## SAVING
    df = pd.DataFrame.from_dict({'Planar':planar,'Macroscopic':macro},orient='index')
    df = df.transpose()
    df.to_csv(output_file)
    return planar, macro

#------------------------------------------------------------------------------

def moving_cube(cube=[1,1,1],vector=[1,1,1],origin=[0,0,0],magnitude = 280,input_file='LOCPOT',output_file='MovingCube.csv',img_file='MovingCube.png'):
    '''
    Calculate the travelling volume average of the electronic potential along a specific vector.

    This function calculates the volume average of the electronic potential as a function of the position along the specified
    vector. The volume average is performed by moving a cube of specified dimensions along the vector from the specified
    origin position. The magnitude parameter determines the distance covered in each direction from the origin. The resulting
    potential values at each position are plotted, and the data is saved to a CSV file.

    Parameters:
        cube (list, optional): The size of the cube used for volume averaging in units of mesh points (NGX/Y/Z).
                               Default is [1, 1, 1].
        vector (list, optional): The vector along which the cube moves for volume averaging. Default is [1, 1, 1].
        origin (list, optional): The starting position of the cube in fractional coordinates. Default is [0, 0, 0].
        magnitude (float, optional): The distance covered by the cube in each direction from the origin along the vector (in Angstroms).
                                     Default is 280.
        input_file (str, optional): The filename of the file containing the electronic potential (e.g., LOCPOT).
                                    Default is 'LOCPOT'.
        output_file (str, optional): Name of the output data file to store the volume-averaged potential data.
                                     Default is 'MovingCube.csv'.
        img_file (str, optional): Name of the output image file for the potential plot. Default is 'MovingCube.png'.

    Returns:
        list: A list containing the volume-averaged potential values at each position along the vector.

    Example:
        >>> cube_size = [2, 2, 2]
        >>> vector = [1, 0, 0]
        >>> origin = [0.5, 0.5, 0.5]
        >>> magnitude = 280
        >>> input_file = 'LOCPOT'
        >>> output_file = 'MovingCube.csv'
        >>> img_file = 'MovingCube.png'
        >>> potential_values = moving_cube(cube_size, vector, origin, magnitude, input_file, output_file, img_file)
        ... (plot generated and potential values saved to 'MovingCube.csv')
    '''
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, vector_2_abscissa, travelling_volume_average
    import matplotlib.pyplot as plt
    import pandas as pd

    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    cubes_potential = travelling_volume_average(grid_pot,cube,origin,vector,NGX,NGY,NGZ,magnitude)
    abscissa = vector_2_abscissa(vector,magnitude,resolution_x,resolution_y,resolution_z)

    ## PLOTTING
    plt.plot(abscissa, cubes_potential)
    plt.xlabel("$z (\AA)$")
    plt.ylabel("Potential (eV)")
    plt.savefig(img_file)

    ##SAVING
    df = pd.DataFrame.from_dict({'Potential':cubes_potential},orient='index')
    df = df.transpose()
    df.to_csv(output_file)
    return cubes_potential

#------------------------------------------------------------------------------

def spherical_average(cube_size,cube_origin,input_file='LOCPOT',print_output=True):
    '''
    Calculate the volume average of the electronic potential within a spherical region.

    This function calculates the volume average of the electronic potential within a spherical region
    defined by a specific size and origin. The size of the spherical region is specified by the cube_size
    parameter, which determines the number of mesh points along each direction (NGX/Y/Z). The origin of the
    sphere is given by the cube_origin parameter, specified in fractional coordinates. The function reads the
    electronic potential data from the specified input file (e.g., LOCPOT) and calculates the potential and variance
    within the spherical region.

    Parameters:
        cube_size (list): The size of the spherical region in units of mesh points (NGX/Y/Z).
        cube_origin (list): The origin of the spherical region in fractional coordinates.
        input_file (str, optional): The filename of the file containing the electronic potential (e.g., LOCPOT).
                                    Default is 'LOCPOT'.
        print_output (bool, optional): If True, the function prints the calculated potential and variance.
                                       Default is True.

    Returns:
        tuple: A tuple containing the volume-averaged potential and the variance within the spherical region.

    Example:
        >>> cube_size = [5, 5, 5]
        >>> cube_origin = [0.5, 0.5, 0.5]
        >>> input_file = 'LOCPOT'
        >>> potential, variance = spherical_average(cube_size, cube_origin, input_file)
        Potential            Variance
        --------------------------------
        1.23456              0.56789
    '''
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, volume_average

    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    cube = cube_size
    origin = cube_origin
    travelled = [0,0,0]
    cube_pot, cube_var = volume_average(origin=cube_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])

    ## PRINTING
    if print_output == True:
        print("Potential            Variance")
        print("--------------------------------")
        print(cube_pot,"   ", cube_var)
    return cube_pot, cube_var

#------------------------------------------------------------------------------

def plot_field_at_point(a_point,b_point,c_point,input_file='LOCPOT'):
    '''
    UNDER DEVELOPMENT. FULL OF BUGS (BEWARE)
    Plot the electric field magnitude and direction on a user-defined plane.

    This function plots the electric field magnitude and direction on a user-defined plane defined by three points: a_point, b_point, and c_point. The function reads the electronic potential data from the specified input file (e.g., LOCPOT) and calculates the electric field (gradient of the potential). The electric field is then visualized by plotting contours of the electric field magnitude and arrows indicating the direction of the electric field on the defined plane.

    Parameters:
        a_point (list): The fractional coordinates of the first point defining the plane.
        b_point (list): The fractional coordinates of the second point defining the plane.
        c_point (list): The fractional coordinates of the third point defining the plane.
        input_file (str, optional): The filename of the file containing the electronic potential (e.g., LOCPOT).
                                    Default is 'LOCPOT'.

    Returns:
        None: This function directly plots the electric field visualization.

    Example:
        >>> a_point = [0.1, 0.2, 0.3]
        >>> b_point = [0.4, 0.5, 0.6]
        >>> c_point = [0.7, 0.8, 0.9]
        >>> input_file = 'LOCPOT'
        >>> plot_field_at_point(a_point, b_point, c_point, input_file)
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors,cm #colour maps; so I can specify cube helix
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    print("Calculating gradients (Electic field, E=-Grad.V )...")
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = points_2_plane(a,b,c)
    ## Calculate magnitude of gradient.
    # Should be able to use numpy.linalg.norm for this, but the Python array indices are causing me grief
    X2 = np.multiply(grad_x,grad_x)
    Y2 = np.multiply(grad_y,grad_y)
    Z2 = np.multiply(grad_z,grad_z)
    grad_mag = np.sqrt(np.add(X2,Y2,Z2))
    # This was my, non working, attempt to use the built in function.
    #grad_mag=np.linalg.norm( [grad_y,grad_y,grad_z], axis=3)

    ## This function in Macrodensity averages Efield ACROSS Z for Slab calculations
    #xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag) #AVG over full volume

    # Here we construct the same xx,yy,grd variables with a SLICE, forming a plane in XY at particular ZSLICE
    xx, yy = np.mgrid[0:NGX,0:NGY]
    ZSLICE= NGZ/2 # Chosses where in the Z axis the XY slice is cut through
    # Slice of magnitude of electric field, for contour plotting
    grd=grad_mag[:,:,ZSLICE]
    # Slices of x and y components for arrow plotting
    grad_x_slice=grad_x[:,:,ZSLICE]
    grad_y_slice=grad_y[:,:,ZSLICE]
    # OK, that's all our data

    # This code tiles the data to (2,2) to re-expand unit cell to a 2x2 supercell in XY
    xx,yy=np.mgrid[0:2*NGX,0:2*NGY]
    grd=np.tile(grd, (2,2))
    grad_x_slice=np.tile(grad_x_slice, (2,2))
    grad_y_slice=np.tile(grad_y_slice, (2,2))
    # End of tiling code

    ## Contours of the above sliced data
    plt.contour(xx,yy,grd,6,cmap=cm.cubehelix)

    # Also generate a set of Efield arrows ('quiver') for this data.
    # Specifying the drawing parameters is quite frustrating - they are very brittle + poorly documented.
    plt.quiver(xx,yy, grad_x_slice, grad_y_slice,
            color='grey',
            units='dots', width=1, headwidth=3, headlength=4
            ) #,
    #        units='xy', scale=10., zorder=3, color='blue',
    #        width=0.007, headwidth=3., headlength=4.)

    plt.axis('equal') #force square aspect ratio; this assuming X and Y are equal.
    plt.show()

#------------------------------------------------------------------------------

def plot_plane_field(a_point,b_point,c_point,input_file='LOCPOT'):
    '''
    UNDER DEVELOPMENT, FULL OF BUGS (BEWARE)
    Plot the electric field on a user-defined plane and display it as a contour plot.

    Parameters:
        a_point (list): Fractional coordinates of the first point that defines the plane.
        b_point (list): Fractional coordinates of the second point that defines the plane.
        c_point (list): Fractional coordinates of the third point that defines the plane.
        input_file (str, optional): The filename of the VASP LOCPOT file containing the electrostatic potential. Default is 'LOCPOT'.

    Returns:
        None

    Note:
        - The function reads the electrostatic potential from the specified VASP LOCPOT file.
        - The plane is defined by three points: a_point, b_point, and c_point.
        - The electric field (gradient of the electrostatic potential) is computed using finite differences.
        - The function creates a contour plot of the electric field on the defined plane.
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = points_2_plane(a,b,c)
    ## Get the gradients
    XY = np.multiply(grad_x,grad_y)
    grad_mag = np.multiply(XY,grad_z)
    ## Create the plane
    #print(NGX,NGY,NGZ,plane_coeff,grad_mag)
    #xx,yy,grd = create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag)
    print(create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag)[10])
    ## Plot the surface
    plt.contour(xx,yy,grd,1)
    plt.show()

#------------------------------------------------------------------------------

def plot_active_plane(cube_size,cube_origin,tolerance=1E-4,input_file='LOCPOT'):
     '''
    UNDER DEVELOPMENT. FULL OF BUGS (BEWARE)
    Plot the active plane with contour and planar average of the electric field and potential.

    Parameters:
        cube_size (list): The size of the cube used for sampling the active plane.
        cube_origin (list): The origin point of the cube in fractional coordinates.
        tolerance (float, optional): The cutoff variance for identifying active and non-active cubes. Default is 1E-4.
        input_file (str, optional): The filename of the VASP LOCPOT file containing the electrostatic potential. Default is 'LOCPOT'.

    Returns:
        None

    Note:
        - The function reads the electrostatic potential from the specified VASP LOCPOT file.
        - The active plane is determined by sampling the potential in a cube around the cube_origin point.
        - The cutoff_varience parameter sets the threshold for distinguishing active and non-active cubes based on their variance in potential.
        - The function creates a contour plot of the electric field on the active plane.
        - It also plots the planar average of the electric field and potential throughout the material.
    '''
    import math
    import numpy as np
    import matplotlib.pyplot as plt
    from macrodensity.density_tools import read_vasp_density, matrix_2_abc, density_2_grid, numbers_2_grid, planar_average, macroscopic_average, volume_average
    from macrodensity.beta_tools import create_plotting_mesh,points_2_plane

    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    cutoff_varience = tolerance
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    ## Convert the fractional points to grid points on the density surface
    a = pot.numbers_2_grid(a_point,NGX,NGY,NGZ)
    b = pot.numbers_2_grid(b_point,NGX,NGY,NGZ)
    c = pot.numbers_2_grid(c_point,NGX,NGY,NGZ)
    plane_coeff = pot.points_2_plane(a,b,c)

    ## Get the gradients
    XY = np.multiply(grad_x,grad_y)
    grad_mag = np.multiply(XY,grad_z)

    ## Create the plane
    xx,yy,grd =  pot.create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
    ## Plot the surface
    plt.contourf(xx,yy,grd,V)
    plt.show()

    ##------------------------------------------------------------------
    ## Plotting a planar average (Field/potential) throughout the material
    ##------------------------------------------------------------------
    ## FIELDS
    planar = pot.planar_average(grad_x,NGX,NGY,NGZ)
    ## POTENTIAL
    planar = pot.planar_average(grid_pot,NGX,NGY,NGZ)
    ## MACROSCOPIC AVERAGE
    macro  = pot.macroscopic_average(planar,4.80,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    plt.show()

#------------------------------------------------------------------------------