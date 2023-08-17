""" 
macrodensity.plotting contains different types of plotting functions 
such as band alignment diagrams and potentials at different grid points.
"""

from __future__ import division, print_function

import ase
from ase.io import vasp, cube
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import pandas as pd
from scipy.interpolate import interp1d

from macrodensity.alpha import points_2_plane, create_plotting_mesh
from macrodensity.density import (
    density_2_grid,
    numbers_2_grid,
    volume_average,
    density_2_grid_gulp,
    macroscopic_average,
    planar_average,
)
from macrodensity.io import read_vasp_density, read_gulp_potential
from macrodensity.utils import matrix_2_abc


def energy_band_alignment_diagram(energies: list, materials:list, limit:float=8., width:float=1.,
                                  cols: list=['#74356C','#efce19'], textsize: int=22,
                                  arrowhead: float=0.7, outfile:str ='BandAlignment',
                                  references: list=[], edge=None) -> None:
   
    """
    Plot an energy band alignment diagram for a list of materials.

    Parameters:
        energies (list): A list of tuples containing the ionization potential (IP) and electron affinity (EA) of each material. The format is [(IP_1, EA_1), ...].

        materials (list): A list of material names corresponding to each set of energies.

        limit (float, optional): The limit for the energy axis (in eV). Default is 8.0.

        width (float, optional): The width of the bars representing IP and EA. Default is 1.0.

        cols (list, optional): A list of colors to use for the bars. Default is ['#74356C','#efce19'].

        textsize (int, optional): The font size for the text in the plot. Default is 22.
        
        arrowhead (float, optional): The size of the arrowhead for the energy arrows. Default is 0.7.

        outfile (str, optional): The base name for the output file (both .eps and .png files will be saved). Default is 'BandAlignment'.

        references (list, optional): A list of reference points (as tuples) to be shown as dashed lines on the plot. Each tuple should be in the format (reference_value, label). Default is an empty list.

        edge (None or str, optional): The edge color for the bars. If None, there will be no edge color. Default is None.

    Returns:
        None: The function generates and displays the energy band alignment diagram.

    Example:
        >>> energies = [(5.2, 2.8), (4.9, 3.1), (5.5, 2.6)]
        >>> materials = ['Material A', 'Material B', 'Material C']
        >>> energy_band_alignment_diagram(energies, materials, limit=8.0, width=0.8,
                                    cols=['#74356C', '#efce19'], textsize=18,
                                    arrowhead=0.5, outfile='BandAlignment',
                                    references=[(3.0, 'Reference 1'), (4.0, 'Reference 2')],
                                    edge='black')
    """
    fig, ax1 = plt.subplots(1, 1, sharex=True)
    fig.set_size_inches(len(energies) * 3, limit * 0.75)
    mpl.rcParams['xtick.labelsize'] = textsize
    mpl.rcParams['ytick.labelsize'] = textsize
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['ytick.major.width'] = 3
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['ytick.minor.size'] = 4
    mpl.rcParams['axes.linewidth'] = 3

    ax1.set_color_cycle(cols)
    ax2 = ax1.twinx()
    ind = np.arange(len(energies))

    ## Bars for the IP and background colour
    for i in ind:
        ax1.bar(i,-limit, width, edgecolor=None)
        ax1.bar(i,-energies[i][1], width, color='w', edgecolor=None)

    ## Reset the colours back to the start and plot the EA
    ax1.set_color_cycle(cols)
    for i in ind:
        ax1.bar(i,-energies[i][0], width, edgecolor=None,alpha=0.8)

    ## Set the limits of the axes
    ax1.set_ylim(-limit,0)
    ax2.set_ylim(-limit,0)
    ax1.set_xlim(-0.5,len(energies)-0.5)

    ## Set the names
    ax1.set_xticks(ind)
    ax1.set_xticklabels(materials,size=textsize)
    ran = [ str(k) for k in np.arange(0,limit+2,2)]
    ax1.set_yticklabels(ran[::-1],size=textsize)
    ran = [ '' for k in np.arange(0,limit+2,2)]
    ran[0] = 'Vacuum Level'
    ax2.set_yticklabels(ran[::-1],size=textsize)
    ax1.set_ylabel('Energy (eV)', size=textsize)


    os1 = 0.15   # Offset of the text 'IP' in the plot
    os2 = 0.2    # Offset of the text 'EA' in the plot

    for i, en in enumerate(energies):
        ax1.arrow(i-0.25,-en[0],0, en[0]-arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i-0.25,0,0, -en[1]+arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        ax1.arrow(i-0.25,0,0, -en[0]+arrowhead, width=0.005,
                  head_length=arrowhead,head_width=0.07, fc='black',ec='None')
        loc_ip = -(en[0] + en[1]) / 2
        ax1.text(i-os1,loc_ip,"IP  %3.1f"%en[1],fontsize=textsize)

        loc_ea = -en[0] / 2
        ax1.text(i-os2,loc_ea,"EA %3.1f"%en[0],fontsize=textsize)

        ax1.minorticks_on()
        # Don't show minor ticks on x-axis
        ax1.tick_params(axis='x',which='minor',bottom='off')
        ax2.minorticks_on()

    for ref in references:
        ax1.hlines(-ref[1], -0.5, len(energies) - 0.5,
                   linestyles='--', colors='r')
        ax1.text(len(energies) - 0.45, -ref[1] - 0.1, ref[0],
                 fontsize=textsize, color='r')

    fig.savefig('%s.eps'%outfile,bbox_inches='tight')
    fig.savefig('%s.png'%outfile,bbox_inches='tight')
    plt.show()
    print("Figure saved as %s.eps and %s.png"%(outfile, outfile))
    plt.close(fig)


def plot_active_space(cube_size: list,cube_origin: list,tolerance: float=1E-4,input_file='LOCPOT',print_output=True, plot_pot= False) -> tuple: 
    '''
    Plot the active space (vacuum and non-vacuum regions) based on potential variations.

    This function analyzes the potential variations within the specified cubes of the given size
    and determines whether each cube belongs to the vacuum or non-vacuum region based on the provided tolerance. 
    This function also plots the cube potentials of vacuum and non vacuum cubes.


    Parameters:
        cube_size (list of int): The size of the cubes in units of mesh points (NGX/Y/Z) for analysis.

        cube_origin (list of float): The starting point (origin) of the cubes in fractional coordinates (range [0, 1]).

        tolerance (float, optional): The cutoff variance value to distinguish vacuum from non-vacuum cubes. Default is 1E-4.

        input_file (str, optional): The file with VASP output for potential. Default is 'LOCPOT'.

        print_output (bool, optional): Whether to print the analysis results. Default is True.

        plot_pot (bool, optional): Whether to plot vacuum and non vacuum potential in a scatter (VOXEL) plot. Default is False.

    Returns:
        tuple: A tuple containing the number of cubes identified as vacuum and non-vacuum regions.

    Note:
        The function calculates the potential variation within each cube and compares it to the tolerance value.
        Cubes with potential variations below the tolerance are considered vacuum regions, while others are non-vacuum regions.

    Example:
        >>> cube_size = [2, 2, 2]
        >>> cube_origin = [0.0, 0.0, 0.0]
        >>> plot_active_space(cube_size, cube_origin, tolerance=1E-5)
    
    '''
    ## GETTING POTENTIAL
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    cutoff_variance = tolerance
    #grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    #travelled = [0,0,0]


    ## DISTNGUISHING VACUUM FROM NON_VACUUM
    vacuum = []
    vac_pot = []
    non_vacuum = []
    nvac_pot = []

    for i in range(0,NGX,cube_size[0]):
        for j in range(0,NGY,cube_size[1]):
            for k in range(0,NGZ,cube_size[2]):
                sub_origin = [float(i)/NGX,float(j)/NGY,float(k)/NGZ]
                cube_pot, cube_var = volume_average(origin=sub_origin,cube=cube_size,grid=grid_pot,nx=NGX,ny=NGY,nz=NGZ,travelled=[0,0,0])
                if cube_var <= cutoff_variance:
                    vacuum.append(sub_origin)
                    vac_pot.append(cube_pot)
                else:
                    non_vacuum.append(sub_origin)
                    nvac_pot.append(cube_pot)

    if print_output == True:
        print("Number of vacuum cubes: ", len(vacuum))
        print("Number of non-vacuum cubes: ", len(non_vacuum))
        print("Percentage of vacuum cubes: ",(float(len(vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
        print("Percentage of non-vacuum cubes: ",(float(len(non_vacuum))/(float(len(vacuum))+float(len(non_vacuum)))*100.))
    

    ## PLOTTING FUNCTIONS (CALYSTA ADDED)
    def plot_cube_potentials(coords, potentials, color_map='viridis', alpha=0.5):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        x, y, z = zip(*coords)

        norm_pots = (potentials - np.min(potentials))/(np.max(potentials) - np.min(potentials))

        cmap = plt.get_cmap(color_map)
        colors = cmap(norm_pots)

        ax.scatter(x, y, z, c= colors, alpha=alpha)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        plt.show()

        return fig
    
    if plot_pot == True:
        plot_cube_potentials(vacuum, vac_pot, color_map= 'viridis')
        plot_cube_potentials(non_vacuum, nvac_pot, color_map= 'viridis')
    else:
        pass 

    return len(vacuum), len(non_vacuum)


def plot_gulp_potential(lattice_vector: float,input_file: str='gulp.out',output_file: str='GulpPotential.csv',img_file: str='GulpPotential.png',new_resolution: int = 3000) -> tuple:
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
def plot_on_site_potential(species: str,sample_cube: list,potential_file: str='LOCPOT',coordinate_file: str='POSCAR',output_file: str='OnSitePotential.csv',img_file: str='OnSitePotential.png') -> list:
    '''
    Plot on-site electrostatic potential for a specific species.

    This function reads the electronic potential from the specified VASP output file (LOCPOT) and the atomic coordinates
    from the POSCAR file. It then calculates the on-site electrostatic potential for the specified species and generates
    a histogram to visualize its distribution across the sample cube.

    Parameters:
        species (str): The chemical symbol of the species whose on-site potential is of interest.

        sample_cube (list of int): The size of the sampling cube in units of mesh points (NGX/Y/Z).

        potential_file (str, optional): The filename of the VASP output file containing the electronic potential (LOCPOT). Default is 'LOCPOT'.

        coordinate_file (str, optional): The filename of the POSCAR file containing atomic coordinates. Default is 'POSCAR'.

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

    ## GETTING POTENTIALS
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(potential_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)
    grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    coords = vasp.read_vasp(coordinate_file)
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


def plot_planar_average(lattice_vector: float,input_file: str='LOCPOT',output_file: str='PlanarAverage.csv',img_file: str='PlanarAverage.png') -> tuple:
    '''
    Plot planar and macroscopic average of the electronic potential.

    This function reads the electronic potential from the specified VASP output file (LOCPOT) and calculates the planar
    and macroscopic averages of the potential along the lattice vector direction. It then generates two plots: one for
    the planar average and another for the macroscopic average. The function also saves the data to a CSV file.

    Parameters:
        lattice_vector (float): The repeating unit over which the potential is averaged to get the macroscopic average (in Angstroms).
        
        input_file (str, optional): The filename of the VASP output file containing the electronic potential (LOCPOT). Default is 'LOCPOT'.

        output_file (str, optional): Name of the output data file to store the planar and macroscopic average data. Default is 'PlanarAverage.csv'.

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
    '''

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


def plot_planar_cube(input_file: str,lattice_vector: float,output_file: str='PlanarCube.csv',img_file: str='PlanarCube.png') -> tuple:
    '''
    Plot planar and macroscopic average of the electronic potential from a cube file.

    This function reads the electronic potential from the specified cube file and calculates the planar and macroscopic
    averages of the potential along the lattice vector direction. It then generates two plots: one for the planar average
    and another for the macroscopic average. The function also saves the data to a CSV file.

    Parameters:
        input_file (str): The filename of the cube file containing the electronic potential.

        lattice_vector (float): The repeating unit over which the potential is averaged to get the macroscopic average (in Angstroms).

        output_file (str, optional): Name of the output data file to store the planar and macroscopic average data. Default is 'PlanarCube.csv'.

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
    '''
   

    # GETTING POTENTIAL
    potential, atoms = cube.read_cube_data(input_file)
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


def plot_field_at_point(a_point: list,b_point: list,c_point: list,input_file: str='LOCPOT', grad_calc: bool=False) -> None:
    '''
    UNDER DEVELOPMENT. FULL OF BUGS (BEWARE)
    Plot the electric field magnitude and direction on a user-defined plane.

    This function plots the electric field magnitude and direction on a user-defined plane defined by three points: a_point, b_point, and c_point. The function reads the electronic potential data from the specified input file (e.g., LOCPOT) and calculates the electric field (gradient of the potential). 
    The electric field is then visualized by plotting contours of the electric field magnitude and arrows indicating the direction of the electric field on the defined plane.

    Parameters:
        a_point (list): The fractional coordinates of the first point defining the plane.

        b_point (list): The fractional coordinates of the second point defining the plane.

        c_point (list): The fractional coordinates of the third point defining the plane.

        input_file (str, optional): The filename of the file containing the electronic potential (e.g., LOCPOT). Default is 'LOCPOT'.

        grad_calc (bool): if True , calculates the gradient of the field. Default is False due to computational expense

    Returns:
        None: This function directly plots the electric field visualization.

    Example:
        >>> a_point = [0.1, 0.2, 0.3]
        >>> b_point = [0.4, 0.5, 0.6]
        >>> c_point = [0.7, 0.8, 0.9]
        >>> input_file = 'LOCPOT'
        >>> plot_field_at_point(a_point, b_point, c_point, input_file)
    '''
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
    if grad_calc == True:
        print("Calculating gradients (Electic field, E=-Grad.V )...")
        grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
    else:
        pass 
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
    grad_mag = np.linalg.norm(np.add(X2,Y2,Z2), axis=-1)
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


def plot_plane_field(a_point: list,b_point: list,c_point: list,input_file: str='LOCPOT') -> None:
    '''
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
    xx,yy,grd = create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag)
    print(create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_mag)[10])
    ## Plot the surface
    plt.contour(xx,yy,grd,1)
    plt.show()


def plot_active_plane(cube_size: list,cube_origin: list,tolerance: float=1E-4,input_file: str='LOCPOT', grad_calc: bool= False) -> None:
    '''
    Plot the active plane with contour and planar average of the electric field and potential.

    Parameters:
        cube_size (list): The size of the cube used for sampling the active plane.

        cube_origin (list): The origin point of the cube in fractional coordinates.

        tolerance (float, optional): The cutoff variance for identifying active and non-active cubes. Default is 1E-4.

        input_file (str, optional): The filename of the VASP LOCPOT file containing the electrostatic potential. Default is 'LOCPOT'.

        grad_calc (bool): if True , calculates the gradient of the field. Default is False due to computational expense

    Returns:
        None

    Note:
        - The function reads the electrostatic potential from the specified VASP LOCPOT file.
        - The active plane is determined by sampling the potential in a cube around the cube_origin point.
        - The cutoff_varience parameter sets the threshold for distinguishing active and non-active cubes based on their variance in potential.
        - The function creates a contour plot of the electric field on the active plane.
        - It also plots the planar average of the electric field and potential throughout the material.
    '''
    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = read_vasp_density(input_file)
    vector_a,vector_b,vector_c,av,bv,cv = matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = density_2_grid(vasp_pot,NGX,NGY,NGZ)

    potential_variance = np.var(grid_pot)
    cutoff_variance = tolerance

    active_cube = potential_variance >= cutoff_variance #using tolerance in input parameter

    #Input section (define the plane with 3 points, fractional coordinates)
    a_point = [0, 0, 0]
    b_point = [1, 0, 1]
    c_point = [0, 1, 0]

    if active_cube:
        ## Get the gradiens (Field), if required.
        ## Comment out if not required, due to compuational expense.

        if grad_calc == True:
            print("Calculating gradients (Electic field, E=-Grad.V )...")
            grad_x,grad_y,grad_z = np.gradient(grid_pot[:,:,:],resolution_x,resolution_y,resolution_z)
        else:
            pass 
        

        ## Convert the fractional points to grid points on the density surface
        a = numbers_2_grid(a_point,NGX,NGY,NGZ)
        b = numbers_2_grid(b_point,NGX,NGY,NGZ)
        c = numbers_2_grid(c_point,NGX,NGY,NGZ)
        plane_coeff = points_2_plane(a,b,c)

        ## Get the gradients
        XY = np.multiply(grad_x,grad_y)
        grad_mag = np.multiply(XY,grad_z)

        ## Create the plane
        xx,yy,grd = create_plotting_mesh(NGX,NGY,NGZ,plane_coeff,grad_x)
        ## Plot the surface
        plt.contourf(xx,yy,grd) #This only plots the surface (no contour for the potentials)
        plt.show()
    else:
        print('The cube is not active (variance is below tolerance set)')


    ##------------------------------------------------------------------
    ## Plotting a planar average (Field/potential) throughout the material
    ##------------------------------------------------------------------
    ## FIELDS
    planar = planar_average(grad_x,NGX,NGY,NGZ)
    ## POTENTIAL
    planar = planar_average(grid_pot,NGX,NGY,NGZ)
    ## MACROSCOPIC AVERAGE
    macro  = macroscopic_average(planar,4.80,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    plt.show()
