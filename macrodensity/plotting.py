"""macrodensity.plotting contains different configs of plotting functions such as band
alignment diagrams and potentials at different grid points."""

from __future__ import division, print_function

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.io import cube, vasp
from matplotlib import cm
from scipy.interpolate import interp1d

from macrodensity.averages import (
    macroscopic_average,
    planar_average,
    travelling_volume_average,
    volume_average,
)
from macrodensity.density import gradient_magnitude
from macrodensity.io import read_gulp_potential, read_vasp_density
from macrodensity.tools import _find_active_space, create_plotting_mesh
from macrodensity.utils import (
    density_2_grid,
    matrix_2_abc,
    numbers_2_grid,
    points_2_plane,
    vector_2_abscissa,
)

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
plt.style.use(f"{MODULE_DIR}/macrodensity.mplstyle")


def energy_band_alignment_diagram(
    energies: dict,
    ylims: list = (-8, 0),
    width: float = 1.0,
    cols: list = ["#74356C", "#efce19"],
    textsize: int = 24,
    arrowhead: float = 0.4,
    outfile: str = "BandAlignment",
    references: dict = {},
    edge=None,
    fig_format: str = "pdf",
) -> plt.figure:
    """Plot an energy band alignment diagram for a list of materials.

    Parameters:
        energies (dict): A dictionary mapping each material name to a tuple
         containing its ionization potential (IP) and electron affinity
         (EA) (e.g. {"CdS": (4.4, 7.7), ...}).

        ylims (tuple, optional): The limits for the energy/y axis (in eV).
            Default is (-8.0, 0.0)

        width (float, optional): The width of the bars representing IP and EA.
            Default is 1.0.

        cols (list, optional): A list of colors to use for the bars.
            Default is ['#74356C','#efce19'].

        textsize (int, optional): The font size for the text in the plot.
            Default is 22.

        arrowhead (float, optional): The size of the arrowhead for the energy arrows.
            Default is 0.7.

        outfile (str, optional): The base name for the output file to be saved.
            Default is 'BandAlignment'.

        references (dict, optional): A dictionary of reference points to be
            shown as dashed lines on the plot (e.g. {label: reference_value}.
            Default is an empty dictionary.

        edge (None or str, optional): The edge color for the bars.
            If None, there will be no edge color. Default is None.

        fig_format (str, optional): The format used to save the image.
            Default is "pdf".

    Returns:
        Figure: A matplotlib figure object containing the energy band alignment
            diagram.

    Example:
        >>> energies = [(5.2, 2.8), (4.9, 3.1), (5.5, 2.6)]
        >>> materials = ['Material A', 'Material B', 'Material C']
        >>> energy_band_alignment_diagram(energies, materials, ylims=(-9.0, 0.0), width=0.8,
                    cols=['#74356C', '#efce19'], textsize=18,
                    arrowhead=0.5, outfile='BandAlignment',
                    references={'Reference 1': 3.0, 'Reference 1', 'Reference 2': 4.0},
                    edge='black')"""
    energies_list = list(energies.values())
    materials = list(energies.keys())

    fig, ax1 = plt.subplots(1, 1, sharex=True)
    fig.set_size_inches(len(energies) * 3, abs(ylims[0]) * 0.75)
    mpl.rcParams["xtick.labelsize"] = textsize
    mpl.rcParams["ytick.labelsize"] = textsize
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams["ytick.major.width"] = 3
    mpl.rcParams["ytick.major.size"] = 7
    mpl.rcParams["ytick.minor.size"] = 4
    mpl.rcParams["axes.linewidth"] = 3
    ax2 = ax1.twinx()
    ind = np.arange(len(energies_list))
    ax1.set_prop_cycle(color=cols)

    # Bars for the IP and background colour
    for i in ind:
        ax1.bar(i, ylims[0], width, edgecolor=None)
        ax1.bar(i, -energies_list[i][1], width, color="w", edgecolor=None)

    # Reset the colours back to the start and plot the EA
    ax1.set_prop_cycle(color=cols)
    for i in ind:
        ax1.bar(i, -energies_list[i][0], width, edgecolor=None, alpha=0.8)

    # Set the limits of the axes
    ax1.set_ylim(ylims[0], ylims[1])
    ax2.set_ylim(ylims[0], ylims[1])
    ax1.set_xlim(-0.5, len(energies_list) - 0.5)

    # Set the names
    ax1.set_xticks(ind)
    ax1.set_xticklabels(materials, size=textsize)
    ran = [float(k) for k in np.arange(0, abs(ylims[0]) + 2, 2)]
    ax1.set_yticks([float(k) for k in np.arange(ylims[0], 2, 2)])
    ax1.set_yticklabels(ran[::-1], size=textsize)
    # Add tick for vacuum level
    ran = ["" for k in np.arange(0, abs(ylims[0]) + 2, 2)]
    ax2.set_yticks([float(k) for k in np.arange(ylims[0], 2, 2)])
    ran[0] = "Vacuum Level"
    ax2.set_yticklabels(ran[::-1], size=textsize)
    ax1.set_ylabel("Energy (eV)", size=textsize)

    # Offsets
    os1 = 0.15  # Offset of the text 'IP' in the plot
    os2 = 0.2  # Offset of the text 'EA' in the plot

    for i, en in enumerate(energies_list):
        ax1.arrow(
            i - 0.25,
            -en[0],
            0,
            en[0] - arrowhead,
            width=0.005,
            head_length=arrowhead,
            head_width=0.07,
            fc="black",
            ec="None",
        )
        ax1.arrow(
            i - 0.25,
            0,
            0,
            -en[1] + arrowhead,
            width=0.005,
            head_length=arrowhead,
            head_width=0.07,
            fc="black",
            ec="None",
        )
        ax1.arrow(
            i - 0.25,
            0,
            0,
            -en[0] + arrowhead,
            width=0.005,
            head_length=arrowhead,
            head_width=0.07,
            fc="black",
            ec="None",
        )
        loc_ip = -(en[0] + en[1]) / 2
        ax1.text(i - os1, loc_ip, "IP  %3.1f" % en[1], fontsize=textsize)

        loc_ea = -en[0] / 2
        ax1.text(i - os2, loc_ea, "EA %3.1f" % en[0], fontsize=textsize)

        # Don't show minor ticks on x-axis
        ax1.tick_params(axis="x", which="minor", bottom="off")
        ax2.tick_params(axis="x", which="minor", bottom="off")

    for label, value in references.items():
        ax1.hlines(
            -value, -0.5, len(energies_list) - 0.5, linestyles="--", colors="r"
        )
        ax1.text(
            len(energies_list) - 0.45,
            -value - 0.1,
            label,
            fontsize=textsize,
            color="r",
        )

    fig.savefig(f"{outfile}.{fig_format}", bbox_inches="tight")
    print(f"Figure saved as {outfile}.{fig_format}")
    plt.close(fig)
    return fig


def plot_active_space(
    cube_size: list,
    cube_origin: list,
    tolerance: float = 1e-4,
    input_file="LOCPOT",
    print_output=True,
) -> tuple:
    """Plot the active space (vacuum and non-vacuum regions) based on potential
    variations.

    This function analyzes the potential variations within the specified
    cubes of the given size and determines whether each cube belongs to
    the vacuum or non-vacuum region based on the provided tolerance.
    This function also plots the cube potentials of vacuum and
    non vacuum cubes.

    Parameters:
        cube_size (list of int): The size of the cubes in units of mesh
        points (NGX/Y/Z) for analysis.

        cube_origin (list of float): The starting point (origin) of the
        cubes in fractional coordinates (range [0, 1]).

        tolerance (float, optional): The cutoff variance value to
        distinguish vacuum from non-vacuum cubes. Default is 1E-4.

        input_file (str, optional): The file with VASP output for potential.
        Default is 'LOCPOT'.

        print_output (bool, optional): Whether to print the analysis results.
        Default is True.

    Returns:
        tuple: A tuple containing the number of cubes identified as
        vacuum and non-vacuum regions.

    Note:
        The function calculates the potential variation within each
        cube and compares it to the tolerance value. Cubes with potential
        variations below the tolerance are considered vacuum regions,
        while others are non-vacuum regions.

    Example:
        >>> cube_size = [2, 2, 2]
        >>> cube_origin = [0.0, 0.0, 0.0]
        >>> plot_active_space(cube_size, cube_origin, tolerance=1E-5)
    """

    def _plot_cube_potentials(
        coords_vac,
        potential_vac,
        coords_non_vac,
        potential_non_vac,
        color_map="viridis",
        alpha=0.5,
    ):
        fig = plt.figure(figsize=plt.figaspect(0.4))
        # First ax
        ax = fig.add_subplot(121, projection="3d")
        x, y, z = zip(*coords_vac)
        norm_pots = (potential_vac - np.min(potential_vac)) / (
            np.max(potential_vac) - np.min(potential_vac)
        )
        cmap = plt.get_cmap(color_map)
        colors = cmap(norm_pots)
        ax.scatter(x, y, z, c=colors, alpha=alpha)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Vacuum region")

        # set up the axes for the second plot
        ax2 = fig.add_subplot(1, 2, 2, projection="3d")
        x, y, z = zip(*coords_non_vac)
        norm_pots = (potential_non_vac - np.min(potential_non_vac)) / (
            np.max(potential_non_vac) - np.min(potential_non_vac)
        )
        cmap = plt.get_cmap(color_map)
        colors = cmap(norm_pots)
        ax2.scatter(x, y, z, c=colors, alpha=alpha)
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z")
        ax2.set_title("Non-vacuum region")

        norm = plt.Normalize(
            vmin=np.min(potential_vac), vmax=np.max(potential_vac)
        )
        fig.colorbar(
            cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=[ax, ax2],
            fraction=0.08,
            pad=0.08,
            shrink=0.5,
        )

        return fig

    dict_output = _find_active_space(
        input_file=input_file,
        cube_size=cube_size,
        cube_origin=cube_origin,
        tolerance=tolerance,
        print_output=True,
    )
    vacuum, non_vacuum, vac_pot, nvac_pot = (
        dict_output["Vacuum"],
        dict_output["Non-vacuum"],
        dict_output["Vacuum Potential"],
        dict_output["Non-vacuum Potential"],
    )

    fig = _plot_cube_potentials(
        coords_vac=vacuum,
        potential_vac=vac_pot,
        coords_non_vac=non_vacuum,
        potential_non_vac=nvac_pot,
        color_map="viridis",
    )

    return fig


def plot_on_site_potential(
    species: str,
    sample_cube: list,
    potential_file: str = "LOCPOT",
    coordinate_file: str = "POSCAR",
    output_file: str = "OnSitePotential.csv",
    img_file: str = "OnSitePotential.png",
) -> tuple:
    """Plot on-site electrostatic potential for a specific species.

    This function reads the electronic potential from the specified VASP
    output file (LOCPOT) and the atomic coordinates from the POSCAR file.
    It then calculates the on-site electrostatic potential for the
    specified species and generates a histogram to visualize
    its distribution across the sample cube.

    Parameters:
        species (str): The chemical symbol of the species whose on-site
        potential is of interest.

        sample_cube (list of int): The size of the sampling cube
        in units of mesh points (NGX/Y/Z).

        potential_file (str, optional): The filename of the VASP output file
        containing the electronic potential (LOCPOT). Default is 'LOCPOT'.

        coordinate_file (str, optional): The filename of the POSCAR file
        containing atomic coordinates. Default is 'POSCAR'.

        output_file (str, optional): Name of the output data file
        to store the on-site potential values.
        Default is 'OnSitePotential.csv'.

        img_file (str, optional): Name of the output image file
        for the histogram plot. Default is 'OnSitePotential.png'.

    Returns:
        tuple: A tuple containing the on-site potential values and the figure object.

    Example:
        >>> species = 'O'
        >>> sample_cube = [2, 2, 2]
        >>> potential_file = 'LOCPOT'
        >>> coordinate_file = 'POSCAR'
        >>> output_file = 'OnSitePotential.csv'
        >>> img_file = 'OnSitePotential.png'
        >>> on_site_potential = plot_on_site_potential(species, sample_cube, potential_file, coordinate_file, output_file, img_file)
        ... (plot generated and on-site potential data saved to 'OnSitePotential.png' and 'OnSitePotential.csv')
    """

    # Get potential
    if (
        "vasp" in potential_file
        or "LOCPOT" in potential_file
        or "CHGCAR" in potential_file
    ):
        vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(potential_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot, electrons = density_2_grid(
            vasp_pot, NGX, NGY, NGZ, config="VASP"
        )
    elif "cube" in potential_file:
        grid_pot, atoms = cube.read_cube_data(potential_file)
        vector_a = np.linalg.norm(atoms.cell[1])
        vector_b = np.linalg.norm(atoms.cell[1])
        vector_c = np.linalg.norm(atoms.cell[2])
        NGX = len(grid_pot)
        NGY = len(grid_pot[0])
        NGZ = len(grid_pot[0][0])
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
    elif "gulp" in potential_file or ".out" in potential_file:
        gulp_pot, NGX, NGY, NGZ, lattice = read_gulp_potential(potential_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot = density_2_grid(gulp_pot, NGX, NGY, NGZ, config="GULP")
    else:
        raise ValueError(f"File {potential_file} not recognised!")

    grad_x, grad_y, grad_z = np.gradient(
        grid_pot[:, :, :], resolution_x, resolution_y, resolution_z
    )
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
        origin = [
            grid_position[0] - 2,
            grid_position[1] - 2,
            grid_position[2] - 1,
        ]
        cube_potential, cube_var = volume_average(
            origin, cube, grid_pot, NGX, NGY, NGZ
        )
        potentials_list.append(cube_potential)

    ## PLOTTING
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(
        potentials_list, num_bins, facecolor="#6400E1", alpha=0.5
    )
    ax.set_xlabel("Hartree potential (V)")
    ax.set_ylabel("% of centres")
    plt.savefig(img_file)

    ## SAVING
    df = pd.DataFrame.from_dict({"Potential": potentials_list}, orient="index")
    df = df.transpose()
    df.to_csv(output_file)
    return df, fig


def plot_planar_average(
    lattice_vector: float,
    input_file: str = "LOCPOT",  # VASP potential file by default
    axis: str = "z",
    output_file: str = "planar_average.csv",
    img_file: str = "planar_average.png",
    new_resolution: int = 3000,
) -> tuple:
    """Calculate planar and macroscopic averages of potential data from different
    fileconfigs like gulp, cube, and vasp.

    Args:
        lattice_vector (float): The lattice vector value.
        input_file (str): Path to the input potential file.
        axis (str, Optional): Axis along which to calculate the
            average ('x', 'y', or 'z'). Default is "z".
        output_file (str): Path to save the output CSV file.
        img_file (str): Path to save the output image file.
        new_resolution (int): New resolution for interpolation.

    Returns:
        tuple: A tuple containing a dataframe with the planar average
        and macroscopic average and a figure object.
    """

    def _plot(planar, macro, img_file):
        fig, ax = plt.subplots()
        ax.set_ylabel("V (V)")
        ax.set_xlabel("Grid Position")
        ax.plot(planar, label="Planar")
        ax.plot(macro, label="Macroscopic")
        ax.legend(frameon=True)
        fig.savefig(img_file)
        return fig

    def _save_df(planar, macro, output_file, interpolated_potential=None):
        if interpolated_potential:
            df = pd.DataFrame.from_dict(
                {
                    "Planar": planar,
                    "Macroscopic": macro,
                    "Interpolated": interpolated_potential,
                },
                orient="index",
            )
        else:
            df = pd.DataFrame.from_dict(
                {"Planar": planar, "Macroscopic": macro}, orient="index"
            )
        df = df.transpose()
        df.to_csv(output_file)
        return df

    # Check axis is valid
    if axis not in ["x", "y", "z"]:
        raise ValueError(
            f'Axis {axis} not recognised! Must be "x", "y", or "z".'
        )

    if "cube" in input_file:
        potential, atoms = cube.read_cube_data(input_file)
        vector_a = np.linalg.norm(atoms.cell[1])
        vector_b = np.linalg.norm(atoms.cell[1])
        vector_c = np.linalg.norm(atoms.cell[2])
        NGX = len(potential)
        NGY = len(potential[0])
        NGZ = len(potential[0][0])
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ

        ## PLANAR AVERAGE
        planar = planar_average(potential, NGX, NGY, NGZ, axis=axis)
        ## MACROSCOPIC AVERAGE
        axis_to_resolution = {
            "x": resolution_x,
            "y": resolution_y,
            "z": resolution_z,
        }
        macro = macroscopic_average(
            planar, lattice_vector, axis_to_resolution[axis]
        )

        ## PLOTTING
        fig = _plot(planar, macro, img_file)
        ## SAVING
        df = _save_df(planar, macro, output_file)

    elif "gulp" in input_file or ".out" in input_file:
        pot, NGX, NGY, NGZ, lattice = read_gulp_potential(input_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        # TODO: Update config parameter in density_2_grid to be consistent with
        # code naming in other functions (e.g. if here we use GULP to refer to GULP,
        # should do the same in other functions)
        # Also use lower case for config variable following python conventions (eg config -> config)
        grid_pot = density_2_grid(pot, NGX, NGY, NGZ, config="GULP")

        ## POTENTIAL PLANAR AVERAGE
        planar = planar_average(grid_pot, NGX, NGY, NGZ, axis=axis)
        np.savetxt(output_file, planar)

        ## MACROSCOPIC AVERAGE
        axis_to_ng = {"x": NGX, "y": NGY, "z": NGZ}
        axis_to_vector = {"x": vector_a, "y": vector_b, "z": vector_c}
        new_abscissa = np.linspace(0, axis_to_ng[axis] - 1, new_resolution)
        f = interp1d(range(axis_to_ng[axis]), planar, kind="cubic")
        interpolated_potential = [f(i) for i in new_abscissa]
        macro = macroscopic_average(
            planar, lattice_vector, axis_to_vector[axis] / new_resolution
        )

        ## PLOTTING
        fig = _plot(planar, macro, img_file)
        ## SAVING
        df = _save_df(planar, macro, output_file, interpolated_potential)

    elif (
        "vasp" in input_file
        or "LOCPOT" in input_file
        or "CHGCAR" in input_file
    ):
        pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
        vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
        resolution_x = vector_a / NGX
        resolution_y = vector_b / NGY
        resolution_z = vector_c / NGZ
        grid_pot, electrons = density_2_grid(pot, NGX, NGY, NGZ, config="VASP")

        ## PLANAR AVERAGE
        planar = planar_average(grid_pot, NGX, NGY, NGZ, axis=axis)
        ## MACROSCOPIC AVERAGE
        axis_to_resolution = {
            "x": resolution_x,
            "y": resolution_y,
            "z": resolution_z,
        }
        macro = macroscopic_average(
            planar, lattice_vector, axis_to_resolution[axis]
        )

        ## PLOTTING
        fig = _plot(planar, macro, img_file)

        ## SAVING
        df = _save_df(planar, macro, output_file)

    else:
        raise ValueError(f"File {input_file} not recognised!")

    return df, fig


def plot_field_at_point(
    a_point: list,
    b_point: list,
    c_point: list,
    input_file: str = "LOCPOT",
    grad_calc: bool = False,
) -> plt.figure:
    """Plot the electric field magnitude and direction on a user-defined plane.

    This function plots the electric field magnitude and direction on a user-defined plane
    defined by three points: a_point, b_point, and c_point.
    The function reads the electronic potential data from the specified
    input file (e.g., LOCPOT) and calculates the electric field (gradient of the potential).
    The electric field is then visualized by plotting contours of the
    electric field magnitude and arrows indicating the direction of the
    electric field on the defined plane.

    Parameters:
        a_point (list): The fractional coordinates of the first point defining the plane.

        b_point (list): The fractional coordinates of the second point defining the plane.

        c_point (list): The fractional coordinates of the third point defining the plane.

        input_file (str, optional): The filename of the file containing the
        electronic potential (e.g., LOCPOT). Default is 'LOCPOT'.

        grad_calc (bool): if True , calculates the gradient of the field.
        Default is False due to computational expense

    Returns:
        Figure: A matplotlib figure object containing the electric field contours and arrows.

    Example:
        >>> a_point = [0.1, 0.2, 0.3]
        >>> b_point = [0.4, 0.5, 0.6]
        >>> c_point = [0.7, 0.8, 0.9]
        >>> input_file = 'LOCPOT'
        >>> plot_field_at_point(a_point, b_point, c_point, input_file)
    """
    # ------------------------------------------------------------------
    # Get the potential
    # ------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    resolution_x = vector_a / NGX
    resolution_y = vector_b / NGY
    resolution_z = vector_c / NGZ
    grid_pot, electrons = density_2_grid(
        vasp_pot, NGX, NGY, NGZ, config="VASP"
    )
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    if grad_calc == True:
        print("Calculating gradients (Electic field, E=-Grad.V )...")
        grad_x, grad_y, grad_z = np.gradient(
            grid_pot[:, :, :], resolution_x, resolution_y, resolution_z
        )
    else:
        pass
    # ------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point, NGX, NGY, NGZ)
    b = numbers_2_grid(b_point, NGX, NGY, NGZ)
    c = numbers_2_grid(c_point, NGX, NGY, NGZ)
    plane_coeff = points_2_plane(a, b, c)
    ## Calculate magnitude of gradient.
    # Should be able to use numpy.linalg.norm for this, but the Python array indices are causing me grief
    X2 = np.multiply(grad_x, grad_x)
    Y2 = np.multiply(grad_y, grad_y)
    Z2 = np.multiply(grad_z, grad_z)
    grad_mag = np.linalg.norm(np.add(X2, Y2, Z2), axis=-1)
    # This was my, non working, attempt to use the built in function.
    # grad_mag=np.linalg.norm( [grad_y,grad_y,grad_z], axis=3)

    ## This function in Macrodensity averages Efield ACROSS Z for Slab calculations
    xx, yy, grd = md.create_plotting_mesh(
        NGX, NGY, NGZ, plane_coeff, grad_mag
    )  # AVG over full volume

    # Here we construct the same xx,yy,grd variables with a SLICE, forming a plane in XY at particular ZSLICE
    xx, yy = np.mgrid[0:NGX, 0:NGY]
    ZSLICE = NGZ / 2  # Chosses where in the Z axis the XY slice is cut through
    # Slice of magnitude of electric field, for contour plotting
    grd = grad_mag[:, :, ZSLICE]
    # Slices of x and y components for arrow plotting
    grad_x_slice = grad_x[:, :, ZSLICE]
    grad_y_slice = grad_y[:, :, ZSLICE]
    # OK, that's all our data

    # This code tiles the data to (2,2) to re-expand unit cell to a 2x2 supercell in XY
    xx, yy = np.mgrid[0 : 2 * NGX, 0 : 2 * NGY]
    grd = np.tile(grd, (2, 2))
    grad_x_slice = np.tile(grad_x_slice, (2, 2))
    grad_y_slice = np.tile(grad_y_slice, (2, 2))
    # End of tiling code

    ## Contours of the above sliced data
    fig, ax = plt.subplots()
    ax.contour(xx, yy, grd, 6, cmap=cm.cubehelix)

    # Also generate a set of Efield arrows ('quiver') for this data.
    # Specifying the drawing parameters is quite frustrating - they are very brittle + poorly documented.
    ax.quiver(
        xx,
        yy,
        grad_x_slice,
        grad_y_slice,
        color="grey",
        units="dots",
        width=1,
        headwidth=3,
        headlength=4,
    )  # ,
    #        units='xy', scale=10., zorder=3, color='blue',
    #        width=0.007, headwidth=3., headlength=4.)

    plt.axis(
        "equal"
    )  # force square aspect ratio; this assuming X and Y are equal.
    plt.show()

    return fig


def plot_plane_field(
    a_point: list, b_point: list, c_point: list, input_file: str = "LOCPOT"
) -> plt.figure:
    """Plot the electric field on a user-defined plane and display it as a contour plot.

    Parameters:
        a_point (list): Fractional coordinates of the first point that defines the plane.

        b_point (list): Fractional coordinates of the second point that defines the plane.

        c_point (list): Fractional coordinates of the third point that defines the plane.

        input_file (str, optional): The filename of the VASP LOCPOT file containing the
        electrostatic potential. Default is 'LOCPOT'.

    Returns:
        Figure: A matplotlib figure object containing the electric field contours.

    Note:
        - The function reads the electrostatic potential from the specified VASP LOCPOT file.
        - The plane is defined by three points: a_point, b_point, and c_point.
        - The electric field (gradient of the electrostatic potential) is computed using finite differences.
        - The function creates a contour plot of the electric field on the defined plane.
    """

    # ------------------------------------------------------------------
    # Get the potential
    # ------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    resolution_x = vector_a / NGX
    resolution_y = vector_b / NGY
    resolution_z = vector_c / NGZ
    grid_pot, electrons = density_2_grid(
        vasp_pot, NGX, NGY, NGZ, config="VASP"
    )
    ## Get the gradiens (Field), if required.
    ## Comment out if not required, due to compuational expense.
    grad_x, grad_y, grad_z = np.gradient(
        grid_pot[:, :, :], resolution_x, resolution_y, resolution_z
    )
    # ------------------------------------------------------------------
    ## Get the equation for the plane
    ## This is the section for plotting on a user defined plane;
    ## uncomment commands if this is the option that you want.
    ##------------------------------------------------------------------
    ## Convert the fractional points to grid points on the density surface
    a = numbers_2_grid(a_point, NGX, NGY, NGZ)
    b = numbers_2_grid(b_point, NGX, NGY, NGZ)
    c = numbers_2_grid(c_point, NGX, NGY, NGZ)
    plane_coeff = points_2_plane(a, b, c)
    ## Get the gradients
    XY = np.multiply(grad_x, grad_y)
    grad_mag = np.multiply(XY, grad_z)
    ## Create the plane
    # print(NGX,NGY,NGZ,plane_coeff,grad_mag)
    xx, yy, grd = create_plotting_mesh(NGX, NGY, NGZ, plane_coeff, grad_mag)
    print(create_plotting_mesh(NGX, NGY, NGZ, plane_coeff, grad_mag)[10])
    ## Plot the surface
    fig, ax = plt.subplots()
    ax.contour(xx, yy, grd, 1)
    plt.show()

    return fig


def plot_active_plane(
    cube_size: list,
    cube_origin: list,
    tolerance: float = 1e-4,
    input_file: str = "LOCPOT",
    grad_calc: bool = False,
) -> plt.figure:
    """Plot the active plane with contour and planar average of the electric field and
    potential.

    Parameters:
        cube_size (list): The size of the cube used for sampling the active plane.

        cube_origin (list): The origin point of the cube in fractional coordinates.

        tolerance (float, optional): The cutoff variance for identifying active and
        non-active cubes. Default is 1E-4.

        input_file (str, optional): The filename of the VASP LOCPOT file containing
        the electrostatic potential. Default is 'LOCPOT'.

        grad_calc (bool): if True , calculates the gradient of the field.
        Default is False due to computational expense

    Returns:
        Figure: A matplotlib figure object containing the electric field contours and planar average.

    Note:
        - The function reads the electrostatic potential from the specified VASP LOCPOT
            file.

        - The active plane is determined by sampling the potential in a cube around the
            cube_origin point.

        - The cutoff_varience parameter sets the threshold for distinguishing active and
            non-active cubes based on their variance in potential.

        - The function creates a contour plot of the electric field on the active plane.

        - It also plots the planar average of the electric field and potential throughout the material.
    """
    # ------------------------------------------------------------------
    # Get the potential
    # ------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    resolution_x = vector_a / NGX
    resolution_y = vector_b / NGY
    resolution_z = vector_c / NGZ
    grid_pot, electrons = density_2_grid(
        vasp_pot, NGX, NGY, NGZ, config="VASP"
    )

    potential_variance = np.var(grid_pot)
    cutoff_variance = tolerance

    active_cube = (
        potential_variance >= cutoff_variance
    )  # using tolerance in input parameter

    # Input section (define the plane with 3 points, fractional coordinates)
    a_point = [0, 0, 0]
    b_point = [1, 0, 1]
    c_point = [0, 1, 0]

    if active_cube:
        ## Get the gradiens (Field), if required.
        ## Comment out if not required, due to compuational expense.

        if grad_calc == True:
            print("Calculating gradients (Electic field, E=-Grad.V )...")
            grad_x, grad_y, grad_z = np.gradient(
                grid_pot[:, :, :], resolution_x, resolution_y, resolution_z
            )
        else:
            pass

        ## Convert the fractional points to grid points on the density surface
        a = numbers_2_grid(a_point, NGX, NGY, NGZ)
        b = numbers_2_grid(b_point, NGX, NGY, NGZ)
        c = numbers_2_grid(c_point, NGX, NGY, NGZ)
        plane_coeff = points_2_plane(a, b, c)

        ## Get the gradients
        XY = np.multiply(grad_x, grad_y)
        grad_mag = np.multiply(XY, grad_z)

        ## Create the plane
        xx, yy, grd = create_plotting_mesh(NGX, NGY, NGZ, plane_coeff, grad_x)
        ## Plot the surface
        plt.contourf(
            xx, yy, grd
        )  # This only plots the surface (no contour for the potentials)
        plt.show()
    else:
        print("The cube is not active (variance is below tolerance set)")

    ##------------------------------------------------------------------
    ## Plotting a planar average (Field/potential) throughout the material
    ##------------------------------------------------------------------
    ## FIELDS
    planar = planar_average(grad_x, NGX, NGY, NGZ)
    ## POTENTIAL
    planar = planar_average(grid_pot, NGX, NGY, NGZ)
    ## MACROSCOPIC AVERAGE
    macro = macroscopic_average(planar, 4.80, resolution_z)

    fig, ax = plt.subplots()
    ax.plot(planar)
    ax.plot(macro)
    plt.savefig("Planar.eps")
    plt.show()

    return fig


def plot_variation_along_vector(
    vector: list = [1, 1, 1],
    cube_size: list = [1, 1, 1],
    origin_point: list = [0, 0, 0],
    vector_magnitude: int = 280,
    input_file: str = "LOCPOT",
    show_electric_field: bool = True,
    img_file: str = "potential_variation.png",
    output_file: str = "potential_variation.csv",
):
    """Plot the potential and field variation along a specified vector.

    This function calculates the volume average of the electronic potential as
    function of the position along the specified vector. The volume average is
    performed by moving a cube of specified dimensions along the vector from the
    specified origin position. The magnitude parameter determines the distance
    covered in each direction from the origin. The resulting potential values
    at each position are plotted, and the data is saved to a CSV file.

    Parameters:
        vector (:obj:`list`, optional): The vector along which the cube moves
            for volume averaging. Default is [1, 1, 1].

        cube_size (:obj:`list`, optional): The size of the cube used for
            volume averaging in units of mesh points (NGX/Y/Z). Default is [1, 1, 1].

        origin_point (:obj:`list`, optional): The starting position of the cube
            in fractional coordinates. Default is [0, 0, 0].

        magnitude (:obj:`float`, optional): The distance covered by the cube in
            each direction from the origin along the vector (in Angstroms).
            Default is 280.

        input_file (:obj:`str`, optional): The filename of the file containing
            the electronic potential (e.g., LOCPOT). Default is 'LOCPOT'.

        output_file (:obj:`str,` optional): Name of the output data file to store
            the volume-averaged potential data. Default is 'MovingCube.csv'.

        img_file (:obj:`str`, optional): Name of the output image file for the
            potential plot. Default is 'MovingCube.png'.

    Returns:
        :obj:`fig`: A figure showing the volume-averaged potential values
        at each position along the vector.
    """
    vasp_pot, NGX, NGY, NGZ, lattice = read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = matrix_2_abc(lattice)
    resolution_x = vector_a / NGX
    resolution_y = vector_b / NGY
    resolution_z = vector_c / NGZ
    grid_pot, electrons = density_2_grid(vasp_pot, NGX, NGY, NGZ)
    grad_x, grad_y, grad_z = np.gradient(
        grid_pot[:, :, :], resolution_x, resolution_y, resolution_z
    )

    # Average along vector
    cubes_potential = travelling_volume_average(
        grid=grid_pot,
        cube=cube_size,
        origin=origin_point,
        vector=vector,
        nx=NGX,
        ny=NGY,
        nz=NGZ,
        magnitude=vector_magnitude,
    )
    abscissa = vector_2_abscissa(
        vector,
        vector_magnitude,
        resolution_x,
        resolution_y,
        resolution_z,
    )

    # Plotting potential
    fig, ax = plt.subplots()
    ax.plot(abscissa, cubes_potential)
    ax.set_xlabel("$z (\AA)$")
    ax.set_ylabel("Potential (V)")
    # Plot field
    if show_electric_field:
        # Getting the field and plotting it
        grad_mag = gradient_magnitude(grad_x, grad_y, grad_z)
        cubes_field = travelling_volume_average(
            grid=grad_mag,
            cube=cube_size,
            origin=origin_point,
            vector=vector,
            nx=NGX,
            ny=NGY,
            nz=NGZ,
            magnitude=vector_magnitude,
        )
        ax.plot(abscissa, cubes_field)
        ax.set_ylabel("Magnitude")
        ax.legend(
            ["Potential (eV)", "Field Magnitude (eV/$\AA$))"], frameon=True
        )

    fig.savefig(img_file)

    # Save dataframe
    df = pd.DataFrame.from_dict({"Potential": cubes_potential}, orient="index")
    df = df.transpose()
    df.to_csv(output_file)

    return fig
