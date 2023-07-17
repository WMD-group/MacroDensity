MacroDensity
====================
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.884521.svg)

A set of Python scripts to read in a electrostatic potentials and electron densities from electronic structure calculations and plot in a number of ways, including:

* Planar average
* Spherical average
* Atom centred average

# Statement of Need

When assessing the potential utility of novel semiconducting devices (p-n juntions, heterostructures, surface terminations) through simulation, an understanding of the variation in the electrostatic potential and electron density across the system is key. However, extraction and useful presentation of this data from the raw output of the simulation can prove cumbersome and often requires the use of visualisation software followed by manual data extraction. This can result in bottlenecks in high throughput screening projects, where the same data extraction procedure is repeatedly applied to large databases of candidate structures.

To address this, the Macrodensity package has been developed as a VASP, FHI-AIMS and GULP post-processing tool. The package contains functions that enable the user to format the data from the VASP LOCPOT and CHGCAR files, the FHI-AIMS *.cube file, and GULP *.out file into physically meaningful quantities, which can then be plotted for user interpretation. The code has been used to rapidly generate data for these publications: [1](http://pubs.acs.org/doi/abs/10.1021/ja4110073),[2](https://aip.scitation.org/doi/10.1063/5.0044866), amongst others. 

# Requirements

* [Python](https://www.python.org)
* [Matplotlib](http://matplotlib.org) (to plot results on the fly)
* [ASE](https://wiki.fysik.dtu.dk/ase/) (for atom centred functionality)
* [Pandas](https://pandas.pydata.org/) (optional - for quicker reading speeds; requires pandas 1.2.0 or newer)
* [Jupyter](https://jupyter.org/) (optional - for `.ipynb` notebooks in the [tutorials](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/tutorials))

# Installation

```
pip install git+git://github.com/WMD-group/MacroDensity.git
```

If you have modified the source code, please run the unit tests with
  ``python setup.py test``.

# Usage

The high level tools of macrodensity can be used as built-in convenience functions and are listed below. For finer input-control using macrodensity's lower level functions, example scripts are given in the [examples](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/examples) folder and can be copied, modified then run via `python <script.py>`. Example input files can be copied from the [tests](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/tests) folder. 

Tutorials in the form of `.ipynb` notebooks are also given in the [tutorials](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/tutorials) folder.

### Planar Average

This is for plotting the planar average and macroscopic average (as defined in [Jackson's Electrodynamics](https://archive.org/details/ClassicalElectrodynamics)) of a potential along a vector. The vector is z by default, edit [examples/PlanarAverage.py](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/examples/PlanarAverage.py) to change the vector).


The function can be run via:

```
## Inputs:
## lattice_vector = Repeating unit over which the potential is averaged to get the macroscopic average (Angstroms)
## input_file = input filename to be read (DEFAULT = 'LOCPOT')
## output file = name of output data file (DEFAULT = 'PlanarAverage.csv')
## img_file = name of output image file (DEFAULT = 'PlanarAverage.png')

## Outputs:
## planar average and macroscopic average
## .csv data file containing: planar average and macroscopic average
## image file plotting .csv data

import macrodensity as md
md.plot_planar_average(lattice_vector=4.75, input_file='LOCPOT', output file='PlanarAverage.csv', img_file='PlanarAverage.png')
```

or by copying, modifying and running [PlanarAverage.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/PlanarAverage.py) separately.

For the best overview of what the lattice_vector setting should be, and how macroscopic averaging in general works, this paper from [Baldereschi](http://iopscience.iop.org/article/10.1088/0022-3727/31/11/002/meta) and the crew can't be beaten.

The output can be plotted as so:

![PA-Heterojunction](/tutorials/HJ-offset.png)

Further analysis of the band offset, deformation potential and volume change is outlined in [/tutorials/HeteroJunction/HeteroJunction.ipynb](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/tutorials/HeteroJunction/HeteroJunction.ipynb). For use in a slab-model style calculation for the Ionisation Potential, see [/tutorials/Slab/SlabCalculation.ipynb](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/tutorials/Slab/SlabCalculation.ipynb).

Additionally, the planar average can be obtained and interpolated from [GULP](http://gulp.curtin.edu.au/gulp/) outputs:

```
## Inputs:
## lattice_vector = 3.0
## input_file = Name of GULP input file (DEFAULT = 'gulp.out')
## output_file = Name of output data file (DEFAULT = 'GulpPotential.csv')
## img_file = Name of output image file (DEFAULT = 'GulpPotential.png')
## new_resolution = Total number of points for the interpolated planar avarage (DEFAULT = 3000)

## Outputs:
## planar average, macroscopic average, interpolated planar average
## .csv data file containing: planar average, macroscopic average, interpolated planar average
## image file plotting .csv data

import macrodensity as md
md.plot_gulp_potential(lattice_vector=3.0, input_file='gulp.out', new_resolution = 3000)
```

as well as `.cube` files:

```
## Inputs:
## input_file = input filename to be read (must be in .cube format)
## lattice_vector = Repeating unit over which the potential is averaged to get the macroscopic average (Angstroms)
## output file = name of output data file (DEFAULT = 'PlanarCube.csv')
## img_file = name of output image file (DEFAULT = 'PlanarCube.png')

## Outputs:
## planar average and macroscopic average
## .csv data file containing: planar average and macroscopic average
## image file plotting .csv data

import macrodensity as md
md.plot_planar_cube(input_file='cube_file.cube', lattice_vector=4.75, output_file='PlanarCube.csv', img_file='PlanarCube.png')
```

See [GulpPotential.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/GulpPotential.py) and [PlanarCube.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/PlanarCube.py)

------------

### Spherical Average

This is for plotting the average potential inside a sphere of given radius. The function can be run via:

```
## Inputs:
## cube_size = size of the cube in units of FFT mesh points (NGX/Y/Z)
## cube_origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
## input_file = VASP LOCPOT input filename to be read (DEFAULT = 'LOCPOT')
## print_output = Print terminal output (DEFAULT = True)

## Outputs:
## cube_potential, cube_variance (Terminal)

import macrodensity as md
md.spherical_average(cube_size=[2,2,2],cube_origin=[0.5,0.5,0.5],input_file='LOCPOT')
```

or by copying, modifying and running [SphericalAverage.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/SphericalAverage.py) separately.

This results in an output of the average potential in the volume, and the variance of the potential:

```
Reading header information...
Reading 3D data using Pandas...
Average of the potential =  -8.597839951107744e-14
Potential            Variance
--------------------------------
7.145660229     2.38371017456777e-05
```

If the variance is too high it means that you are not sampling a plateau in the electrostatic potential. Typically values below 10e-4 are acceptable, but you can also use [MovingCube.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/MovingCube.py) to verify this (See below). Further examples can be found in [/tutorials/Porous/Porous.ipynb](https://github.com/WMD-group/MacroDensity/tree/master/tutorials/Porous)

------------

### Moving Cube

This example takes the same approach as the spherical average above, but moves the sample volume
along a defined vector. This allows you to create a 1D profile of the travelling average of the
cube. 

The function can be run via:

```
## Inputs:
## cube = size of the cube in units of FFT mesh points (NGX/Y/Z)
## origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
## vector = vector across which the unit cell is traversed (hkl convention)
## magnitude = length travelled along the selected vector in units of FFT mesh points (NGX/Y/Z)
## input_file = VASP LOCPOT input filename to be read (DEFAULT = 'LOCPOT')
## output file = name of output data file (DEFAULT = 'MovingCube.csv')
## img_file = name of output image file (DEFAULT = 'MovingCube.png')

## Outputs:
## averaged electrostatic potential for the set cube size (list)
## .csv file containing the above data
## .png file presenting the above data

import macrodensity as md
md.moving_cube(cube=[1,1,1],vector=[1,1,1],origin=[0.17,0.17,0.17],magnitude=85,input_file='LOCPOT')
```

or by copying, modifying and running [MovingCube.py](https://github.com/WMD-group/MacroDensity/blob/V3.1.0/examples/MovingCube.py) separately.

The output for a ZnS (100x100x100 FFT) unit cell is given below. The electrostatic potential at the interstices (0.5 Å and 2.1 Å) can be extracted using [SphericalAverage.py](https://github.com/WMD-group/MacroDensity/blob/master/examples/SphericalAverage.py) or [bulk_interstital_alignment](https://github.com/WMD-group/MacroDensity/blob/1bc91530d2badbfbf1843fe614ca379ef31a122c/macrodensity/alpha_tools.py#L15).

![MovingCube](/tutorials/moving_cube.png)

------------

### Bulk Interstitial Alignment

A convenience function that combines multiple instances of SphericalAverage.py, averages them and prints the valence band and conduction band positions relative ot this average. This function is based on an analysis akin to that of [Frensley and Kroemer's](https://avs.scitation.org/doi/pdf/10.1116/1.568995) method of band alignment. An example of its use for the ZnS (zincnlende polymorph) interstices is given below:

The function can be run via:

```
## Inputs:
## intersices = Positions of the pores/interstices ([[interstice1],[interstice2],...])
## outcar = VASP OUTCAR input filename (DEFAULT = OUTCAR)
## locpot = VASP LOCPOT nput filename (DEFAULT = LOCPOT)
## cube_size = a cube defined by LOCPOT FFT mesh points (DEFAULT = [2,2,2])

## Outputs:
## Aligned Valence Band, Aligned Conduction Band, Interstitial variances

import macrodensity as md
md.bulk_interstitial_alignment(interstices=([0.5,0.5,0.5],[0.25,0.25,0.25]),outcar="OUTCAR",locpot="LOCPOT",cube_size=[2,2,2])
```

Data is presented in a similar format to [SphericalAverage.py](https://github.com/WMD-group/MacroDensity/blob/master/examples/SphericalAverage.py):

```
Reading header information...
Reading 3D data using Pandas...
Reading band edges from file: OUTCAR
Reading potential from file: LOCPOT
Interstital variances: [2.38371017456777e-05, 2.4181806506113166e-05]
VB_aligned      CB_aligned
--------------------------------
-4.95           -2.95
```
Further data gathered on zincblende type structures using this function can be found [here](https://aip.scitation.org/doi/10.1063/5.0044866)

------------

### On Site Potential

This example is for calculating the potential at the sites of a certain atomic nucleus, for example the O nuclei in an oxide. This on-site potential calculated this way is analagous to a Madelung potential and can be useful for predicting electron energy levels (see [this publication](http://pubs.acs.org/doi/abs/10.1021/ar400115x) for details).

The function can be run via:

```
## Inputs:
## potential_file = The file with VASP output for potential (DEFAULT = 'LOCPOT')
## coordinate_file = The coordinates file (DEFAULT = 'POSCAR')
## species = The species whose on-site potential you are interested in (string)
## sample_cube = The size of the sampling cube in units of mesh points (NGX/Y/Z)
## output file = name of output data file (DEFAULT = 'OnSitePotential.csv')
## img_file = name of output image file (DEFAULT = 'OnSitePotential.png')

## Outputs:
## .png histogram output
## .csv data output

import macrodensity as md
md.plot_on_site_potential(species='O',sample_cube=[5,5,5],potential_file='LOCPOT',coordinate_file='POSCAR')
```

The cube parameter determines the size of the sample area, the units are mesh points, the magnitude of the mesh point is calculated by dividing the appropriate lattice vector by the number of points (NGX/Y/Z in [OUTCAR](https://www.vasp.at/wiki/index.php/OUTCAR)).

------------

### Active Space

This checks for plateaus in the electrostatic potential within a volume for a set `tolerance` and prints the result. `cube_size` defines the size of the cube in units of mesh points. `cube_origin` defines the bottom left point of the cube the "0,0,0" point in fractional coordinates.

The function can be run via:

```
## Inputs:
## cube_size = size of the cube in units of FFT mesh points (NGX/Y/Z)
## cube_origin = real-space positioning of the bottom left point of the sampling cube (fractional coordinates of the unit cell).
## tolerance = threshold below which the electrostatic potential is considered to be plateaued (DEFAULT = 1E-4).
## input_file = VASP LOCPOT input filename to be read (DEFAULT = 'LOCPOT')
## print_output = Print terminal output (DEFAULT = True)

## Outputs:
## Percentage of vaccum vs non-vacuum cubes

import macrodensity as md
md.plot_active_space(cube_size=[2,2,2],cube_origin=[0.5,0.5,0.5],tolerance=1E-4,input_file='LOCPOT')
```

Which will give a terminal output that looks like this:

```
Reading header information...
Reading 3D data using Pandas...
Average of the potential =  -4.0478731477833207e-13
Number of vacuum cubes:  17
Number of non-vacuum cubes:  4079
Percentage of vacuum cubes:  0.4150390625
Percentage of non-vacuum cubes:  99.5849609375
```

------------

### Plotting

Macrodensity also allows the facile generation of band alignment diagrams across a range of structures. A more aesthetic alternative is [bapt](https://github.com/utf/bapt).

![BandAlignment](/tutorials/Plotting/BandAlignment05.png)

See [tutorials/Plotting/Plotting.ipynb](https://github.com/WMD-group/MacroDensity/tree/master/tutorials/Plotting) for further instructions.
