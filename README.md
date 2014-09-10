MacroDensity
====================

A python script to read in a VASP LOCPOT file and plots the electrostatic potential in a number of ways.
Plots:
Planar average.
Spherical average.
Atom centred averages (requires Atomistic Simulations Environment)

Author
------------
Keith Butler, Aron Walsh, Chris Hendon.


Requirements
------------
Python
ASE (for atom centred functionality)


Installation
------------

Download the file NewPotentialModule.py. Include this file in your pyhtonpath, or in the directory where you are performing the analysis.

There are several example files included with the package, they are documented below.

PlanarAverage.py
------------
This example is for plotting the planar average of a potential along a vector. In the example we plot along the z vector.
The only variables which need to be set are in the first three lines.
```
input_file = 'LOCPOT.slab'
lattice_vector = 4.75
output_file = 'planar.dat'
```

The variable lattice vector refers to the lattice vector of the bulk crystal structure in the direction of the plotting. It is used to get the macroscopic average, as defined in [Jackson's book](https://archive.org/details/ClassicalElectrodynamics).

The code is then simply executed as:
```
python PlanarAverage.py
```
This results in a plot of the planar average and an output of the potential called planar.dat.

SphericalAverage.py
------------

This example is for plotting the average potential inside a sphere of given radius. This is the method used in or J.A.C.S. paper http://pubs.acs.org/doi/abs/10.1021/ja4110073

A full tutorial of the methods applie in the paper is available here: http://people.bath.ac.uk/chh34/

The lines which need to be edited for this are:
```
input_file = 'LOCPOT.MiL'
cube_size = [2,2,2]    # This size is in units of mesh points
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates
cube_origin = [0,0,0]
```
To run the code simply type:
```
python SphericalAverage.py
```
This results in some outpus telling you the average potential in the volume and the variance of the potential. If the variance is too high it means that you are not sampling  a well converged plateau in the potential; typically values below 10e-4 are acceptable.

OnSitePotential.py
------------

This is for calculating the potentials at the sites of a certain atomic nucleus, for example the O nucleii in an oxide. This on site potential calculated this way is equivalent to a Madelung potential and can be useful fo predicting electron energy levels (see http://pubs.acs.org/doi/abs/10.1021/ar400115x for details).

The input lines to edit are :
```
potential_file = 'LOCPOT' # The file with VASP output for potential
coordinate_file = 'POSCAR' # The coordinates file NOTE NOTE This must be in vasp 4 format 
species = "O"  # The species whose on-site potential you are interested in 
sample_cube = [5,5,5] # The size of the sampling cube in units of mesh points (NGX/Y/Z)
```

The sample cube parameter determines the size of the sample area, the units are mesh points, the magnitude of the mesh point is calculated by dividing the appropriate lattice vector by the number of points (NGX/Y/Z in OUTCAR).

To run the code simply type:
```
python OnSitePotentail.py
```
The result is a histogram plot using Matplotlib. If you prefer data ouput simply edit the final lines of the script.

PlaneField.py
------------
This plots the countour lines of the iso-surface of an electric field in an arbitrary plane as defined in the preamble part of the file.
```
a_point = [0, 0, 0]
b_point = [1, 0, 1]
c_point = [0, 1, 0]

input_file = 'LOCPOT.slab'
```
The execution is simply:
```
python PlaneField.py
```
This creates a contour plot of the field lines.


To-do
------------
- * Add a full description of the input file format.
- * Create input modules for other electronic structure codes.


Disclaimer
----------
This file is not affiliated with *VASP*. Feel free to use and modify, but do so at your own risk.
