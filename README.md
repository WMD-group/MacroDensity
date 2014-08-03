MacroDensity
====================

A python script to read in a VASP LOCPOT file and plots the electrostatic potential in a number of ways.
Plots:
Planar average.
Spherical average.
Atom centred averages (requires Atomistic Simulations Environment)

Author
------------
Keith Butler


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


Description
------------


To-do
------------
- Add a full description of the input file format.

Execution
------------
python InputControl.py

Disclaimer
----------
This file is not affiliated with *VASP*. Feel free to use and modify, but do so at your own risk.
