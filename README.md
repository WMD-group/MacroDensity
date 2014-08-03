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
The code is then simply executed as:
```
python PlanarAverage.py
```
This results in a plot of the planar average and an output of the potential called planar.dat.


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
