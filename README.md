MacroDensity
====================
IMPORTANT NOTE::
Please use only the NewPotentialModule.py and NewInput.py files from now on. If you use older files you will probably get what you deserve. 



A python script to read in a VASP LOCPOT file and plots the electrostatic potential in a number of ways. It can plot planar average, macroscopic average or spherical average for any of the x, y or z axes. Can do the spherical averaging either along a line (S) or at specified points (Po).

Author
------------
Keith Butler


Requirements
------------
Python


Description
------------

The code can read in VASP CHGCAR/LOCPOT files and plot the potential/charge and the gradients of either. There are a number of ways in which these can be plotted and the file InputControl.py is arranged for the different modes.

In all cases if you require the gradient (field), uncomment the line beginning with "grad_x,grad_y,grad_z = ".

(i) Plotting in a plane.
    The plane is defined by three points.
    The plotting uses the contourf function.

(ii) Plotting a planar average throughout a sample.
    This uses the function planar_average, taking the quantity to be plotted & NGX/Y/Z as input.
    You may also wish to use macroscopic averaging to smear out the lattice oscillations, this  uses the macroscopic_average function, taking the quantity to be averaged, the period of smearing and the resolution_z(calculated earlier, do not worry about this) as input. For more info on macroscopic averaging see Jackson "Classical Electrodynamics".
    This plots with the plot function.

(iii) Getting the average inside a cube.
    This is useful for obtaining the potential inside a region of a pore, as described in J. Am. Chem. Soc., 2014, 136 (7), pp 2703–2706.
    You must define the cube size and the point where you want to place the origin of that cube, the fuction cube_potential takes these, plus the value to be sampled (grid_pot) and NGX/Y/Z. It returns the average value and the varience within the sample space. 

(iv) Plotting the average ithin a cube travelling along a vector.
    This performs the same task as (iii), but the cube is propogated along a vector with the unit vector defined by "vector" and a magnitude defined by "magnitude"

To-do
------------
- Add a full description of the input file format.

Execution
------------
python InputControl.py

Disclaimer
----------
This file is not affiliated with *VASP*. Feel free to use and modify, but do so at your own risk.
