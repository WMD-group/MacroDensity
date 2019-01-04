.. _tutorials:

Tutorials
=========

Planar averaging of a potential
-------------------------------

This example is for plotting the planar average of a potential along a vector (here it is z). The only variables which need to be set are in the first three lines. 
Note ``LOCPOT`` file is just a regular ``VASP`` grid file.

Here ``lattice_vector`` is the length of a repeat unit in the direction
of the averaging, it is used to calculate the macroscopic average.

.. code:: python

    import macrodensity as md
    import math
    import numpy as np
    import matplotlib.pyplot as plt

    input_file = 'LOCPOT'
    lattice_vector = 4.75
    planar_file = 'planar.dat'
    macroscopic_file = 'planar.dat'

Nex we read in the potential and convert it to a grid.

.. code:: python

    vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density(input_file)
    vector_a, vector_b, vector_c, av, bv, cv = md.matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)

Finally we get planar and macroscopic averages and save and plot these.

.. code:: python

    planar = md.planar_average(grid_pot, NGX, NGY, NGZ)
    macro  = md.macroscopic_average(planar,lattice_vector,resolution_z)
    plt.plot(planar)
    plt.plot(macro)
    plt.savefig('Planar.eps')
    plt.show()
    np.savetxt(planar_file, planar)
    np.savetxt(macroscopic_file, macro)








