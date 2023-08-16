.. MacroDensity documentation master file, created by
   sphinx-quickstart on Tue Aug  1 15:20:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MacroDensity
=======================

Summary
=======
``MacroDensity`` is a set of Python scripts to read in a electrostatic potentials 
and electron densities from electronic structure calculations and plot in a number of ways, including:

1. Planar average
2. Spherical average
3. Atom centred average

Statement of Need
-----------------
When assessing the potential utility of novel semiconducting devices (p-n juntions, heterostructures, surface terminations) through simulation, 
an understanding of the variation in the electrostatic potential and electron density across the system is key. 
However, extraction and useful presentation of this data from the raw output of the simulation can prove cumbersome and often requires the use of visualisation software followed by manual data extraction. 
This can result in bottlenecks in high throughput screening projects, where the same data extraction procedure is repeatedly applied to large databases of candidate structures.

To address this, the Macrodensity package has been developed as a ``VASP``, ``FHI-AIMS`` and ``GULP`` post-processing tool. 
The package contains functions that enable the user to format the data from the ``VASP`` LOCPOT and CHGCAR files, the FHI-AIMS *.cube file, and GULP *.out file into physically meaningful quantities, which can then be plotted for user interpretation.

Requirements
------------
To use ``MacroDensity`` you will need the following Python packages:

- `Matplotlib <https://matplotlib.org/>`_
- `Pandas <https://pandas.pydata.org/>`_
- `Ase <https://wiki.fysik.dtu.dk/ase/>`_
- `Jupyter <https://jupyter.org/>`_

Installation 
================

User installation:
------------------
``MacroDensity`` can be installed using ``pip``:

.. code:: bash

  pip install macrodensity

.. Alternatively if needed, it can also be installed from ``conda`` with:

.. .. code:: bash

..   conda install macrodensity

Developer installation
----------------------

For development of ``MacroDensity``, you can install a copy of the package from the source directory:

1. Download ``MacroDensity`` source code using the command:

.. code:: bash

    git clone https://github.com/WMD-group/MacroDensity.git
    cd MacroDensity
    pip install -e .


Literature
----------
- General Approach: Butler, K. T., Hendon, C. H., & Walsh, A. `Electronic chemical potentials of porous Metal–Organic frameworks.`_ *Journal of the American Chemical Society*, 136(7), 2703–2706, **2014**
- Theoretical Background: Politzer, P., & Murray, J. S. `The fundamental nature and role of the electrostatic potential in atoms and molecules.`_ *Theoretical Chemistry Accounts*, 108(3), 134–142, **2002**

.. _Electronic chemical potentials of porous Metal–Organic frameworks.: https://doi.org/10.1021/ja4110073
.. _The fundamental nature and role of the electrostatic potential in atoms and molecules.: https://doi.org/10.1007/s00214-002-0360-0


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. toctree::
   :maxdepth: 3
   :hidden:
   .. :caption: Links

   modules
   tutorials
   studies