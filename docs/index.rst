.. MacroDensity documentation master file, created by
   sphinx-quickstart on Tue Aug  1 15:20:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MacroDensity
=======================

``MacroDensity`` is a Python package to post-process electrostatic potential and
electron density files from electronic structure calculations and plot them a number of ways, including
planar, sperical and atom centred averages.

Statement of Need
-----------------
To assess the potential utility of novel semiconducting devices (like p-n junctions, heterostructures,
surface terminations), it is key to understand how the electrostatic potential and electron density
change across the system. However, extraction and post-proccessing of this data from the raw output
of the simulation can prove cumbersome and often requires the use of visualisation software followed
by manual data extraction. This can result in bottlenecks in high throughput screening projects,
where the same data extraction procedure is repeatedly applied to large databases of candidate structures.

To address this, the Macrodensity package has been developed to post-process the output of ab-initio codes like ``VASP``, ``FHI-AIMS`` and ``GULP``.
The package contains functions to format the data from the ``VASP`` ``LOCPOT`` and ``CHGCAR`` files, the ``FHI-AIMS`` ``*.cube`` file,
and ``GULP`` ``*.out`` file into physically meaningful quantities, which can then be plotted for user interpretation.


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

.. code:: bash

    git clone https://github.com/WMD-group/MacroDensity.git
    cd MacroDensity
    pip install -e .


Literature
========================
For more information on the theory behind the package, please see the following references:

- General Approach: Butler, K. T., Hendon, C. H., & Walsh, A. `Electronic chemical potentials of porous Metal–Organic frameworks. <https://doi.org/10.1021/ja4110073>`_ *Journal of the American Chemical Society*, 136(7), 2703–2706, **2014**
- Theoretical Background:
   * Politzer, P., & Murray, J. S. `The fundamental nature and role of the electrostatic potential in atoms and molecules. <https://link.springer.com/article/10.1007/s00214-002-0363-9>`_ *Theoretical Chemistry Accounts*, 108(3), 134–142, **2002**
   * Peressi, M., Binggeli, N. & Baldereschi, A. `Band engineering at interfaces: theory and numerical experiments. <https://iopscience.iop.org/article/10.1088/0022-3727/31/11/002/meta>`_ *Journal of Physics D: Applied Physics*,31(11), 1273, **1998**
- Application to MOFs:
   * Butler, K. T., Hendon, C. H. & Aron Walsh, A. `Electronic Chemical Potentials of Porous Metal–Organic Frameworks. <https://doi.org/10.1021/ja4110073>`_ *Journal of the American Chemical Society*, 136(7), 2703–2706, **2014**


Contributing
========================

Bugs reports, feature requests and questions
----------------------------------------------

Please use the `Issue Tracker <https://github.com/WMD-group/MacroDensity/issues>`_
to report bugs or request new features.

Contributions to extend this package are very welcome! Please use the
`"Fork and Pull" <https://docs.github.com/en/get-started/quickstart/contributing-to-projects>`_
workflow to do so and follow the `PEP8 <https://peps.python.org/pep-0008/>`_ style guidelines.

Tests
----------------------------------------------

Unit tests are in the ``tests`` directory and can be run from the top directory using
`unittest <https://docs.python.org/3/library/unittest.html>`_.
Automatic testing is run on the master and develop branches using Github Actions. Please
run tests and add new tests for any new features whenever submitting pull requests.


License and citation
=======================
``MacroDensity`` is made available under the MIT License.

If you use it in your research, please cite:

* Method: Harnett-Caulfield, L., & Walsh, A. `Assessment of interstitial potentials for rapid prediction of absolute band energies in crystals.`_ *Journal of Chemical Physics*, 155(2). **2021**
* Code: JOSS paper

.. _Assessment of interstitial potentials for rapid prediction of absolute band energies in crystals.: https://doi.org/10.1063/5.0056141


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. toctree::
   :hidden:
   :maxdepth: 4

   installation
   Python API <modules>
   tutorials
   studies
