Getting Started
===============

Requirements
------------

MacroDensity is intended to be relatively dependency free. However, some of the features
require other packages. 

* Python 3
* `Matplotlib <https://matplotlib.org/>`_ (only required for interactive plotting features).
* `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`_ (only required for reading Gaussian cube files and for atom centred averaging).
* `Pandas <https://pandas.pydata.org/>`_ (can speed up I/O of files significantly).

Installation
------------

MacroDensity can be installed directly using ``pip``.

.. code:: bash

    pip install git+git://github.com/WMD-group/MacroDensity.git

You should now be able to import MacroDensity to your own Python 
scripts.

.. code:: python

    import macrodensity as md

For examples of how to use MacroDensity for common tasks see our 
tutorials.
