MacroDensity
====================
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.884521.svg)

A set of Python scripts to read in a electrostatic potentials and electron densities from electronic structure calculations and plot in a number of ways, including:

* Planar average
* Spherical average
* Atom centred average

# Statement of Need

To assess the potential utility of novel semiconducting devices (like p-n junctions, heterostructures, surface terminations), it is key to understand how the electrostatic potential and electron density change across the system. However, extraction and useful presentation of this data from the raw output of the simulation can prove cumbersome and often requires the use of visualisation software followed by manual data extraction. This can result in bottlenecks in high throughput screening projects, where the same data extraction procedure is repeatedly applied to large databases of candidate structures.

To address this, the Macrodensity package has been developed to post-process the output of ab-initio codes like `VASP`, `FHI-AIMS` and `GULP`. The package contains functions that enable the user to format the data from the `VASP LOCPOT` and `CHGCAR` files, the `FHI-AIMS *.cube` file, and `GULP *.out` file into physically meaningful quantities, which can then be plotted for user interpretation.

# Requirements

* [Python](https://www.python.org)
* [Matplotlib](http://matplotlib.org) (to plot results on the fly)
* [ASE](https://wiki.fysik.dtu.dk/ase/) (for atom centred functionality)
* [Pandas](https://pandas.pydata.org/) (optional - for quicker reading speeds; requires pandas 1.2.0 or newer)
* [Jupyter](https://jupyter.org/) (optional - for `.ipynb` notebooks in the [tutorials](https://github.com/WMD-group/MacroDensity/tree/V3.1.0/tutorials))

# Installation

## User installation

`Macrodensity` can be installed using `pip`:

```
pip install macrodensity
```

## Developer installation
For development of MacroDensity, you can install a copy of the package from the source directory:

```
git clone https://github.com/WMD-group/MacroDensity.git
cd MacroDensity
pip install -e .
```


Literature
========================
For more information on the theory behind the package, please see the following references:

- General Approach: Butler, K. T., Hendon, C. H., & Walsh, A. [Electronic chemical potentials of porous Metal–Organic frameworks](https://doi.org/10.1021/ja4110073). *Journal of the American Chemical Society*, 136(7), 2703–2706, 2014
- Theoretical Background: 
   * Politzer, P., & Murray, J. S. [The fundamental nature and role of the electrostatic potential in atoms and molecules](https://link.springer.com/article/10.1007/s00214-002-0363-9). *Theoretical Chemistry Accounts*, 108(3), 134–142, 2002
   * Peressi, M., Binggeli, N. & Baldereschi, A. [Band engineering at interfaces: theory and numerical experiments](https://iopscience.iop.org/article/10.1088/0022-3727/31/11/002/meta). *Journal of Physics D: Applied Physics*,31(11), 1273, 1998
- Application to MOFs:
   * Butler, K. T., Hendon, C. H. & Aron Walsh, A. [Electronic Chemical Potentials of Porous Metal–Organic Frameworks](https://doi.org/10.1021/ja4110073). *Journal of the American Chemical Society*, 136(7), 2703–2706, 2014
  


# Contributing

## Bugs reports, feature requests and questions

Please use the [Issue Tracker](https://github.com/WMD-group/MacroDensity/issues)
to report bugs or request new features.

Contributions to extend this package are very welcome! Please use the
[Fork and Pull](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)
workflow to do so and follow the [PEP8](https://peps.python.org/pep-0008/) style guidelines.

## Tests

Unit tests are in the ``tests`` directory and can be run from the top directory using
[unittest](https://docs.python.org/3/library/unittest.html).
Automatic testing is run on the master and develop branches using Github Actions. Please
run tests and add new tests for any new features whenever submitting pull requests.


# License and citation

MacroDensity is made available under the MIT License.

If you use it in your research, please cite:

* Method: Harnett-Caulfield, L., & Walsh, A.  [*Journal of Chemical Physics*](https://doi.org/10.1063/5.0044866) 2021

