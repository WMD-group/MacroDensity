# TODO's

- [x] Checking function docstrings (e.g. should describe what they do, describe inputs/outputs & include type's of input & output variables) 
- [x] Testing with GitHub actions
- [x] Formatting & cleanup of files
    * Remove unused imports (e.g. if package/function is imported but not used etc)
    * Maybe some reorganization and renaming of input files to make it more user-friendly (alpha_tools.py etc)?
- [x] Generate documentation (with `sphinx`)
    * Quickstart guide on using `sphinx` [here](https://www.sphinx-doc.org/en/master/usage/quickstart.html)
    * The docs should include 3 main pages:
        1. Intro to what the package can do (index.rst)
        2. Python API for main functions (main_functions.rst)
        3. A good tutorial (as Jupyter Notebook) exemplifying the main functions & workflow
            * The notebook should combine `Markdown` and `Code` cells to explain the workflow, showing the outputs/plots etc
- [x] Publish docs in `Readthedocs`
- [x] Installation
  - [x] Enable `pip` install: this requires adding package to [pypi](https://pypi.org/)
- [x] Update tutorial 1 to rewrite long code snippets as functions
- [x] Checking code & docstrings formatting with linters
    * Using the packages `isort`, `black`, `pycodestyle`, `flake8`
- [ ] Write up `JOSS` paper
  - Some examples here: 
    - [Example1](https://joss.theoj.org/papers/10.21105/joss.04817)
    - [Example2](https://joss.theoj.org/papers/10.21105/joss.03171)
- [x] Update `tools.spherical_average` and `plotting.plot_on_site_potential` to handle different output files rather than just `LOCPOT`. 
- [x] Format line length of docstrings to follow PEP conventions
- [x] Automatically publish version to pypi when a new release is made on GitHub
- [ ] Check which functions are not currently tested in `unit_tests.py` (for instance
 `plot_field_at_point`, `plot_plane_field`, `energy_band_alignment_diagram`, `get_band_extrema` are not currently tested). 
 Can use matplotlib to compare figures as done [here](https://github.com/SMTG-UCL/ShakeNBreak/blob/develop/tests/test_plotting.py#L1038). Add tests for these functions
- [ ] Add [codecov](https://docs.codecov.com/docs#step-5-get-coverage-analysis-from-codecov) badge to README to see how many lines of code are tested
- [ ] Refactor `plot_planar_average` to show coordinate position in x axis rather than grid position
- [ ] Add linting workflow (similar to [this](https://github.com/SMTG-UCL/ShakeNBreak/blob/develop/.github/workflows/lint.yml))
- [ ] Update conda installation instructions (need to update current conda package, https://anaconda.org/hcc/macrodensity):
    * Could look to add the package to the conda-forge channel following the instructions here: https://github.com/conda-forge/staged-recipes
    * If you want to update the package on the HCC channel, then the instructions are here though it appears there isn't a recipe for Macrodensity in this repo (probably due to age): https://github.com/unlhcc/hcc-conda-recipes
- [ ] Automatically get lattice vector from `POSCAR`` file, rather than taking it from user input
