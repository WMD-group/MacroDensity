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
- [ ] Update `tools.spherical_average` and `plotting.plot_on_site_potential` to handle different output files rather than just `LOCPOT`. 
- [ ] Update conda installation instructions
- [ ] Automatically publish version to pypi when a new release is made on GitHub
- [ ] Improve aesthetics of plots
- [ ] Automatically get lattice vector from POSCAR file, rather than taking it from user input
- [ ] update plot_onsite_potential function to fix histogram function
