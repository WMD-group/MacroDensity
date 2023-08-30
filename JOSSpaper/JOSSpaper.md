---
title: 'MacroDensity: Understanding The Electrostatic Potential and Electron Density Landscapes within Systems of Quantum Mechanical Simulations'
tags:
  - python
  - density functional theory
  - DFT
  - electrostatic potentials
  - electron density 
  - ab initio
  - vasp
authors:
# @aron, could you check list of authors is ok?
  - name: Calysta A. Tesiman
    orcid: 0009-0008-7784-4320
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Keith T. Butler
    corresponding: False 
    affiliation: 1
  - name: Jarvist M. Frost
    corresponding: False 
    affiliation: 1
  - name: Aron Walsh
    corresponding: true 
    affiliation: 1
affiliations:
 - name: Thomas Young Centre & Department of Materials, Imperial College London, London SW7 2AZ, UK 
   index: 1
date: 29 August 2023
bibliography: JOSSpaper.bib

---

# Summary
 
When assessing the potential utility of novel semiconducting devices (p-n juntions, 
heterostructures, surface terminations) through simulation, an understanding of the 
variation in the electrostatic potential and electron density across the system is key [@Politzer:2002]. 
However, extraction and useful presentation of this data from the raw output of the 
simulation can prove cumbersome and often requires the use of visualisation software followed 
by manual data extraction. This can result in bottlenecks in high throughout screening projects, 
where the same data extraction procedure is repeatedly applied to large databases of candidate structures.


# Statement of need

The general approach of electrostatic potential and electron density data as well as its translation to a grid mesh 
has been discussed in @Bulter:2014. This approach uses the Kohn-Sham density functional theory framework. This 
approach then samples the spherical averages over points within the system onto a matrix, where our raw data is generated
To proceess this data appropriately, ``MacroDensity`` was developed to simplify the data extraction and visualisation processes. 
By defining planes or vectors along the landscape of electrostatic potentials and electronic density matrix,
it becomes trivial to produce meaningful analysis and plots for visualisation across a user defined area. 

# MacroDensity

``MacroDensity`` is a set of Python modules developed to read and analyse electrostatic potentials and electron 
density data from electronic structure calculations derived from Density Functional Theory (DFT) [@Kohn:1996]. The package 
allows users to read from VASP LOCPOT [@vasp], CHGCAR files, FHI-AIMS [@fhi_aims] , *.cube file, and GULP *.out 
files and format the data into physically meaningful quantities, which can then be plotted for user interpretation.

The package formats datasets containing information about a system's lattice parameters electron density, and 
electrostatic potentials. ``MacroDensity`` contains some high-level tools and functions to calculate and plot
the planar and macroscopic average (as defined in Jackson's Electrodynamics [@Jackson:2003]). The determination of 
the lattice vector settings and how the macroscopic averaging is calculated in this package is best described from the work of
@Peressi:1998. 

``MacroDensity`` can also calculate and plot the localised potential around a certain atomic nucleus of a system. The approach
to calculating this on site (Hartree) potential is similar to calculating the Madelung potential. this is useful for 
electron energy level predictions [@Walsh:2013]. 

In addition, the spherical average around a user defined point within the system can be calculated using the package. Calculations 
and averaging of this average at different points in space can be used to quantify the valence band and conduction band positions relative ot this average.
This is a convenience functions which is included within the package, which calculates the bulk interstitial alignment similarly to that
from @Frensley:1976. 

``MacroDensity`` also contains other functions including the Moving Average, which calculates the average and macroscopic potentials within
a specified volume of cube which moves along a plane of the system's lattice. The approach to this calculation is similar to the spherical average 
function.

``MacroDensity`` has been used to rapidly generate data for the publications @Bulter:2014 and @Walsh:2013 amongst others. 

<!-- # Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->

<!-- # Acknowledgements -->

# References