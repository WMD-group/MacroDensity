---
title: 'MacroDensity: Analysis of electrostatic potential and electron density landscapes of crystals'  
tags:
  - python
  - density functional theory
  - DFT
  - electrostatic potentials
  - electron density 
  - ab initio
  - vasp
authors:
  - name: Calysta A. Tesiman
    orcid: 0009-0008-7784-4320
    equal-contrib: false
    affiliation: "1" # (Multiple affiliations must be quoted)
authors:
  - name: Liam Harnett-Caulfield 
    equal-contrib: false
    affiliation: "1" # (Multiple affiliations must be quoted)
authors:
  - name: Irea Mosquera-Lois
    equal-contrib: false
    affiliation: "1" # (Multiple affiliations must be quoted)
authors:
  - name: Keith T. Butler
    equal-contrib: false
    corresponding: true
	affiliation: "2" # (Multiple affiliations must be quoted)
authors:
  - name: Aron Walsh 
    orcid: 0000-0001-5460-7033
    equal-contrib: false
    corresponding: true
	affiliation: "1" # (Multiple affiliations must be quoted)
	
affiliations:
 - name: Department of Materials, Imperial College London, London, United Kingdom
   index: 1
 - name: Department of Chemistry, University College London, London, United Kingdom
   index: 2
   
date: 31 August 2023
bibliography: JOSSpaper.bib

---

# Summary

We report a Python package to simplify the analysis of electrostatic potentials and electron density of crystals. Macrodensity can read volumetric output files from the first-principles materials modelling codes VASP (LOCPOT format) and FHI-AIMS (cube format), as well as the classical electrostatic potentials from GULP (standard output). The code consists of functions that calculate and plot planar macroscopic averages, and spherical average, as well as calculating conduction and valence band alignments over a user-defined vector or plane. As a result, this code has been used aid the data analysis and generation for a few publications. 

# Statement of need

When assessing the novel semiconducting materials and devices (p-n junctions, 
heterostructures, surface terminations) through simulation, an understanding of the 
variation in the electrostatic potential and electron density across the system is key [@Politzer:2002]. However, extraction and useful presentation of this data from the raw output of the simulation can prove cumbersome and often require the use of visualisation software followed 
by manual data extraction. This can result in bottlenecks in high throughput screening projects, 
where the same data extraction procedure is repeatedly applied to large databases of candidate structures.

The general approach of electrostatic potential and electron density data as well as its translation to a grid mesh is discussed in @Bulter:2014. This approach uses the Kohn-Sham density functional theory framework. This approach then samples the spherical averages over points within the system onto a matrix, where our raw data is generated. To process this data appropriately, ``MacroDensity`` was developed to simplify the data extraction and visualisation processes.  By defining planes or vectors along the landscape of electrostatic potentials and electronic density matrix, it becomes straightforward to produce meaningful analysis and visualisation plots across a user-defined area. 

# MacroDensity

``MacroDensity`` is a set of Python modules developed to read and analyse electrostatic potentials and electron density data from electronic structure calculations derived from Density Functional Theory (DFT) [@Kohn:1996]. The package allows users to read from VASP LOCPOT [@vasp], CHGCAR files, FHI-AIMS [@fhi_aims] (cube file), and GULP (standard output) files and format the data into physically meaningful quantities, which can then be plotted for user interpretation.

The package formats datasets containing information about a system's lattice parameters electron density, and electrostatic potentials. ``MacroDensity`` contains some high-level tools and functions to calculate and plot the planar and macroscopic average as defined in Jackson's Electrodynamics [@Jackson:2003] (Figure 1a). The determination of the lattice vector settings and how the macroscopic averaging is calculated in this package is best described from the work of @Peressi:1998. 

\begin{equation}
\label{eq:Planar-average}
  \overline{f}\left( z\right) =\dfrac{1}{S}\int _{S}f\left( x,y,z\right) dxdy
\end{equation}

``MacroDensity`` can also calculate and plot the localised potential around a certain atomic nucleus of a system. The approach to calculating this on-site (Hartree) potential is similar to calculating the Madelung potential (Figure 1b). This is useful for electron energy level predictions [@Walsh:2013]. In addition, the spherical average around a user-defined point within the system can be calculated using the package using an approach akin to the planar average function (integrating over a sphere instead of a plane). Calculations and averaging at different points in space can be used to quantify the valence band and conduction band positions relative to this average. This is a convenience function that is included within the package, which calculates the bulk interstitial alignment similarly to that from @Frensley:1976. 

\begin{equation}
\label{eq:Madelung-potential}
  \V_M = \sum_{i,j} \frac{(-1)^{n_i + n_j}}{4\pi\varepsilon_0 r_{ij}}
\end{equation}

``MacroDensity`` also calculates the average and macroscopic potentials within a specified volume of cube which moves along a plane of the system's lattice. The approach to this calculation is similar to the spherical average function (Figure 1c).

``MacroDensity`` has been used to rapidly generate data for the publications @Bulter:2014 and @Walsh:2013 amongst others. 

![Example analysis done with the package for AlAs, CsPbI<sub>3</sub>, and MgO: a) plots of the planar (blue) and macroscopic (orange) averages of the potential, b) plots of the mean potential along the [111] vector, c) onsite (Hartree) potentials of the constituent atoms of the compounds analysed. \label{fig1}](figure.png){ width=20% }

# Acknowledgements

We acknowledge input from Adam J. Jackson and Jarvist M. Frost in the early stages of the project. The work received financial support by Samsung Advanced Institute of Technology.     