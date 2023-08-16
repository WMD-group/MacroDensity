"""
MacroDensity is a package to read, process and plot electrostatic potential and electron density 
files from electronic structure calculations.
"""
import math

import numpy as np
from scipy import interpolate

from macrodensity.vasp import *
from macrodensity.density import *
from macrodensity.alpha import *
from macrodensity.plotting import *


