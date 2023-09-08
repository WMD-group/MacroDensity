#!/usr/bin/env python

__author__ = "Keith T. Butler"
__copyright__ = "Copyright Keith T. Butler (2013)"
__version__ = "3.1.0"
__maintainer__ = "Keith T. Butler"
__email__ = "k.t.butler@ucl.ac.uk"
__date__ = "17 July 2023"

import os
import unittest

from setuptools import setup

module_dir = os.path.dirname(os.path.abspath(__file__))

def unit_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='unit_tests.py')
    return test_suite


if __name__ == "__main__":
    setup(
        name='MacroDensity',
        version='3.1.0',
        description='Manipulation of electron density',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/WMD-group/MacroDensity',
        author='Keith T. Butler',
        author_email='k.t.butler@ucl.ac.uk',
        license='MIT License',
        packages=['macrodensity'],
        zip_safe=False,
        install_requires=['scipy', 'numpy', 'spglib', 'ase', "pandas"],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Scientific/Engineering'],
        test_suite='setup.unit_tests'
    )
