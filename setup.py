
#!/usr/bin/env python

__author__ = "Keith T. Butler"
__copyright__ = "Copyright Keith T. Butler (2013)"
__version__ = "1.0"
__maintainer__ = "Keith T. Butler"
__email__ = "keith.butler@stfc.ac.uk"
__date__ = "Jan 11 2019"

from setuptools import setup
import os
import unittest

module_dir = os.path.dirname(os.path.abspath(__file__))

def unit_tests():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='unit_tests.py')
    return test_suite

if __name__ == "__main__":
    setup(
        name='MacroDensity',
        version='1.0',
        description='Manipulation of electron density',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/WMD-group/MacroDensity',
        author='Keith T. Butler',
        author_email='keith.butler@stfc.ac.uk',
        license='GNU General Public License (GPL) v3',
        packages=['macrodensity'],
        zip_safe=False,
        install_requires=['scipy', 'numpy', 'spglib', 'ase'],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Scientific/Engineering'],
        test_suite='setup.unit_tests'
    )
