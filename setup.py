
#!/usr/bin/env python

__author__ = "Keith T. Butler"
__copyright__ = "Copyright Keith T. Butler (2013)"
__version__ = "1.0"
__maintainer__ = "Keith T. Butler"
__email__ = "k.t.butler@bath.ac.uk"
__date__ = "Aug 24 2017"

from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='MacroDensity',
        version='1.0',
        description='Manipulation of electron density',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/WMD-group/MacroDensity',
        author='Keith T. Butler',
        author_email='k.t.butler@bath.ac.uk',
        license='GNU General Public License (GPL) v3',
        packages=['macrodensity'],
        zip_safe=False,
        install_requires=['scipy','numpy','spglib','ase'],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Scientific/Engineering']
    )
