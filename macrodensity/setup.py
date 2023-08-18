from setuptools import find_packages, setup

setup(
    name='macrodensity',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'spglib', 'ase', 'pandas', 'matplotlib']
    )

