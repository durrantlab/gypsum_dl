#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="gypsum_dl",
    version="1.2.1",
    description="Gypsum-DL is a free, open-source program that converts 1D and 2D small-molecule representations (SMILES strings or flat SDF files) into 3D models. It outputs models with alternate ionization, tautomeric, chiral, cis/trans isomeric, and ring-conformational states.",
    author="Durrant Lab",
    url="https://github.com/durrantlab/gypsum_dl",
    install_requires=["rdkit", "numpy", "scipy", "mpi4py", "setuptools"],
    packages=find_packages(),
)
