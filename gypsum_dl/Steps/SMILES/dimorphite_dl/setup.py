#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="dimorphite_dl",
    version="1.2.4",
    description="Dimorphite-DL adds hydrogen atoms to molecular representations, as appropriate for a user-specified pH range. It is a fast, accurate, accessible, and modular open-source program for enumerating small-molecule ionization states.",
    author="Durrant Lab",
    url="https://github.com/durrantlab/dimorphite_dl/",
    install_requires=["rdkit", "numpy", "scipy", "mpi4py", "setuptools"],
    packages=find_packages(),
)
