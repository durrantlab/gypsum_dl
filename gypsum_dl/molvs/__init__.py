# -*- coding: utf-8 -*-
"""
MolVS - Molecule Validation and Standardization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolVS is a python tool built on top of RDKit that performs validation and standardization of chemical structures.

"""

from __future__ import division, print_function, unicode_literals

import logging

from .errors import MolVSError, StandardizeError, ValidateError
from .standardize import (
    Standardizer,
    canonicalize_tautomer_smiles,
    enumerate_tautomers_smiles,
    standardize_smiles,
)
from .validate import Validator, validate_smiles

__title__ = "MolVS"
__version__ = "0.1.1"
__author__ = "Matt Swain"
__email__ = "m.swain@me.com"
__license__ = "MIT"
__copyright__ = "Copyright 2019 Matt Swain"


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())
