# -*- coding: utf-8 -*-
"""
molvs.errors
~~~~~~~~~~~~

This module contains exceptions that are raised by MolVS.

"""

from __future__ import division, print_function, unicode_literals


class MolVSError(Exception):
    pass


class StandardizeError(MolVSError):
    pass


class ValidateError(MolVSError):
    pass


class StopValidateError(ValidateError):
    """Called by Validations to stop any further validations from being performed."""

    pass
