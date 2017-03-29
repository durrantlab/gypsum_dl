from __future__ import absolute_import
import unittest
from . import MultiprocessorTests as mt

class UnitTests(object):
    """
    Unit testing object for scoria.
    """
    def __init__(self):
        """
        Initalizes the unit tests.
        """
        self._suite = unittest.TestSuite()
        self._runner = unittest.TextTestRunner()

    # Running Suite

    def run(self):
        """
        Runs the currently queued suite of tests.
        """
        self._runner.run(self._suite)

    def run_all(self):
        """
        Quickly runs all unit tests.
        """
        self.add_all_tests()
        self.run()

    # Add specific module tests

    def add_all_tests(self):
        """
        Adds all available tests to the suite.
        """
        self.add_multiprocessor_tests()

    def add_multiprocessor_tests(self):
        """
        Adds the information tests.
        """
        information_tests = unittest.makeSuite(mt.MultiprocessorTests)
        self._suite.addTests(information_tests)
