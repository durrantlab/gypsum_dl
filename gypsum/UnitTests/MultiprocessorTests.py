from __future__ import absolute_import
import unittest
import time
from gypsum import mp_queue

def function_testing(thing, state, secs):
    time.sleep(secs)
    this = "The " + thing + " is " + state + "."
    return this

def single_input(input):
    return "Happiness is " + input + "."

class MultiprocessorTests(unittest.TestCase):
    """
    Base Test Suite
    """

    # Initialization and destruction for each test.

    def setUp(self):
        """
        Setting up the test molecule.
        """
        pass

    def tearDown(self):
        """
        Cleans up variables for the next test.
        """
        pass

    ### Tests
    # Loaders

    def test_one_job_one_processor_zero_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0)]

        this = mp_queue.MultiThreading(tasks, 1, function_testing)

        self.assertEqual(this[0], "The bun is done.")

    def test_three_jobs_one_processor_zero_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", 0),
                 ("scone", "burned", 0)
                ]

        this = mp_queue.MultiThreading(tasks, 1, function_testing)

        self.assertEqual(this[0], "The bun is done.")
        self.assertEqual(this[1], "The cake is baked.")
        self.assertEqual(this[2], "The scone is burned.")


    def test_three_jobs_one_processor_variable_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 1),
                 ("cake", "baked", 2),
                 ("scone", "burned", 3)
                ]

        this = mp_queue.MultiThreading(tasks, 1, function_testing)

        self.assertEqual(this[0], "The bun is done.")
        self.assertEqual(this[1], "The cake is baked.")
        self.assertEqual(this[2], "The scone is burned.")


    def test_three_jobs_three_processors_zero_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", 0),
                 ("scone", "burned", 0)
                ]

        this = mp_queue.MultiThreading(tasks, 3, function_testing)

        self.assertIn("The bun is done.", this)
        self.assertIn("The cake is baked.", this)
        self.assertIn("The scone is burned.", this)

    def test_three_jobs_three_processors_variable_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", .5),
                 ("scone", "burned", .1)
                ]

        this = mp_queue.MultiThreading(tasks, 3, function_testing)

        self.assertIn("The bun is done.", this)
        self.assertIn("The cake is baked.", this)
        self.assertIn("The scone is burned.", this)

    def test_five_jobs_one_processors_zero_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", 0),
                 ("scone", "burned", 0),
                 ("cupcake", "frosted", 0),
                 ("bread", "risen", 0)
                ]

        this = mp_queue.MultiThreading(tasks, 1, function_testing)

        self.assertIn("The bun is done.", this)
        self.assertIn("The cake is baked.", this)
        self.assertIn("The scone is burned.", this)
        self.assertIn("The cupcake is frosted.", this)
        self.assertIn("The bread is risen.", this)

    def test_five_jobs_three_processors_zero_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", 0),
                 ("scone", "burned", 0),
                 ("cupcake", "frosted", 0),
                 ("bread", "risen", 0)
                ]

        this = mp_queue.MultiThreading(tasks, 3, function_testing)

        self.assertIn("The bun is done.", this)
        self.assertIn("The cake is baked.", this)
        self.assertIn("The scone is burned.", this)
        self.assertIn("The cupcake is frosted.", this)
        self.assertIn("The bread is risen.", this)

    def test_five_jobs_one_processors_variable_sleep(self):
        """
        Tests a basic run of the MultiThreading function.
        """
        tasks = [("bun", "done", 0),
                 ("cake", "baked", 0.1),
                 ("scone", "burned", 0.5),
                 ("cupcake", "frosted", 0.2),
                 ("bread", "risen", 0.2)
                ]

        this = mp_queue.MultiThreading(tasks, 1, function_testing)

        self.assertIn("The bun is done.", this)
        self.assertIn("The cake is baked.", this)
        self.assertIn("The scone is burned.", this)
        self.assertIn("The cupcake is frosted.", this)
        self.assertIn("The bread is risen.", this)

    def test_flexible_inputs(self):
        """
        Testing basic tuple wrapping on individual inputs.
        """

        tasks = ["cats", "tea", "giant robots"]

        this = mp_queue.MultiThreading(tasks, 1, single_input)

        self.assertIn("Happiness is cats.", this)
        self.assertIn("Happiness is tea.", this)
        self.assertIn("Happiness is giant robots.", this)
