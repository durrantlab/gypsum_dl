import os

import pytest

from gypsum_dl import enable_logging

TEST_DIR = os.path.dirname(__file__)


@pytest.fixture
def test_dir():
    return os.path.abspath(TEST_DIR)


@pytest.fixture(scope="session", autouse=True)
def turn_on_logging():
    enable_logging(10)
