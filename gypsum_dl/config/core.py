from typing import Any

from pydantic import BaseModel, Field

from gypsum_dl.config.utils import Render, YamlIO


class GypsumConfig(BaseModel, YamlIO, Render):
    """Configuration for gypsum."""
