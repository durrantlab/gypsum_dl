from typing import Any

from abc import ABC

import yaml


class YamlIO(ABC):
    """Handles YAML inputs and outputs."""

    def update(self, data: dict[str, Any]) -> None:
        """Iteratively update pydantic model.

        Args:
            data: Key-value mapping to update attributes with.

        Notes:
            Many of our pydantic models have fields that need to specified by
            instantiating other objects or models. In order to instantiate these
            objects, you can use the `import` root key to specify which class to use.
            For example, if you need to specify the
            [`AmberTopoGen`][simulation.amber.topo.AmberTopoGen] as the
            simlify configuration `engine`, you can add this to your YAML file.

            ```yaml
            import:
                engine: simlify.simulation.amber.topo.AmberTopoGen
            ```

            This will call [`get_obj_from_string`][utils.get_obj_from_string] and
            set the [`engine`][config.SimlifyConfig.engine] attribute to
            [`AmberTopoGen`][simulation.amber.topo.AmberTopoGen].
            We handle these imports before any other field is handled.
        """
        for key, value in data.items():
            if key in self.model_fields:  # type: ignore # pylint: disable=no-member
                setattr(self, key, value)

    def from_yaml(self, yaml_paths: str | list[str]) -> None:
        """Update the instance's attributes from one or more YAML files.

        Args:
            yaml_paths: A sequence of YAML file paths or a single YAML file path.

        Raises:
            FileNotFoundError: If any of the YAML files cannot be found.
        """
        if isinstance(yaml_paths, str):
            yaml_paths = [yaml_paths]
        for yaml_path in yaml_paths:
            with open(yaml_path, "r", encoding="utf-8") as f:
                yaml_data = yaml.safe_load(f)
            self.update(yaml_data)

    def to_yaml(self, file_path: str) -> None:
        """Serialize a Pydantic BaseModel instance to a YAML file.

        Args:
            file_path: Path to the YAML file to write the serialized data to.

        Raises:
            YamlIOError: If the file cannot be written to.
        """
        config_dict = self.model_dump()  # type: ignore # pylint: disable=no-member
        with open(file_path, "w", encoding="utf-8") as f:
            yaml.dump(config_dict, f, default_flow_style=False)
