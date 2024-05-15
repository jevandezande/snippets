"""Configuration reader."""

from importlib import resources
from importlib.resources.abc import Traversable
from pathlib import Path
from typing import Any

import tomllib

PROJECT_NAME = "snippets"


def read_configs(section: str | None = None) -> dict[str, Any]:
    """
    Read configurations.

    Reads:
        PROJECT_ROOT/PROJECT_NAME/default_config.toml
        ~/.config/PROJECT_NAME/config.toml

    :param section: section to select
    :return: merged configurations
    """
    config: dict[str, Any] = {}

    for reader in (read_default_config, read_user_config):
        try:
            config |= reader(section)
        except (FileNotFoundError, KeyError):
            pass

    return config


def read_config(file: Path | Traversable, section: str | None = None) -> dict[str, Any]:
    """
    Read specified configuration.

    :param file: file to read
    :param section: section to select
    :return: configurations
    """
    with file.open("rb") as f:
        if section:
            return tomllib.load(f)[section]  # type:ignore[no-any-return]
        return tomllib.load(f)


def read_default_config(section: str | None = None) -> dict[str, Any]:
    """
    Read configurations from PROJECT_ROOT/PROJECT_NAME/default_config.toml.

    :param section: section to select
    :return: configurations
    """
    file = resources.files(PROJECT_NAME).joinpath("default_config.toml")
    return read_config(file, section)


def read_user_config(section: str | None = None) -> dict[str, Any]:
    """
    Read configurations from ~/.config/PROJECT_NAME/config.toml.

    :param section: section to select
    :return: configurations
    """
    return read_config(Path(f"~/.config/{PROJECT_NAME}/config.toml").expanduser(), section)
