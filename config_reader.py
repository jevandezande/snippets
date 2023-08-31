import os
from importlib import resources
from typing import Any

import tomllib

PROJECT_NAME = ""


def read_configs(section: str | None = None) -> dict[str, Any]:
    """
    Reads configurations from PROJECT_ROOT/PROJECT_NAME/default_config.toml and
    ~/.config/PROJECT_NAME/config.toml if available

    :param section: section to select
    """
    defaults = read_default_config(section)

    try:
        defaults |= read_user_config(section)
    except (FileNotFoundError, KeyError):
        pass

    return defaults


def read_default_config(section: str | None = None) -> dict[str, Any]:
    """
    Reads configurations from PROJECT_ROOT/PROJECT_NAME/default_config.toml

    :param section: section to select
    """
    file = resources.files(f"{PROJECT_NAME}").joinpath("default_config.toml")
    with file.open("rb") as f:
        config = tomllib.load(f)

    if section:
        return config[section]  # type:ignore[no-any-return]

    return config


def read_user_config(section: str | None = None) -> dict[str, Any]:
    """
    Reads configurations from ~/.config/PROJECT_NAME/config.toml

    :param section: section to select
    """
    config_file = os.path.join(os.path.expanduser("~"), ".config", f"{PROJECT_NAME}", "config.toml")
    if not os.path.exists(config_file):
        raise FileNotFoundError("Could not find user {config_file=}")

    with open(config_file, "rb") as f:
        config = tomllib.load(f)

    if section:
        return config[section]  # type:ignore[no-any-return]

    return config
