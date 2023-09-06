from pathlib import Path

from pytest import raises

from snippets.config_reader import (
    read_config,
    read_configs,
    read_default_config,
    read_user_config,
)


def test_read_configs() -> None:
    assert read_configs() == {"test_section": {"spam": "eggs"}}
    assert read_configs("test_section") == {"spam": "eggs"}
    assert read_configs("spam") == {}


def test_read_config() -> None:
    with raises(FileNotFoundError):
        read_config(Path("spam/eggs.conf"))


def test_read_default_config() -> None:
    assert read_default_config() == {"test_section": {"spam": "eggs"}}


def test_read_user_config() -> None:
    with raises(FileNotFoundError):
        read_user_config()
