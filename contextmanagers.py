import os
from contextlib import contextmanager
from typing import Iterator


@contextmanager
def cd(
    directory: str, makedir: bool = False, exist_ok: bool = False, verbose: bool = False
) -> Iterator[None]:
    """
    Manage entering and exiting directories.

    :param directory: directory to enter
    :param makedir: make the directory
    :param exist_ok: if making a directory, don't error if it already exists
    :param verbose: print when entering and exiting a directory
    """
    if makedir:
        os.makedirs(directory, exist_ok=exist_ok)

    previous_dir = os.getcwd()

    os.chdir(os.path.expanduser(directory))

    try:
        if verbose:
            print(f"Switched to {directory=}")

        yield
    finally:
        os.chdir(previous_dir)

        if verbose:
            print(f"Switched to {previous_dir=}")
