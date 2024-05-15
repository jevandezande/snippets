"""Context managers for common tasks."""

import os
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator


@contextmanager
def cd(
    directory: str | Path,
    mkdir: bool = False,
    exist_ok: bool = False,
    copy_along: Iterable[str | Path] = (),
    verbose: bool = False,
) -> Iterator[None]:
    """
    Manage entering and exiting directories.

    :param directory: directory to enter
    :param mkdir: make the directory
    :param exist_ok: if making a directory, don't error if it already exists
    :param copy_along: copy files and directory along to directory
    :param verbose: print when entering and exiting a directory
    :yield: None
    """
    directory = Path(directory)

    if mkdir:
        directory.mkdir(parents=True, exist_ok=exist_ok)

    previous_dir = Path.cwd()

    for item in copy_along:
        if verbose:
            print(f"Copying: {item}")
        shutil.copy(item, directory)

    os.chdir(directory)

    try:
        if verbose:
            print(f"Switched to {directory=}")

        yield
    finally:
        os.chdir(previous_dir)

        if verbose:
            print(f"Switched to {previous_dir=}")
