import shutil
import os
from contextlib import contextmanager
from typing import Iterator, Iterable
from pathlib import Path


@contextmanager
def cd(
    directory: str | Path,
    makedir: bool = False,
    exist_ok: bool = False,
    copy_along: Iterable[str] = (),
    verbose: bool = False,
) -> Iterator[None]:
    """
    Manage entering and exiting directories.

    :param directory: directory to enter
    :param makedir: make the directory
    :param exist_ok: if making a directory, don't error if it already exists
    :param copy_along: copy files and directory along to directory
    :param verbose: print when entering and exiting a directory
    """
    directory = Path(directory)

    if makedir:
        directory.mkdir(parents=True, exist_ok=True)

    previous_dir = Path.cwd()

    for item in copy_along:
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
