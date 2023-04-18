#!/usr/bin/env python3

import argparse
import os
import platform
import subprocess
from contextlib import contextmanager
from typing import Iterator
from os.path import expanduser

import requests
from bs4 import BeautifulSoup
from natsort import natsorted

OB_URL = "https://build-download.schrodinger.com/OB"
NAME_FORM = "Schrodinger_Suites_{version}_Advanced_{op_sys}-x86_64"


def cli() -> None:
    parser = argparse.ArgumentParser(description="CLI for Schrodinger Suites downloader/installer")
    parser.add_argument(
        "-v",
        "--version",
        help="Version to install [%(default)s]",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-b", "--build", help="Build to install [%(default)s]", type=str, default=None
    )
    parser.add_argument(
        "-o", "--os", help="Operating system [%(default)s]", type=str, default="{guess}"
    )
    parser.add_argument(
        "-D",
        "--skip_download",
        help="Skip download [%(default)s]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-I",
        "--skip_install",
        help="Skip install [%(default)s]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-Q",
        "--query",
        help="Query available versions and builds [%(default)s]",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-m",
        "--multi-install",
        nargs="*",
        help="Install multiple versions and builds (e.g. 2023-2/094 2023-3/013; defaults to last 4)",
    )
    parser.add_argument(
        "--verbose",
        help="Print information during process [%(default)s]",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    op_sys = platform.system() if args.os == "{guess}" else args.os

    if args.query:
        query(op_sys, args.version, args.build, verbose=True)
    elif args.multi_install is not None:
        if args.version or args.build:
            raise ValueError("Cannot use version or build flags with multi-install")

        multi_install = args.multi_install
        if len(args.multi_install) == 0:
            version, _ = query(op_sys, verbose=args.verbose)
            multi_install = [change_version(version, i) for i in range(-3, 1)]

        if args.verbose:
            print(f"Multi-install: {multi_install}")

        for v in multi_install:
            version, *b = v.split("/")
            build = b[0] if b else None
            install(
                op_sys,
                version,
                build,
                args.skip_download,
                args.skip_install,
                verbose=args.verbose,
            )
    else:
        install(
            op_sys,
            args.version,
            args.build,
            args.skip_download,
            args.skip_install,
            verbose=args.verbose,
        )


def install(
    op_sys: str,
    version: str | None = None,
    build: str | None = None,
    skip_download: bool = False,
    skip_install: bool = False,
    install_dir: str | None = None,
    scratch_dir: str | None = None,
    thirdparty_dir: str | None = None,
    verbose: bool = False,
) -> None:
    """
    Install the specified version/build, if unspecified, install the latest
    """
    if not skip_download:
        version, build = get_installer(op_sys, version, build, verbose)
    else:
        assert version

    name = NAME_FORM.format(op_sys=op_sys, version=version)
    install_dir = expanduser(install_dir or f"~/progs/schrodinger/{version}/{build}")
    scratch_dir = expanduser(scratch_dir or "/scr/")
    thirdparty_dir = expanduser(thirdparty_dir or f"{install_dir}/thirdparty")

    if not skip_install:
        if op_sys == "Windows":
            raise NotImplementedError("Windows is not yet supported")

        with cd(name):
            subprocess.check_call(
                [
                    "./INSTALL",
                    "-s",
                    install_dir,  # install directory
                    "-k",
                    scratch_dir,  # scratch directory
                    "-t",
                    thirdparty_dir,  # thirdparty modules directory
                    # "-b",  # batchmode
                ]
            )


def query_versions(
    op_sys: str,
    verbose: bool = False,
) -> list[str]:
    """
    Find all versions
    """
    ob_page = requests.get(OB_URL)
    ob_soup = BeautifulSoup(ob_page.text, features="html.parser")
    versions = natsorted(li.text for li in ob_soup.find_all(attrs={"class": "build"}))

    if verbose:
        print(f"Versions:\n{versions}\n")

    return versions


def query_builds(
    op_sys: str,
    version: str,
    verbose: bool = False,
) -> list[str]:
    """
    Find all builds corresponding to the given version
    >>> query_builds("linux", "2023-2")
    ['014', '028', '042', '055', '069', '081', '094']
    """
    version_url = f"{OB_URL}/{version}"
    version_page = requests.get(version_url)
    version_soup = BeautifulSoup(version_page.text, features="html.parser")

    builds = natsorted(li.text[-3:] for li in version_soup.find_all(attrs={"class": "build"}))
    if verbose:
        print(f"Builds:\n{builds}\n")

    return builds


def query(
    op_sys: str,
    version: str | None = None,
    build: str | None = None,
    verbose: bool = False,
) -> tuple[str, str]:
    """
    Find the latest version

    :return: version, build
    """
    version = version or query_versions(op_sys, verbose)[-1]
    build = build or query_builds(op_sys, version, verbose)[-1]

    if verbose:
        print(f"Selected: {version}/{build}")

    return version, build


def get_installer(
    op_sys: str,
    version: str | None = None,
    build: str | None = None,
    verbose: bool = False,
) -> tuple[str, str]:
    """
    Download the installer for the specified version

    :return: version, build
    """
    version, build = query(op_sys, version, build, verbose)
    name = NAME_FORM.format(op_sys=op_sys, version=version)
    tar_file = f"{name}.tar"
    url = f"{OB_URL}/{version}/build-{build}/{tar_file}"

    try:
        subprocess.check_call(["wget", url])
    except subprocess.CalledProcessError as e:
        e.add_note(
            f"Could not download {version}/{build}\n"
            + f"Perhaps try version: {change_version(version, -1)}"
        )
        raise

    subprocess.check_call(["tar", "xf", tar_file])

    return version, build


@contextmanager
def cd(
    directory: str,
    makedir: bool = False,
    exist_ok: bool = False,
    verbose: bool = False,
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

    os.chdir(expanduser(directory))

    try:
        if verbose:
            print(f"Switched to {directory=}")

        yield
    finally:
        os.chdir(previous_dir)

        if verbose:
            print(f"Switched to {previous_dir=}")


def change_version(version: str, increment: int) -> str:
    """
    Increment/decrement the version (assumes 4 releases in a year)
    >>> change_version("2023-2", 3)
    '2024-1'
    >>> change_version("2023-2", -2)
    '2022-4'
    """
    year, release = map(int, version.split("-"))

    year += (release + increment - 1) // 4
    release = (release + increment - 1) % 4 + 1

    return f"{year}-{release}"


if __name__ == "__main__":
    cli()
