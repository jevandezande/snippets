"""Tests for the contextmanagers."""

from pathlib import Path

from pytest import raises

from snippets.contextmanagers import cd


def test_cd(tmpdir: Path) -> None:
    """Test cd."""

    def cwd() -> Path:
        return Path(".").cwd()

    start_path = cwd()

    with cd(tmpdir):
        tmpdir_path = cwd()
        assert tmpdir_path != start_path

        with cd("new_dir", mkdir=True):
            assert tmpdir_path / "new_dir" == cwd()
        assert tmpdir_path == cwd()

        with cd("new_dir", mkdir=True, exist_ok=True):
            pass

        with raises(FileExistsError):
            with cd("new_dir", mkdir=True):
                pass

        tmpfile = "file.tmp"
        with open(tmpfile, "w") as f:
            f.write("Temporary data")

        with cd(tmpdir_path / "new_dir", copy_along=[tmpfile], verbose=True):
            assert (cwd() / tmpfile).is_file()

    assert start_path == cwd()
