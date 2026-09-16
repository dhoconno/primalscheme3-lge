"""Independent immutable-origin copies, including real APFS clone semantics."""

import ctypes
import errno
import importlib.util
import sys

import pytest


def api():
    name = "primalscheme3.panel.immutable_copy"
    assert importlib.util.find_spec(name) is not None, (
        "independent copy helper is missing"
    )
    return __import__(name, fromlist=["copy_immutable_file"])


@pytest.mark.parametrize("fallback", [False, True])
def test_copies_are_independent_on_edit_delete_and_relocation(
    tmp_path, monkeypatch, fallback
):
    module = api()
    if fallback:
        monkeypatch.setattr(module, "_clonefile_api", lambda: None)
    source = tmp_path / "source"
    original = bytes(range(256)) * 8192
    source.write_bytes(original)
    destination = tmp_path / "destination"
    method = module.copy_immutable_file(source, destination)
    assert method in ("clonefile", "copyfile")
    if fallback:
        assert method == "copyfile"
    assert destination.read_bytes() == original
    assert source.stat().st_ino != destination.stat().st_ino
    assert source.stat().st_nlink == destination.stat().st_nlink == 1
    with destination.open("r+b") as handle:
        handle.write(b"DEST")
    assert source.read_bytes() == original
    with source.open("r+b") as handle:
        handle.seek(100)
        handle.write(b"SOURCE")
    expected = b"DEST" + original[4:]
    assert destination.read_bytes() == expected
    relocated = tmp_path / "relocated"
    destination.rename(relocated)
    source.unlink()
    assert relocated.read_bytes() == expected
    another = tmp_path / "another"
    module.copy_immutable_file(relocated, another)
    another.unlink()
    assert relocated.read_bytes() == expected


def test_real_darwin_clonefile_path_when_supported(tmp_path):
    if sys.platform != "darwin":
        pytest.skip("Darwin clonefile is not available on this platform")
    module = api()
    source = tmp_path / "source"
    source.write_bytes(b"clone me" * 8192)
    destination = tmp_path / "clone"
    if not module._try_clonefile(source, destination):
        pytest.skip("filesystem does not support clonefile")
    assert source.read_bytes() == destination.read_bytes()
    assert source.stat().st_ino != destination.stat().st_ino
    destination.write_bytes(b"independent")
    assert source.read_bytes() == b"clone me" * 8192


@pytest.mark.parametrize(
    "error", [errno.ENOTSUP, errno.EXDEV, errno.ENOSYS, errno.EINVAL]
)
def test_unsupported_clone_falls_back_portably(tmp_path, monkeypatch, error):
    module = api()

    def unsupported(*args):
        ctypes.set_errno(error)
        return -1

    monkeypatch.setattr(module, "_clonefile_api", lambda: unsupported)
    source, destination = tmp_path / "source", tmp_path / "destination"
    source.write_bytes(b"portable bytes")
    assert module.copy_immutable_file(source, destination) == "copyfile"
    assert destination.read_bytes() == source.read_bytes()
    assert source.stat().st_ino != destination.stat().st_ino


def test_existing_destination_and_clone_io_errors_never_overwrite(
    tmp_path, monkeypatch
):
    module = api()
    source, destination = tmp_path / "source", tmp_path / "destination"
    source.write_bytes(b"source")
    destination.write_bytes(b"keep")
    with pytest.raises(FileExistsError):
        module.copy_immutable_file(source, destination)
    assert destination.read_bytes() == b"keep"
    destination.unlink()

    def failed(*args):
        ctypes.set_errno(errno.EIO)
        return -1

    monkeypatch.setattr(module, "_clonefile_api", lambda: failed)
    with pytest.raises(OSError) as error:
        module.copy_immutable_file(source, destination)
    assert error.value.errno == errno.EIO
    assert not destination.exists()
    assert source.read_bytes() == b"source"


def test_fallback_exclusive_creation_and_partial_cleanup(tmp_path, monkeypatch):
    module = api()
    source, destination = tmp_path / "source", tmp_path / "destination"
    source.write_bytes(b"source")

    def raced(*args):
        destination.write_bytes(b"other writer")
        return False

    monkeypatch.setattr(module, "_try_clonefile", raced)
    with pytest.raises(FileExistsError):
        module.copy_immutable_file(source, destination)
    assert destination.read_bytes() == b"other writer"
    destination.unlink()
    monkeypatch.setattr(module, "_try_clonefile", lambda *args: False)

    def broken(reader, writer, **kwargs):
        writer.write(b"partial")
        raise OSError(errno.ENOSPC, "disk full")

    monkeypatch.setattr(module.shutil, "copyfileobj", broken)
    with pytest.raises(OSError):
        module.copy_immutable_file(source, destination)
    assert not destination.exists()
    assert source.read_bytes() == b"source"


def test_helper_rejects_symlink_source_and_destination(tmp_path):
    module = api()
    source = tmp_path / "source"
    source.write_bytes(b"source")
    link = tmp_path / "link"
    link.symlink_to(source)
    with pytest.raises(ValueError, match="regular file"):
        module.copy_immutable_file(link, tmp_path / "new")
    with pytest.raises(FileExistsError):
        module.copy_immutable_file(source, link)
    assert source.read_bytes() == b"source"
