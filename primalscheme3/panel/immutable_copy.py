"""Independent file copies, with Darwin copy-on-write cloning when available.

This is a storage optimization only. Callers remain responsible for scientific
receipt/hash verification. No path creates a hard link or replaces a destination.
"""

from __future__ import annotations

import ctypes
import errno
import os
import shutil
import stat
import sys
from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=1)
def _clonefile_api():
    if sys.platform != "darwin":
        return None
    try:
        function = ctypes.CDLL(None, use_errno=True).clonefile
    except (AttributeError, OSError):
        return None
    function.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int)
    function.restype = ctypes.c_int
    return function


def _try_clonefile(source: Path, destination: Path) -> bool:
    function = _clonefile_api()
    if function is None:
        return False
    paths = tuple(os.fsencode(path) for path in (source, destination))
    if any(b"\0" in path for path in paths):
        raise ValueError("embedded null byte in copy path")
    # clonefile(2), flags=0: distinct inode/attributes, shared data blocks until
    # either file is written. Atomic success or no newly created destination.
    if function(*paths, 0) == 0:
        return True
    error = ctypes.get_errno()
    if error in (errno.ENOTSUP, errno.EXDEV, errno.ENOSYS, errno.EINVAL):
        return False
    raise OSError(error, os.strerror(error), str(destination))


def copy_immutable_file(source, destination) -> str:
    """Copy one regular file to a fresh path; return clonefile or copyfile.

    clonefile support depends on platform and filesystem. Unsupported/cross-
    volume clones fall back to a bounded-memory byte copy. Existing destinations,
    permission failures and I/O failures are errors, not overwrite requests.
    """
    source, destination = Path(source), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(errno.EEXIST, os.strerror(errno.EEXIST), str(destination))
    if not stat.S_ISREG(source.lstat().st_mode):
        raise ValueError("immutable copy source must be a regular file")
    if _try_clonefile(source, destination):
        return "clonefile"
    created = False
    try:
        with source.open("rb") as reader, destination.open("xb") as writer:
            created = True
            shutil.copyfileobj(reader, writer, length=1024 * 1024)
    except BaseException:
        if created:
            destination.unlink(missing_ok=True)
        raise
    return "copyfile"
