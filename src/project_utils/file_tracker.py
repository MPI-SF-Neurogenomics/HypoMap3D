import hashlib
import os
import getpass
from datetime import datetime
from pathlib import Path

import toml

# ---------------------------------------------------------------------------- #
#                            Exposed Class(es)                                 #
# ---------------------------------------------------------------------------- #

class FileTracker:
    """
    Manager for tracking files and logging their metadata and checksums in a TOML file.

    FileTracker provides a simple interface to track files created or modified
    by your program. It automatically computes MD5 checksums, records file size
    and modification time, and allows optional descriptions for each file.
    All metadata is stored in a TOML log file.

    Features
    --------
    - Tracks files relative to the log file location.
    - Computes MD5 checksum for integrity verification.
    - Records file size (bytes and human-readable), modification timestamp, and user.
    - Supports optional file descriptions via `_TrackedFile.set_description()` or `.description`.
    - Returns a context manager for safe file writing and automatic logging.

    Attributes
    ----------
    log_file : pathlib.Path
        Absolute path to the TOML log file.

    Methods
    -------
    register(path)
        Returns a `_TrackedFile` context manager for writing and tracking a file.

    Examples
    --------
    >>> tracker = FileTracker("tracked_files.toml")
    >>> with tracker.register("some_file.txt") as tf:
    ...     tf.set_description("Example data file")
    ...     with open(tf.filepath, 'w') as f:
    ...         f.write("hello world")
    """

    def __init__(self, log_file="tracked_files.toml"):
        self.log_file = Path(log_file).absolute()
        
        if not self.log_file.parent.exists():
            raise FileNotFoundError(f'No such directory `{self.log_file.parent}`.')

    def register(self, path):
        """
        Return a context manager for writing to a file at `path` while tracking it.

        The returned `_TrackedFile` object handles computing checksum and logging
        metadata upon exiting the context.

        Parameters
        ----------
        path : str or pathlib.Path
            Path of the file to track.

        Returns
        -------
        _TrackedFile
            A context manager for the specified file.
        """
        return _TrackedFile(path, self.log_file)
    
    def __repr__(self):
        return f"<FileTracker log_file='{self.log_file}'>"

# ---------------------------------------------------------------------------- #
#                               Helper Class(es)                               #
# ---------------------------------------------------------------------------- #

class _TrackedFile:
    """
    Context manager for tracking a single file within a FileTracker.

    Handles computing MD5 checksum, logging metadata to a TOML file, and optional
    user descriptions. Designed to be returned by `FileTracker.register()`.

    Features
    --------
    - Computes MD5 checksum of file content on exit.
    - Logs metadata: size (bytes and human-readable), modification time, user, and description.
    - Stores file paths relative to the log file location.
    - Provides `.filepath` property for convenient access to the file path.
    - Supports `.description` attribute or `.set_description()` method.

    Attributes
    ----------
    path : pathlib.Path
        Absolute path of the file being tracked.
    log_file : pathlib.Path
        Path to the TOML log file where metadata is stored.
    checksum : str or None
        MD5 checksum of the file content; computed on exit.
    description : str
        Optional user-provided description; defaults to 'No description provided.'

    Properties
    ----------
    filepath : str
        Returns the file path as a string, suitable for use with `open()`.

    Methods
    -------
    set_description(descr)
        Set or update the description for this file.

    Context Management
    ------------------
    __enter__()
        Returns self to allow configuration (e.g., setting description).
    __exit__(exc_type, exc_val, exc_tb)
        Computes checksum, gathers metadata, and writes the record to the TOML log.
        Skips logging if an exception occurred during file writing.

    Usage Example
    -------------
    >>> tracker = FileTracker("tracked_files.toml")
    >>> with tracker.register("data/bounding_box.json") as tf:
    ...     tf.set_description("Bounding boxes for spatial dataset")
    ...     with open(tf.filepath, 'w') as out:
    ...         json.dump({"boxes": []}, out)
    ...     print(tf.checksum)
    ...     print(tf.description)
    """

    def __init__(self, path, log_file):
        self.path = Path(path)
        self.log_file = log_file
        self.checksum = None
        self.description = 'No description provided.'

    def set_description(self, descr):
        """Set or update the description for this file."""
        self.description = descr
    
    @property
    def filepath(self):
        """Return the file path as a string."""
        return str(self.path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None or not self.path.exists():
            return  # Skip logging if write failed

        # Compute checksum
        data = self.path.read_bytes()
        self.checksum = hashlib.md5(data).hexdigest()
        size = self.path.stat().st_size
        mod_time = datetime.fromtimestamp(self.path.stat().st_mtime).isoformat()
        user = getpass.getuser()

        # Relative path to log
        log_dir = self.log_file.parent
        rel_path = os.path.relpath(self.path, log_dir)

        # Load existing log or initialize
        log_data = {}
        if self.log_file.exists():
            log_data = toml.loads(self.log_file.read_text() or "{}")

        record = {
            "description": self.description,
            "relative_path": rel_path,
            "md5checksum": self.checksum,
            "size_bytes": size,
            "size_human": _human_readable_size(size),
            "modified": mod_time,
            "user": user,
        }

        log_data[rel_path] = record
        sorted_log = dict(sorted(log_data.items(), key=lambda kv: kv[0]))
        self.log_file.write_text(toml.dumps(sorted_log))

    def __repr__(self):
        return (f"<TrackedFile path='{self.path}' "
                f"checksum='{self.checksum or 'pending'}'>")

# ---------------------------------------------------------------------------- #
#                             Helper Functions(s)                              #
# ---------------------------------------------------------------------------- #

def _human_readable_size(size_bytes):
    """
    Convert a file size in bytes to a human-readable string.

    Parameters
    ----------
    size_bytes : int
        Size of the file in bytes.

    Returns
    -------
    str
        Human-readable string, e.g., '1.2 KB', '42 B', '3.4 MB'.
    """
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
