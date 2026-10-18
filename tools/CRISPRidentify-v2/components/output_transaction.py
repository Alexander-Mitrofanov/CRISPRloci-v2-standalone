"""Rollback-capable publication for multi-file CRISPRidentify outputs."""

import errno
import os
from pathlib import Path
import tempfile


_UNSUPPORTED_DIRECTORY_FSYNC = frozenset((
    errno.EBADF,
    errno.EINVAL,
    getattr(errno, "ENOTSUP", errno.EINVAL),
))


def _absolute_path(path):
    return Path(os.path.abspath(os.fspath(path)))


def _lexists(path):
    return os.path.lexists(os.fspath(path))


def _unlink_if_present(path):
    try:
        os.unlink(os.fspath(path))
    except FileNotFoundError:
        pass


def _fsync_directory(directory):
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(os.fspath(directory), flags)
    except OSError as error:
        if error.errno in _UNSUPPORTED_DIRECTORY_FSYNC:
            return
        raise
    try:
        try:
            os.fsync(descriptor)
        except OSError as error:
            if error.errno not in _UNSUPPORTED_DIRECTORY_FSYNC:
                raise
    finally:
        os.close(descriptor)


def _fsync_directories(paths):
    for directory in sorted(
            {path.parent for path in paths
             if path.parent.exists()},
            key=lambda item: os.fspath(item)):
        _fsync_directory(directory)


def remove_completion_marker(path):
    """Remove and durably retire a stale successful-run marker."""

    target = _absolute_path(path)
    try:
        os.unlink(os.fspath(target))
    except FileNotFoundError:
        return
    _fsync_directory(target.parent)


class OutputTransaction:
    """Stage a set of files and publish it with rollback on failure.

    A completion target, when supplied to publish(), is removed before
    publication and replaced last. It is never restored during rollback:
    its presence therefore means every preceding artifact replacement and
    directory sync completed successfully.
    """

    def __init__(self):
        self._artifacts = {}
        self._committed = False

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception, traceback):
        if not self._committed:
            self.discard()
        return False

    def _register(self, target, staged):
        if self._committed:
            raise RuntimeError("output transaction is already committed")
        target = _absolute_path(target)
        if target in self._artifacts:
            raise ValueError(
                "output target was staged more than once: {}".format(target)
            )
        self._artifacts[target] = staged
        return target

    def stage_writer(
            self, target, writer, *, encoding="utf-8", newline=None):
        """Run writer(stream) into a durable same-directory stage file."""

        return self._stage_writer(target, writer, encoding=encoding,
                                  newline=newline, binary=False)

    def stage_binary_writer(self, target, writer):
        """Stage a binary artifact with the same rollback/durability rules."""

        return self._stage_writer(target, writer, binary=True)

    def _stage_writer(self, target, writer, *, binary,
                      encoding="utf-8", newline=None):

        target = _absolute_path(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        descriptor, staged_name = tempfile.mkstemp(
            dir=os.fspath(target.parent),
            prefix=".crispridentify-stage-",
        )
        staged = Path(staged_name)
        try:
            stream = (os.fdopen(descriptor, mode="wb") if binary else
                      os.fdopen(descriptor, mode="w", encoding=encoding,
                                newline=newline))
            with stream:
                writer(stream)
                stream.flush()
                os.fsync(stream.fileno())
        except BaseException:
            _unlink_if_present(staged)
            raise
        try:
            self._register(target, staged)
        except BaseException:
            _unlink_if_present(staged)
            raise
        return staged

    def stage_text(self, target, content):
        if not isinstance(content, str):
            raise TypeError("staged text content must be a string")
        return self.stage_writer(target, lambda stream: stream.write(content))

    def delete(self, target):
        """Stage removal of an optional artifact as part of this transaction."""

        self._register(target, None)

    def discard(self):
        """Remove unpublished stage files without changing output targets."""

        for staged in self._artifacts.values():
            if staged is not None:
                _unlink_if_present(staged)
        self._artifacts.clear()

    @property
    def targets(self):
        return tuple(self._artifacts)

    def _validated_order(self, completion_target):
        if not self._artifacts:
            raise ValueError("cannot publish an empty output transaction")
        completion = (
            None if completion_target is None
            else _absolute_path(completion_target)
        )
        if completion is not None:
            if completion not in self._artifacts:
                raise ValueError("completion target was not staged")
            if self._artifacts[completion] is None:
                raise ValueError("completion target cannot be a deletion")

        ordered = [
            target for target in self._artifacts
            if target != completion
        ]
        if completion is not None:
            ordered.append(completion)
        for target in ordered:
            staged = self._artifacts[target]
            if staged is not None:
                if not staged.is_file():
                    raise FileNotFoundError(
                        "staged output disappeared: {}".format(staged)
                    )
                if staged.parent != target.parent:
                    raise ValueError(
                        "staged output must share the target directory"
                    )
        return ordered, completion

    @staticmethod
    def _backup_path(target):
        descriptor, name = tempfile.mkstemp(
            dir=os.fspath(target.parent),
            prefix=".crispridentify-backup-",
        )
        os.close(descriptor)
        return Path(name)

    def publish(self, completion_target=None):
        """Replace every target, with the completion target strictly last."""

        if self._committed:
            raise RuntimeError("output transaction is already committed")
        ordered, completion = self._validated_order(completion_target)
        if completion is not None:
            remove_completion_marker(completion)

        non_completion = [
            target for target in ordered if target != completion
        ]
        backups = {}
        backed_up = set()
        published = set()
        try:
            for target in non_completion:
                if not _lexists(target):
                    backups[target] = None
                    continue
                backup = self._backup_path(target)
                backups[target] = backup
                os.replace(os.fspath(target), os.fspath(backup))
                backed_up.add(target)

            for target in non_completion:
                staged = self._artifacts[target]
                if staged is None:
                    continue
                os.replace(os.fspath(staged), os.fspath(target))
                published.add(target)

            _fsync_directories(non_completion)

            if completion is not None:
                staged = self._artifacts[completion]
                os.replace(os.fspath(staged), os.fspath(completion))
                _fsync_directory(completion.parent)
        except BaseException:
            # The completion marker is removed first and never restored.
            if completion is not None:
                _unlink_if_present(completion)
            restoration_errors = []
            for target in reversed(non_completion):
                try:
                    backup = backups.get(target)
                    if target in backed_up and backup is not None:
                        _unlink_if_present(target)
                        os.replace(os.fspath(backup), os.fspath(target))
                        backed_up.remove(target)
                    elif target in published:
                        _unlink_if_present(target)
                except BaseException as error:
                    restoration_errors.append(error)
            try:
                _fsync_directories(non_completion)
            except BaseException as error:
                restoration_errors.append(error)
            self.discard()
            for target, backup in backups.items():
                if backup is not None and target not in backed_up:
                    _unlink_if_present(backup)
            if restoration_errors:
                raise RuntimeError(
                    "output publication failed and rollback was incomplete; "
                    "recoverable backups were retained: {}".format(
                        restoration_errors[0]
                    )
                )
            raise

        self._committed = True
        self._artifacts.clear()
        for backup in backups.values():
            if backup is not None:
                try:
                    _unlink_if_present(backup)
                except OSError:
                    # Publication is complete; a hidden backup is safer than
                    # invalidating the successful completion marker.
                    pass
        try:
            _fsync_directories(non_completion)
        except OSError:
            # All replacements and the completion sync already succeeded.
            pass
