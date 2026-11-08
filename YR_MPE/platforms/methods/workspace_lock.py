"""Single-writer lock for a YR-MPE workspace.

A workspace is a directory shared by an index, imported sources, run manifests
and temporary artifacts.  Two workbenches must not write into the same
directory at the same time.  The lock is advisory and intentionally small:
an atomically-created marker in the workspace root, containing the owner PID
and host so that a stale marker left by a crashed process can be recovered.
"""

from __future__ import annotations

import ctypes
import json
import os
import socket
import time
import uuid
from typing import Optional, Tuple

LOCK_NAME = ".yrmpe.lock"


def pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    if os.name == "nt":
        # os.kill(pid, 0) is not a harmless probe on Windows.  OpenProcess
        # with query-only access is the safe equivalent.
        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(0x1000, False, pid)  # QUERY_LIMITED
        if handle:
            kernel32.CloseHandle(handle)
            return True
        return kernel32.GetLastError() == 5  # access denied means it exists
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _read_owner(path: str) -> Optional[dict]:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            value = json.load(handle)
        return value if isinstance(value, dict) else None
    except (OSError, ValueError):
        return None


class WorkspaceLock:
    """Advisory single-writer lock for one workspace path."""

    def __init__(self):
        self.path = ""
        self.owner: Optional[dict] = None
        self.token = uuid.uuid4().hex

    def acquire(self, workspace_path: str) -> Tuple[bool, str]:
        """Acquire ``workspace_path/.yrmpe.lock``.

        Returns ``(True, "")`` on success and ``(False, reason)`` when another
        live process owns the workspace.  A lock left by a dead process is
        removed and retried once.
        """
        self.release()
        workspace_path = os.path.abspath(str(workspace_path))
        try:
            os.makedirs(workspace_path, exist_ok=True)
        except OSError as exc:
            return False, "Cannot create workspace directory: %s" % exc

        lock_path = os.path.join(workspace_path, LOCK_NAME)
        payload = {
            "pid": os.getpid(),
            "host": socket.gethostname(),
            "started_at": time.time(),
            "application": "YPhylo Workbench",
            "token": self.token,
        }
        encoded = json.dumps(payload, ensure_ascii=False).encode("utf-8")

        for attempt in range(2):
            try:
                fd = os.open(
                    lock_path,
                    os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                    0o600,
                )
            except FileExistsError:
                current = _read_owner(lock_path)
                if current and current.get("token") == self.token:
                    # Re-entering the same workspace through the same lock.
                    self.path = workspace_path
                    self.owner = current
                    return True, ""
                if current and pid_is_alive(int(current.get("pid", 0))):
                    who = "pid %s on %s" % (
                        current.get("pid", "?"),
                        current.get("host", "another host"),
                    )
                    return False, (
                        "This workspace is already open (%s). "
                        "Close the other Workbench or choose another directory."
                        % who
                    )
                # Stale or unreadable marker: remove once and retry.
                if attempt == 0:
                    try:
                        os.unlink(lock_path)
                    except OSError:
                        pass
                    continue
                return False, "Workspace lock exists and could not be replaced."
            except OSError as exc:
                return False, "Cannot create workspace lock: %s" % exc
            else:
                with os.fdopen(fd, "wb") as handle:
                    handle.write(encoded)
                    handle.flush()
                    os.fsync(handle.fileno())
                self.path = workspace_path
                self.owner = payload
                return True, ""

        return False, "Workspace lock could not be acquired."

    def release(self) -> None:
        if not self.path:
            return
        lock_path = os.path.join(self.path, LOCK_NAME)
        current = _read_owner(lock_path)
        if current and current.get("token") == self.token:
            try:
                os.unlink(lock_path)
            except OSError:
                pass
        self.path = ""
        self.owner = None

    def is_held(self) -> bool:
        return bool(self.path and os.path.isfile(
            os.path.join(self.path, LOCK_NAME)
        ))


__all__ = ["WorkspaceLock", "LOCK_NAME", "pid_is_alive"]
