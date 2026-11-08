"""Helpers for starting console programs from a GUI process.

When YRTools is launched with ``pythonw.exe`` there is no attached console.
Windows then creates a console window for every console-subsystem child
process (IQ-TREE, MAFFT, MrBayes, ...).  ``CREATE_NO_WINDOW`` prevents that
allocation; ``STARTF_USESHOWWINDOW``/``SW_HIDE`` is kept as a fallback for
older Python versions.
"""

from __future__ import annotations

import os
import subprocess


def hidden_process_kwargs() -> dict:
    """Return Popen/run keyword arguments that suppress a console window."""
    if os.name != "nt":
        return {}
    kwargs = {}
    flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    if flags:
        kwargs["creationflags"] = flags
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = subprocess.SW_HIDE
    kwargs["startupinfo"] = startupinfo
    return kwargs


def _merge_hidden(kwargs: dict) -> dict:
    hidden = hidden_process_kwargs()
    if "creationflags" in kwargs:
        hidden["creationflags"] = kwargs["creationflags"] | hidden.get("creationflags", 0)
    kwargs.update(hidden)
    return kwargs


def run_hidden(cmd, **kwargs):
    """``subprocess.run`` with the console window suppressed on Windows."""
    return subprocess.run(cmd, **_merge_hidden(kwargs))


def popen_hidden(cmd, **kwargs):
    """``subprocess.Popen`` with the console window suppressed on Windows."""
    return subprocess.Popen(cmd, **_merge_hidden(kwargs))
