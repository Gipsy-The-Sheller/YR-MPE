"""Tests for the Windows console-window suppression helpers."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.runtime.subprocess_utils import (  # noqa: E402
    hidden_process_kwargs,
    popen_hidden,
    run_hidden,
)


class SubprocessUtilsTests(unittest.TestCase):
    def test_hidden_kwargs(self):
        kwargs = hidden_process_kwargs()
        if os.name == "nt":
            self.assertTrue(kwargs["creationflags"] & subprocess.CREATE_NO_WINDOW)
            self.assertTrue(kwargs["startupinfo"].dwFlags & subprocess.STARTF_USESHOWWINDOW)
        else:
            self.assertEqual(kwargs, {})

    def test_run_hidden(self):
        result = run_hidden(
            [sys.executable, "-c", "print('hidden-ok')"],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("hidden-ok", result.stdout)

    def test_popen_hidden(self):
        process = popen_hidden(
            [sys.executable, "-c", "print('popen-ok')"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True,
        )
        stdout, _ = process.communicate(timeout=30)
        self.assertEqual(process.returncode, 0)
        self.assertIn("popen-ok", stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
