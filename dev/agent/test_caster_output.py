"""Regression test for CASTER-site output handling.

ASTER writes the tree to ``--output``. The plugin used to reopen that file
with mode ``"w"`` and write an empty stdout string over it, discarding a valid
tree. This test uses a tiny Python stub in place of the binary.

Run from the repository root:

    python dev/agent/test_caster_output.py
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from YR_MPE.plugins.caster_site_plugin import CasterSiteThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class CasterOutputTests(unittest.TestCase):
    def test_output_written_by_tool_is_not_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            alignment = root / "input.fas"
            alignment.write_text(">a\nATGC\n>b\nATGC\n", encoding="utf-8")
            output = root / "tree.tre"
            code = (
                "import pathlib, sys\n"
                "i = sys.argv.index('--output')\n"
                "pathlib.Path(sys.argv[i + 1]).write_text('(a,b);\\n')\n"
            )
            thread = CasterSiteThread(
                sys.executable,
                [str(alignment)],
                ["-c", code, "--output", str(output)],
                output_file=str(output),
                workdir=tmp,
            )
            thread.is_running = True
            thread.execute_commands()
            self.assertTrue(output.is_file())
            self.assertEqual(output.read_text(encoding="utf-8").strip(), "(a,b);")


if __name__ == "__main__":
    unittest.main(verbosity=2)
