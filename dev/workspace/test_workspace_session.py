"""Tests for temporary-session cleanup and the workspace status strip."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class WorkspaceSessionTests(unittest.TestCase):
    def test_temporary_session_is_removed_on_close(self):
        host = YR_MPEA_Widget()
        session = host._temporary_workdir
        self.assertTrue(os.path.isdir(session))
        self.assertIn("Temporary workspace", host.workspace_status_label.text())
        self.assertEqual(host.get_workdir(), session)
        host.close()
        APP.processEvents()
        self.assertFalse(os.path.exists(session))


if __name__ == "__main__":
    unittest.main(verbosity=2)
