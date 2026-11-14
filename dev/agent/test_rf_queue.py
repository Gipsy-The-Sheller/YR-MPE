"""Regression test: the RF dialog is non-modal and queue-attached.

Run from the repository root:

    python dev/agent/test_rf_queue.py
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

from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402
from YR_MPE.plugins.rf_plugin import RFComparisonPlugin  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class RFQueueTests(unittest.TestCase):
    def test_rf_dialog_is_non_modal_and_attached(self):
        with tempfile.TemporaryDirectory() as tmp:
            trees = Path(tmp) / "trees.nwk"
            trees.write_text("(a,b);\n(a,c);\n", encoding="utf-8")
            host = YR_MPEA_Widget()
            host.open_rf_wrapper(import_data=str(trees))
            dialog = host.rf_dialogs[-1]
            self.addCleanup(dialog.close)
            plugin = dialog.findChild(RFComparisonPlugin)
            self.assertIsNotNone(plugin)
            self.assertFalse(dialog.isModal())
            self.assertTrue(dialog.isVisible())
            self.assertIsNotNone(plugin.task_queue)


if __name__ == "__main__":
    unittest.main(verbosity=2)
