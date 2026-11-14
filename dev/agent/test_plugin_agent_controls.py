"""Regression tests for the generic plugin agent controls.

These cover two bridge bugs found during the open-task benchmark:

- setting an output path used to overwrite ``imported_files`` because every
  path control was treated as an input;
- ``plugin.run`` returned no ``result`` payload even when a job was created.

Run from the repository root:

    python dev/agent/test_plugin_agent_controls.py
"""

from __future__ import annotations

import os
import sys
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QLineEdit

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from YR_MPE.runtime.agent_bridge import AgentBridge  # noqa: E402
from YR_MPE.runtime.plugin_agent import agent_invoke, apply_control  # noqa: E402
from YR_MPE.templates.base_process_thread import BaseProcessThread  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class DummyPlugin:
    def __init__(self):
        self.plugin_name = "Dummy"
        self.tool_name = "dummy"
        self.tool_path = ""
        self.is_running = False
        self.imported_files = ["/data/in.fna"]
        self.input_edit = QLineEdit()
        self.output_edit = QLineEdit()
        self.run_button = None
        self.stop_button = None
        self.queue_job_id = ""

    def help_doc(self):
        return {"controls": [
            {"id": "input", "widget": "input_edit",
             "label": "Input alignment", "kind": "path", "cli": []},
            {"id": "output", "widget": "output_edit",
             "label": "Output tree", "kind": "path", "cli": []},
        ]}

    def run_analysis(self):
        self.queue_job_id = "job-123"
        return "job-123"


class PluginAgentControlTests(unittest.TestCase):
    def test_output_path_does_not_replace_inputs(self):
        plugin = DummyPlugin()
        self.assertTrue(apply_control(plugin, "output", "/tmp/out.nwk"))
        self.assertEqual(plugin.imported_files, ["/data/in.fna"])

        self.assertTrue(apply_control(plugin, "input", "/data/new.fna"))
        self.assertEqual(plugin.imported_files, ["/data/new.fna"])

    def test_dotted_advanced_dialog_widgets_resolve(self):
        from PyQt5.QtWidgets import QSpinBox

        class Dialog:
            def __init__(self):
                self.seed_spin = QSpinBox()

        class Plugin(DummyPlugin):
            def __init__(self):
                super().__init__()
                self.advanced_dialog = Dialog()

            def help_doc(self):
                return {"controls": [
                    {"id": "seed", "widget": "advanced_dialog.seed_spin",
                     "label": "Seed", "kind": "int", "cli": ["--seed", "{value}"]},
                ]}

        plugin = Plugin()
        self.assertTrue(apply_control(plugin, "seed", 42))
        self.assertEqual(plugin.advanced_dialog.seed_spin.value(), 42)

    def test_run_returns_result_payload(self):
        plugin = DummyPlugin()
        response = agent_invoke(plugin, "plugin.run")
        self.assertTrue(response["ok"])
        self.assertEqual(response["job_id"], "job-123")
        self.assertEqual(response["result"]["job_id"], "job-123")
        self.assertIn("status", response["result"])

    def test_no_skill_mode_strips_markdown(self):
        plugin = DummyPlugin()
        plugin.help_doc = lambda: {
            "markdown": "secret skill text", "summary": "summary",
            "controls": plugin.control_bindings() if hasattr(
                plugin, "control_bindings"
            ) else [],
        }
        previous = os.environ.get("YR_MPE_DISABLE_SKILLS")
        try:
            os.environ["YR_MPE_DISABLE_SKILLS"] = "1"
            response = agent_invoke(plugin, "plugin.help_doc")
            self.assertTrue(response["result"]["skills_disabled"])
            self.assertNotIn("markdown", response["result"])
            self.assertNotIn("summary", response["result"])
        finally:
            if previous is None:
                os.environ.pop("YR_MPE_DISABLE_SKILLS", None)
            else:
                os.environ["YR_MPE_DISABLE_SKILLS"] = previous

    def test_session_temp_dir_is_used_by_process_threads(self):
        import tempfile
        from pathlib import Path

        previous = os.environ.get("YR_MPE_SESSION_TMP")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                os.environ["YR_MPE_SESSION_TMP"] = tmp
                thread = BaseProcessThread("dummy-tool", [], [])
                self.assertEqual(thread.workdir, tmp)
                self.assertTrue(Path(tmp).is_dir())
        finally:
            if previous is None:
                os.environ.pop("YR_MPE_SESSION_TMP", None)
            else:
                os.environ["YR_MPE_SESSION_TMP"] = previous

    def test_bridge_creates_session_temp_dir(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            bridge = AgentBridge(session_dir=tmp, allow_anonymous=True)
            self.assertTrue((Path(tmp) / "tmp").is_dir())
            bridge.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
