"""Start an agent bridge from About, as a YRTools-hosted window would."""

from __future__ import annotations

import json
import os
import sys
import time
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QPlainTextEdit, QPushButton

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.platforms.about_dialog import AboutDialog  # noqa: E402
from YR_MPE.platforms.mpea import YR_MPEA_Widget  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class AgentSessionStartTests(unittest.TestCase):
    def test_about_starts_a_token_protected_session(self):
        host = YR_MPEA_Widget()
        try:
            dialog = AboutDialog(
                host,
                bridge=None,
                start_session=host.start_agent_session,
            )
            start_buttons = [
                button for button in dialog.findChildren(QPushButton)
                if button.text() == "Start agent session"
            ]
            self.assertEqual(len(start_buttons), 1)
            self.assertTrue(start_buttons[0].isEnabled())

            start_buttons[0].click()

            bridge = dialog.bridge
            self.assertIsNotNone(bridge)
            self.assertTrue(bridge.token)
            self.assertTrue((bridge.session_dir / "READY").exists())
            prompts = dialog.findChildren(QPlainTextEdit)
            self.assertEqual(len(prompts), 1)
            self.assertIn("Agent token: %s" % bridge.token,
                          prompts[0].toPlainText())
            self.assertIn("workbench.invoke_plugin",
                          prompts[0].toPlainText())
            action_ids = {
                action["id"] for action in bridge.registry.list_actions()
            }
            self.assertIn("workbench.info", action_ids)
            self.assertIn("workbench.list_windows", action_ids)
            self.assertIn("workbench.invoke_plugin", action_ids)

            info = bridge.registry.invoke("workbench.info", {})
            self.assertTrue(info["ok"])
            self.assertEqual(info["result"]["application"],
                             "YPhylo Workbench")

            command = {
                "schema": "yrmpe.bridge.command.v1",
                "id": "test-info",
                "client": "unit-test",
                "op": "workbench.info",
                "params": {},
                "token": bridge.token,
            }
            with bridge.commands_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(command) + "\n")
            response = None
            for _ in range(50):
                host._poll_agent_bridge()
                if bridge.responses_path.exists():
                    lines = bridge.responses_path.read_text(
                        encoding="utf-8"
                    ).splitlines()
                    for line in lines:
                        candidate = json.loads(line)
                        if candidate.get("id") == "test-info":
                            response = candidate
                            break
                if response is not None:
                    break
                time.sleep(0.01)
            self.assertIsNotNone(response)
            self.assertTrue(response["ok"])
            self.assertEqual(response["result"]["application"],
                             "YPhylo Workbench")

            host.stop_agent_session()
            self.assertIsNone(host.agent_bridge)
            self.assertFalse((bridge.session_dir / "READY").exists())
        finally:
            host.stop_agent_session()
            host.close()
            host.deleteLater()

    def test_second_start_reuses_the_running_session(self):
        host = YR_MPEA_Widget()
        try:
            first = host.start_agent_session()
            second = host.start_agent_session()
            self.assertIs(first, second)
            self.assertIs(host.agent_bridge, first)
        finally:
            host.stop_agent_session()
            host.close()
            host.deleteLater()


if __name__ == "__main__":
    unittest.main(verbosity=2)
