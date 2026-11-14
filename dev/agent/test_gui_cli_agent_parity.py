"""Check that the agent controls and the CLI preview tell the same story.

The same control value is applied through the agent API, exposed in the live
help document, and checked against the command the plugin would execute. This
is a regression test for the shared semantic layer, not a benchmark.

Run from the repository root:

    python dev/agent/test_gui_cli_agent_parity.py
"""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from YR_MPE.runtime.plugin_agent import apply_control  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)


class GuiCliAgentParityTests(unittest.TestCase):
    def test_mafft_control_and_preview_agree(self):
        from YR_MPE.plugins.mafft_plugin import MAFFTPlugin

        plugin = MAFFTPlugin(headless=True)
        self.assertTrue(apply_control(plugin, "input", "/tmp/in.fna"))
        self.assertTrue(apply_control(plugin, "threads", 4))
        self.assertTrue(apply_control(plugin, "mode", "Auto mode (--auto)"))
        self.assertTrue(apply_control(
            plugin, "direction",
            "Adjust direction accurately (--adjustdirectionaccurately)",
        ))

        rows = {row["id"]: row for row in plugin.help_doc()["controls"]}
        self.assertEqual(rows["threads"]["cli_args"], ["--thread", "4"])
        self.assertEqual(rows["mode"]["cli_args"], ["--auto"])
        self.assertEqual(
            rows["direction"]["cli_args"], ["--adjustdirectionaccurately"]
        )

        preview = plugin.cli_preview()
        for token in (
            "--auto", "--thread", "4", "--adjustdirectionaccurately",
            "/tmp/in.fna",
        ):
            self.assertIn(token, preview)

    def test_iqtree_control_and_preview_agree(self):
        from YR_MPE.plugins.iqtree_plugin import IQTreePlugin

        plugin = IQTreePlugin(headless=True)
        self.assertTrue(apply_control(plugin, "input", "/tmp/in.fna"))
        self.assertTrue(apply_control(plugin, "seq_type", "DNA"))
        self.assertTrue(apply_control(plugin, "threads", 4))

        rows = {row["id"]: row for row in plugin.help_doc()["controls"]}
        self.assertEqual(rows["seq_type"]["cli_args"], ["-st", "DNA"])
        self.assertEqual(rows["threads"]["cli_args"], ["-nt", "4"])

        preview = plugin.cli_preview()
        for token in ("-s", "/tmp/in.fna", "-st", "DNA", "-nt", "4"):
            self.assertIn(token, preview)

    def test_astral_advanced_controls_are_reachable(self):
        from YR_MPE.plugins.astral_plugin import AstralPlugin

        plugin = AstralPlugin(headless=True)
        self.assertTrue(apply_control(plugin, "seed", 999))
        self.assertTrue(apply_control(plugin, "length", "CULength"))
        self.assertEqual(plugin.advanced_dialog.seed_spin.value(), 999)
        self.assertEqual(
            plugin.advanced_dialog.length_combo.currentText(), "CULength"
        )
        preview = plugin.cli_preview()
        self.assertIn("astral4", preview)
        self.assertIn("--seed", preview)
        self.assertIn("--length", preview)

    def test_iqtree_partition_action(self):
        import tempfile
        from pathlib import Path
        from YR_MPE.plugins.iqtree_plugin import IQTreePlugin

        with tempfile.TemporaryDirectory() as tmp:
            alignment = Path(tmp) / "super.fas"
            alignment.write_text(
                ">a\nATGCATGCATGCATGCATGC\n>b\nATGCATGCATGCATGCATGC\n",
                encoding="utf-8",
            )
            plugin = IQTreePlugin(headless=True)
            response = plugin.agent_invoke(
                "plugin.apply_partitions",
                {
                    "mode": "TUL",
                    "alignment_file": str(alignment),
                    "partitions": [
                        {"name": "p1", "range": "1-10", "seq_type": "DNA"},
                        {"name": "p2", "range": "11-20", "seq_type": "DNA"},
                    ],
                },
            )
            self.assertTrue(response["ok"])
            inspect = plugin.agent_invoke("plugin.inspect_partitions")
            self.assertEqual(inspect["result"]["mode"], "TUL")
            self.assertEqual(len(inspect["result"]["partitions"]), 2)

    def test_iqtree_partition_outputs_resolve_in_session_workdir(self):
        import tempfile
        from YR_MPE.plugins.iqtree_partition_thread import IQTreePartitionThread
        from YR_MPE.plugins.partition_mode import PartitionMode

        previous = os.environ.get("YR_MPE_SESSION_TMP")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                os.environ["YR_MPE_SESSION_TMP"] = tmp
                thread = IQTreePartitionThread(
                    "iqtree3", "/tmp/input.fas", [], PartitionMode.EL, {}
                )
                thread.output_prefix = "supermatrix"
                self.assertEqual(
                    Path(thread._output_path(".treefile")),
                    Path(tmp) / "supermatrix.treefile",
                )
        finally:
            if previous is None:
                os.environ.pop("YR_MPE_SESSION_TMP", None)
            else:
                os.environ["YR_MPE_SESSION_TMP"] = previous

    def test_rf_result_spec_exports_matrix(self):
        import tempfile
        from YR_MPE.plugins.rf_plugin import RFComparisonPlugin

        plugin = RFComparisonPlugin(headless=True)
        plugin._out_dir = tempfile.mkdtemp()
        plugin._payload = {
            "labels": ["a", "b", "c"],
            "matrix": [[0, 1, 2], [1, 0, 3], [2, 3, 0]],
            "ntaxa": 10,
            "groups": [],
            "rfdist": "",
        }
        result = plugin._build_result_spec()
        self.assertEqual(result.validate(), [])
        item = result.datasets[0].items[0]
        self.assertEqual(item.item_type, "table")
        self.assertTrue(Path(item.file_path).is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
