"""Tests for the BEST agent surface and block round-trip."""

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

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.plugins.best_plugin import (  # noqa: E402
    BestConfig, BestDialog, build_mrbayes_block, parse_mrbayes_block,
)

APP = QApplication.instance() or QApplication(sys.argv)


class BestAgentTests(unittest.TestCase):
    def setUp(self):
        self.dialog = BestDialog(workdir=tempfile.mkdtemp(prefix="best-agent-"))

    def tearDown(self):
        self.dialog.close()

    def test_inspect_options_and_help(self):
        result = self.dialog.agent_invoke("best.inspect_options", {})
        self.assertTrue(result["ok"])
        self.assertIn("block", result["result"])
        doc = self.dialog.agent_invoke("best.help_doc", {})
        self.assertTrue(doc["ok"])
        self.assertEqual(doc["result"]["id"], "best-mrbayes")
        self.assertTrue(self.dialog.agent_contract()["actions"])

    def test_block_roundtrip_returns_canonical_and_unparsed(self):
        self.dialog.apply_options({
            "species_partition": "SpA: 1-2, SpB: 3-4",
            "ngen": 2500,
        })
        block = build_mrbayes_block(self.dialog._config())
        result = self.dialog.agent_invoke("best.apply_block", {"block": block})
        self.assertTrue(result["ok"])
        payload = result["result"]
        self.assertTrue(payload["roundtrip_equal"])
        self.assertEqual(payload["unparsed_lines"], [])
        self.assertIn("ngen=2500", payload["canonical_block"])

    def test_unparsed_lines_are_reported(self):
        block = build_mrbayes_block(BestConfig(ngen=1000))
        block = block.replace(
            "    mcmc ngen=1000",
            "    prset custompr=custom(1);\n    mcmc ngen=1000",
        )
        result = self.dialog.agent_invoke("best.apply_block", {"block": block})
        self.assertTrue(result["ok"])
        self.assertTrue(any("custompr" in line for line in
                            result["result"]["unparsed_lines"]))

    def test_modify_options_updates_form(self):
        result = self.dialog.agent_invoke(
            "best.modify_options", {"options": {"ngen": 1234, "ploidy": "haploid"}}
        )
        self.assertTrue(result["ok"])
        self.assertEqual(self.dialog.ngen_spin.value(), 1234)
        self.assertEqual(self.dialog.ploidy_combo.currentText(), "haploid")


if __name__ == "__main__":
    unittest.main(verbosity=2)
