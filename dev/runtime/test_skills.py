"""Tests for skill documents and GUI-to-CLI control bindings."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QSpinBox, QWidget

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.runtime.skills import (  # noqa: E402
    SkillLibrary, _resolve_widget, build_plugin_help_doc,
)

APP = QApplication.instance() or QApplication(sys.argv)


class ManifestDocumentTests(unittest.TestCase):
    def test_every_catalog_entry_has_a_document(self):
        library = SkillLibrary()
        missing = [entry["id"] for entry in library.entries()
                   if not entry.get("file")
                   or not (library.root / entry["file"]).is_file()]
        self.assertEqual(missing, [])


class DottedWidgetTests(unittest.TestCase):
    def test_resolves_nested_attribute_path(self):
        class Panel(QWidget):
            def __init__(self):
                super().__init__()
                self.spin = QSpinBox()

        class Plugin(QWidget):
            def __init__(self):
                super().__init__()
                self.panel = Panel()

        plugin = Plugin()
        self.assertIs(_resolve_widget(plugin, "panel.spin"), plugin.panel.spin)
        self.assertIsNone(_resolve_widget(plugin, "panel.missing.spin"))


class LiveSkillDocumentTests(unittest.TestCase):
    def _check(self, plugin, skill_id, expected_flag):
        try:
            doc = build_plugin_help_doc(plugin)
            self.assertEqual(doc["id"], skill_id)
            self.assertTrue(doc["markdown"])
            self.assertIn(expected_flag, doc["markdown"])
            missing = [row["id"] for row in doc["controls"]
                       if not row["widget_found"]]
            self.assertEqual(missing, [])
        finally:
            plugin.close()

    def test_astral_advanced_controls_resolve(self):
        from YR_MPE.plugins.astral_plugin import AstralPlugin

        plugin = AstralPlugin()
        doc = build_plugin_help_doc(plugin)
        controls = {row["id"]: row for row in doc["controls"]}
        self.assertTrue(controls["rounds"]["widget_found"])
        self.assertEqual(controls["version"]["cli_args"], ["astral4"])
        self.assertTrue(controls["support"]["cli_args"][0] == "-u")
        plugin.close()

    def test_caster_site_controls_resolve(self):
        from YR_MPE.plugins.caster_site_plugin import CasterSitePlugin

        self._check(CasterSitePlugin(), "caster-site", "--support")

    def test_model_parameter_controls_resolve(self):
        from YR_MPE.plugins.model_parameter_plugin import ModelParameterPlugin

        self._check(ModelParameterPlugin(), "model-parameters", "-m MODEL")

    def test_rf_controls_resolve(self):
        from YR_MPE.plugins.rf_plugin import RFComparisonPlugin

        self._check(RFComparisonPlugin(), "rf", "--tree-dist-all")

    def test_extended_skill_set_has_no_broken_bindings(self):
        from YR_MPE.plugins.clustal_omega_plugin import ClustalOmegaPlugin
        from YR_MPE.plugins.decenttree_plugin import DecentTreePlugin
        from YR_MPE.plugins.gblocks_plugin import GBlocksPlugin
        from YR_MPE.plugins.macse_plugin import MACSEPlugin
        from YR_MPE.plugins.mpboot_plugin import MPBootPlugin
        from YR_MPE.plugins.mrbayes_plugin import MrBayesPlugin
        from YR_MPE.plugins.muscle5_plugin import Muscle5Plugin
        from YR_MPE.plugins.phylobayes_plugin import PhyloBayesPlugin
        from YR_MPE.plugins.treerecs_plugin import TreeRecsPlugin
        from YR_MPE.plugins.trimal_plugin import TrimAlPlugin

        cases = [
            (ClustalOmegaPlugin, "clustal-omega"),
            (DecentTreePlugin, "decenttree"),
            (GBlocksPlugin, "gblocks"),
            (MACSEPlugin, "macse"),
            (MPBootPlugin, "mpboot"),
            (MrBayesPlugin, "mrbayes"),
            (Muscle5Plugin, "muscle5"),
            (PhyloBayesPlugin, "phylobayes"),
            (TreeRecsPlugin, "treerecs"),
            (TrimAlPlugin, "trimal"),
        ]
        for plugin_class, skill_id in cases:
            plugin = plugin_class()
            try:
                doc = build_plugin_help_doc(plugin)
                self.assertEqual(doc["id"], skill_id)
                self.assertTrue(doc["markdown"])
                missing = [row["id"] for row in doc["controls"]
                           if not row["widget_found"]]
                self.assertEqual(missing, [], skill_id)
            finally:
                plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
