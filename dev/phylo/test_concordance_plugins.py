"""Tests for gCF/sCF/MSC-diagnostics helpers and table results."""

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

from YR_MPE.plugins import concordance_common as common  # noqa: E402
from YR_MPE.plugins.gcf_plugin import GeneConcordancePlugin  # noqa: E402
from YR_MPE.plugins.msc_diagnostics_plugin import MSCDiagnosticsPlugin  # noqa: E402
from YR_MPE.plugins.scf_plugin import SiteConcordancePlugin  # noqa: E402
from YR_MPE.plugins.site_diagnostics_plugin import SiteDiagnosticsPlugin  # noqa: E402
from YR_MPE.runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec  # noqa: E402

APP = QApplication.instance() or QApplication(sys.argv)

# Repo config.json stores tool paths relative to YR_MPE/ with a leading
# slash; on Windows os.path.isabs() treats that as absolute, so the base
# plugin would block on its configuration-guide dialog. These tests draw
# widgets only and never execute a tool.
from YR_MPE.templates.base_plugin_ui import BasePlugin  # noqa: E402

BasePlugin.show_config_guide = lambda self: None


class ConcordanceCommonTests(unittest.TestCase):
    def test_smaller_side_is_canonical(self):
        taxa = {"A", "B", "C", "D", "E"}
        self.assertEqual(common.canonical_split({"A", "B"}, taxa), ("A", "B"))
        self.assertEqual(
            common.canonical_split({"C", "D", "E"}, taxa), ("A", "B")
        )
        self.assertEqual(common.canonical_split({"D", "E"}, taxa), ("D", "E"))

    def test_x_from_gcf(self):
        self.assertEqual(common.x_from_gcf("100"), "")
        self.assertEqual(common.x_from_gcf("33.333"), "")
        self.assertAlmostEqual(float(common.x_from_gcf("60")), 0.510826, places=5)

    def test_wastral_label_check(self):
        root = Path(tempfile.mkdtemp())
        good = root / "good.trees"
        good.write_text("((A:1,B:1)0.9:1,(C:1,(D:1,E:1)0.8:1)0.9:1);\n")
        ok, reason, detail = common.check_wastral_labels(str(good))
        self.assertTrue(ok, reason)
        self.assertEqual(detail["missing"], 0)

        bad = root / "bad.trees"
        bad.write_text("((A:1,B:1)/1/100:1,(C:1,(D:1,E:1)0.8:1)0.9:1);\n")
        ok, reason, detail = common.check_wastral_labels(str(bad))
        self.assertFalse(ok)
        self.assertIn("composite", reason)
        self.assertGreater(detail["composite"], 0)

        missing = root / "missing.trees"
        missing.write_text("((A:1,B:1):1,(C:1,(D:1,E:1):1):1);\n")
        ok, reason, detail = common.check_wastral_labels(str(missing))
        self.assertFalse(ok)
        self.assertIn("no numeric label", reason)

    def test_cf_parser_groups_root_sides_once(self):
        root = Path(tempfile.mkdtemp())
        stat = root / "x.cf.stat"
        branch = root / "x.cf.branch"
        stat.write_text(
            "ID\tgCF\tgCF_N\tgDF1\tgDF1_N\tgDF2\tgDF2_N\tgDFP\tgDFP_N\tgN\tLabel\tLength\n"
            "7\t60\t3\t0\t0\t0\t0\t40\t2\t5\t\t-1\n"
            "8\t60\t3\t0\t0\t0\t0\t40\t2\t5\t\t-1\n"
            "9\t80\t4\t0\t0\t20\t1\t0\t0\t5\t\t-1\n"
        )
        branch.write_text("((A,B)7,(C,(D,E)9)8)6:0.0;\n")
        _, rows, _ = common.cf_rows(str(root / "x"), "gcf")
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["split"], "A,B")

    def test_astral_annotation_parser(self):
        root = Path(tempfile.mkdtemp())
        tree = root / "a.tre"
        tree.write_text(
            "((A,B)'[CULength=0.5;q1=0.6;q2=0.2;q3=0.2;localPP=0.9]',C,D);\n"
        )
        rows = common.parse_astral_tree(str(tree))
        self.assertTrue(rows)
        self.assertEqual(next(iter(rows.values()))["CULength"], "0.5")

    def test_caster_nan_length_keeps_quartet_labels(self):
        root = Path(tempfile.mkdtemp())
        tree = root / "caster.tre"
        tree.write_text(
            "((A,(B,C)'[CULength=nan;q1=0.6;q2=-0.1;q3=0.5]':nan),D);\n"
        )
        rows = common.parse_astral_tree(str(tree))
        self.assertTrue(rows)
        annotations = next(iter(rows.values()))
        self.assertEqual(annotations["q1"], "0.6")
        self.assertEqual(annotations["q2"], "-0.1")

    def test_tree_taxa_keeps_full_labels(self):
        root = Path(tempfile.mkdtemp())
        tree = root / "tree.nwk"
        tree.write_text(
            "(('Taxon_A_very_long_name':1,Taxon_B:1):1,Taxon_C:1);\n"
        )
        self.assertEqual(
            common.tree_taxa(str(tree)),
            ["Taxon_A_very_long_name", "Taxon_B", "Taxon_C"],
        )

    def test_site_rows_merge_and_negative_flag(self):
        columns = common.columns_for("site_diagnostics")
        cf_rows = [
            {
                "split": "A,B", "split_size": "2", "branch_id": "5",
                "sCF": "70", "sDF1": "20", "sDF2": "10", "sN": "10",
            }
        ]
        caster_rows = {
            ("A", "B"): {"q1": "0.6", "q2": "-0.1", "q3": "0.5",
                         "CULength": "0.8"}
        }
        cols, rows, summary = common.merge_site_diagnostics_rows(
            columns, cf_rows, caster_rows
        )
        self.assertEqual(cols, columns)
        self.assertEqual(rows[0]["sCF"], "70")
        self.assertEqual(rows[0]["q2_caster"], "-0.1")
        self.assertEqual(rows[0]["CU_caster"], "0.8")
        self.assertEqual(rows[0]["q_negative"], "yes")
        self.assertEqual(summary["caster_branches"], 1)

    def test_table_result_contract(self):
        root = Path(tempfile.mkdtemp())
        table = root / "cf.tsv"
        table.write_text("split\tgCF\nA,B\t60\n")
        spec = DatasetItemSpec(
            loci_name="cf_table",
            item_type="table",
            file_path=str(table),
            format="tsv",
            data={"columns": ["split", "gCF"], "rows": [["A,B", "60"]]},
        )
        result = ResultSpec(
            source="test",
            datasets=[DatasetSpec(name="cf", items=[spec], source="test")],
        )
        self.assertEqual(result.validate(), [])
        self.assertEqual(spec.item_type, "table")


class ConcordancePluginTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def _plugin(self, cls):
        return cls(workdir=str(self.root))

    def test_plugin_metadata(self):
        gcf = self._plugin(GeneConcordancePlugin)
        scf = self._plugin(SiteConcordancePlugin)
        diag = self._plugin(MSCDiagnosticsPlugin)
        try:
            self.assertEqual(gcf.MODE, "gcf")
            self.assertEqual(scf.MODE, "scf")
            self.assertEqual(diag.MODE, "diagnostics")
            self.assertIn("minh2020gcf", gcf.citation_keys)
            self.assertIn("mo2023scf", scf.citation_keys)
            self.assertIn("sayyari2016astral", diag.citation_keys)
        finally:
            for plugin in (gcf, scf, diag):
                plugin.close()

    def test_diagnostics_preflight_label(self):
        plugin = self._plugin(MSCDiagnosticsPlugin)
        try:
            good = self.root / "good.trees"
            good.write_text("((A:1,B:1)0.9:1,(C:1,(D:1,E:1)0.8:1)0.9:1);\n")
            plugin.gene_trees_edit.setText(str(good))
            plugin._update_wastral_check()
            self.assertIn("wASTRAL is enabled", plugin.wastral_check_label.text())

            bad = self.root / "bad.trees"
            bad.write_text("((A:1,B:1)/1/100:1,(C:1,(D:1,E:1)0.8:1)0.9:1);\n")
            plugin.gene_trees_edit.setText(str(bad))
            plugin._update_wastral_check()
            self.assertIn("will be rejected", plugin.wastral_check_label.text())
        finally:
            plugin.close()

    def test_diagnostics_overview_has_bar_panels(self):
        plugin = self._plugin(MSCDiagnosticsPlugin)
        try:
            plugin._draw_overview([
                {"split": "A,B", "gCF": "71.43", "gDF1": "0",
                 "gDF2": "0", "gDFP": "28.57", "CU_astral": "0.86",
                 "CU_wastral": "0.75", "x_gcf": "0.85",
                 "qCF": "0.82", "qDF1": "0.14", "qDF2": "0.04"},
                {"split": "A,B,C", "gCF": "30", "gDF1": "50",
                 "gDF2": "20", "gDFP": "0", "CU_astral": "0.98",
                 "CU_wastral": "0.48", "x_gcf": "1.54",
                 "qCF": "0.30", "qDF1": "0.50", "qDF2": "0.20"},
                {"split": "A,B,C,D", "gCF": "60", "gDF1": "20",
                 "gDF2": "20", "gDFP": "0", "CU_astral": "0.55",
                 "CU_wastral": "0.50", "x_gcf": "0.51",
                 "qCF": "0.30", "qDF1": "0.40", "qDF2": "0.30"},
                {"split": "B,C,D", "gCF": "20", "gDF1": "10",
                 "gDF2": "10", "gDFP": "60", "CU_astral": "0.20",
                 "CU_wastral": "0.18", "x_gcf": "",
                 "qCF": "0.20", "qDF1": "0.10", "qDF2": "0.10"},
            ])
            self.assertEqual(len(plugin.overview_figure.axes), 2)
            gcf_axis, qcf_axis = plugin.overview_figure.axes
            self.assertEqual(
                [patch.get_height() for patch in gcf_axis.patches], [3, 1]
            )
            self.assertEqual(
                [patch.get_height() for patch in qcf_axis.patches], [2, 2]
            )
        finally:
            plugin.close()

    def test_bipartition_info_matrix_uses_full_taxa(self):
        plugin = self._plugin(GeneConcordancePlugin)
        try:
            payload = {
                "table_path": str(self.root / "x.tsv"),
                "columns": ["split", "gCF"],
                "rows": [
                    {"split": "Taxon_A_very_long_name,Taxon_B", "gCF": "70"},
                    {"split": "Taxon_B,Taxon_C", "gCF": "40"},
                ],
                "taxa": [
                    "Taxon_A_very_long_name", "Taxon_B", "Taxon_C",
                    "Taxon_D_not_in_any_split",
                ],
                "summary": {"branches": 2},
            }
            plugin._display_payload(payload)
            table = plugin.bipartition_table
            self.assertEqual(plugin.view_tabs.count(), 4)
            self.assertEqual(table.rowCount(), 4)
            self.assertEqual(table.columnCount(), 3)
            self.assertEqual(table.item(0, 0).text(), "Taxon_A_very_long_name")
            self.assertEqual(table.horizontalHeaderItem(1).text(), "S1")
            self.assertEqual(
                table.item(0, 1).background().color().name(), "#4e699a"
            )
            self.assertNotEqual(
                table.item(3, 1).background().color().name(), "#4e699a"
            )
            self.assertIn("S1 (2 taxa)", plugin.bipartition_key.text())
            self.assertFalse(plugin.save_graph_button.isEnabled())
            plugin.view_tabs.setCurrentIndex(2)
            self.assertTrue(plugin.save_graph_button.isEnabled())
        finally:
            plugin.close()

    def test_site_diagnostics_metadata(self):
        plugin = self._plugin(SiteDiagnosticsPlugin)
        try:
            self.assertEqual(plugin.MODE, "site_diagnostics")
            self.assertTrue(plugin.RUNS_CASTER)
            self.assertFalse(plugin.NEEDS_GENE_TREES)
            self.assertTrue(plugin.NEEDS_ALIGNMENT)
            self.assertIn("mo2023scf", plugin.citation_keys)
            self.assertIn("zhang2025caster", plugin.citation_keys)
        finally:
            plugin.close()

    def test_site_overview_has_bar_panels(self):
        plugin = self._plugin(SiteDiagnosticsPlugin)
        try:
            plugin._draw_overview([
                {"split": "A,B", "sCF": "70", "sDF1": "20", "sDF2": "10",
                 "q1_caster": "0.6", "q2_caster": "-0.1", "q3_caster": "0.5"},
                {"split": "A,B,C", "sCF": "30", "sDF1": "50", "sDF2": "20",
                 "q1_caster": "0.2", "q2_caster": "0.7", "q3_caster": "0.1"},
            ])
            self.assertEqual(len(plugin.overview_figure.axes), 2)
            scf_axis, qcf_axis = plugin.overview_figure.axes
            self.assertEqual(
                [patch.get_height() for patch in scf_axis.patches], [1, 1]
            )
            self.assertEqual(
                [patch.get_height() for patch in qcf_axis.patches], [1, 1]
            )
        finally:
            plugin.close()

    def test_simplex_panels_for_site_and_gcf(self):
        site = self._plugin(SiteDiagnosticsPlugin)
        try:
            site._draw_simplex([
                {"split": "A,B", "sCF": "70", "sDF1": "20", "sDF2": "10",
                 "q1_caster": "0.6", "q2_caster": "0.2", "q3_caster": "0.2",
                 "CU_caster": "0.8"},
                {"split": "A,B,C", "sCF": "30", "sDF1": "50", "sDF2": "20",
                 "q1_caster": "0.2", "q2_caster": "0.7", "q3_caster": "0.1",
                 "CU_caster": "0.4"},
            ])
            # Two independent figures: sCF and CASTER qCF.
            self.assertEqual(set(site.simplex_views), {"scf", "caster"})
            self.assertEqual(len(site.simplex_views["scf"][0].axes), 1)
            # CASTER figure carries a CULength colourbar.
            self.assertEqual(len(site.simplex_views["caster"][0].axes), 2)
        finally:
            site.close()

        gcf = self._plugin(GeneConcordancePlugin)
        try:
            gcf._draw_simplex([
                {"split": "A,B", "gCF": "70", "gDF1": "20", "gDF2": "10",
                 "gDFP": "0"},
            ])
            self.assertEqual(len(gcf.simplex_figure.axes), 1)
        finally:
            gcf.close()

        diag = self._plugin(MSCDiagnosticsPlugin)
        try:
            diag._draw_simplex([
                {"split": "A,B", "gCF": "70", "gDF1": "20", "gDF2": "10",
                 "gDFP": "10", "qCF": "0.6", "qDF1": "0.25",
                 "qDF2": "0.15", "CU_astral": "0.9"},
            ])
            # One Simplex tab with two panels; each panel has its own
            # colourbar (gDFP and CULength).
            self.assertEqual(set(diag.simplex_views), {"combined"})
            self.assertEqual(len(diag.simplex_figure.axes), 4)
        finally:
            diag.close()

    def test_table_result_spec(self):
        plugin = self._plugin(GeneConcordancePlugin)
        try:
            table = self.root / "table.tsv"
            table.write_text("split\tgCF\nA,B\t60\n")
            plugin._payload = {
                "table_path": str(table),
                "columns": ["split", "gCF"],
                "rows": [{"split": "A,B", "gCF": "60"}],
                "summary": {"branches": 1},
                "artifacts": [str(table)],
            }
            result = plugin._build_result_spec()
            self.assertEqual(result.validate(), [])
            self.assertEqual(result.datasets[0].items[0].item_type, "table")
        finally:
            plugin.close()

    def test_site_table_result_spec(self):
        plugin = self._plugin(SiteDiagnosticsPlugin)
        try:
            table = self.root / "site.tsv"
            table.write_text("split\tsCF\tq1_caster\nA,B\t70\t0.6\n")
            plugin._payload = {
                "table_path": str(table),
                "columns": ["split", "sCF", "q1_caster"],
                "rows": [{"split": "A,B", "sCF": "70", "q1_caster": "0.6"}],
                "summary": {"branches": 1},
                "artifacts": [str(table)],
            }
            result = plugin._build_result_spec()
            self.assertEqual(result.validate(), [])
            self.assertEqual(result.datasets[0].items[0].item_type, "table")
        finally:
            plugin.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
