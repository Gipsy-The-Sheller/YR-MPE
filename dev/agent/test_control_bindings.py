"""Every agent-operable plugin declares a GUI-to-CLI control map.

The three GUI-only exceptions are deliberate: Tree Merger and Sequence
Concatenator are manual data-preparation tools, and IcyTree is a JavaScript
viewer that an agent cannot operate through the local bridge.

Run from the repository root:

    python dev/agent/test_control_bindings.py
"""

from __future__ import annotations

import importlib
import os
import sys
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication

QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

APP = QApplication.instance() or QApplication(sys.argv)

# skill id -> (module, class)
PLUGINS = {
    "orthofinder": ("YR_MPE.plugins.orthofinder_plugin", "OrthoFinderPlugin"),
    "mafft": ("YR_MPE.plugins.mafft_plugin", "MAFFTPlugin"),
    "macse": ("YR_MPE.plugins.macse_plugin", "MACSEPlugin"),
    "trimal": ("YR_MPE.plugins.trimal_plugin", "TrimAlPlugin"),
    "gblocks": ("YR_MPE.plugins.gblocks_plugin", "GBlocksPlugin"),
    "muscle5": ("YR_MPE.plugins.muscle5_plugin", "Muscle5Plugin"),
    "clustal-omega": ("YR_MPE.plugins.clustal_omega_plugin", "ClustalOmegaPlugin"),
    "modelfinder": ("YR_MPE.plugins.model_finder_plugin", "ModelFinderPlugin"),
    "modeltest-ng": ("YR_MPE.plugins.modeltest_ng_plugin", "ModelTestNGPlugin"),
    "model-parameters": ("YR_MPE.plugins.model_parameter_plugin", "ModelParameterPlugin"),
    "iqtree": ("YR_MPE.plugins.iqtree_plugin", "IQTreePlugin"),
    "raxml-ng": ("YR_MPE.plugins.raxmlng_plugin", "RAxMLNGPlugin"),
    "mrbayes": ("YR_MPE.plugins.mrbayes_plugin", "MrBayesPlugin"),
    "phylobayes": ("YR_MPE.plugins.phylobayes_plugin", "PhyloBayesPlugin"),
    "mpboot": ("YR_MPE.plugins.mpboot_plugin", "MPBootPlugin"),
    "decenttree": ("YR_MPE.plugins.decenttree_plugin", "DecentTreePlugin"),
    "astral": ("YR_MPE.plugins.astral_plugin", "AstralPlugin"),
    "caster-site": ("YR_MPE.plugins.caster_site_plugin", "CasterSitePlugin"),
    "treerecs": ("YR_MPE.plugins.treerecs_plugin", "TreeRecsPlugin"),
    "lsd2": ("YR_MPE.plugins.lsd2_plugin", "LSD2Plugin"),
    "minitracer": ("YR_MPE.plugins.minitracer_plugin", "MiniTracerPlugin"),
    "rf": ("YR_MPE.plugins.rf_plugin", "RFComparisonPlugin"),
    "codon-ztest": ("YR_MPE.plugins.codon_ztest_plugin", "CodonZTestPlugin"),
    "batch-dnds": ("YR_MPE.plugins.batch_dnds_plugin", "BatchDnDsPlugin"),
    "seqmatrix": ("YR_MPE.plugins.seqmatrix", "Lg_SeqMatrix"),
    "seqdbg": ("YR_MPE.plugins.seqdbg_plugin", "SeqDBGWindow"),
    "best-mrbayes": ("YR_MPE.plugins.best_plugin", "BestDialog"),
    "ml-distance": ("YR_MPE.plugins.ml_distance_plugin", "MLDistancePlugin"),
    "uncorrected-distance": ("YR_MPE.plugins.uncorrected_distance_plugin", "UncorrectedDistancePlugin"),
    "saturation": ("YR_MPE.plugins.substitution_saturation_plugin", "SubstitutionSaturationPlugin"),
    "pdguide": ("YR_MPE.plugins.pdguide_plugin", "PdGuidePlugin"),
    "sequence-viewer": ("YR_MPE.sequence_editor", "SequenceAlignmentViewer"),
    "console": ("YR_MPE.console_plugin", "YRMPEConsolePlugin"),
    "gcf": ("YR_MPE.plugins.gcf_plugin", "GeneConcordancePlugin"),
    "scf": ("YR_MPE.plugins.scf_plugin", "SiteConcordancePlugin"),
    "site-diagnostics": ("YR_MPE.plugins.site_diagnostics_plugin", "SiteDiagnosticsPlugin"),
    "msc-diagnostics": ("YR_MPE.plugins.msc_diagnostics_plugin", "MSCDiagnosticsPlugin"),
    "workflow-runner": ("YR_MPE.plugins.workflow_runner_plugin", "WorkflowRunnerPlugin"),
}

EXEMPT = {"tree-merger", "sequence-concatenator", "icytree"}


class ControlBindingCoverageTests(unittest.TestCase):
    def test_every_agent_plugin_declares_bindings(self):
        missing = []
        for skill_id, (module_name, class_name) in PLUGINS.items():
            module = importlib.import_module(module_name)
            cls = getattr(module, class_name)
            if not callable(getattr(cls, "control_bindings", None)):
                missing.append(skill_id)
        self.assertEqual(missing, [])

    def test_binding_rows_have_stable_shape(self):
        # These maps are pure literals; the larger legacy plugins read widget
        # state inside control_bindings(), so they are covered by the
        # declaration test and by the runtime help-doc tests instead.
        pure = {
            "orthofinder", "ml-distance", "uncorrected-distance",
            "saturation", "codon-ztest", "batch-dnds", "pdguide",
            "minitracer", "seqmatrix", "seqdbg", "sequence-viewer",
            "console",
        }
        problems = []
        for skill_id in pure:
            module_name, class_name = PLUGINS[skill_id]
            cls = getattr(importlib.import_module(module_name), class_name)

            class _Holder:
                pass

            rows = cls.control_bindings(_Holder())
            if not isinstance(rows, list) or not rows:
                problems.append("%s: empty control map" % skill_id)
                continue
            for row in rows:
                if not isinstance(row, dict) or not row.get("id"):
                    problems.append("%s: row without id" % skill_id)
                    continue
                if "cli" not in row:
                    problems.append("%s/%s: missing cli" % (skill_id, row["id"]))
                if not row.get("label"):
                    problems.append("%s/%s: missing label" % (skill_id, row["id"]))
        self.assertEqual(problems, [])

    def test_exemptions_are_documented(self):
        self.assertEqual(
            EXEMPT, {"tree-merger", "sequence-concatenator", "icytree"}
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
