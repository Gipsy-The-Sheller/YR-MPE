"""sCF and CASTER-site q1/q2/q3 diagnostics on a fixed species tree.

This plugin is deliberately separate from the gene-tree based MSC
diagnostics. CASTER-site q1-q3 are site-based quartet scores and can carry
negative weights; sCF is sensitive to multiple hits. They are reported as
descriptive concordance diagnostics only, not as inputs to MSC parameter
estimation.
"""

from __future__ import annotations

from .concordance_plugin_base import ConcordancePluginBase


class SiteDiagnosticsPlugin(ConcordancePluginBase):
    MODE = "site_diagnostics"
    TITLE = "Site concordance diagnostics (sCF + CASTER q1-q3)"
    NEEDS_GENE_TREES = False
    NEEDS_ALIGNMENT = True
    RUNS_ASTRAL = False
    RUNS_CASTER = True

    def __init__(self, import_from=None, import_data=None, **kwargs):
        super().__init__(
            import_from=import_from, import_data=import_data, **kwargs
        )
        if hasattr(self, "run_button"):
            self.run_button.setText("Run sCF + CASTER")

    def init_plugin_info(self):
        super().init_plugin_info()
        self.plugin_name = self.TITLE
        self.citation_keys = [
            "minh2020gcf", "mo2023scf", "zhang2025aster", "zhang2025caster"
        ]
        self.input_types = {
            "Newick": ["nwk", "newick", "tre", "tree"],
            "Alignment": ["fas", "fasta", "phy", "nex", "nexus"],
        }
        self.output_types = {"TSV table": ".tsv"}


class SiteDiagnosticsPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return SiteDiagnosticsPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
