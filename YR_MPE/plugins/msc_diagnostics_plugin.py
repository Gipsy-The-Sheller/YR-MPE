"""gCF + ASTRAL/wASTRAL branch diagnostics on a fixed species tree."""

from __future__ import annotations

from .concordance_plugin_base import ConcordancePluginBase


class MSCDiagnosticsPlugin(ConcordancePluginBase):
    MODE = "diagnostics"
    TITLE = "MSC branch diagnostics (gCF + ASTER CU)"
    NEEDS_GENE_TREES = True
    NEEDS_ALIGNMENT = False
    RUNS_ASTRAL = True

    def init_plugin_info(self):
        super().init_plugin_info()
        self.plugin_name = self.TITLE
        self.citation_keys = ["minh2020gcf", "sayyari2016astral", "zhang2022wastral"]
        self.input_types = {"Newick": ["nwk", "newick", "tre", "tree"]}
        self.output_types = {"TSV table": ".tsv"}


class MSCDiagnosticsPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return MSCDiagnosticsPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
