"""Gene concordance factor (gCF) plugin."""

from __future__ import annotations

from .concordance_plugin_base import ConcordancePluginBase


class GeneConcordancePlugin(ConcordancePluginBase):
    MODE = "gcf"
    TITLE = "Gene concordance factors (gCF)"
    NEEDS_GENE_TREES = True
    NEEDS_ALIGNMENT = False
    RUNS_ASTRAL = False

    def init_plugin_info(self):
        super().init_plugin_info()
        self.plugin_name = self.TITLE
        self.input_types = {"Newick": ["nwk", "newick", "tre", "tree"]}
        self.output_types = {"TSV table": ".tsv"}


class GeneConcordancePluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return GeneConcordancePlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
