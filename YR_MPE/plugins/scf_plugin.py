"""Site concordance factor (sCF) plugin."""

from __future__ import annotations

from .concordance_plugin_base import ConcordancePluginBase


class SiteConcordancePlugin(ConcordancePluginBase):
    MODE = "scf"
    TITLE = "Site concordance factors (sCF)"
    NEEDS_GENE_TREES = False
    NEEDS_ALIGNMENT = True
    RUNS_ASTRAL = False

    def init_plugin_info(self):
        super().init_plugin_info()
        self.plugin_name = self.TITLE
        self.citation_keys = ["minh2020gcf", "mo2023scf"]
        self.input_types = {
            "Newick": ["nwk", "newick", "tre", "tree"],
            "Alignment": ["fas", "fasta", "phy", "nex", "nexus"],
        }
        self.output_types = {"TSV table": ".tsv"}


class SiteConcordancePluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return SiteConcordancePlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
