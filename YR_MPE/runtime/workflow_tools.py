"""Map workflow step ops to plugin entry classes for configuration windows."""

from __future__ import annotations

import importlib
from typing import Optional, Tuple

ENTRIES = {
    "mafft": ("YR_MPE.plugins.mafft_plugin", "MAFFTPluginEntry"),
    "macse": ("YR_MPE.plugins.macse_plugin", "MACSEPluginEntry"),
    "trimal": ("YR_MPE.plugins.trimal_plugin", "TrimAlPluginEntry"),
    "gblocks": ("YR_MPE.plugins.gblocks_plugin", "GBlocksPluginEntry"),
    "muscle5": ("YR_MPE.plugins.muscle5_plugin", "Muscle5PluginEntry"),
    "clustal_omega": ("YR_MPE.plugins.clustal_omega_plugin", "ClustalOmegaPluginEntry"),
    "modelfinder": ("YR_MPE.plugins.model_finder_plugin", "ModelFinderPluginEntry"),
    "model_finder": ("YR_MPE.plugins.model_finder_plugin", "ModelFinderPluginEntry"),
    "modeltest": ("YR_MPE.plugins.modeltest_ng_plugin", "ModelTestNGPluginEntry"),
    "modeltest_ng": ("YR_MPE.plugins.modeltest_ng_plugin", "ModelTestNGPluginEntry"),
    "iqtree": ("YR_MPE.plugins.iqtree_plugin", "IQTreePluginEntry"),
    "raxmlng": ("YR_MPE.plugins.raxmlng_plugin", "RAxMLNGPluginEntry"),
    "raxml_ng": ("YR_MPE.plugins.raxmlng_plugin", "RAxMLNGPluginEntry"),
    "mrbayes": ("YR_MPE.plugins.mrbayes_plugin", "MrBayesPluginEntry"),
    "phylobayes": ("YR_MPE.plugins.phylobayes_plugin", "PhyloBayesPluginEntry"),
    "mpboot": ("YR_MPE.plugins.mpboot_plugin", "MPBootPluginEntry"),
    "decenttree": ("YR_MPE.plugins.decenttree_plugin", "DecentTreePluginEntry"),
    "astral": ("YR_MPE.plugins.astral_plugin", "AstralPluginEntry"),
    "caster": ("YR_MPE.plugins.caster_site_plugin", "CasterSitePluginEntry"),
    "caster_site": ("YR_MPE.plugins.caster_site_plugin", "CasterSitePluginEntry"),
    "treerecs": ("YR_MPE.plugins.treerecs_plugin", "TreeRecsPluginEntry"),
    "lsd2": ("YR_MPE.plugins.lsd2_plugin", "LSD2PluginEntry"),
    "rf": ("YR_MPE.plugins.rf_plugin", "RFComparisonPluginEntry"),
    "model_parameter": ("YR_MPE.plugins.model_parameter_plugin", "ModelParameterPluginEntry"),
}


def normalize_op(op: str) -> str:
    token = str(op or "").strip().lower()
    token = token.split(":", 1)[0].split("/", 1)[0]
    if "." in token:
        token = token.split(".", 1)[0]
    return token.replace("-", "_")


def entry_for_op(op: str) -> Optional[Tuple[str, str]]:
    token = normalize_op(op)
    if token in ENTRIES:
        return ENTRIES[token]
    for key, value in ENTRIES.items():
        if token.startswith(key + "_"):
            return value
    return None


def create_plugin(op: str, workdir=None):
    entry = entry_for_op(op)
    if entry is None:
        raise ValueError("no plugin entry for workflow op %r" % op)
    module_name, class_name = entry
    module = importlib.import_module(module_name)
    plugin_class = getattr(module, class_name)
    instance = plugin_class()
    if hasattr(instance, "run"):
        return instance.run(workdir=workdir)
    return instance
