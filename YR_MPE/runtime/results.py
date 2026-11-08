"""Standard result objects exchanged between plugins and the host.

A GUI plugin that produces a dataset should build a :class:`ResultSpec` and
emit it through ``BasePlugin.export_result``. The host (YR-MPEA) converts the
spec into its own ``DatasetInfo``/``DatasetItem`` model. Plugins therefore do
not import the platform's dataset manager, and the host does not need a
bespoke signal/handler for every plugin.

The classes are stdlib-only; ``sequences`` may contain any object with a
``.seq`` attribute (Biopython ``SeqRecord`` in practice).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

RESULT_SCHEMA = "yrmpe.result.v1"

ITEM_TYPE_SEQUENCE = "sequence"
ITEM_TYPE_ALIGNMENT = "alignment"
ITEM_TYPE_MODEL = "model"
ITEM_TYPE_DISTANCE = "distance"
ITEM_TYPE_PHYLOGENY = "phylogeny"
ITEM_TYPE_VARIANT = "variant"
ITEM_TYPE_COALESCENT = "coalescent"
ITEM_TYPE_CLOCK = "clock"
ITEM_TYPE_CHAIN = "chain"
ITEM_TYPE_TREESET = "treeset"
ITEM_TYPE_TABLE = "table"

ALL_ITEM_TYPES = {
    ITEM_TYPE_SEQUENCE,
    ITEM_TYPE_ALIGNMENT,
    ITEM_TYPE_MODEL,
    ITEM_TYPE_DISTANCE,
    ITEM_TYPE_PHYLOGENY,
    ITEM_TYPE_VARIANT,
    ITEM_TYPE_COALESCENT,
    ITEM_TYPE_CLOCK,
    ITEM_TYPE_CHAIN,
    ITEM_TYPE_TREESET,
    ITEM_TYPE_TABLE,
}


def _jsonable(value: Any) -> Any:
    """Best-effort conversion for JSON persistence (paths, tuples, etc.)."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


@dataclass
class DatasetItemSpec:
    """One locus/result item inside a dataset."""

    loci_name: str
    item_type: str = ITEM_TYPE_SEQUENCE
    sequences: List[Any] = field(default_factory=list)
    file_path: str = ""
    format: str = "fasta"
    is_aligned: bool = False
    length: int = 0
    sequence_count: int = 0
    data: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    parent_index: Optional[int] = None
    selected: bool = False
    import_id: Optional[str] = None
    run_id: Optional[str] = None

    def __post_init__(self):
        if self.sequences and not self.sequence_count:
            self.sequence_count = len(self.sequences)
        if self.sequences and not self.length:
            try:
                self.length = max(len(str(record.seq)) for record in self.sequences)
            except Exception:
                self.length = 0

    def validate(self) -> List[str]:
        problems: List[str] = []
        if not self.loci_name:
            problems.append("item has no loci_name")
        if self.item_type not in ALL_ITEM_TYPES:
            problems.append("unknown item_type %r" % self.item_type)
        if not self.sequences and not self.file_path:
            problems.append("%s has neither sequences nor file_path" % (self.loci_name or "<unnamed>"))
        if self.sequences and self.sequence_count != len(self.sequences):
            problems.append(
                "%s declares %d sequences but carries %d"
                % (self.loci_name, self.sequence_count, len(self.sequences))
            )
        return problems

    def to_settings_dict(self) -> Dict[str, Any]:
        """Serializable form stored in DatasetInfo.settings['dataset_items']."""
        return {
            "loci_name": self.loci_name,
            "item_type": self.item_type,
            "file_path": self.file_path,
            "format": self.format,
            "is_aligned": self.is_aligned,
            "length": self.length,
            "sequence_count": self.sequence_count,
            "data": _jsonable(self.data),
            "metadata": _jsonable(self.metadata),
            "selected": self.selected,
            "import_id": self.import_id,
            "run_id": self.run_id,
        }


@dataclass
class DatasetSpec:
    """A dataset produced by a plugin."""

    name: str
    items: List[DatasetItemSpec] = field(default_factory=list)
    description: str = ""
    source: str = ""
    partition_scheme: Optional[str] = None
    summary: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    dataset_id: Optional[str] = None

    def validate(self) -> List[str]:
        problems: List[str] = []
        if not self.name:
            problems.append("dataset has no name")
        if not self.items:
            problems.append("dataset %r has no items" % self.name)
        for index, item in enumerate(self.items):
            for problem in item.validate():
                problems.append("item %d: %s" % (index, problem))
        return problems


@dataclass
class ResultSpec:
    """Top-level result emitted by a plugin."""

    source: str
    datasets: List[DatasetSpec] = field(default_factory=list)
    artifacts: List[str] = field(default_factory=list)
    summary: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    schema: str = RESULT_SCHEMA

    def validate(self) -> List[str]:
        problems: List[str] = []
        if self.schema != RESULT_SCHEMA:
            problems.append("unsupported result schema %r" % self.schema)
        if not self.source:
            problems.append("result has no source")
        if not self.datasets:
            problems.append("result has no datasets")
        for index, dataset in enumerate(self.datasets):
            for problem in dataset.validate():
                problems.append("dataset %d: %s" % (index, problem))
        return problems

    def to_settings_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "source": self.source,
            "summary": _jsonable(self.summary),
            "warnings": _jsonable(self.warnings),
            "artifacts": _jsonable(self.artifacts),
        }
