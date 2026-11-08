"""Shared helpers for the alignment plugins.

The aligners differ in their CLI and parameter widgets, but their input
normalisation, output parsing and ResultSpec shape are the same. Keeping this
logic in one module makes the MAFFT/Clustal/MUSCLE/MACSE migration consistent.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from Bio import AlignIO, SeqIO
from PyQt5.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout

from ...runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec


def records_temp(plugin, records, label: str = "import") -> str:
    """Write SeqRecords to a plugin temp file and remember it for cleanup."""
    temp_file = plugin.create_temp_file(suffix="_%s.fas" % label)
    with open(temp_file, "w", encoding="utf-8") as handle:
        SeqIO.write(list(records), handle, "fasta")
    plugin.temp_files.append(temp_file)
    return temp_file


def items_to_files(plugin, items) -> list[str]:
    """Convert dataset-item-like objects to existing paths or temp FASTA."""
    files = []
    for item in items:
        path = getattr(item, "file_path", "") or ""
        if path and os.path.exists(path):
            files.append(path)
            continue
        records = list(getattr(item, "sequences", []) or [])
        if records:
            files.append(records_temp(plugin, records, "dataset_item"))
    return files


def prepare_import(plugin, import_data):
    """Return ``(files, batch_mode)`` for the common import payload shapes."""
    files = []
    batch_mode = False

    if isinstance(import_data, dict):
        items = import_data.get("dataset_items")
        if items is not None:
            files = items_to_files(plugin, items)
            batch_mode = len(files) > 1
        elif import_data.get("file_path"):
            files = [str(import_data["file_path"])]
        elif import_data.get("files"):
            files = [str(path) for path in import_data["files"]]
        elif import_data.get("sequences"):
            files = [records_temp(plugin, import_data["sequences"], "import")]
    elif isinstance(import_data, str):
        files = [import_data]
    elif isinstance(import_data, (list, tuple)):
        if import_data and isinstance(import_data[0], (list, tuple)):
            batch_mode = True
            for index, records in enumerate(import_data):
                files.append(records_temp(plugin, records, "batch_%d" % index))
        elif import_data and isinstance(import_data[0], dict):
            files = items_to_files(plugin, import_data)
            batch_mode = len(files) > 1
        elif import_data and hasattr(import_data[0], "seq"):
            files = [records_temp(plugin, import_data, "import")]

    return [path for path in files if path and os.path.exists(path)], batch_mode


def collect_partition_sequences(items):
    """Return ``(taxon_ids, records_by_item)`` across all partitions.

    A supermatrix must use the union of taxa from every partition. Using the
    first partition's taxon list silently drops taxa that are absent there but
    present in other loci (common with missing data).
    """
    records_by_item = []
    taxon_ids = []
    seen = set()
    for item in items:
        records = list(getattr(item, "sequences", []) or [])
        if not records:
            path = getattr(item, "file_path", "") or ""
            if path and os.path.exists(path):
                try:
                    records = list(SeqIO.parse(path, "fasta"))
                except Exception:
                    records = []
        records_by_item.append(records)
        for record in records:
            taxon = str(getattr(record, "id", "") or "")
            if taxon and taxon not in seen:
                seen.add(taxon)
                taxon_ids.append(taxon)
    return taxon_ids, records_by_item


def read_alignment(path: str):
    """Read FASTA or Clustal output and return ``(records, format)``."""
    records = list(SeqIO.parse(str(path), "fasta"))
    if records:
        return records, "fasta"
    try:
        alignment = AlignIO.read(str(path), "clustal")
    except Exception as exc:
        raise ValueError(
            "Cannot parse alignment %s as FASTA or Clustal: %s" % (path, exc)
        )
    return list(alignment), "clustal"


def build_result(tool: str, source: str, output_files, input_files,
                 parameters, mode: str, batch: bool,
                 suffix: str = "", dataset_name: Optional[str] = None):
    """Build the standard alignment ResultSpec and validate rectangularity."""
    if not output_files:
        raise ValueError("No alignment output is available.")

    items = []
    total_sequences = 0
    max_length = 0
    for output_file in output_files:
        records, file_format = read_alignment(output_file)
        if len(records) < 2:
            raise ValueError("%s contains fewer than two aligned sequences." % output_file)
        lengths = {len(str(record.seq)) for record in records}
        if len(lengths) > 1:
            raise ValueError("%s is not a rectangular alignment." % output_file)
        total_sequences += len(records)
        max_length = max(max_length, max(lengths))
        path = Path(output_file)
        loci_name = path.stem
        if suffix and loci_name.endswith(suffix):
            loci_name = loci_name[: -len(suffix)]
        items.append(DatasetItemSpec(
            loci_name=loci_name,
            item_type="alignment",
            sequences=records,
            file_path=str(path),
            format=file_format,
            is_aligned=True,
            length=max(lengths),
            sequence_count=len(records),
            data={"tool": tool, "parameters": list(parameters or [])},
            metadata={"source": tool},
        ))

    default_name = "%s alignment" % tool if len(items) == 1 else "%s alignments" % tool
    inputs = list(input_files or [])
    dataset = DatasetSpec(
        name=dataset_name or default_name,
        items=items,
        description="%d alignment(s), %d sequences" % (len(items), total_sequences),
        source=source,
        summary={
            "alignments": len(items),
            "sequences": total_sequences,
            "length": max_length,
        },
        metadata={
            "input": inputs[0] if len(inputs) == 1 else "",
            "input_files": inputs,
            "parameters": list(parameters or []),
            "mode": mode,
            "batch": bool(batch),
        },
    )
    return ResultSpec(source=source, datasets=[dataset])


def setup_preview(plugin, placeholder: str):
    """Install the lightweight FASTA/Clustal text preview on a plugin."""
    layout = QVBoxLayout(plugin.output_tab)
    layout.setContentsMargins(12, 12, 12, 12)
    plugin.summary_label = QLabel("No alignment yet.")
    plugin.summary_label.setWordWrap(True)
    plugin.summary_label.setStyleSheet("font-weight: 600;")
    layout.addWidget(plugin.summary_label)
    plugin.preview_text = QPlainTextEdit()
    plugin.preview_text.setReadOnly(True)
    plugin.preview_text.setPlaceholderText(placeholder)
    layout.addWidget(plugin.preview_text, 1)


def show_preview(plugin, output_files):
    if not output_files:
        plugin.summary_label.setText("No alignment output.")
        plugin.preview_text.clear()
        return
    first = Path(output_files[0])
    plugin.summary_label.setText(
        "%d alignment file(s) · %s" % (len(output_files), first.name)
    )
    try:
        text = first.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        plugin.preview_text.setPlainText("Cannot read %s: %s" % (first, exc))
        return
    if len(text) > 200_000:
        text = text[:200_000] + "\n... (preview truncated)"
    plugin.preview_text.setPlainText(text)
