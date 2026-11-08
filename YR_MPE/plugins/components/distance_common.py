"""Shared distance-matrix parsing and ResultSpec construction."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from ...runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec


def parse_phylip(path: str):
    """Parse the PHYLIP-style matrix written by IQ-TREE (.mldist)."""
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if not lines:
        raise ValueError("Distance matrix is empty: %s" % path)
    first = lines[0].split()
    try:
        count = int(first[0]) if first else len(lines)
    except ValueError:
        count = len(lines)
    if count < 1 or len(lines) < count + 1:
        raise ValueError("Malformed distance matrix: %s" % path)
    labels = []
    matrix = []
    for line in lines[1:count + 1]:
        parts = line.split()
        if len(parts) < count + 1:
            raise ValueError("Malformed distance row in %s" % path)
        labels.append(parts[0])
        matrix.append([float(value) for value in parts[1:count + 1]])
    return labels, matrix


def build_result(tool: str, source: str, files: Sequence[str],
                 parameters=None, mode: str = ""):
    """Build a ResultSpec with one distance item per matrix file."""
    if not files:
        raise ValueError("No distance matrix output is available.")
    items = []
    total = 0
    for path in files:
        labels, matrix = parse_phylip(path)
        total = max(total, len(labels))
        items.append(DatasetItemSpec(
            loci_name=Path(path).stem,
            item_type="distance",
            file_path=str(path),
            format="phylip",
            data={"labels": labels, "matrix": matrix},
            metadata={"source": tool, "tool": tool},
        ))
    dataset = DatasetSpec(
        name="%s distance" % tool if len(items) == 1 else "%s distances" % tool,
        items=items,
        description="%d distance matrix(es)" % len(items),
        source=source,
        summary={"matrices": len(items), "sequences": total},
        metadata={
            "parameters": list(parameters or []),
            "mode": mode,
        },
    )
    return ResultSpec(
        source=source,
        datasets=[dataset],
        artifacts=[str(path) for path in files],
    )
