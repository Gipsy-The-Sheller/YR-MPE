"""Emit a compact machine-readable summary on stdout.

Terminal analysis plugins (batch dN/dS, codon Z-test, RF comparison) do not
produce a dataset, but an external agent still needs their results. Each run
prints one JSON line prefixed with ``YMPE_SUMMARY`` so a human or an agent can
grep the process output without parsing the GUI tables.
"""

from __future__ import annotations

import json
import sys
from typing import Any

SUMMARY_PREFIX = "YMPE_SUMMARY "


def emit_summary(tool: str, **fields: Any) -> str:
    """Print one JSON summary line; return the line (without the prefix)."""
    record = {"schema": "yrmpe.summary.v1", "tool": tool}
    record.update(fields)
    line = json.dumps(record, ensure_ascii=False, default=str)
    if sys.stdout is not None:
        print(SUMMARY_PREFIX + line, flush=True)
    return line
