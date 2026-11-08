"""Shared helpers for gCF, sCF and ASTER concordance diagnostics.

The three plugins use the same canonical split representation: each bipartition
is named by the smaller side of the split (ties are resolved lexicographically),
so a table row is unambiguous even when the two sides are written in a
different order by different programs.
"""

from __future__ import annotations

import json
import math
import os
import re
from io import StringIO
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from Bio import Phylo

Split = Tuple[str, ...]
Row = Dict[str, str]


# ----------------------------------------------------------------------
# configuration and small utilities
# ----------------------------------------------------------------------
def resolve_tool(plugin_path: str, name: str) -> str:
    """Resolve a tool path from YR_MPE/config.json."""
    config_path = Path(plugin_path) / "config.json"
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    wanted = name.lower()
    for record in payload:
        if str(record.get("name", "")).lower() != wanted:
            continue
        raw = str(record.get("path", ""))
        if not raw:
            return ""
        path = Path(raw)
        if not path.is_absolute():
            path = Path(plugin_path) / raw.lstrip("/\\")
        return str(path.resolve()) if path.is_file() else ""
    return ""


def newick_lines(path: str) -> List[str]:
    return [
        line.strip() for line in Path(path).read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


_NONFINITE_BRANCH_LENGTH = re.compile(
    r":\s*(?:nan|[-+]?inf)(?=[,);])", re.IGNORECASE
)


def parse_newick(text: str):
    """Parse Newick, tolerating ASTER/CASTER ``:nan`` branch lengths.

    Biopython drops an internal node label when the following branch length
    is ``nan`` and stores ``nan`` in ``clade.confidence`` instead. The branch
    length is irrelevant to split/support parsing here, so it is replaced by
    zero before parsing.
    """
    text = _NONFINITE_BRANCH_LENGTH.sub(":0", text)
    return Phylo.read(StringIO(text), "newick")


def tree_taxa(path: str) -> List[str]:
    """Return all tip labels of a Newick tree, sorted for stable display."""
    tree = parse_newick(Path(path).read_text(encoding="utf-8", errors="replace"))
    return sorted({tip.name for tip in tree.get_terminals() if tip.name})


def canonical_split(side: Iterable[str], all_taxa: Iterable[str]) -> Split:
    """Return the smaller side of a bipartition, with a stable tie-break."""
    left = set(side)
    right = set(all_taxa) - left
    if not left or not right:
        return ()
    left_key = (len(left), tuple(sorted(left)))
    right_key = (len(right), tuple(sorted(right)))
    return tuple(sorted(left if left_key <= right_key else right))


def split_label(split: Split) -> str:
    return ",".join(split)


def _tree_splits_with_annotations(tree) -> Dict[Split, Dict[str, Any]]:
    taxa = {tip.name for tip in tree.get_terminals() if tip.name}
    out: Dict[Split, Dict[str, Any]] = {}
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        side = {tip.name for tip in clade.get_terminals() if tip.name}
        if not (1 < len(side) < len(taxa) - 1):
            continue
        split = canonical_split(side, taxa)
        if not split:
            continue
        label = clade.name
        if label is None and clade.confidence is not None:
            label = str(clade.confidence)
        out[split] = {
            "label": str(label) if label is not None else "",
            "branch_length": clade.branch_length,
        }
    return out


def _node_label(clade) -> str:
    value = clade.name
    if value is None and clade.confidence is not None:
        value = clade.confidence
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def parse_cf_branch(path: str) -> Dict[str, Dict[str, Any]]:
    """Map IQ-TREE branch IDs (``.cf.branch`` labels) to canonical splits."""
    tree = parse_newick(Path(path).read_text(encoding="utf-8", errors="replace"))
    taxa = {tip.name for tip in tree.get_terminals() if tip.name}
    out: Dict[str, Dict[str, Any]] = {}
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        label = _node_label(clade)
        side = {tip.name for tip in clade.get_terminals() if tip.name}
        if not label or not (1 < len(side) < len(taxa) - 1):
            continue
        split = canonical_split(side, taxa)
        if split:
            out[label] = {
                "split": split,
                "branch_length": clade.branch_length,
            }
    return out


def parse_cf_stat(path: str) -> List[Dict[str, Any]]:
    """Parse IQ-TREE's ``.cf.stat`` table (comment header tolerated)."""
    lines = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        lines.append(line)
    if not lines:
        raise ValueError("Concordance-factor statistics file is empty.")
    header = [part.strip() for part in lines[0].split("\t")]
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        values = {
            header[index]: (parts[index] if index < len(parts) else "")
            for index in range(len(header))
        }
        for key, value in list(values.items()):
            if key in {"gCF", "gDF1", "gDF2", "gDFP", "sCF", "sDF1", "sDF2"}:
                values[key] = _number(value)
            elif key in {"gCF_N", "gDF1_N", "gDF2_N", "gDFP_N", "gN",
                         "sCF_N", "sDF1_N", "sDF2_N", "sN"}:
                values[key] = _integer(value)
            elif key == "Length":
                values[key] = _number(value)
        rows.append(values)
    return rows


def _number(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return ""
    try:
        return "%.6g" % float(text)
    except ValueError:
        return text


def _integer(value: Any) -> str:
    text = str(value).strip()
    try:
        return str(int(float(text)))
    except ValueError:
        return text


def cf_rows(prefix: str, mode: str) -> tuple[List[str], List[Row], Dict[str, Any]]:
    """Build gCF/sCF rows from an IQ-TREE prefix."""
    stat_path = str(prefix) + ".cf.stat"
    branch_path = str(prefix) + ".cf.branch"
    if not Path(stat_path).is_file():
        raise ValueError("IQ-TREE did not write %s" % stat_path)
    if not Path(branch_path).is_file():
        raise ValueError("IQ-TREE did not write %s" % branch_path)
    branch_map = parse_cf_branch(branch_path)
    stat_rows = parse_cf_stat(stat_path)
    rows: List[Row] = []
    seen_splits = set()
    for stat in stat_rows:
        branch_id = str(stat.get("ID", "")).strip()
        branch = branch_map.get(branch_id)
        if not branch:
            continue
        split = branch["split"]
        if split in seen_splits:
            # A rooted Newick tree exposes both sides of the root edge as
            # separate branch IDs; in the unrooted table they are one split.
            continue
        seen_splits.add(split)
        row: Row = {
            "split": split_label(split),
            "split_size": str(len(split)),
            "branch_id": branch_id,
            "branch_length": str(stat.get("Length", "")),
        }
        if mode == "gcf":
            row.update({
                "gCF": stat.get("gCF", ""),
                "gDF1": stat.get("gDF1", ""),
                "gDF2": stat.get("gDF2", ""),
                "gDFP": stat.get("gDFP", ""),
                "gN": stat.get("gN", ""),
            })
        else:
            row.update({
                "sCF": stat.get("sCF", ""),
                "sDF1": stat.get("sDF1", ""),
                "sDF2": stat.get("sDF2", ""),
                "sN": stat.get("sN", ""),
            })
        rows.append(row)
    rows.sort(key=lambda row: (int(row["split_size"]), row["split"]))
    summary = {
        "branches": len(rows),
        "mode": mode,
        "source": "IQ-TREE 3",
    }
    return columns_for(mode), rows, summary


def columns_for(mode: str) -> List[str]:
    if mode == "gcf":
        return ["split", "split_size", "branch_id", "gCF", "gDF1", "gDF2",
                "gDFP", "gN", "branch_length"]
    if mode == "scf":
        return ["split", "split_size", "branch_id", "sCF", "sDF1", "sDF2",
                "sN", "branch_length"]
    if mode == "site_diagnostics":
        return ["split", "split_size", "branch_id", "sCF", "sDF1", "sDF2",
                "sN", "q1_caster", "q2_caster", "q3_caster", "q_negative"]
    return ["split", "split_size", "branch_id", "gCF", "gDF1", "gDF2",
            "gDFP", "gN", "qCF", "qDF1", "qDF2", "CU_astral", "CU_wastral",
            "x_gcf", "gDFP_flag", "wastral_status"]


# ----------------------------------------------------------------------
# ASTER output
# ----------------------------------------------------------------------
def parse_astral_tree(path: str) -> Dict[Split, Dict[str, str]]:
    """Parse ASTER/CASTER annotated Newick labels into split annotations."""
    tree = parse_newick(Path(path).read_text(encoding="utf-8", errors="replace"))
    taxa = {tip.name for tip in tree.get_terminals() if tip.name}
    out: Dict[Split, Dict[str, str]] = {}
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        side = {tip.name for tip in clade.get_terminals() if tip.name}
        if not (1 < len(side) < len(taxa) - 1):
            continue
        split = canonical_split(side, taxa)
        if not split:
            continue
        label = _node_label(clade)
        annotations = {}
        for key in ("CULength", "SULength", "q1", "q2", "q3", "localPP"):
            match = re.search(r"\b" + key + r"=([^;\]]+)", label)
            if match:
                value = match.group(1).strip()
                try:
                    if not math.isfinite(float(value)):
                        value = ""
                except ValueError:
                    pass
                annotations[key] = value
        out[split] = annotations
    return out


def check_wastral_labels(path: str) -> Tuple[bool, str, Dict[str, Any]]:
    """Check whether every internal gene-tree branch has a numeric label.

    wASTRAL's weighted parser needs one numeric support/length value per
    internal branch. Missing labels, composite labels such as ``/1/100`` and
    non-numeric labels are rejected instead of being silently coerced.
    """
    lines = newick_lines(path)
    if not lines:
        return False, "gene-tree file is empty", {"trees": 0}
    total = 0
    missing = 0
    nonnumeric = 0
    composite = 0
    for line in lines:
        try:
            tree = parse_newick(line)
        except Exception:
            return False, "gene-tree file contains an unparseable Newick tree", {
                "trees": len(lines), "missing": missing,
                "nonnumeric": nonnumeric, "composite": composite,
            }
        total += 1
        for clade in tree.get_nonterminals():
            if clade is tree.root:
                continue
            value = clade.name
            if value is None and clade.confidence is not None:
                value = str(clade.confidence)
            text = str(value).strip() if value is not None else ""
            if not text:
                missing += 1
                continue
            if "/" in text:
                composite += 1
                continue
            try:
                number = float(text)
            except ValueError:
                nonnumeric += 1
                continue
            if not math.isfinite(number):
                nonnumeric += 1
    detail = {
        "trees": total,
        "missing": missing,
        "composite": composite,
        "nonnumeric": nonnumeric,
    }
    if composite:
        return False, "composite branch labels (for example /1/100)", detail
    if missing:
        return False, "%d internal branches have no numeric label" % missing, detail
    if nonnumeric:
        return False, "%d internal branches have non-numeric labels" % nonnumeric, detail
    return True, "", detail


def x_from_gcf(value: Any) -> str:
    """Return the gCF-derived CU length, or blank when p <= 1/3."""
    try:
        p = float(value) / 100.0
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(p) or p <= 1.0 / 3.0 or p >= 1.0:
        return ""
    return "%.6g" % (-math.log(1.5 * (1.0 - p)))


def normalize_cf(value: Any) -> str:
    try:
        p = float(value) / 100.0
    except (TypeError, ValueError):
        return ""
    return "%.6g" % p if math.isfinite(p) else ""


def merge_diagnostics_rows(
    cf_columns: List[str],
    cf_rows_: List[Row],
    astral_rows: Dict[Split, Dict[str, str]],
    wastral_rows: Optional[Dict[Split, Dict[str, str]]],
    wastral_status: str,
) -> tuple[List[str], List[Row], Dict[str, Any]]:
    rows = []
    for cf in cf_rows_:
        split = tuple(cf["split"].split(",")) if cf.get("split") else ()
        astral = astral_rows.get(split, {})
        weighted = (wastral_rows or {}).get(split, {})
        gcf = cf.get("gCF", "")
        gdfp = cf.get("gDFP", "0")
        try:
            has_polyphyly = float(gdfp or 0) > 0
        except ValueError:
            has_polyphyly = False
        row: Row = {
            "split": cf["split"],
            "split_size": cf["split_size"],
            "branch_id": cf["branch_id"],
            "gCF": gcf,
            "gDF1": cf.get("gDF1", ""),
            "gDF2": cf.get("gDF2", ""),
            "gDFP": gdfp,
            "gN": cf.get("gN", ""),
            "qCF": astral.get("q1", ""),
            "qDF1": astral.get("q2", ""),
            "qDF2": astral.get("q3", ""),
            "CU_astral": astral.get("CULength", ""),
            "CU_wastral": weighted.get("CULength", ""),
            "x_gcf": x_from_gcf(gcf),
            "gDFP_flag": "yes" if has_polyphyly else "",
            "wastral_status": wastral_status if wastral_status else "available",
        }
        rows.append(row)
    rows.sort(key=lambda row: (int(row["split_size"]), row["split"]))
    summary = {
        "branches": len(rows),
        "astral_branches": len(astral_rows),
        "wastral_branches": len(wastral_rows or {}),
        "wastral_status": wastral_status or "available",
    }
    return columns_for("diagnostics"), rows, summary


def merge_site_diagnostics_rows(
    cf_columns: List[str],
    cf_rows_: List[Row],
    caster_rows: Dict[Split, Dict[str, str]],
) -> tuple[List[str], List[Row], Dict[str, Any]]:
    """Merge IQ-TREE sCF with CASTER-site q1/q2/q3 on the same splits."""
    rows: List[Row] = []
    seen = set()
    for cf in cf_rows_:
        split = tuple(cf["split"].split(",")) if cf.get("split") else ()
        caster = caster_rows.get(split, {})
        q1 = caster.get("q1", "")
        q2 = caster.get("q2", "")
        q3 = caster.get("q3", "")
        negative = False
        for value in (q1, q2, q3):
            try:
                if float(value) < 0:
                    negative = True
            except (TypeError, ValueError):
                pass
        rows.append({
            "split": cf["split"],
            "split_size": cf["split_size"],
            "branch_id": cf["branch_id"],
            "sCF": cf.get("sCF", ""),
            "sDF1": cf.get("sDF1", ""),
            "sDF2": cf.get("sDF2", ""),
            "sN": cf.get("sN", ""),
            "q1_caster": q1,
            "q2_caster": q2,
            "q3_caster": q3,
            "q_negative": "yes" if negative else "",
        })
        seen.add(split)
    for split, caster in caster_rows.items():
        if split in seen:
            continue
        q1 = caster.get("q1", "")
        q2 = caster.get("q2", "")
        q3 = caster.get("q3", "")
        negative = False
        for value in (q1, q2, q3):
            try:
                if float(value) < 0:
                    negative = True
            except (TypeError, ValueError):
                pass
        rows.append({
            "split": split_label(split),
            "split_size": str(len(split)),
            "branch_id": "",
            "sCF": "",
            "sDF1": "",
            "sDF2": "",
            "sN": "",
            "q1_caster": q1,
            "q2_caster": q2,
            "q3_caster": q3,
            "q_negative": "yes" if negative else "",
        })
    rows.sort(key=lambda row: (int(row["split_size"]), row["split"]))
    summary = {
        "branches": len(rows),
        "scf_branches": len(cf_rows_),
        "caster_branches": len(caster_rows),
        "source": "IQ-TREE 3 + CASTER-site",
    }
    return columns_for("site_diagnostics"), rows, summary


# ----------------------------------------------------------------------
# table output
# ----------------------------------------------------------------------
def write_table(path: str, columns: Sequence[str], rows: Sequence[Row]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(column, "")) for column in columns) + "\n")


def table_payload(path: str, columns: Sequence[str], rows: Sequence[Row],
                  table_type: str, **extra: Any) -> Dict[str, Any]:
    payload = {
        "table_type": table_type,
        "columns": list(columns),
        "rows": [dict(row) for row in rows],
        "table_path": path,
    }
    payload.update(extra)
    return payload
