"""Parse OrthoFinder outputs and map single-copy orthogroups back to CDS."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_cds_map(path: Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["protein_id"]] = row
    return mapping


def find_orthogroups_dir(root: Path) -> Path:
    candidates = sorted(Path(root).rglob("Orthogroups.tsv"))
    if not candidates:
        raise FileNotFoundError(f"Orthogroups.tsv not found under {root}")
    return candidates[0].parent


def read_orthogroups(path: Path) -> dict[str, dict[str, list[str]]]:
    """Return {orthogroup: {species: [gene_id, ...]}}."""
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0].lower() != "orthogroup":
            raise ValueError(f"Unexpected Orthogroups.tsv header: {reader.fieldnames}")
        species_names = reader.fieldnames[1:]
        result: dict[str, dict[str, list[str]]] = {}
        for row in reader:
            og = row["Orthogroup"]
            result[og] = {}
            for species in species_names:
                raw = (row.get(species) or "").strip()
                genes = [g.strip() for g in raw.split(",") if g.strip()]
                result[og][species] = genes
        return result


def read_single_copy_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split()[0].rstrip(":"))
    return ids


def _resolve_gene(gene: str, cds_map: dict[str, dict[str, str]]) -> dict[str, str] | None:
    if gene in cds_map:
        return cds_map[gene]
    # OrthoFinder sometimes adds a species prefix or a suffix; try the
    # longest matching protein ID as a fallback.
    candidates = [key for key in cds_map if key.endswith(gene) or gene.endswith(key)]
    if len(candidates) == 1:
        return cds_map[candidates[0]]
    return None


def write_single_copy_ogs(
    orthogroups: dict[str, dict[str, list[str]]],
    single_copy_ids: set[str],
    cds_map: dict[str, dict[str, str]],
    output_dir: Path,
) -> dict[str, int]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary: list[dict[str, str]] = []
    written = 0

    for og in sorted(single_copy_ids):
        if og not in orthogroups:
            continue
        by_species = orthogroups[og]
        if not by_species or any(len(genes) != 1 for genes in by_species.values()):
            continue

        protein_records: list[SeqRecord] = []
        cds_records: list[SeqRecord] = []
        map_rows: list[dict[str, str]] = []
        ok = True
        for species, genes in sorted(by_species.items()):
            entry = _resolve_gene(genes[0], cds_map)
            if entry is None:
                ok = False
                break
            protein_id = entry["protein_id"]
            # Use the species label as the sequence id.  This is the key that
            # the concatenator and partitioned-NEXUS exporter use to match
            # orthologues across single-copy OG files; the original protein
            # id remains in the map table and in dbxrefs.
            protein_records.append(
                SeqRecord(Seq(entry["aa_seq"]), id=species,
                          description=species,
                          dbxrefs=["protein_id:" + protein_id])
            )
            cds_records.append(
                SeqRecord(Seq(entry["cds_seq"]), id=species,
                          description=species,
                          dbxrefs=["protein_id:" + protein_id])
            )
            map_rows.append({
                "orthogroup": og,
                "species": species,
                "protein_id": protein_id,
                "cds_id": protein_id,
                "gene": entry.get("gene", ""),
                "product": entry.get("product", ""),
                "record_accession": entry.get("record_accession", ""),
                "start": entry.get("start", ""),
                "end": entry.get("end", ""),
                "strand": entry.get("strand", ""),
                "transl_table": entry.get("transl_table", ""),
            })
        if not ok or not protein_records:
            continue

        SeqIO.write(protein_records, output_dir / f"{og}.faa", "fasta")
        SeqIO.write(cds_records, output_dir / f"{og}.cds.fna", "fasta")
        with (output_dir / f"{og}.map.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(map_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(map_rows)

        summary.append({
            "orthogroup": og,
            "species_count": str(len(map_rows)),
            "protein_count": str(len(protein_records)),
            "status": "ok",
        })
        written += 1

    with (output_dir / "single_copy_og_summary.tsv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["orthogroup", "species_count", "protein_count", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    return {"orthogroups_written": written, "summary_rows": len(summary)}


def derive_single_copy_ids(orthogroups: dict[str, dict[str, list[str]]]) -> set[str]:
    return {
        og for og, by_species in orthogroups.items()
        if by_species and all(len(genes) == 1 for genes in by_species.values())
    }
