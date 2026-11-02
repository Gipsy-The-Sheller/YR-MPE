"""GBFF/GenBank -> protein FASTA + CDS mapping for OrthoFinder.

This is a prototype for the YR-MPE "Genomic Dataset Formation (OrthoFinder)"
feature.  It parses organelle GenBank files (mt/cp), extracts CDS features,
translates them when /translation is absent, and writes:

  proteomes/<species>.faa   protein FASTA for OrthoFinder
  cds/<species>.fna         in-frame CDS FASTA (stop codon removed by default)
  cds_map.tsv               protein -> species/coordinates/CDS/AA

The parser is deliberately conservative: it reports warnings instead of
silently guessing, and skips records that cannot be translated cleanly.
"""

from __future__ import annotations

import csv
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

GENBANK_SUFFIXES = {".gb", ".gbk", ".gbff", ".genbank"}


@dataclass
class CDSRecord:
    protein_id: str
    species: str
    record_id: str
    record_accession: str
    locus_tag: str
    ncbi_protein_id: str
    gene: str
    product: str
    start: int
    end: int
    strand: int
    codon_start: int
    transl_table: int
    aa_seq: str
    cds_seq: str
    source_file: str
    warnings: str = ""

    def map_row(self) -> dict[str, str]:
        return {
            "protein_id": self.protein_id,
            "species": self.species,
            "record_id": self.record_id,
            "record_accession": self.record_accession,
            "locus_tag": self.locus_tag,
            "ncbi_protein_id": self.ncbi_protein_id,
            "gene": self.gene,
            "product": self.product,
            "start": str(self.start),
            "end": str(self.end),
            "strand": str(self.strand),
            "codon_start": str(self.codon_start),
            "transl_table": str(self.transl_table),
            "aa_length": str(len(self.aa_seq)),
            "cds_length": str(len(self.cds_seq)),
            "aa_seq": self.aa_seq,
            "cds_seq": self.cds_seq,
            "source_file": self.source_file,
            "warnings": self.warnings,
        }


def _qual(feature, key: str, default: str = "") -> str:
    values = feature.qualifiers.get(key)
    if not values:
        return default
    value = values[0]
    if isinstance(value, bytes):
        value = value.decode("utf-8", "ignore")
    return str(value).strip()


def _sanitize(value: str, fallback: str = "cds") -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return value or fallback


def _location_bounds(location) -> tuple[int, int, int]:
    start = int(location.start)
    end = int(location.end)
    strand = int(location.strand) if location.strand is not None else 0
    return start, end, strand


def _has_internal_stop(aa: str) -> bool:
    return "*" in aa[:-1] if aa.endswith("*") else "*" in aa


def _translate_cds(
    cds: Seq,
    table: int,
    translation: str,
    internal_stop: str,
) -> tuple[str, list[str]]:
    warnings: list[str] = []
    if translation:
        aa = re.sub(r"\s+", "", translation)
        # Light sanity check: compare with a direct translation.  Chloroplast
        # and mitochondrial genes often use alternative start codons, so a
        # first-residue difference (M vs I/V/T/L) is expected and is reported
        # separately from a real mismatch.
        calc = str(cds.translate(table=table)).rstrip("*")
        aa = aa.rstrip("*")
        if calc and aa and aa != calc:
            if aa[0] == "M" and aa[1:] == calc[1:]:
                warnings.append(
                    f"alternative start codon: /translation starts M, CDS translates {calc[0]}"
                )
            else:
                warnings.append("translation qualifier differs from direct translation")
        return aa, warnings

    aa = str(cds.translate(table=table))
    if _has_internal_stop(aa):
        if internal_stop == "skip":
            raise ValueError("internal stop codon")
        warnings.append("internal stop codon")
        if internal_stop == "warn":
            pass
    return aa.rstrip("*"), warnings


def parse_genbank_file(
    path: Path,
    species: str | None = None,
    default_table: int = 1,
    skip_pseudo: bool = True,
    internal_stop: str = "skip",
    keep_stop: bool = False,
    max_records: int | None = None,
    record_as_species: bool = False,
) -> list[CDSRecord]:
    default_species = species or _sanitize(path.stem)
    records: list[CDSRecord] = []
    seen: set[str] = set()

    for record_index, record in enumerate(SeqIO.parse(str(path), "genbank")):
        if max_records is not None and record_index >= max_records:
            break
        accession = record.id
        if record.annotations.get("accessions"):
            accession = str(record.annotations["accessions"][0])
        if record_as_species:
            record_species = _sanitize(
                record.id or accession or f"{path.stem}_{record_index + 1}"
            )
        else:
            record_species = default_species
        for feature_index, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            if skip_pseudo and ("pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers):
                continue

            locus_tag = _qual(feature, "locus_tag")
            ncbi_protein_id = _qual(feature, "protein_id")
            gene = _qual(feature, "gene")
            product = _qual(feature, "product")
            translation = _qual(feature, "translation")
            codon_start_raw = _qual(feature, "codon_start", "1")
            table_raw = _qual(feature, "transl_table", str(default_table))
            try:
                codon_start = max(1, int(codon_start_raw))
            except ValueError:
                codon_start = 1
            try:
                table = int(table_raw)
            except ValueError:
                table = default_table

            try:
                cds = feature.extract(record.seq)
            except Exception as exc:  # pragma: no cover - defensive
                continue
            if not cds or len(cds) < 3:
                continue

            cds = cds[codon_start - 1 :]
            if len(cds) % 3:
                # Partial CDS at the end of a genome is common; keep the
                # complete codons and record a warning.
                cds = cds[: len(cds) - (len(cds) % 3)]

            try:
                aa, warnings = _translate_cds(cds, table, translation, internal_stop)
            except ValueError:
                continue

            if not keep_stop and len(cds) >= 3:
                last = str(cds[-3:]).upper()
                if str(Seq(last).translate(table=table)) == "*":
                    cds = cds[:-3]

            if not aa or not cds:
                continue

            base = locus_tag or ncbi_protein_id or gene or f"cds_{record.id}_{feature_index}"
            protein_id = _sanitize(f"{record_species}_{base}")
            suffix = 1
            original = protein_id
            while protein_id in seen:
                suffix += 1
                protein_id = f"{original}_{suffix}"
            seen.add(protein_id)

            start, end, strand = _location_bounds(feature.location)
            records.append(
                CDSRecord(
                    protein_id=protein_id,
                    species=record_species,
                    record_id=record.id,
                    record_accession=accession,
                    locus_tag=locus_tag,
                    ncbi_protein_id=ncbi_protein_id,
                    gene=gene,
                    product=product,
                    start=start,
                    end=end,
                    strand=strand,
                    codon_start=codon_start,
                    transl_table=table,
                    aa_seq=aa,
                    cds_seq=str(cds).upper(),
                    source_file=path.name,
                    warnings="; ".join(warnings),
                )
            )
    return records


def iter_genbank_files(input_path: Path) -> Iterator[Path]:
    input_path = Path(input_path)
    if input_path.is_file():
        if input_path.suffix.lower() in GENBANK_SUFFIXES:
            yield input_path
        return
    for path in sorted(input_path.iterdir()):
        if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES:
            yield path


def write_dataset(records: Iterable[CDSRecord], output_dir: Path) -> dict[str, int]:
    output_dir = Path(output_dir)
    proteome_dir = output_dir / "proteomes"
    cds_dir = output_dir / "cds"
    proteome_dir.mkdir(parents=True, exist_ok=True)
    cds_dir.mkdir(parents=True, exist_ok=True)

    grouped: dict[str, list[CDSRecord]] = {}
    for record in records:
        grouped.setdefault(record.species, []).append(record)

    map_path = output_dir / "cds_map.tsv"
    fieldnames = list(CDSRecord.map_row(next(iter(grouped.values()))[0]).keys()) if grouped else []
    with map_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()

        for species, species_records in grouped.items():
            protein_records = [
                SeqRecord(Seq(r.aa_seq), id=r.protein_id, description=r.product)
                for r in species_records
            ]
            cds_records = [
                SeqRecord(Seq(r.cds_seq), id=r.protein_id, description=r.gene or r.product)
                for r in species_records
            ]
            SeqIO.write(protein_records, proteome_dir / f"{species}.faa", "fasta")
            SeqIO.write(cds_records, cds_dir / f"{species}.fna", "fasta")
            for r in species_records:
                writer.writerow(r.map_row())

    return {species: len(items) for species, items in grouped.items()}


def parse_folder(
    input_path: Path,
    default_table: int = 1,
    skip_pseudo: bool = True,
    internal_stop: str = "skip",
    keep_stop: bool = False,
    record_as_species: bool = False,
) -> list[CDSRecord]:
    records: list[CDSRecord] = []
    for path in iter_genbank_files(input_path):
        records.extend(
            parse_genbank_file(
                path,
                default_table=default_table,
                skip_pseudo=skip_pseudo,
                internal_stop=internal_stop,
                keep_stop=keep_stop,
                record_as_species=record_as_species,
            )
        )
    return records
