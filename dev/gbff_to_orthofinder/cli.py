#!/usr/bin/env python3
"""Prototype CLI: GBFF -> OrthoFinder -> single-copy orthogroups -> CDS."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from gbff import parse_folder, write_dataset
from orthogroups import (
    derive_single_copy_ids,
    find_orthogroups_dir,
    read_cds_map,
    read_orthogroups,
    read_single_copy_ids,
    write_single_copy_ogs,
)


def run_orthofinder(
    output_dir: Path,
    orthofinder: Path,
    search_dirs: list[Path],
    threads: int,
) -> int:
    proteomes = (output_dir / "proteomes").resolve()
    results_root = (output_dir / "orthofinder").resolve()
    if results_root.exists():
        shutil.rmtree(results_root)

    env = os.environ.copy()
    path_parts = [str(p) for p in search_dirs if p]
    path_parts.append(str(orthofinder.parent))
    new_path = os.pathsep.join(path_parts + [env.get("PATH", "")])
    env["PATH"] = new_path
    # On Windows, CreateProcess resolves the child executable (blastp) using
    # the parent process PATH, not only the env= block.  Set it temporarily so
    # OrthoFinder can find blastp/makeblastdb/mcl.
    old_path = os.environ.get("PATH")
    os.environ["PATH"] = new_path

    # This OrthoFinder 2.5.5 build accepts only '-S blast' (its help text
    # mentions diamond, but the bundled options are blast-only).
    cmd = [
        str(orthofinder),
        "-f", str(proteomes),
        "-o", str(results_root),
        "-t", str(threads),
        "-S", "blast",
        "-og",
        "-X",
    ]
    print("running:", " ".join(cmd))
    try:
        result = subprocess.run(
            cmd, env=env, capture_output=True, text=True, encoding="utf-8", errors="replace"
        )
    finally:
        if old_path is None:
            os.environ.pop("PATH", None)
        else:
            os.environ["PATH"] = old_path
    if result.stdout:
        print(result.stdout[-4000:])
    if result.stderr:
        print(result.stderr[-4000:], file=sys.stderr)
    return result.returncode


def postprocess(output_dir: Path, orthogroups_root: Path | None = None) -> dict[str, int]:
    root = orthogroups_root or output_dir
    og_dir = find_orthogroups_dir(root)
    cds_map = read_cds_map(output_dir / "cds_map.tsv")
    orthogroups = read_orthogroups(og_dir / "Orthogroups.tsv")
    single_path = og_dir / "Orthogroups_SingleCopyOrthologues.txt"
    if single_path.exists():
        single_ids = read_single_copy_ids(single_path)
    else:
        single_ids = derive_single_copy_ids(orthogroups)
    return write_single_copy_ogs(
        orthogroups, single_ids, cds_map, output_dir / "single_copy_ogs"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="folder of .gb/.gbff/.gbk files")
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    parser.add_argument("--default-table", type=int, default=1,
                        help="fallback NCBI genetic code when /transl_table is absent")
    parser.add_argument("--skip-pseudo", action="store_true", default=True)
    parser.add_argument("--keep-pseudo", dest="skip_pseudo", action="store_false")
    parser.add_argument("--internal-stop", choices=["skip", "warn", "keep"], default="skip")
    parser.add_argument("--keep-stop", action="store_true",
                        help="keep the terminal stop codon in the CDS FASTA")
    parser.add_argument("--record-as-taxon", action="store_true",
                        help="treat each GenBank record as a separate taxon (multi-record .gb file)")
    parser.add_argument("--run-orthofinder", action="store_true")
    parser.add_argument("--orthofinder", type=Path,
                        default=Path("E:/Codes/gf_mito_b/orthofinder.exe"))
    parser.add_argument("--blast-dir", type=Path,
                        default=Path("E:/Codes/gf_mito_b/tools/blast_tb"),
                        help="directory containing blastp/makeblastdb")
    parser.add_argument("--diamond-dir", type=Path, default=None,
                        help="optional DIAMOND directory (not used by this OrthoFinder build)")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--postprocess", action="store_true")
    parser.add_argument("--orthogroups-root", type=Path,
                        help="existing OrthoFinder results root for --postprocess")
    args = parser.parse_args(argv)

    records = parse_folder(
        args.input,
        default_table=args.default_table,
        skip_pseudo=args.skip_pseudo,
        internal_stop=args.internal_stop,
        keep_stop=args.keep_stop,
        record_as_species=args.record_as_taxon,
    )
    if not records:
        print("no CDS records parsed", file=sys.stderr)
        return 1
    summary = write_dataset(records, args.out)
    print("parsed species/CDS:", summary)
    print("warnings:", sum(1 for r in records if r.warnings))

    if args.run_orthofinder:
        search_dirs = [args.blast_dir]
        if args.diamond_dir:
            search_dirs.append(args.diamond_dir)
        code = run_orthofinder(args.out, args.orthofinder, search_dirs, args.threads)
        if code != 0:
            print(f"OrthoFinder exited with code {code}", file=sys.stderr)
            return code

    if args.postprocess:
        result = postprocess(args.out, args.orthogroups_root)
        print("postprocess:", result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
