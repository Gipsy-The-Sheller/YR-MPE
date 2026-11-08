"""Smoke-test every configured YR-MPE binary and a small functional subset.

Usage:
    python dev/integration/smoke_all_tools.py
    python dev/integration/smoke_all_tools.py --functional
    python dev/integration/smoke_all_tools.py --json out.json

The version pass only checks that the binary starts and prints a version/help
message. The functional pass runs tiny, local analyses (four taxa) for tools
that have a non-interactive smoke command. GUI-only plugins and long MCMC
programs are reported as skipped rather than executed.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "YR_MPE" / "config.json"
FIXTURE = ROOT / "dev" / "model" / "testdata" / "four_taxa.fas"


def resolve_tool(path: str) -> Path:
    # config.json uses POSIX-style rooted paths ("/bin/...") that are
    # relative to the YR_MPE package, not absolute Windows paths.
    if path.startswith(("/", "\\")):
        path = path.lstrip("/\\")
    elif os.path.isabs(path):
        return Path(path)
    return ROOT / "YR_MPE" / path


def command_for(path: Path, args: List[str]) -> List[str]:
    if path.suffix.lower() in (".bat", ".cmd"):
        return ["cmd", "/c", str(path), *args]
    return [str(path), *args]


def run_command(path: Path, args: List[str], timeout: int = 30):
    try:
        result = subprocess.run(
            command_for(path, args),
            capture_output=True,
            text=True,
            errors="replace",
            timeout=timeout,
            cwd=str(path.parent),
        )
        output = ((result.stdout or "") + (result.stderr or "")).strip()
        return {
            "ok": result.returncode == 0 or bool(output),
            "returncode": result.returncode,
            "output": output.splitlines()[:8],
        }
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": "timeout after %ss" % timeout}
    except Exception as exc:
        return {"ok": False, "error": str(exc)}


VERSION_ARGS = {
    "GBlocks": ["-h"],
    "astral": ["astral4", "--version"],
    "MrBayes-MPI-BEAGLE": ["-h"],
    "LSD2": [],
    "DecentTree": [],
    "MPBoot": ["-h"],
    "TrimAl": [],
    "PhyloBayes-bpcomp": ["-h"],
    "MPIRun": ["-help"],
}


def version_check(tool: Dict[str, Any]) -> Dict[str, Any]:
    path = resolve_tool(tool["path"])
    if not path.exists():
        return {"name": tool["name"], "path": str(path), "ok": False,
                "error": "missing binary"}
    args = VERSION_ARGS.get(tool["name"], ["--version"])
    result = run_command(path, args, timeout=12)
    result.update({"name": tool["name"], "path": str(path),
                   "args": args})
    return result


def functional_check(name: str, path: Path, args: List[str],
                     output_candidates: List[Path],
                     timeout: int = 120) -> Dict[str, Any]:
    try:
        result = subprocess.run(
            command_for(path, args),
            capture_output=True,
            text=True,
            errors="replace",
            timeout=timeout,
            cwd=str(path.parent),
        )
        output = ((result.stdout or "") + (result.stderr or "")).strip()
        produced = [str(item) for item in output_candidates if item.exists()
                    and item.stat().st_size > 0]
        ok = result.returncode == 0 and bool(produced)
        return {
            "name": name,
            "ok": ok,
            "returncode": result.returncode,
            "output": output.splitlines()[:5],
            "artifacts": produced,
        }
    except subprocess.TimeoutExpired:
        return {"name": name, "ok": False,
                "error": "timeout after %ss" % timeout}
    except Exception as exc:
        return {"name": name, "ok": False, "error": str(exc)}


def run_functional() -> List[Dict[str, Any]]:
    if not FIXTURE.is_file():
        return [{"name": "functional", "ok": False,
                 "error": "fixture not found: %s" % FIXTURE}]
    work = Path(tempfile.mkdtemp(prefix="yrmpe-smoke-"))
    source = work / "four_taxa.fas"
    source.write_text(FIXTURE.read_text(encoding="utf-8"), encoding="utf-8")
    tools = {tool["name"]: resolve_tool(tool["path"])
             for tool in json.loads(CONFIG.read_text(encoding="utf-8"))}
    results: List[Dict[str, Any]] = []

    # MAFFT writes an alignment to stdout.
    mafft_out = work / "mafft.fas"
    mafft = tools["MAFFT"]
    try:
        proc = subprocess.run(
            command_for(mafft, ["--auto", str(source)]),
            capture_output=True, text=True, errors="replace",
            timeout=120, cwd=str(work),
        )
        mafft_out.write_text(proc.stdout or "", encoding="utf-8")
        results.append({
            "name": "MAFFT", "ok": proc.returncode == 0
            and mafft_out.stat().st_size > 0,
            "returncode": proc.returncode,
            "artifacts": [str(mafft_out)],
            "output": ((proc.stderr or "") + (proc.stdout or "")).splitlines()[:4],
        })
    except Exception as exc:
        results.append({"name": "MAFFT", "ok": False, "error": str(exc)})

    muscle_out = work / "muscle.fas"
    results.append(functional_check(
        "MUSCLE5", tools["Muscle5"],
        ["-align", str(source), "-output", str(muscle_out)],
        [muscle_out],
    ))

    clustal_out = work / "clustalo.fas"
    results.append(functional_check(
        "Clustal Omega", tools["Clustal Omega"],
        ["-i", str(source), "-o", str(clustal_out),
         "--force", "--threads=1"],
        [clustal_out],
    ))

    trimal_out = work / "trimal.fas"
    results.append(functional_check(
        "TrimAl", tools["TrimAl"],
        ["-in", str(source), "-out", str(trimal_out), "-automated1"],
        [trimal_out],
    ))

    iqtree_prefix = work / "iqtree"
    results.append(functional_check(
        "IQ-TREE 3", tools["IQ-TREE 3"],
        ["-s", str(source), "-m", "GTR+G4", "-nt", "1",
         "-pre", str(iqtree_prefix), "-redo"],
        [Path(str(iqtree_prefix) + ".treefile")],
        timeout=180,
    ))

    modeltest_prefix = work / "modeltest"
    results.append(functional_check(
        "ModelTest-NG", tools["ModelTest-NG"],
        ["-i", str(source), "-d", "nt", "-t", "ml", "-p", "1",
         "-o", str(modeltest_prefix), "--force", "--disable-checkpoint"],
        [Path(str(modeltest_prefix) + ".out")],
        timeout=180,
    ))

    raxml_prefix = work / "raxml"
    results.append(functional_check(
        "RAxML-NG", tools["RAxML-NG"],
        ["--search1", "--msa", str(source), "--model", "GTR+G4",
         "--threads", "1", "--prefix", str(raxml_prefix), "--redo"],
        [Path(str(raxml_prefix) + ".raxml.bestTree")],
        timeout=180,
    ))

    decenttree_out = work / "decenttree.nwk"
    results.append(functional_check(
        "DecentTree", tools["DecentTree"],
        ["-fasta", str(source), "-out", str(decenttree_out)],
        [decenttree_out],
    ))

    mpboot_prefix = work / "mpboot"
    results.append(functional_check(
        "MPBoot", tools["MPBoot"],
        ["-s", str(source), "-pre", str(mpboot_prefix), "-n", "10"],
        [Path(str(mpboot_prefix) + ".treefile")],
        timeout=180,
    ))

    # GBlocks is interactive unless given at least a second parameter; leave
    # it out of the automated functional pass.
    results.append({
        "name": "GBlocks", "ok": True, "skipped": True,
        "reason": "interactive CLI; only version/help smoke-tested",
    })
    for name in ("MrBayes-MPI-BEAGLE", "PhyloBayes-mpi", "PhyloBayes-bpcomp",
                 "astral", "caster-site", "MACSE2", "LSD2", "TreeRecs",
                 "OrthoFinder"):
        results.append({
            "name": name, "ok": True, "skipped": True,
            "reason": "long-running/interactive; version smoke-test only",
        })
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--functional", action="store_true",
                        help="also run tiny functional analyses")
    parser.add_argument("--json", default="",
                        help="write the report as JSON")
    args = parser.parse_args()

    tools = json.loads(CONFIG.read_text(encoding="utf-8"))
    report = {
        "config": str(CONFIG),
        "version_checks": [version_check(tool) for tool in tools],
    }
    if args.functional:
        report["functional_checks"] = run_functional()
    failed = [row for row in report["version_checks"] if not row.get("ok")]
    if args.functional:
        failed += [row for row in report["functional_checks"]
                   if not row.get("ok") and not row.get("skipped")]
    report["ok"] = not failed

    for row in report["version_checks"]:
        status = "OK" if row.get("ok") else "FAIL"
        print("%-28s %s %s" % (row["name"], status,
                               (row.get("output") or [row.get("error", "")])[0]))
    if args.functional:
        print("\nFunctional checks:")
        for row in report["functional_checks"]:
            status = "SKIP" if row.get("skipped") else (
                "OK" if row.get("ok") else "FAIL")
            detail = row.get("reason") or (row.get("output") or [""])[0]
            print("%-28s %s %s" % (row["name"], status, detail))

    if args.json:
        Path(args.json).write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
