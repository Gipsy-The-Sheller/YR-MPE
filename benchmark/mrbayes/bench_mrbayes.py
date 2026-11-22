"""Benchmark the CYAN MrBayes 3.2.8 + BEAGLE 4 package.

For each case the script runs the same MrBayes binary twice:

* ``serial``  - one process, BEAGLE disabled;
* ``mpi``     - ``mpiexec -n N`` with CPU BEAGLE threading enabled.

The alignments are simulated under a Jukes-Cantor process on a random tree,
so the comparison is about likelihood cost, not about a particular dataset.
Wall-clock time and Robinson-Foulds distance between the two consensus trees
are written to JSON.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import subprocess
import time
from pathlib import Path

from Bio import Phylo

BASES = "ACGT"
CASES = [
    ("20OTU_1000bp", 20, 1000),
    ("20OTU_100bp", 20, 100),
    ("50OTU_1000bp", 50, 1000),
    ("50OTU_100bp", 50, 100),
]


def random_tree(ntaxa: int, rng: random.Random):
    """Return a Newick string for a random binary tree with JC branch lengths."""
    nodes = [f"t{i + 1}" for i in range(ntaxa)]
    lengths = {name: rng.expovariate(20.0) + 0.01 for name in nodes}
    next_id = 0
    while len(nodes) > 2:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        label = f"i{next_id}"
        next_id += 1
        subtree = f"({a}:{lengths[a]:.6f},{b}:{lengths[b]:.6f}){label}"
        nodes.append(subtree)
        lengths[subtree] = rng.expovariate(20.0) + 0.01
    a, b = nodes
    return f"({a}:{lengths[a]:.6f},{b}:{lengths[b]:.6f});"


def parse_newick(newick: str):
    """Return parent/child structure with branch lengths."""
    text = newick.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        if text[pos] == "(":
            pos += 1
            children = [node()]
            while text[pos] == ",":
                pos += 1
                children.append(node())
            assert text[pos] == ")"
            pos += 1
            while pos < len(text) and text[pos] not in ",():":
                pos += 1
            length = 0.0
            if pos < len(text) and text[pos] == ":":
                pos += 1
                start = pos
                while pos < len(text) and text[pos] not in ",()":
                    pos += 1
                length = float(text[start:pos])
            return {"children": children, "length": length, "name": ""}
        start = pos
        while pos < len(text) and text[pos] not in ",():":
            pos += 1
        name = text[start:pos]
        length = 0.0
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            length = float(text[start:pos])
        return {"children": [], "length": length, "name": name}

    return node()


def simulate(ntaxa: int, length: int, seed: int):
    rng = random.Random(seed)
    tree = parse_newick(random_tree(ntaxa, rng))
    sequences = {}

    def evolve(node, sequence):
        for child in node["children"]:
            child_sequence = []
            distance = child["length"]
            change = 1.0 - math.exp(-distance * 0.75)
            for base in sequence:
                if rng.random() < change:
                    child_sequence.append(rng.choice(
                        [b for b in BASES if b != base]))
                else:
                    child_sequence.append(base)
            if child["children"]:
                evolve(child, child_sequence)
            else:
                sequences[child["name"]] = "".join(child_sequence)

    evolve(tree, [rng.choice(BASES) for _ in range(length)])
    return sequences


def write_nexus(path: Path, sequences: dict, mpi: bool, ngen: int,
                mpi_procs: int, beagle_threads: int, nruns: int = 2) -> None:
    names = list(sequences)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#NEXUS\nbegin data;\n")
        handle.write(" dimensions ntax=%d nchar=%d;\n" % (
            len(names), len(next(iter(sequences.values())))))
        handle.write(" format datatype=dna missing=? gap=-;\n matrix\n")
        for name in names:
            handle.write("  %s %s\n" % (name, sequences[name]))
        handle.write(" ;\nend;\nbegin mrbayes;\n")
        if mpi:
            handle.write(
                " set usebeagle=yes beagledevice=cpu beagleprecision=double "
                "beaglescaling=dynamic beaglethreads=yes "
                "beaglethreadcount=%d;\n" % beagle_threads)
        # Serial runs leave BEAGLE off by default; MrBayes 3.2.8 aborts at
        # exit when `set usebeagle=no` is combined with autoclose.
        handle.write(" lset nst=6 rates=gamma ngammacat=4;\n")
        handle.write(
            " prset statefreqpr=dirichlet(1,1,1,1) "
            "revmatpr=dirichlet(1,1,1,1,1,1) shapepr=exponential(1.0);\n")
        handle.write(
            " mcmcp ngen=%d samplefreq=100 nchains=4 nruns=%d "
            "printfreq=1000 savebrlens=yes;\n" % (ngen, nruns))
        handle.write(" mcmc;\n")
        burnin = max(1, ngen // 4)
        handle.write(" sump burnin=%d;\n" % burnin)
        handle.write(
            " sumt burnin=%d contype=majorityrule conformat=simple;\n" % burnin)
        handle.write("end;\n")


def load_tree(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"tree\s+\S+\s*=\s*(?:\[&[UR]\]\s*)?(\(.*?;)", text,
        re.S | re.I,
    )
    newick = match.group(1) if match else text.strip()
    newick = re.sub(r"\[&[UR]\]", "", newick)
    tmp = path.with_suffix(".clean.tre")
    tmp.write_text(newick, encoding="utf-8")
    return Phylo.read(str(tmp), "newick")


def splits(tree):
    taxa = {terminal.name for terminal in tree.get_terminals()}
    result = set()
    for clade in tree.get_nonterminals():
        tips = frozenset(t.name for t in clade.get_terminals())
        if not tips or tips == taxa:
            continue
        other = frozenset(taxa - tips)
        result.add(tips if len(tips) <= len(other) else other)
    return result


def rf_distance(first, second) -> int:
    return len(splits(first) ^ splits(second))


def run(cmd, cwd, env, timeout):
    started = time.perf_counter()
    result = subprocess.run(
        cmd, cwd=str(cwd), env=env, timeout=timeout,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(
            "%s failed (%s):\n%s" % (" ".join(cmd), result.returncode,
                                    result.stdout[-3000:]))
    return elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mb-dir", required=True,
                        help="extracted CYAN mrbayes package directory")
    parser.add_argument("--serial-mb", default="",
                        help="optional separate serial binary (for example "
                             "the official MrBayes build); default: <mb-dir>/mb.exe")
    parser.add_argument("--out", default="mrbayes_beagle.json")
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--ngen", type=int, default=50000)
    parser.add_argument("--mpi", type=int, default=8)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--nruns", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--serial-only", action="store_true",
                        help="validate the harness without MPI")
    parser.add_argument("--case", default="",
                        help="run one case name (for example 50OTU_1000bp)")
    args = parser.parse_args()

    mb_dir = Path(args.mb_dir).resolve()
    mb = mb_dir / "mb.exe"
    serial_mb = Path(args.serial_mb).resolve() if args.serial_mb else mb
    mpiexec = mb_dir / "mpiexec.exe"
    if not mb.is_file():
        raise SystemExit("mb.exe not found in %s" % mb_dir)
    if not serial_mb.is_file():
        raise SystemExit("serial MrBayes not found in %s" % serial_mb)
    if not args.serial_only and not mpiexec.is_file():
        raise SystemExit("mpiexec.exe not found in %s" % mb_dir)

    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(mb_dir), env.get("PATH", "")])

    root = Path("mrbayes_benchmark").resolve()
    root.mkdir(exist_ok=True)
    rows = []
    for case_name, ntaxa, length in CASES:
        if args.case and case_name != args.case:
            continue
        case_dir = root / case_name
        case_dir.mkdir(exist_ok=True)
        sequences = simulate(ntaxa, length, seed=ntaxa * 10000 + length)
        serial_nex = case_dir / "serial.nex"
        mpi_nex = case_dir / "mpi.nex"
        write_nexus(serial_nex, sequences, False, args.ngen, args.mpi,
                    args.threads, args.nruns)
        write_nexus(mpi_nex, sequences, True, args.ngen, args.mpi,
                    args.threads, args.nruns)

        serial_times, mpi_times, rf_values = [], [], []
        for rep in range(1, args.reps + 1):
            # MrBayes resumes from .ckp files; a timed replicate must start
            # from a clean directory.
            for prefix in ("serial.nex", "mpi.nex"):
                for old in case_dir.glob(prefix + ".*"):
                    old.unlink()
            serial_time = run([str(serial_mb), str(serial_nex)], case_dir,
                              env, args.timeout)
            serial_tree = case_dir / "serial.nex.con.tre"
            serial_times.append(serial_time)
            if args.serial_only:
                print("%s rep %d: serial %.1fs" % (
                    case_name, rep, serial_time), flush=True)
                continue
            for old in case_dir.glob("*.con.tre"):
                if old.name != "serial.nex.con.tre":
                    old.unlink()
            mpi_time = run(
                [str(mpiexec), "-n", str(args.mpi), str(mb),
                 str(mpi_nex)],
                case_dir, env, args.timeout)
            mpi_tree = case_dir / "mpi.nex.con.tre"
            mpi_times.append(mpi_time)
            if serial_tree.is_file() and mpi_tree.is_file():
                rf_values.append(rf_distance(
                    load_tree(serial_tree), load_tree(mpi_tree)))
            print("%s rep %d: serial %.1fs mpi %.1fs" % (
                case_name, rep, serial_time, mpi_time), flush=True)

        serial_mean = sum(serial_times) / len(serial_times)
        mpi_mean = sum(mpi_times) / len(mpi_times) if mpi_times else None
        rows.append({
            "case": case_name,
            "ntaxa": ntaxa,
            "length": length,
            "serial_times": serial_times,
            "mpi_times": mpi_times,
            "time_ratio": (mpi_mean / serial_mean) if mpi_mean else None,
            "speedup": (serial_mean / mpi_mean) if mpi_mean else None,
            "rf": rf_values,
            "ngen": args.ngen,
            "mpi_processes": args.mpi,
            "nruns": args.nruns,
            "beagle_threads": args.threads,
            "serial_binary": serial_mb.name,
        })

    Path(args.out).write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
