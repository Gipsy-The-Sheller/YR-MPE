"""Assemble a YR-MPE binary bundle from CYAN release assets.

The CYAN release publishes one archive per tool and target:
``<tool>-<version>-<Platform>-<arch>-<tag>.zip|tar.gz``.  This script selects
the assets for one target, extracts them under ``bin/`` (or ``src/`` for the
platform-neutral IcyTree site) and writes a single archive that extracts at
the ``YR_MPE/`` level.
"""

from __future__ import annotations

import argparse
import io
import shutil
import tarfile
import tempfile
import zipfile
from pathlib import Path

# Tools whose archive root belongs under YR_MPE/bin.
BIN_TOOLS = [
    "aster",
    "aster2",
    "bmge",
    "clustalo",
    "decenttree",
    "hyphy",
    "iqtree3",
    "lsd2",
    "mafft",
    "modeltest-ng",
    "mpboot",
    "mrbayes",
    "muscle5",
    "orthofinder3",
    "paml",
    "pbmpi",
    "raxml-ng",
    "treerecs",
    "trimal",
    "ygblocks",
]

# Network-independent web assets live under YR_MPE/src.
SRC_TOOLS = ["icytree"]

# One package serves every target (Java jar / static web site), but CYAN
# builds it in the Windows x64 job.
NEUTRAL_TOOLS = {"bmge", "icytree"}


def find_asset(assets: list[Path], tool: str, platform: str, arch: str):
    """Return the asset for ``tool`` on one target, or None."""
    if tool in NEUTRAL_TOOLS:
        platform, arch = "Windows", "x64"
    marker = f"-{platform}-{arch}-"
    for path in assets:
        name = path.name
        if not name.startswith(tool + "-") or marker not in name:
            continue
        if name.endswith(".zip") or name.endswith(".tar.gz"):
            return path
    return None


def _safe_member(name: str) -> bool:
    path = Path(name)
    return not path.is_absolute() and ".." not in path.parts


def extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as bundle:
            for member in bundle.infolist():
                if not _safe_member(member.filename):
                    raise ValueError(f"Unsafe path in {archive}: {member.filename}")
            bundle.extractall(destination)
    elif tarfile.is_tarfile(archive):
        with tarfile.open(archive) as bundle:
            for member in bundle.getmembers():
                if not _safe_member(member.name):
                    raise ValueError(f"Unsafe path in {archive}: {member.name}")
            bundle.extractall(destination)
    else:
        raise ValueError(f"Unsupported archive: {archive}")


def write_archive(staging: Path, output: Path, fmt: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if fmt == "zip":
        with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(staging.rglob("*")):
                if path.is_file():
                    bundle.write(path, path.relative_to(staging).as_posix())
        return
    with tarfile.open(output, "w:gz", compresslevel=9) as bundle:
        for path in sorted(staging.rglob("*")):
            if not path.is_file():
                continue
            info = bundle.gettarinfo(
                str(path), arcname=path.relative_to(staging).as_posix())
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            with path.open("rb") as handle:
                bundle.addfile(info, handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assets-dir", required=True)
    parser.add_argument("--platform", required=True,
                        choices=["Windows", "Linux", "macOS"])
    parser.add_argument("--arch", required=True, choices=["x64", "arm64"])
    parser.add_argument("--format", required=True, choices=["zip", "tar.gz"])
    parser.add_argument("--output", required=True)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()

    assets = sorted(Path(args.assets_dir).iterdir())
    staging = Path(tempfile.mkdtemp(prefix="yrmpe_bundle_"))
    found, missing = [], []
    try:
        for target, tools in (("bin", BIN_TOOLS), ("src", SRC_TOOLS)):
            for tool in tools:
                archive = find_asset(assets, tool, args.platform, args.arch)
                if archive is None:
                    missing.append(tool)
                    continue
                extract(archive, staging / target)
                found.append(f"{tool}: {archive.name}")

        (staging / "BUNDLE_INFO.txt").write_text(
            "YR-MPE binary bundle\n"
            f"CYAN release: {args.tag or 'unknown'}\n"
            f"Target: {args.platform}-{args.arch}\n"
            f"Format: {args.format}\n\n"
            "Included:\n" + "".join(f"  {row}\n" for row in found) +
            "\nMissing on this target (install from PATH or skip):\n" +
            "".join(f"  {row}\n" for row in missing),
            encoding="utf-8",
        )
        write_archive(staging, Path(args.output), args.format)
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    print(f"wrote {args.output}")
    print(f"included {len(found)}; missing {len(missing)}")
    for row in missing:
        print(f"  missing: {row}")


if __name__ == "__main__":
    main()
