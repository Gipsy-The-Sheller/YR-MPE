"""Build the versioned source archive from tracked and new source files."""

from __future__ import annotations

import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERSION = "1.2.1"
PREFIX = "YR-MPE-%s" % VERSION
OUT = ROOT / "dist" / ("%s-src.tar.gz" % PREFIX)

EXCLUDE = {
    "settings.json",  # user workspace history; ship settings.template.json instead
    "YR_MPE/config.json",  # local executable paths; the Initializer writes it
}

EXTRA = [
    "YR_MPE/plugins/workflow_runner_plugin.py",
    "YR_MPE/plugins/components/masked_alignment_report.py",
    "YR_MPE/skills/workflow-runner.md",
    "dev/workflow/test_workflow_runner.py",
    "dev/workspace/test_workspace_reopen.py",
    "dev/workspace/test_provenance.py",
    "dev/agent/test_control_bindings.py",
    "dev/release/build_platform_bundle.py",
    "dev/release/test_build_platform_bundle.py",
    ".github/workflows/build.yml",
    "CYAN_RELEASE",
    "cyan_tools.json",
    "benchmark/mrbayes/README.md",
    "benchmark/mrbayes/bench_mrbayes.py",
    "benchmark/mrbayes/summarize_mrbayes_bench.py",
    "benchmark/mrbayes/redraw_figS2.py",
    "benchmark/mrbayes/mrbayes_beagle.json",
    "benchmark/mrbayes/summary.txt",
]


def main() -> None:
    files = subprocess.check_output(
        ["git", "ls-files"], cwd=ROOT, text=True
    ).splitlines()
    files.extend(EXTRA)
    seen = set()
    selected = []
    for rel in files:
        rel = rel.strip().replace("\\", "/")
        if not rel or rel in seen or rel in EXCLUDE:
            continue
        seen.add(rel)
        path = ROOT / rel
        if path.is_file():
            selected.append((rel, path))

    OUT.parent.mkdir(exist_ok=True)
    with tarfile.open(OUT, "w:gz", compresslevel=9) as tar:
        for rel, path in selected:
            info = tar.gettarinfo(str(path), arcname="%s/%s" % (PREFIX, rel))
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            with path.open("rb") as handle:
                tar.addfile(info, handle)
    print(OUT, OUT.stat().st_size, "files", len(selected))


if __name__ == "__main__":
    main()
