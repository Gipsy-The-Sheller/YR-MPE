"""Audit skill-manifest coverage against settings.json.

Usage (from the repository root):

    python dev/skills/audit_skills.py

Exit status is 1 when a manifest document file is missing or a configured
analysis plugin has no catalog entry.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SETTINGS = REPO / "settings.json"
MANIFEST = REPO / "YR_MPE" / "skills" / "manifest.json"
SKILLS = MANIFEST.parent

# Application shells and the template plugin are not analysis skills.
EXCLUDED_PREFIXES = ("YPhylo Workbench", "About YPhylo Workbench", "Example Plugin")
EXCLUDED_NAMES = {"YPhylo Workbench"}


def is_excluded(name: str) -> bool:
    return name in EXCLUDED_NAMES or name.startswith(EXCLUDED_PREFIXES)


def norm(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").lower())


def class_stem(entry_point: str) -> str:
    name = str(entry_point or "").rsplit(":", 1)[-1]
    for suffix in ("PluginEntry", "Plugin", "Entry"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return norm(name)


def main() -> int:
    settings = json.loads(SETTINGS.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    entries = manifest.get("skills", [])

    catalogued = set()
    for entry in entries:
        for key in ("id", "title", "plugin_names", "tool_names", "class_names"):
            value = entry.get(key, [])
            if isinstance(value, str):
                value = [value]
            catalogued.update(norm(item) for item in value)

    problems = []
    for entry in entries:
        file_name = entry.get("file", "")
        if not file_name:
            problems.append("NO DOC  %-28s %s" % (entry["id"], entry.get("title", "")))
        elif not (SKILLS / file_name).is_file():
            problems.append("MISSING %-28s %s" % (entry["id"], file_name))

    listed = set()
    for plugin in settings.get("plugins", []):
        name = plugin.get("name", "")
        if not name or is_excluded(name):
            continue
        listed.add(name)
        candidates = {norm(name), class_stem(plugin.get("entry_point", ""))}
        if not any(
            candidate and (
                candidate in catalog_key or catalog_key in candidate
            )
            for candidate in candidates
            for catalog_key in catalogued
        ):
            problems.append("NO ENTRY %-27s %s" % (name, plugin.get("entry_point", "")))

    if problems:
        print("\n".join(problems))
        print("\n%d issue(s)" % len(problems))
        return 1
    print("skills: %d entries, %d configured plugins, no gaps" %
          (len(entries), len(listed)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
