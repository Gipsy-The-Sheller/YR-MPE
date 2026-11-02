"""Check the YR-MPE citation contract.

Usage:
    python dev/citations/check_citations.py
    python dev/citations/check_citations.py --list

Exit code 1 means at least one config/plugin key is missing, a plugin still
hardcodes citation text, or the citation files did not parse cleanly.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "YR_MPE" / "runtime" / "citations.py"
CONFIG_PATH = REPO / "YR_MPE" / "config.json"
CITATION_DIR = REPO / "YR_MPE" / "citations"
PLUGIN_DIR = REPO / "YR_MPE" / "plugins"


def load_citations():
    spec = importlib.util.spec_from_file_location("yrmpe_citations", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def tool_keys(entry):
    keys = entry.get("citations", entry.get("citation", []))
    if isinstance(keys, str):
        return [keys]
    return list(keys)


def scan_plugins(library):
    hardcoded = []
    missing = []
    for path in sorted(PLUGIN_DIR.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            for target in node.targets:
                if not (
                    isinstance(target, ast.Attribute)
                    and isinstance(target.value, ast.Name)
                    and target.value.id == "self"
                ):
                    continue
                if target.attr == "citation" and isinstance(node.value, ast.List) and node.value.elts:
                    hardcoded.append("%s:%s" % (path.name, node.lineno))
                if target.attr == "citation_keys" and isinstance(node.value, ast.List):
                    for element in node.value.elts:
                        if isinstance(element, ast.Constant) and isinstance(element.value, str):
                            if element.value not in library.entries:
                                missing.append("%s:%s" % (path.name, element.value))
    return hardcoded, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list tools and citation keys")
    args = parser.parse_args()

    citations = load_citations()
    library = citations.CitationLibrary(str(CITATION_DIR))
    with open(CONFIG_PATH, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    unknown = []
    uncited = []
    exempt = {"MPIRun"}
    for entry in config:
        name = entry.get("name", "?")
        keys = tool_keys(entry)
        if args.list:
            print("%-22s %s" % (name, ", ".join(keys) if keys else "-"))
        if name not in exempt and not keys:
            uncited.append(name)
        for key in keys:
            if key not in library.entries:
                unknown.append((name, key))

    hardcoded, plugin_missing = scan_plugins(library)

    problems = 0
    print("Citation library: %s" % CITATION_DIR)
    print("  entries: %d" % len(library.entries))
    if library.errors:
        problems += 1
        print("  parse errors:")
        for error in library.errors:
            print("    - %s" % error)
    if library.duplicates:
        problems += 1
        print("  duplicate keys:")
        for duplicate in library.duplicates:
            print("    - %s" % duplicate)

    if unknown:
        problems += 1
        print("Unknown keys referenced by config.json:")
        for name, key in unknown:
            print("  - %s -> %s" % (name, key))

    if uncited:
        problems += 1
        print("Tools without citations:")
        for name in uncited:
            print("  - %s" % name)

    if hardcoded:
        problems += 1
        print("Plugins with hardcoded citation text:")
        for location in hardcoded:
            print("  - %s" % location)

    if plugin_missing:
        problems += 1
        print("Unknown keys in plugin citation_keys:")
        for location in plugin_missing:
            print("  - %s" % location)

    if not problems:
        print("OK: config, plugins and citation files agree.")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
