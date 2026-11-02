"""Tests for the YR-MPE citation library and its config contract.

Run from the repository root:

    python dev/citations/test_citations.py

The parser tests load ``YR_MPE/runtime/citations.py`` directly so they do not
need PyQt. The contract tests inspect config.json and the plugin sources.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "YR_MPE" / "runtime" / "citations.py"
CONFIG_PATH = REPO / "YR_MPE" / "config.json"
CITATION_DIR = REPO / "YR_MPE" / "citations"
PLUGIN_DIR = REPO / "YR_MPE" / "plugins"


def _load_module():
    spec = importlib.util.spec_from_file_location("yrmpe_citations", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


citations = _load_module()


SAMPLE_BIB = r"""
@article{alpha2020,
  author = {Katoh, Kazutaka and Standley, Daron M.},
  title = {A {Nested {Brace}} Title with \& Ampersand},
  journal = {Molecular Biology and Evolution},
  year = {2020},
  volume = {30},
  number = {4},
  pages = {772-780},
  doi = {10.1000/example}
}

@book{beta2021,
  author = {S{\"o}ding, Johannes},
  title = {A Book},
  publisher = {Example Press},
  year = {2021}
}

@phdthesis{gamma2000,
  author = {van Dongen, Stijn},
  title = {Graph Clustering},
  school = {University of Utrecht},
  year = {2000}
}
"""


class BibTeXTests(unittest.TestCase):
    def test_parse_and_render(self):
        entries, errors = citations.parse_bibtex(SAMPLE_BIB, "sample.bib")
        self.assertEqual(errors, [])
        self.assertEqual([entry.key for entry in entries], ["alpha2020", "beta2021", "gamma2000"])

        alpha = entries[0]
        self.assertEqual(alpha.authors, "Katoh, K. & Standley, D. M.")
        text = alpha.format_text()
        self.assertIn("A Nested Brace Title with & Ampersand", text)
        self.assertIn("Molecular Biology and Evolution 30(4): 772–780", text)
        self.assertIn("https://doi.org/10.1000/example", text)

        html = alpha.format_html()
        self.assertIn("&amp;", html)
        self.assertIn("<i>Molecular Biology and Evolution</i>", html)
        self.assertIn("<b>30</b>(4)", html)

    def test_accents(self):
        entries, _ = citations.parse_bibtex(SAMPLE_BIB)
        self.assertIn("Söding", entries[1].authors)
        self.assertEqual(citations.clean_latex(r"S{\'e}bastien"), "Sébastien")

    def test_book_and_thesis(self):
        entries, _ = citations.parse_bibtex(SAMPLE_BIB)
        self.assertEqual(
            entries[1].format_text(),
            "Söding, J. (2021). A Book. Example Press.",
        )
        self.assertIn("PhD thesis, University of Utrecht", entries[2].format_text())

    def test_duplicate_and_missing(self):
        text = SAMPLE_BIB + SAMPLE_BIB
        library_dir = tempfile.mkdtemp()
        path = os.path.join(library_dir, "dup.bib")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        library = citations.CitationLibrary(library_dir)
        self.assertEqual(len(library.entries), 3)
        self.assertEqual(len(library.duplicates), 3)
        found, missing = library.get_many(["alpha2020", "nope"])
        self.assertEqual([entry.key for entry in found], ["alpha2020"])
        self.assertEqual(missing, ["nope"])


class RISTests(unittest.TestCase):
    def test_parse_ris(self):
        text = """TY  - JOUR
AU  - Katoh, Kazutaka
AU  - Standley, Daron M.
TI  - MAFFT version 7
JO  - Molecular Biology and Evolution
PY  - 2013
VL  - 30
IS  - 4
SP  - 772
EP  - 780
DO  - 10.1093/molbev/mst010
ER  -
"""
        entries, errors = citations.parse_ris(text, "sample.ris")
        self.assertEqual(errors, [])
        self.assertEqual(len(entries), 1)
        entry = entries[0]
        self.assertEqual(entry.entry_type, "article")
        self.assertEqual(entry.authors, "Katoh, K. & Standley, D. M.")
        self.assertIn("30(4): 772–780", entry.format_text())
        self.assertIn("https://doi.org/10.1093/molbev/mst010", entry.format_text())


class LibraryTests(unittest.TestCase):
    def test_reload_on_change(self):
        directory = tempfile.mkdtemp()
        path = os.path.join(directory, "software.bib")
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(SAMPLE_BIB)
        library = citations.CitationLibrary(directory)
        self.assertIn("alpha2020", library.entries)

        with open(path, "w", encoding="utf-8") as handle:
            handle.write(SAMPLE_BIB.replace("A {Nested {Brace}} Title", "A New Title"))
        future = time.time_ns() + 2_000_000_000
        os.utime(path, ns=(future, future))

        self.assertTrue(library.reload_if_changed())
        self.assertIn("A New Title", library.get("alpha2020").format_text())

    def test_missing_directory(self):
        library = citations.CitationLibrary(os.path.join(tempfile.mkdtemp(), "absent"))
        self.assertEqual(library.entries, {})
        self.assertEqual(library.errors, [])


class ConfigContractTests(unittest.TestCase):
    def setUp(self):
        self.library = citations.CitationLibrary(str(CITATION_DIR))
        with open(CONFIG_PATH, "r", encoding="utf-8") as handle:
            self.config = json.load(handle)

    def test_library_is_clean(self):
        self.assertEqual(self.library.errors, [])
        self.assertEqual(self.library.duplicates, [])

    def test_every_config_key_resolves(self):
        missing = []
        for tool in self.config:
            keys = tool.get("citations", tool.get("citation", []))
            if isinstance(keys, str):
                keys = [keys]
            for key in keys:
                if key not in self.library.entries:
                    missing.append((tool.get("name"), key))
        self.assertEqual(missing, [], "config.json references unknown citation keys")

    def test_external_tools_have_citations(self):
        exempt = {"MPIRun"}  # launcher, not a citable tool
        missing = [
            tool.get("name")
            for tool in self.config
            if tool.get("name") not in exempt and not tool.get("citations")
        ]
        self.assertEqual(missing, [], "external tools without citations")

    def test_plugins_have_no_hardcoded_citations(self):
        offenders = []
        for path in sorted(PLUGIN_DIR.glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Assign):
                    continue
                for target in node.targets:
                    if (
                        isinstance(target, ast.Attribute)
                        and target.attr == "citation"
                        and isinstance(target.value, ast.Name)
                        and target.value.id == "self"
                        and isinstance(node.value, ast.List)
                        and node.value.elts
                    ):
                        offenders.append("%s:%s" % (path.name, node.lineno))
        self.assertEqual(offenders, [], "citation text must live in citations/*.bib")

    def test_plugin_citation_keys_resolve(self):
        missing = []
        for path in sorted(PLUGIN_DIR.glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.Assign):
                    continue
                for target in node.targets:
                    if not (
                        isinstance(target, ast.Attribute)
                        and target.attr == "citation_keys"
                        and isinstance(target.value, ast.Name)
                        and target.value.id == "self"
                    ):
                        continue
                    if not isinstance(node.value, ast.List):
                        continue
                    for element in node.value.elts:
                        if isinstance(element, ast.Constant) and isinstance(element.value, str):
                            if element.value not in self.library.entries:
                                missing.append("%s:%s" % (path.name, element.value))
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
