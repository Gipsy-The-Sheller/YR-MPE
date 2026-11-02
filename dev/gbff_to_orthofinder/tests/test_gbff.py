import csv
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

from gbff import parse_genbank_file, parse_folder, write_dataset  # noqa: E402
from orthogroups import (  # noqa: E402
    derive_single_copy_ids,
    read_cds_map,
    read_orthogroups,
    write_single_copy_ogs,
)

GB_DIR = Path("E:/Codes/gf_mito_b/analysis/broad_pmc/gb_cache")
GB_FILES = ["HM347959.gb", "JN637766.gb", "KF150644.gb"]


class GBFFParserTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not all((GB_DIR / name).exists() for name in GB_FILES):
            raise unittest.SkipTest(f"test GB files not found under {GB_DIR}")

    def test_parse_single_file(self):
        records = parse_genbank_file(GB_DIR / GB_FILES[0], default_table=11)
        self.assertGreater(len(records), 50)
        self.assertEqual(len({r.protein_id for r in records}), len(records))
        for record in records:
            self.assertTrue(record.aa_seq)
            self.assertTrue(record.cds_seq)
            self.assertEqual(len(record.cds_seq) % 3, 0)
            self.assertEqual(record.transl_table, 11)

    def _parse_test_files(self):
        with tempfile.TemporaryDirectory() as tmp_in:
            input_dir = Path(tmp_in)
            for name in GB_FILES:
                shutil.copy(GB_DIR / name, input_dir / name)
            return parse_folder(input_dir, default_table=11)

    def test_write_dataset(self):
        records = self._parse_test_files()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            summary = write_dataset(records, out)
            self.assertEqual(len(summary), len(GB_FILES))
            self.assertTrue((out / "proteomes").is_dir())
            self.assertTrue((out / "cds").is_dir())
            self.assertTrue((out / "cds_map.tsv").exists())
            mapping = read_cds_map(out / "cds_map.tsv")
            self.assertEqual(len(mapping), sum(summary.values()))

    def test_multirecord_file_as_taxa(self):
        with tempfile.TemporaryDirectory() as tmp:
            multi = Path(tmp) / "multi.gb"
            with multi.open("wb") as handle:
                for name in GB_FILES[:2]:
                    handle.write((GB_DIR / name).read_bytes())
            records = parse_genbank_file(multi, default_table=11, record_as_species=True)
            species = {record.species for record in records}
            self.assertEqual(len(species), 2)
            self.assertGreater(len(records), 100)
            # parse_folder must also accept a file path, not only a directory.
            records2 = parse_folder(multi, default_table=11, record_as_species=True)
            self.assertEqual(len(records2), len(records))

    def test_single_copy_backmap(self):
        records = self._parse_test_files()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_dataset(records, out)
            mapping = read_cds_map(out / "cds_map.tsv")
            by_species: dict[str, list[str]] = {}
            for record in mapping.values():
                by_species.setdefault(record["species"], []).append(record["protein_id"])
            species = sorted(by_species)
            self.assertGreaterEqual(len(species), 3)

            og_dir = out / "orthofinder" / "Orthogroups"
            og_dir.mkdir(parents=True)
            with (og_dir / "Orthogroups.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["Orthogroup"] + species)
                writer.writerow(["OG0000001"] + [by_species[s][0] for s in species])
            (og_dir / "Orthogroups_SingleCopyOrthologues.txt").write_text(
                "OG0000001\n", encoding="utf-8"
            )

            orthogroups = read_orthogroups(og_dir / "Orthogroups.tsv")
            self.assertEqual(len(derive_single_copy_ids(orthogroups)), 1)
            result = write_single_copy_ogs(
                orthogroups,
                {"OG0000001"},
                mapping,
                out / "single_copy_ogs",
            )
            self.assertEqual(result["orthogroups_written"], 1)
            self.assertTrue((out / "single_copy_ogs" / "OG0000001.faa").exists())
            self.assertTrue((out / "single_copy_ogs" / "OG0000001.cds.fna").exists())
            self.assertTrue((out / "single_copy_ogs" / "OG0000001.map.tsv").exists())


if __name__ == "__main__":
    unittest.main()
