"""Tests for the CYAN-to-YR-MPE platform bundle assembler."""

from __future__ import annotations

import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "dev" / "release" / "build_platform_bundle.py"


def make_tar(path: Path, root: str, binary: str) -> None:
    payload = path.parent / binary
    payload.write_bytes(b"binary")
    with tarfile.open(path, "w:gz") as archive:
        archive.add(payload, arcname="%s/%s" % (root, binary))
    payload.unlink()


def make_zip(path: Path, root: str, binary: str) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("%s/%s" % (root, binary), b"binary")


class PlatformBundleTests(unittest.TestCase):
    def test_selects_target_and_neutral_assets(self):
        with tempfile.TemporaryDirectory() as tmp:
            assets = Path(tmp) / "assets"
            assets.mkdir()
            make_tar(assets / "diamond-2.2.9-Linux-x64-v1.tar.gz",
                     "diamond", "diamond")
            make_tar(assets / "diamond-2.2.9-Linux-arm64-v1.tar.gz",
                     "diamond", "diamond")
            make_zip(assets / "icytree-abc-Windows-x64-v1.zip",
                     "icytree", "index.html")
            output = Path(tmp) / "bundle.tar.gz"
            subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--assets-dir", str(assets),
                 "--platform", "Linux", "--arch", "x64",
                 "--format", "tar.gz", "--output", str(output),
                 "--tag", "v1"],
                check=True,
                stdout=subprocess.PIPE,
                text=True,
            )
            with tarfile.open(output) as archive:
                names = archive.getnames()
            self.assertIn("bin/diamond/diamond", names)
            self.assertIn("src/icytree/index.html", names)
            self.assertNotIn("bin/diamond/diamond.exe", names)

    def test_rejects_path_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            assets = Path(tmp) / "assets"
            assets.mkdir()
            bad = assets / "diamond-2.2.9-Linux-x64-v1.tar.gz"
            with tarfile.open(bad, "w:gz") as archive:
                info = tarfile.TarInfo("../escape")
                info.size = 0
                archive.addfile(info, None)
            result = subprocess.run(
                [sys.executable, str(SCRIPT),
                 "--assets-dir", str(assets),
                 "--platform", "Linux", "--arch", "x64",
                 "--format", "tar.gz", "--output", str(Path(tmp) / "x.tar.gz")],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Unsafe path", result.stderr + result.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
