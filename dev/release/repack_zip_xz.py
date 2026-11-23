"""Convert a Windows .zip bundle to .tar.xz for the Gitee mirror."""

from __future__ import annotations

import sys
import tarfile
import time
import zipfile

SKIP = (
    "bin/diamond/", "bin/famsa/", "bin/mcl/",
    "bin/FastTree/", "bin/fastme/",
)


def main(source: str, target: str) -> None:
    with zipfile.ZipFile(source) as src, tarfile.open(
            target, "w:xz", preset=9) as dst:
        for info in src.infolist():
            if info.is_dir() or info.filename.startswith(SKIP):
                continue
            data = src.read(info)
            member = tarfile.TarInfo(info.filename)
            member.size = len(data)
            member.mtime = int(time.time())
            member.mode = 0o755 if info.filename.lower().endswith(
                (".exe", ".bat", ".dll")) else 0o644
            import io
            dst.addfile(member, io.BytesIO(data))
    print(target)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
