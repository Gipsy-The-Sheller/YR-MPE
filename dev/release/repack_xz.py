"""Repack a platform bundle as .tar.xz, dropping duplicate OrthoFinder deps.

OrthoFinder ships its own diamond/famsa/mcl under ``orthofinder/bin``; the
top-level copies are not referenced by any plugin and make the archive exceed
Gitee's 100 MB attachment limit.
"""

from __future__ import annotations

import sys
import tarfile

SKIP = (
    "bin/diamond/", "bin/famsa/", "bin/mcl/",
    "bin/FastTree/", "bin/fastme/",
)


def main(source: str, target: str) -> None:
    with tarfile.open(source, "r:*") as src, tarfile.open(
            target, "w:xz", preset=9) as dst:
        for member in src:
            if member.name.startswith(SKIP):
                continue
            handle = src.extractfile(member) if member.isfile() else None
            dst.addfile(member, handle)
    print(target)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
