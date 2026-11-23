"""Mirror a GitHub release to Gitee from a Gitee Go pipeline.

The Gitee runner is Linux; it downloads the assets built by GitHub Actions,
converts the large bundles to .tar.xz (Gitee's 100 MB attachment limit) and
publishes them to the matching Gitee release.  Requires ``GITEE_TOKEN`` in
the repository's Gitee Go environment variables.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import requests

GITHUB_REPO = "Gipsy-The-Sheller/YR-MPE"
GITEE_REPO = os.environ.get("GITEE_REPO", "ZJXMolls/YR-MPE")
API = "https://gitee.com/api/v5"
TOKEN = os.environ.get("GITEE_TOKEN", "").strip()

sys.path.insert(0, str(Path(__file__).resolve().parent))
from repack_xz import main as repack_tar_xz  # noqa: E402
from repack_zip_xz import main as repack_zip_xz  # noqa: E402


def fail(message: str) -> None:
    print("ERROR:", message, file=sys.stderr)
    raise SystemExit(1)


def tag_name() -> str:
    for key in ("GITEE_TAG", "GITEE_REF_NAME", "GITEE_REF"):
        value = os.environ.get(key, "").strip()
        if value.startswith("refs/tags/"):
            return value.split("/", 2)[-1]
        if value.startswith("v"):
            return value
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--exact-match"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except subprocess.CalledProcessError:
        response = requests.get(
            "https://api.github.com/repos/%s/releases/latest" % GITHUB_REPO,
            timeout=60, headers={"Accept": "application/vnd.github+json"})
        if response.status_code == 200:
            return response.json()["tag_name"]
        fail("cannot determine the tag; run the pipeline on a tag")


def gitee(method: str, path: str, **kwargs):
    headers = kwargs.pop("headers", {})
    headers["Authorization"] = "token %s" % TOKEN
    response = requests.request(
        method, API + path, headers=headers, timeout=120, **kwargs)
    if response.status_code >= 400:
        fail("%s %s -> %s %s" % (
            method, path, response.status_code, response.text[:500]))
    return response


def ensure_release(tag: str, body: str) -> int:
    response = requests.get(
        API + "/repos/%s/releases/tags/%s" % (GITEE_REPO, tag),
        headers={"Authorization": "token %s" % TOKEN}, timeout=60)
    data = response.json() if response.status_code == 200 else None
    if isinstance(data, dict) and data.get("id"):
        return int(data["id"])
    payload = {
        "tag_name": tag,
        "target_commitish": os.environ.get("GITEE_TARGET", "master"),
        "name": "YPhylo Workbench %s" % tag.lstrip("v"),
        "body": body or "Mirror of the GitHub release.",
        "prerelease": False,
    }
    created = gitee("POST", "/repos/%s/releases" % GITEE_REPO, json=payload)
    return int(created.json()["id"])


def existing_names(release_id: int) -> set[str]:
    response = gitee(
        "GET", "/repos/%s/releases/%d/attach_files" % (GITEE_REPO, release_id),
        params={"per_page": 100})
    return {row["name"] for row in response.json()}


def upload(release_id: int, path: Path) -> None:
    with path.open("rb") as handle:
        gitee(
            "POST",
            "/repos/%s/releases/%d/attach_files" % (GITEE_REPO, release_id),
            files={"file": (path.name, handle)},
        )
    print("uploaded", path.name)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if not TOKEN:
        fail("GITEE_TOKEN is not configured in the Gitee Go environment")
    tag = tag_name()
    print("mirroring", tag)

    github = requests.get(
        "https://api.github.com/repos/%s/releases/tags/%s"
        % (GITHUB_REPO, tag), timeout=60,
        headers={"Accept": "application/vnd.github+json"})
    if github.status_code != 200:
        fail("GitHub release %s not found: %s" % (tag, github.status_code))
    release = github.json()
    release_id = ensure_release(tag, release.get("body", ""))
    present = existing_names(release_id)

    work = Path(tempfile.mkdtemp(prefix="gitee_mirror_"))
    published = []
    for asset in release.get("assets", []):
        name = asset["name"]
        if name == "SHA256SUMS":
            continue
        target_name = name
        if name.endswith(".zip") and "-src." not in name:
            target_name = name[:-4] + ".tar.xz"
        elif name.endswith(".tar.gz") and "-src." not in name:
            target_name = name[:-7] + ".tar.xz"
        local = work / name
        print("downloading", name)
        with requests.get(asset["browser_download_url"], stream=True,
                          timeout=300) as response:
            response.raise_for_status()
            with local.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    handle.write(chunk)
        if target_name != name:
            converted = work / target_name
            if name.endswith(".zip"):
                repack_zip_xz(str(local), str(converted))
            else:
                repack_tar_xz(str(local), str(converted))
            local.unlink()
            local = converted
        if target_name in present:
            print("exists, skipping", target_name)
        else:
            upload(release_id, local)
        published.append(local)

    sums = work / "SHA256SUMS"
    sums.write_text(
        "".join("%s *%s\n" % (sha256(path), path.name) for path in published),
        encoding="utf-8")
    if "SHA256SUMS" in present:
        rows = gitee(
            "GET",
            "/repos/%s/releases/%d/attach_files" % (GITEE_REPO, release_id),
            params={"per_page": 100}).json()
        for row in rows:
            if row["name"] == "SHA256SUMS":
                gitee(
                    "DELETE",
                    "/repos/%s/releases/%d/attach_files/%s"
                    % (GITEE_REPO, release_id, row["id"]))
    upload(release_id, sums)
    print("done")


if __name__ == "__main__":
    main()
