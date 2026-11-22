"""YR-MPE first-run initializer.

The initializer has three jobs:

1. install the bundled platform binaries from a GitHub or Gitee mirror (or a
   user-supplied URL), or a single tool from the pinned CYAN release;
2. detect bundled and PATH executables and write ``YR_MPE/config.json``;
3. write ``settings.json`` for the YRTools host.

The binary archive is expected to contain a ``bin`` directory with the
subdirectories used by ``YR_MPE/bin`` (``mafft``, ``iqtree3``, ``mb``,
``aster`` ...).  Both a flat archive and one wrapped in ``YR_MPE/bin`` are
accepted.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

from PyQt5.QtCore import QThread, Qt, pyqtSignal
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

try:
    import requests
except ImportError:  # pragma: no cover - requests is a declared dependency
    requests = None

CYAN_REPO = "Gipsy-The-Sheller/CYAN"

# CYAN asset prefix -> menu label.  The archive root is copied as-is, so the
# layout under YR_MPE/bin matches config.template.json.
CYAN_TOOLS = [
    ("mafft", "MAFFT"),
    ("iqtree3", "IQ-TREE 3"),
    ("modeltest-ng", "ModelTest-NG"),
    ("raxml-ng", "RAxML-NG"),
    ("mrbayes", "MrBayes"),
    ("pbmpi", "PhyloBayes-MPI"),
    ("aster", "ASTER v1"),
    ("aster2", "ASTER2"),
    ("lsd2", "LSD2"),
    ("trimal", "TrimAl"),
    ("ygblocks", "YGBlocks"),
    ("clustalo", "Clustal Omega"),
    ("muscle5", "MUSCLE 5"),
    ("mpboot", "MPBoot"),
    ("decenttree", "DecentTree"),
    ("fasttree", "FastTree"),
    ("fastme", "FastME"),
    ("paml", "PAML"),
    ("hyphy", "HyPhy"),
    ("diamond", "DIAMOND"),
    ("mcl", "MCL"),
    ("famsa", "FAMSA"),
    ("orthofinder3", "OrthoFinder 3"),
    ("treerecs", "TreeRecs"),
    ("bmge", "BMGE"),
    ("icytree", "IcyTree"),
]

CYAN_SRC_TOOLS = {"icytree"}


def cyan_target():
    """Return the CYAN ``(platform, arch)`` label for this machine."""
    system = platform.system()
    platform_name = {"Windows": "Windows", "Linux": "Linux",
                     "Darwin": "macOS"}.get(system, system)
    machine = platform.machine().lower()
    arch = "arm64" if machine in ("arm64", "aarch64") else "x64"
    if platform_name == "Windows" and arch == "arm64":
        # PyQt5 has no Windows-arm64 wheel; the x64 build runs under emulation.
        arch = "x64"
    return platform_name, arch


def cyan_release_tag():
    """Read the release pinned in CYAN_RELEASE, else ask GitHub."""
    pin = Path(__file__).resolve().parent / "CYAN_RELEASE"
    if pin.is_file():
        tag = pin.read_text(encoding="utf-8").strip()
        if tag:
            return tag
    if requests is None:
        raise RuntimeError("requests is required to resolve the CYAN release")
    response = requests.get(
        "https://api.github.com/repos/%s/releases/latest" % CYAN_REPO,
        timeout=30,
        headers={"Accept": "application/vnd.github+json"},
    )
    response.raise_for_status()
    return response.json()["tag_name"]


def load_cyan_manifest():
    path = Path(__file__).resolve().parent / "cyan_tools.json"
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_cyan_asset(tool, platform_name, arch, tag=None):
    """Return ``(asset_name, download_url)`` from the pinned CYAN manifest.

    The manifest records the tool versions published in the CYAN release, so
    the Initializer does not need the GitHub API (which is often unreachable
    from mainland China).  IcyTree is one platform-neutral package built in
    the Windows job.
    """
    manifest = load_cyan_manifest()
    release = tag or manifest.get("release") or cyan_release_tag()
    entry = manifest.get("tools", {}).get(tool)
    if not entry:
        raise RuntimeError("CYAN manifest has no entry for %s" % tool)
    target = "Windows-x64" if tool in CYAN_SRC_TOOLS else "%s-%s" % (
        platform_name, arch)
    if target not in entry.get("targets", []):
        raise RuntimeError(
            "CYAN %s has no package for %s on %s"
            % (release, tool, target)
        )
    extension = "zip" if target.startswith("Windows") else "tar.gz"
    name = "%s-%s-%s-%s.%s" % (
        tool, entry["version"], target, release, extension)
    url = "https://github.com/%s/releases/download/%s/%s" % (
        CYAN_REPO, release, name)
    return name, url


class LicenseDialog(QDialog):
    """Display the GPL-3.0 notice."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("License Information")
        self.resize(620, 420)
        layout = QVBoxLayout()

        title = QLabel("GNU General Public License v3.0")
        font = QFont()
        font.setPointSize(16)
        font.setBold(True)
        title.setFont(font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        text = QTextEdit()
        text.setReadOnly(True)
        text.setText(
            "This program is free software: you can redistribute it and/or "
            "modify it under the terms of the GNU General Public License as "
            "published by the Free Software Foundation, either version 3 of "
            "the License, or (at your option) any later version.\n\n"
            "This program is distributed in the hope that it will be useful, "
            "but WITHOUT ANY WARRANTY; without even the implied warranty of "
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the "
            "GNU General Public License for more details.\n\n"
            "You should have received a copy of the GNU General Public "
            "License along with this program. If not, see "
            "https://www.gnu.org/licenses/."
        )
        layout.addWidget(text)

        close = QPushButton("Close")
        close.clicked.connect(self.close)
        layout.addWidget(close)
        self.setLayout(layout)


class ManualConfigDialog(QDialog):
    """Manual executable-path editor."""

    def __init__(self, tools_data, parent=None):
        super().__init__(parent)
        self.tools_data = tools_data
        self.tool_paths = {}
        self.setWindowTitle("Manual Configuration")
        self.resize(760, 460)

        layout = QVBoxLayout()
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Software Name", "Path", "Select"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Fixed)
        self.table.setColumnWidth(2, 100)

        self.table.setRowCount(len(self.tools_data))
        for row, tool in enumerate(self.tools_data):
            name_item = QTableWidgetItem(tool["name"])
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(row, 0, name_item)
            path_item = QTableWidgetItem(tool.get("path", ""))
            self.table.setItem(row, 1, path_item)
            browse = QPushButton("Browse")
            browse.clicked.connect(lambda checked, r=row: self.browse_path(r))
            self.table.setCellWidget(row, 2, browse)
            self.tool_paths[row] = path_item

        layout.addWidget(self.table)
        buttons = QHBoxLayout()
        ok = QPushButton("OK")
        ok.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        buttons.addWidget(ok)
        buttons.addWidget(cancel)
        layout.addLayout(buttons)
        self.setLayout(layout)

    def browse_path(self, row):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select %s executable" % self.tools_data[row]["name"]
        )
        if path:
            self.tool_paths[row].setText(path)

    def get_configured_tools(self):
        configured = []
        for row, tool in enumerate(self.tools_data):
            new_tool = tool.copy()
            new_tool["path"] = self.tool_paths[row].text().strip()
            configured.append(new_tool)
        return configured


class BinaryDownloadThread(QThread):
    """Download one archive without blocking the Qt event loop."""

    progress = pyqtSignal(int, int)
    status = pyqtSignal(str)
    finished = pyqtSignal(str)
    error = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(self, url, destination, parent=None):
        super().__init__(parent)
        self.url = url
        self.destination = destination
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        try:
            if self.url.startswith("file://"):
                source = self.url[7:]
                shutil.copyfile(source, self.destination)
                self.finished.emit(self.destination)
                return
            if os.path.exists(self.url):
                shutil.copyfile(self.url, self.destination)
                self.finished.emit(self.destination)
                return
            if requests is None:
                raise RuntimeError(
                    "The 'requests' package is required to download binaries."
                )

            self.status.emit("Connecting...")
            headers = {"User-Agent": "YR-MPE-Initializer"}
            with requests.get(
                self.url,
                stream=True,
                timeout=(20, 120),
                headers=headers,
                allow_redirects=True,
            ) as response:
                response.raise_for_status()
                total = int(response.headers.get("Content-Length", 0) or 0)
                downloaded = 0
                with open(self.destination, "wb") as handle:
                    for chunk in response.iter_content(chunk_size=256 * 1024):
                        if self._cancelled:
                            self.status.emit("Download cancelled.")
                            self.cancelled.emit()
                            return
                        if not chunk:
                            continue
                        handle.write(chunk)
                        downloaded += len(chunk)
                        self.progress.emit(downloaded, total)
            self.finished.emit(self.destination)
        except Exception as exc:  # noqa: BLE001 - reported to the GUI
            self.error.emit(str(exc))


def _safe_extract_zip(archive, destination):
    destination = Path(destination).resolve()
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            target = (destination / member.filename).resolve()
            if destination != target and destination not in target.parents:
                raise RuntimeError("Unsafe archive member: %s" % member.filename)
        bundle.extractall(destination)


def _safe_extract_tar(archive, destination):
    destination = Path(destination).resolve()
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            target = (destination / member.name).resolve()
            if destination != target and destination not in target.parents:
                raise RuntimeError("Unsafe archive member: %s" % member.name)
        bundle.extractall(destination)


def _find_binary_root(root: Path):
    """Locate the directory that should be copied into ``YR_MPE/bin``."""
    root = Path(root)
    candidates = []
    for name in ("bin", "bins", "binaries"):
        candidates.extend(path for path in root.rglob(name) if path.is_dir())
    if (root / "mafft").is_dir() and (root / "iqtree3").is_dir():
        candidates.append(root)

    def score(path: Path) -> int:
        markers = ("mafft", "iqtree3", "mb", "aster", "trimAl", "pbmpi-1.9")
        return sum(1 for marker in markers if (path / marker).exists())

    if not candidates:
        return None
    return max(candidates, key=score)


class Initializer(QWidget):
    """YR-MPE first-run wizard."""

    def __init__(self):
        super().__init__()
        self.plugin_path = Path(__file__).resolve().parent
        self.tools_data = []
        self.configured_tools = []
        self.binary_sources = []
        platform_name, arch = cyan_target()
        self.binary_archive = "YR-MPE-bins-%s-%s.%s" % (
            platform_name, arch, "zip" if platform_name == "Windows" else "tar.gz")
        self.download_thread = None
        self._download_path = None
        self._cyan_mode = False
        self._cyan_tool = ""
        self.init_ui()
        self.load_manifest()
        self.load_templates()

    # ------------------------------------------------------------------ UI
    def init_ui(self):
        self.setWindowTitle("YR-MPE Initializer")
        self.setMinimumSize(720, 720)

        layout = QVBoxLayout()
        title = QLabel("Welcome to the YR-MPE Initializer")
        font = QFont()
        font.setBold(True)
        title.setFont(font)
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        license_group = QGroupBox("License")
        license_layout = QHBoxLayout()
        license_layout.addWidget(QLabel("View the GPL-3.0 license:"))
        self.license_button = QPushButton("View License")
        self.license_button.clicked.connect(self.show_license)
        license_layout.addWidget(self.license_button)
        license_group.setLayout(license_layout)
        layout.addWidget(license_group)

        path_group = QGroupBox("Installation Path")
        path_layout = QHBoxLayout()
        path_layout.addWidget(QLabel("YR-MPE Path:"))
        self.path_edit = QLineEdit("YR-MPE")
        self.path_button = QPushButton("Browse")
        self.path_button.clicked.connect(self.browse_path)
        path_layout.addWidget(self.path_edit)
        path_layout.addWidget(self.path_button)
        path_group.setLayout(path_layout)
        layout.addWidget(path_group)

        binary_group = QGroupBox("Bundled platform binaries")
        binary_layout = QVBoxLayout()

        source_row = QHBoxLayout()
        source_row.addWidget(QLabel("Source:"))
        self.source_combo = QComboBox()
        self.source_combo.currentIndexChanged.connect(self.on_source_changed)
        source_row.addWidget(self.source_combo, 1)
        self.download_button = QPushButton("Download and install binaries")
        self.download_button.clicked.connect(self.start_download)
        source_row.addWidget(self.download_button)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel_download)
        source_row.addWidget(self.cancel_button)
        binary_layout.addLayout(source_row)

        url_row = QHBoxLayout()
        url_row.addWidget(QLabel("Archive URL:"))
        self.url_edit = QLineEdit()
        self.url_edit.setPlaceholderText("https://.../YR-MPE-bins-win64.zip")
        url_row.addWidget(self.url_edit, 1)
        binary_layout.addLayout(url_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        binary_layout.addWidget(self.progress_bar)
        self.status_label = QLabel(
            "Choose a mirror, or paste a custom URL for the platform binary archive."
        )
        self.status_label.setWordWrap(True)
        binary_layout.addWidget(self.status_label)
        binary_group.setLayout(binary_layout)
        layout.addWidget(binary_group)

        cyan_group = QGroupBox("Single tool from CYAN")
        cyan_layout = QVBoxLayout()
        cyan_row = QHBoxLayout()
        cyan_row.addWidget(QLabel("Tool:"))
        self.cyan_tool_combo = QComboBox()
        for key, label in CYAN_TOOLS:
            self.cyan_tool_combo.addItem(label, key)
        cyan_row.addWidget(self.cyan_tool_combo, 1)
        self.cyan_download_button = QPushButton("Download tool")
        self.cyan_download_button.clicked.connect(self.start_cyan_download)
        cyan_row.addWidget(self.cyan_download_button)
        cyan_layout.addLayout(cyan_row)
        try:
            cyan_target_label = "%s-%s" % cyan_target()
            cyan_release_label = cyan_release_tag()
        except Exception as exc:
            cyan_target_label = "unknown"
            cyan_release_label = "unresolved (%s)" % exc
        cyan_note = QLabel(
            "Target %s · release %s. A single tool can be installed without "
            "the full platform bundle."
            % (cyan_target_label, cyan_release_label)
        )
        cyan_note.setWordWrap(True)
        cyan_note.setStyleSheet("color: palette(mid);")
        cyan_layout.addWidget(cyan_note)
        cyan_group.setLayout(cyan_layout)
        layout.addWidget(cyan_group)

        config_group = QGroupBox("Tool Configuration")
        config_layout = QVBoxLayout()
        auto_row = QHBoxLayout()
        auto_row.addWidget(QLabel("Detect bundled tools and tools in PATH:"))
        self.auto_config_button = QPushButton("Auto Config")
        self.auto_config_button.clicked.connect(self.auto_config)
        auto_row.addWidget(self.auto_config_button)
        config_layout.addLayout(auto_row)

        manual_row = QHBoxLayout()
        manual_row.addWidget(QLabel("Set executable paths manually:"))
        self.manual_config_button = QPushButton("Manual Config")
        self.manual_config_button.clicked.connect(self.manual_config)
        manual_row.addWidget(self.manual_config_button)
        config_layout.addLayout(manual_row)
        config_group.setLayout(config_layout)
        layout.addWidget(config_group)

        self.finish_button = QPushButton("Finish && Install")
        self.finish_button.clicked.connect(self.finish_install)
        layout.addWidget(self.finish_button)
        self.setLayout(layout)

    # -------------------------------------------------------------- loading
    def load_manifest(self):
        manifest_path = self.plugin_path / "binary_sources.json"
        sources = []
        archive = self.binary_archive
        if manifest_path.exists():
            try:
                data = json.loads(manifest_path.read_text(encoding="utf-8"))
                archives = data.get("archives") or {}
                if archives:
                    archive = archives.get(
                        "%s-%s" % cyan_target(), archive)
                else:
                    archive = data.get("archive", archive)
                sources = data.get("sources", [])
            except Exception as exc:
                QMessageBox.warning(
                    self, "Binary sources", "Cannot read binary_sources.json:\n%s" % exc
                )
        if not sources:
            sources = [
                {
                    "name": "GitHub",
                    "url": "https://github.com/Gipsy-The-Sheller/YR-MPE/"
                    "releases/download/v{version}/{archive}",
                },
                {
                    "name": "Gitee",
                    "url": "https://gitee.com/ZJXMolls/YR-MPE/"
                    "releases/download/v{version}/{archive}",
                },
            ]
        self.binary_archive = archive
        self.binary_sources = sources
        self.source_combo.blockSignals(True)
        self.source_combo.clear()
        for source in sources:
            self.source_combo.addItem(source.get("name", "Mirror"), source)
        self.source_combo.addItem("Custom URL", None)
        self.source_combo.blockSignals(False)
        self.on_source_changed(0)

    def load_templates(self):
        template = self.plugin_path / "config.template.json"
        if not template.exists():
            return
        with open(template, "r", encoding="utf-8") as handle:
            self.tools_data = json.load(handle)

        existing = {}
        config_path = self.plugin_path / "YR_MPE" / "config.json"
        if config_path.exists():
            try:
                for tool in json.loads(config_path.read_text(encoding="utf-8")):
                    existing[tool.get("name", "")] = tool.get("path", "")
            except Exception:
                existing = {}
        for tool in self.tools_data:
            if tool.get("name") in existing and existing[tool["name"]]:
                tool["path"] = existing[tool["name"]]

    def _version(self):
        try:
            from YR_MPE import __version__

            return __version__
        except Exception:
            return "1.2.1"

    def _resolved_url(self, template):
        return (template
                .replace("{version}", self._version())
                .replace("{archive}", self.binary_archive))

    # ------------------------------------------------------------- events
    def show_license(self):
        LicenseDialog(self).exec_()

    def browse_path(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select Installation Path", str(self.plugin_path)
        )
        if path:
            self.path_edit.setText(os.path.relpath(path, str(self.plugin_path)))

    def on_source_changed(self, index):
        source = self.source_combo.itemData(index)
        if source is None:
            self.url_edit.setReadOnly(False)
            self.url_edit.setPlaceholderText(
                "Paste a direct URL to a .zip or .tar.gz binary archive"
            )
            return
        self.url_edit.setReadOnly(False)
        self.url_edit.setText(self._resolved_url(source.get("url", "")))
        self.url_edit.setCursorPosition(0)

    # ----------------------------------------------------------- download
    def start_cyan_download(self):
        if self.download_thread is not None and self.download_thread.isRunning():
            return
        tool = self.cyan_tool_combo.currentData()
        platform_name, arch = cyan_target()
        try:
            tag = cyan_release_tag()
            name, url = resolve_cyan_asset(tool, platform_name, arch, tag)
        except Exception as exc:
            QMessageBox.critical(
                self, "CYAN download",
                "Cannot resolve the CYAN package:\n%s" % exc,
            )
            return

        suffix = ".zip" if name.endswith(".zip") else ".tar.gz"
        handle = tempfile.NamedTemporaryFile(
            prefix="yrmpe_cyan_", suffix=suffix, delete=False
        )
        self._download_path = handle.name
        handle.close()

        self._cyan_mode = True
        self._cyan_tool = tool
        self.download_button.setEnabled(False)
        self.cyan_download_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.status_label.setText("Downloading %s (%s)" % (name, tag))
        self.download_thread = BinaryDownloadThread(url, self._download_path, self)
        self.download_thread.progress.connect(self.on_download_progress)
        self.download_thread.status.connect(self.status_label.setText)
        self.download_thread.finished.connect(self.on_download_finished)
        self.download_thread.error.connect(self.on_download_error)
        self.download_thread.cancelled.connect(self.on_download_cancelled)
        self.download_thread.start()

    def start_download(self):
        if self.download_thread is not None and self.download_thread.isRunning():
            return
        url = self.url_edit.text().strip()
        if not url:
            QMessageBox.warning(self, "Download", "Enter an archive URL first.")
            return
        suffix = ".zip"
        lower = url.lower().split("?")[0]
        for candidate in (".tar.gz", ".tgz", ".zip"):
            if lower.endswith(candidate):
                suffix = candidate
                break
        handle = tempfile.NamedTemporaryFile(
            prefix="yrmpe_bins_", suffix=suffix, delete=False
        )
        self._download_path = handle.name
        handle.close()

        self.download_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.status_label.setText("Downloading %s" % url)
        self.download_thread = BinaryDownloadThread(url, self._download_path, self)
        self.download_thread.progress.connect(self.on_download_progress)
        self.download_thread.status.connect(self.status_label.setText)
        self.download_thread.finished.connect(self.on_download_finished)
        self.download_thread.error.connect(self.on_download_error)
        self.download_thread.cancelled.connect(self.on_download_cancelled)
        self.download_thread.start()

    def cancel_download(self):
        if self.download_thread is not None:
            self.download_thread.cancel()
            self.status_label.setText("Cancelling...")

    def on_download_progress(self, downloaded, total):
        if total > 0:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(downloaded)
        else:
            self.progress_bar.setRange(0, 0)
        if total > 0:
            self.status_label.setText(
                "Downloaded %.1f / %.1f MB"
                % (downloaded / 1048576.0, total / 1048576.0)
            )
        else:
            self.status_label.setText("Downloaded %.1f MB" % (downloaded / 1048576.0))

    def on_download_finished(self, archive_path):
        self.download_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.status_label.setText("Extracting binaries...")
        QApplication.processEvents()
        try:
            if self._cyan_mode:
                self.install_cyan_tool(archive_path, self._cyan_tool)
                self.status_label.setText(
                    "CYAN package installed. Detecting tools..."
                )
                self.auto_config()
                QMessageBox.information(
                    self,
                    "Tool installed",
                    "The CYAN package for %s was installed into YR_MPE."
                    % self._cyan_tool,
                )
            else:
                self.install_bundle(archive_path)
                self.status_label.setText("Binaries installed. Detecting tools...")
                self.auto_config()
                QMessageBox.information(
                    self,
                    "Binaries installed",
                    "The platform binaries were installed into YR_MPE/bin.",
                )
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Extraction failed",
                "The archive downloaded, but could not be installed:\n%s" % exc,
            )
            self.status_label.setText("Extraction failed: %s" % exc)
        finally:
            self.progress_bar.setVisible(False)
            self._cyan_mode = False
            self.cyan_download_button.setEnabled(True)
            if archive_path and os.path.exists(archive_path):
                try:
                    os.remove(archive_path)
                except OSError:
                    pass
            self.download_thread = None

    def on_download_cancelled(self):
        self.download_button.setEnabled(True)
        self.cyan_download_button.setEnabled(True)
        self._cyan_mode = False
        self.cancel_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.status_label.setText("Download cancelled.")
        if self._download_path and os.path.exists(self._download_path):
            try:
                os.remove(self._download_path)
            except OSError:
                pass
        self.download_thread = None

    def on_download_error(self, message):
        self.download_button.setEnabled(True)
        self.cyan_download_button.setEnabled(True)
        self._cyan_mode = False
        self.cancel_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.status_label.setText("Download failed.")
        QMessageBox.critical(
            self,
            "Download failed",
            "Could not download the binary archive:\n%s\n\n"
            "Try the other mirror, or paste a direct URL. On Windows, a "
            "certificate-revocation error usually means the system CA store "
            "is offline; use a mirror or a current CA bundle." % message,
        )
        if self._download_path and os.path.exists(self._download_path):
            try:
                os.remove(self._download_path)
            except OSError:
                pass
        self.download_thread = None

    def install_bundle(self, archive_path):
        target = self.plugin_path / "YR_MPE" / "bin"
        target.mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="yrmpe_extract_"))
        try:
            if zipfile.is_zipfile(archive_path):
                _safe_extract_zip(archive_path, work)
            elif tarfile.is_tarfile(archive_path):
                _safe_extract_tar(archive_path, work)
            else:
                raise RuntimeError("Unsupported archive format; use .zip or .tar.gz.")
            binary_root = _find_binary_root(work)
            if binary_root is None:
                names = sorted(path.name for path in work.iterdir())
                raise RuntimeError(
                    "No bin/ directory found in the archive. Top-level entries: %s"
                    % ", ".join(names[:12])
                )
            shutil.copytree(binary_root, target, dirs_exist_ok=True)
            # Platform bundles also carry the IcyTree static site under src/.
            for src_root in (binary_root / "src", binary_root.parent / "src"):
                if not src_root.is_dir():
                    continue
                src_target = self.plugin_path / "YR_MPE" / "src"
                src_target.mkdir(parents=True, exist_ok=True)
                shutil.copytree(src_root, src_target, dirs_exist_ok=True)
                break
        finally:
            shutil.rmtree(work, ignore_errors=True)

    def install_cyan_tool(self, archive_path, tool):
        """Install one CYAN package under YR_MPE/bin (or YR_MPE/src)."""
        target = self.plugin_path / "YR_MPE" / (
            "src" if tool in CYAN_SRC_TOOLS else "bin")
        target.mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="yrmpe_cyan_extract_"))
        try:
            if zipfile.is_zipfile(archive_path):
                _safe_extract_zip(archive_path, work)
            elif tarfile.is_tarfile(archive_path):
                _safe_extract_tar(archive_path, work)
            else:
                raise RuntimeError("Unsupported archive format; use .zip or .tar.gz.")
            roots = [path for path in work.iterdir() if path.is_dir()]
            if not roots:
                shutil.copytree(work, target, dirs_exist_ok=True)
                return
            for root in roots:
                shutil.copytree(root, target / root.name, dirs_exist_ok=True)
        finally:
            shutil.rmtree(work, ignore_errors=True)

    # ---------------------------------------------------------- configure
    def _bundled_path(self, tool):
        bundled = tool.get("bundled", "")
        if not bundled:
            return None
        candidate = self.plugin_path / "YR_MPE" / "bin" / bundled
        candidates = [candidate]
        if bundled.lower().endswith(".exe"):
            candidates.append(candidate.with_suffix(""))
        else:
            candidates.append(candidate.with_name(candidate.name + ".exe"))
        for path in candidates:
            if path.exists():
                return path
        return None

    def auto_config(self):
        if not self.tools_data:
            QMessageBox.warning(self, "Warning", "No tools data loaded.")
            return
        self.configured_tools = []
        not_found = []
        for tool in self.tools_data:
            config_entry = tool.copy()
            path = ""
            bundled = self._bundled_path(tool)
            if bundled is not None:
                relative = bundled.relative_to(
                    self.plugin_path / "YR_MPE" / "bin")
                path = "/" + relative.as_posix()
            else:
                for command in tool.get("cmd", []):
                    found = shutil.which(command)
                    if found:
                        path = found
                        break
            config_entry["path"] = path
            self.configured_tools.append(config_entry)
            if not path:
                not_found.append(tool)

        if not_found:
            message = "The following tools were not found:\n"
            message += "".join("- %s\n" % tool["name"] for tool in not_found)
            message += "\nUse Manual Config to point at them, or install the bundled binaries."
            QMessageBox.information(self, "Auto Config Result", message)
        else:
            QMessageBox.information(
                self, "Auto Config Result", "All tools configured successfully!"
            )

    def manual_config(self):
        if not self.tools_data:
            QMessageBox.warning(self, "Warning", "No tools data loaded.")
            return
        tools = self.configured_tools if self.configured_tools else self.tools_data
        dialog = ManualConfigDialog(tools, self)
        if dialog.exec_() == QDialog.Accepted:
            self.configured_tools = dialog.get_configured_tools()
            QMessageBox.information(self, "Manual Config", "Tool paths updated.")

    # ------------------------------------------------------------- finish
    def finish_install(self):
        if not self.configured_tools:
            reply = QMessageBox.question(
                self,
                "Warning",
                "Tools have not been configured. Run Auto Config first?",
                QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel,
            )
            if reply == QMessageBox.Cancel:
                return
            if reply == QMessageBox.Yes:
                self.auto_config()
                if not self.configured_tools:
                    return

        mpe_path = self.path_edit.text().strip() or "YR-MPE"

        settings_template = self.plugin_path / "settings.template.json"
        if settings_template.exists():
            content = settings_template.read_text(encoding="utf-8")
            content = content.replace("%version%", self._version())
            content = content.replace("%MPEPATH%", mpe_path)
            (self.plugin_path / "settings.json").write_text(
                content, encoding="utf-8"
            )

        config_data = []
        for tool in self.configured_tools:
            config_data.append(
                {"name": tool["name"], "path": tool.get("path", "")}
            )
        config_path = self.plugin_path / "YR_MPE" / "config.json"
        config_path.parent.mkdir(parents=True, exist_ok=True)
        with open(config_path, "w", encoding="utf-8") as handle:
            json.dump(config_data, handle, indent=4, ensure_ascii=False)

        QMessageBox.information(
            self,
            "Installation Complete",
            "YR-MPE has been configured.\nRestart YR-MPE to apply changes.",
        )
        self.close()


class InitializerEntry:
    def run(self):
        return Initializer()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    initializer = Initializer()
    initializer.show()
    sys.exit(app.exec_())
