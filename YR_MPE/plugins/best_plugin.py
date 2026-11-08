"""BEST-style configuration for MrBayes' multi-species coalescent.

BEST is the full-likelihood multispecies-coalescent analysis implemented by
MrBayes 3. This dialog writes the MRBAYES block; it does not invent a second
MCMC engine. The layout follows the BEAUTi / StarBEAST3 idea: partitions,
taxa/ploidy, priors, and MCMC are separate configuration pages.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio.Nexus import Nexus
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit,
    QPushButton, QSpinBox, QTabWidget, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from ..templates.base_plugin_ui import BasePlugin
from ..templates.base_process_thread import BaseProcessThread
from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec


@dataclass
class BestConfig:
    partitions: List[Tuple[str, str, int]] = field(default_factory=list)
    species_partition: str = ""
    ploidy: str = "diploid"
    species_prior: str = "yule"
    clock_model: str = "strict"
    clock_rate_prior: str = "fixed(1.0)"
    partition_rate_prior: str = "fixed"
    gene_rate_prior: str = "fixed"
    speciation_prior: str = "exponential(10.0)"
    extinction_prior: str = "beta(1,1)"
    sample_prob: str = "0.5"
    popsize_prior: str = "gamma(1,100)"
    ngen: int = 100000
    samplefreq: int = 100
    printfreq: int = 1000
    nruns: int = 2
    nchains: int = 4
    extra_lines: List[str] = field(default_factory=list)


def parse_nexus_taxa(path: str) -> List[str]:
    """Return matrix taxon labels from a NEXUS data file."""
    try:
        nexus = Nexus.Nexus(str(path))
        labels = [str(label) for label in nexus.taxlabels if str(label).strip()]
        if labels:
            return labels
    except Exception:
        pass
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    match = re.search(r"\bmatrix\b(.*?);", text, re.I | re.S)
    if not match:
        return []
    labels = []
    for line in match.group(1).splitlines():
        line = line.strip()
        if line:
            labels.append(line.split(None, 1)[0])
    return labels


def build_mrbayes_block(config: BestConfig) -> str:
    """Render a MrBayes MRBAYES block for the multi-species coalescent."""
    lines = ["begin mrbayes;", "    set autoclose=yes nowarn=yes;"]
    if config.partitions:
        for name, range_text, _group in config.partitions:
            if str(range_text).strip():
                lines.append("    charset %s = %s;"
                             % (name, str(range_text).strip().rstrip(";")))
        names = ",".join(name for name, _, _ in config.partitions)
        lines.append("    partition genes = %d:%s;" % (len(config.partitions), names))
        lines.append("    set partition=genes;")
        lines.append("    unlink topology=(all);")
        for parameter in ("revmat", "statefreq", "shape", "pinvar"):
            lines.append("    unlink %s=(all);" % parameter)
        groups = {}
        for index, (_, _, group) in enumerate(config.partitions, start=1):
            groups.setdefault(group, []).append(index)
        for members in groups.values():
            if len(members) > 1:
                member_text = ",".join(str(index) for index in members)
                lines.append(
                    "    link revmat=(%s) statefreq=(%s) shape=(%s) "
                    "pinvar=(%s);"
                    % (member_text, member_text, member_text, member_text)
                )
    if config.species_partition.strip():
        species_text = config.species_partition.strip().rstrip(";")
        species_text = ", ".join(
            line.strip() for line in species_text.splitlines()
            if line.strip()
        )
        lines.append("    speciespartition species = %s;"
                     % species_text)
        lines.append("    set speciespartition=species;")
    lines.append("    lset ploidy=%s;" % config.ploidy)
    lines.append("    prset topologypr=speciestree;")
    lines.append("    prset brlenspr=clock:speciestree;")
    lines.append("    prset clockvarpr=%s;" % config.clock_model)
    lines.append("    prset clockratepr=%s;" % config.clock_rate_prior)
    # MrBayes does not give every gene tree its own StarBEAST3-style clock
    # model: partitions within a gene use Ratepr, genes use Generatepr.
    if config.partition_rate_prior != "fixed":
        lines.append("    prset ratepr=%s;" % config.partition_rate_prior)
    if config.gene_rate_prior != "fixed":
        lines.append("    prset generatepr=%s;" % config.gene_rate_prior)
    lines.append("    prset popvarpr=variable;")
    lines.append("    prset popsizepr=%s;" % config.popsize_prior)
    speciation = (config.speciation_prior or "exponential(10.0)").strip()
    # MrBayes rejects gamma() as a speciation prior; keep the configuration
    # renderable when an imported block contains an older invalid value.
    if "gamma" in speciation.lower():
        speciation = "exponential(10.0)"
    lines.append("    prset speciationpr=%s;" % speciation)
    if config.species_prior == "yule":
        # MrBayes rejects an explicit fixed extinction prior for this
        # clock tree; its default Yule setting is already extinction-free.
        pass
    else:
        extinction = (config.extinction_prior or "beta(1,1)").strip()
        # A fixed zero extinction is outside MrBayes' allowed (0, 1) range.
        if extinction.replace(" ", "").lower() == "fixed(0)":
            extinction = "beta(1,1)"
        lines.append("    prset extinctionpr=%s;" % extinction)
        lines.append("    prset sampleprob=%s;" % config.sample_prob)
    if config.extra_lines:
        lines.extend(config.extra_lines)
    lines.append("    mcmc ngen=%d samplefreq=%d printfreq=%d nruns=%d nchains=%d;"
                 % (config.ngen, config.samplefreq, config.printfreq,
                    config.nruns, config.nchains))
    lines.append("    sump;")
    lines.append("    sumt;")
    lines.append("end;")
    return "\n".join(lines) + "\n"


def parse_nexus_charsets(path: str) -> List[Tuple[str, str, int]]:
    """Read charset name/range pairs from a NEXUS file."""
    partitions = []
    pattern = re.compile(r"^\s*charset\s+([A-Za-z0-9_.-]+)\s*=\s*([^;]+);", re.I)
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = pattern.match(line)
            if match:
                partitions.append((match.group(1), match.group(2).strip(), 1))
    return partitions


def parse_mrbayes_block(block: str) -> BestConfig:
    """Parse the subset of MRBAYES used by this configurator.

    Unknown executable lines are retained in ``extra_lines`` so an agent can
    round-trip custom MrBayes commands without silently losing them.
    """
    config = BestConfig()
    text = re.sub(r"\[[^\]]*\]", " ", block or "", flags=re.S)
    commands = [item.strip() for item in text.split(";") if item.strip()]
    groups: Dict[int, int] = {}
    charset_ranges: Dict[str, str] = {}
    for command in commands:
        lower = command.lower()
        if lower.startswith("link "):
            match = re.search(r"revmat\s*=\s*\(([^)]+)\)", command, re.I)
            if match:
                for index, token in enumerate(
                    part.strip() for part in match.group(1).split(",")
                ):
                    if token.isdigit():
                        groups[int(token)] = index + 1
    for command in commands:
        command = " ".join(command.split())
        lower = command.lower()
        if lower.startswith("begin") or lower.startswith("end"):
            continue
        if lower.startswith("charset "):
            match = re.match(r"charset\s+(\S+)\s*=\s*(.+)", command, re.I)
            if match:
                charset_ranges[match.group(1)] = match.group(2).strip()
            continue
        if lower.startswith("partition genes"):
            match = re.search(
                r"partition\s+genes\s*=\s*\d+\s*:\s*(.+)", command, re.I
            )
            if match:
                names = [name.strip() for name in match.group(1).split(",")
                         if name.strip()]
                config.partitions = [
                    (name, charset_ranges.get(name, ""),
                     groups.get(index + 1, 1))
                    for index, name in enumerate(names)
                ]
            continue
        if lower.startswith("speciespartition species"):
            match = re.search(r"speciespartition\s+species\s*=\s*(.+)",
                              command, re.I)
            if match:
                config.species_partition = match.group(1).strip()
            continue
        if lower.startswith("link ") or lower.startswith("unlink "):
            continue
        if lower.startswith("lset"):
            match = re.search(r"ploidy\s*=\s*([^\s]+)", command, re.I)
            if match:
                config.ploidy = match.group(1).strip()
            continue
        if lower.startswith("prset"):
            changed = False
            for key in (
                "clockvarpr", "clockratepr", "ratepr", "generatepr",
                "popsizepr", "speciationpr", "extinctionpr", "sampleprob",
                "popvarpr",
            ):
                match = re.search(r"%s\s*=\s*([^\s]+)" % key, command, re.I)
                if not match:
                    continue
                value = match.group(1).strip()
                changed = True
                if key == "clockvarpr":
                    config.clock_model = value
                elif key == "clockratepr":
                    config.clock_rate_prior = value
                elif key == "ratepr":
                    config.partition_rate_prior = value
                elif key == "generatepr":
                    config.gene_rate_prior = value
                elif key == "popsizepr":
                    config.popsize_prior = value
                elif key == "speciationpr":
                    config.speciation_prior = value
                elif key == "extinctionpr":
                    config.extinction_prior = value
                    config.species_prior = (
                        "yule" if value.replace(" ", "").lower() == "fixed(0)"
                        else "birth-death"
                    )
                elif key == "sampleprob":
                    config.sample_prob = value
            if lower.startswith("prset topologypr") or lower.startswith(
                "prset brlenspr"
            ):
                changed = True
            if not changed:
                config.extra_lines.append(command + ";")
            continue
        if lower.startswith("mcmc"):
            for key in ("ngen", "samplefreq", "printfreq", "nruns", "nchains"):
                match = re.search(r"%s\s*=\s*(\d+)" % key, command, re.I)
                if match:
                    setattr(config, key, int(match.group(1)))
            continue
        if lower.startswith("set "):
            if lower.startswith(("set autoclose", "set partition",
                                 "set speciespartition")):
                continue
            config.extra_lines.append(command + ";")
            continue
        if lower.startswith("sump") or lower.startswith("sumt"):
            continue
        config.extra_lines.append(command + ";")
    return config


def find_mrbayes_tool() -> str:
    """Resolve MrBayes, allowing an explicit environment override."""
    override = os.environ.get("YRMPE_MRBAYES", "").strip()
    if override and os.path.isfile(override):
        return override
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(root, "config.json")
    try:
        with open(config_path, "r", encoding="utf-8") as handle:
            tools = json.load(handle)
        for tool in tools:
            if str(tool.get("name", "")).lower() == "mrbayes-mpi-beagle":
                path = str(tool.get("path", ""))
                return os.path.join(root, path.lstrip("/\\"))
    except (OSError, ValueError):
        pass
    return ""


class BestThread(BaseProcessThread):
    """Run MCMC and summarize the species tree with a second MrBayes pass.

    MrBayes writes BEST species trees under ``.tree<N>.run<R>.t`` names, while
    ``sumt`` expects ``.run<R>.t``. It also expects the taxa in those trees to
    match the data block. The second pass therefore copies the last tree
    partition (the species tree) to the conventional names and supplies a
    species-only dummy data block before calling ``sumt``.
    """

    result_ready = pyqtSignal(object)

    def __init__(self, tool_path: str, input_file: str, block: str,
                 workdir: Optional[str] = None):
        super().__init__(tool_path, [input_file], [], [input_file],
                         workdir=workdir)
        self.input_file = input_file
        self.block = block

    def get_tool_name(self) -> str:
        return "BEST / MrBayes MSC"

    @staticmethod
    def _species_names(block: str) -> List[str]:
        match = re.search(
            r"speciespartition\s+\w+\s*=\s*(.+?);", block or "",
            re.I | re.S,
        )
        if not match:
            return []
        names = []
        for item in match.group(1).split(","):
            name = item.split(":", 1)[0].strip()
            if name and name not in names:
                names.append(name)
        return names

    @staticmethod
    def _strip_summary(block: str) -> str:
        block = re.sub(r"(?im)^\s*sump\s*;\s*$", "", block or "")
        return re.sub(r"(?im)^\s*sumt\s*;\s*$", "", block)

    def execute_commands(self):
        if not os.path.isfile(self.input_file):
            raise RuntimeError("Input NEXUS file not found: %s" % self.input_file)
        root = self.workdir or tempfile.mkdtemp(prefix="best_")
        os.makedirs(root, exist_ok=True)
        run_dir = tempfile.mkdtemp(prefix="best_run_", dir=root)
        stem = os.path.splitext(os.path.basename(self.input_file))[0]
        nexus_path = os.path.join(run_dir, stem + "_best.nex")
        try:
            text = open(self.input_file, "r", encoding="utf-8",
                        errors="replace").read()
        except OSError as exc:
            raise RuntimeError("Cannot read input NEXUS: %s" % exc)
        text = re.sub(r"begin\s+mrbayes\s*;.*?end\s*;", "", text,
                      flags=re.I | re.S)

        # Stage 1: MCMC without the final summaries.
        mcmc_block = self._strip_summary(self.block)
        with open(nexus_path, "w", encoding="utf-8") as handle:
            handle.write(text.rstrip() + "\n\n" + mcmc_block.rstrip() + "\n")
        self.progress.emit("Running BEST/MrBayes MCMC...")
        result = self.execute_command([self.tool_path, nexus_path],
                                      cwd=run_dir)
        if result.returncode != 0:
            raise RuntimeError(
                "MrBayes MCMC failed: %s"
                % ((result.stderr or "").strip() or "non-zero exit code")
            )

        prefix = os.path.basename(nexus_path)
        files = sorted(
            os.path.join(run_dir, name)
            for name in os.listdir(run_dir)
            if name.startswith(prefix)
            and os.path.isfile(os.path.join(run_dir, name))
        )

        # Stage 2: copy the species-tree sample and run standard sumt.
        tree_files = [path for path in files if path.endswith(".t")]
        species = self._species_names(self.block)
        if not tree_files or not species:
            raise RuntimeError("BEST did not leave a species-tree sample for sumt.")
        tree_numbers = [
            int(match.group(1))
            for path in tree_files
            for match in [re.search(r"\.tree(\d+)\.run\d+\.t$", path)]
            if match
        ]
        if not tree_numbers:
            raise RuntimeError("Could not identify the BEST tree partitions.")
        max_tree = max(tree_numbers)
        copied = []
        for path in tree_files:
            match = re.search(r"\.tree(\d+)\.run(\d+)\.t$", path)
            if not match or int(match.group(1)) != max_tree:
                continue
            run_number = int(match.group(2))
            target = os.path.join(run_dir, "%s.run%d.t" % (prefix, run_number))
            shutil.copy2(path, target)
            copied.append(target)
        if not copied:
            raise RuntimeError("Could not identify the BEST species-tree files.")
        rows = "\n".join(" %s A" % name for name in species)
        stage2 = (
            "#NEXUS\nbegin data;\n"
            " dimensions ntax=%d nchar=1;\n"
            " format datatype=dna;\n"
            " matrix\n%s\n ;\nend;\n\n"
            "begin mrbayes;\n"
            " set autoclose=yes nowarn=yes;\n"
            " sumt;\nend;\n"
            % (len(species), rows)
        )
        with open(nexus_path, "w", encoding="utf-8") as handle:
            handle.write(stage2)
        self.progress.emit("Summarizing BEST species trees with sumt...")
        summary = self.execute_command([self.tool_path, nexus_path],
                                       cwd=run_dir)
        if summary.returncode != 0:
            raise RuntimeError(
                "MrBayes sumt failed: %s"
                % ((summary.stderr or "").strip() or "non-zero exit code")
            )
        summary_files = sorted(
            os.path.join(run_dir, name)
            for name in os.listdir(run_dir)
            if name.startswith(prefix)
            and os.path.isfile(os.path.join(run_dir, name))
            and name.endswith(".con.tre")
        )
        if not summary_files:
            raise RuntimeError("MrBayes sumt produced no consensus tree.")
        output_files = files + summary_files
        self.result_ready.emit({"output_files": output_files})
        self.finished.emit(output_files, [])


class BestDialog(QDialog):
    """BEAUTi-like configuration page for MrBayes MSC / BEST.

    The dialog is agent-operable through two complementary paths:
    ``best.apply_block`` accepts a complete MRBAYES block and reverse-parses
    it into the form, while ``best.modify_options`` addresses the common
    controls directly.
    """

    result_ready = pyqtSignal(object)

    def __init__(self, parent=None, workdir=None, **kwargs):
        super().__init__(parent)
        self.plugin_name = "BEST / MrBayes MSC"
        self.tool_name = "MrBayes-MPI-BEAGLE"
        self.workdir = workdir
        self.tool_path = find_mrbayes_tool()
        self.is_running = False
        self.analysis_thread = None
        self.current_output_files = []
        self.task_queue = None
        self.task_window = None
        self.queue_job_id = ""
        self._console = []
        self._last_error = ""
        self._extra_lines: List[str] = []
        self.setWindowTitle("BEST / MrBayes Multi-Species Coalescent")
        self.resize(860, 700)
        self._build_ui()
        self.refresh_preview()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)

        input_row = QHBoxLayout()
        input_row.addWidget(QLabel("NEXUS with charset partitions:"))
        self.input_edit = QLineEdit()
        input_row.addWidget(self.input_edit, 1)
        browse = QPushButton("Browse")
        browse.clicked.connect(self._browse_input)
        input_row.addWidget(browse)
        load = QPushButton("Load partitions")
        load.clicked.connect(self._load_partitions)
        input_row.addWidget(load)
        layout.addLayout(input_row)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, 1)

        self.partition_table = QTableWidget(0, 3)
        self.partition_table.setHorizontalHeaderLabels(
            ["Partition", "Range", "Link group"]
        )
        self.partition_table.horizontalHeader().setStretchLastSection(True)
        self.tabs.addTab(self.partition_table, "Partitions")

        taxa = QWidget()
        taxa_layout = QVBoxLayout(taxa)
        self.species_edit = QPlainTextEdit()
        self.species_edit.setPlaceholderText(
            "Sp1: 1, Sp2: 2, ...\n\nUse taxon names or matrix row numbers."
        )
        species_header = QHBoxLayout()
        species_header.addWidget(QLabel(
            "Species partition (MrBayes speciespartition syntax; e.g. "
            "SpA: 1-2, SpB: 3-4):"
        ))
        self.auto_species_button = QPushButton("One individual = one species")
        self.auto_species_button.setToolTip(
            "Fill the species partition using one matrix row per species."
        )
        self.auto_species_button.clicked.connect(self._auto_species_partition)
        species_header.addWidget(self.auto_species_button)
        species_header.addStretch(1)
        taxa_layout.addLayout(species_header)
        taxa_layout.addWidget(self.species_edit, 1)
        self.tabs.addTab(taxa, "Taxa")

        priors = QWidget()
        pform = QFormLayout(priors)
        self.ploidy_combo = QComboBox()
        self.ploidy_combo.addItems(["haploid", "diploid", "zlinked"])
        self.ploidy_combo.setCurrentText("diploid")
        pform.addRow("Ploidy:", self.ploidy_combo)
        self.prior_combo = QComboBox()
        self.prior_combo.addItems(["Yule (extinction = 0)", "Birth-death"])
        pform.addRow("Species-tree prior:", self.prior_combo)
        self.clock_combo = QComboBox()
        self.clock_combo.addItems(["strict"])
        self.clock_combo.setEnabled(False)
        pform.addRow("Species-tree clock (clockvarpr):", self.clock_combo)
        clock_note = QLabel(
            "MrBayes BEST/MSC is documented with a strict clock. Relaxed "
            "gene/species clocks belong to StarBEAST3, not this entry."
        )
        clock_note.setWordWrap(True)
        clock_note.setStyleSheet("color: palette(mid);")
        pform.addRow("", clock_note)
        self.clock_rate_edit = QLineEdit("fixed(1.0)")
        pform.addRow("Species-tree clock rate prior:", self.clock_rate_edit)
        self.partition_rate_combo = QComboBox()
        self.partition_rate_combo.addItems(["fixed", "variable"])
        pform.addRow("Partition rate (Ratepr):", self.partition_rate_combo)
        self.gene_rate_combo = QComboBox()
        self.gene_rate_combo.addItems(["fixed", "variable"])
        pform.addRow("Gene rate (Generatepr):", self.gene_rate_combo)
        self.speciation_edit = QLineEdit("exponential(10.0)")
        pform.addRow("Speciation prior:", self.speciation_edit)
        self.extinction_edit = QLineEdit("beta(1,1)")
        pform.addRow("Extinction prior (BD):", self.extinction_edit)
        self.sample_edit = QLineEdit("0.5")
        pform.addRow("Sampling probability:", self.sample_edit)
        self.popsize_edit = QLineEdit("gamma(1,100)")
        pform.addRow("Population-size prior:", self.popsize_edit)
        self.tabs.addTab(priors, "Priors")

        mcmc = QWidget()
        mform = QFormLayout(mcmc)
        self.ngen_spin = QSpinBox(); self.ngen_spin.setRange(1000, 1000000000)
        self.ngen_spin.setValue(100000)
        self.sample_spin = QSpinBox(); self.sample_spin.setRange(1, 1000000)
        self.sample_spin.setValue(100)
        self.printfreq_spin = QSpinBox(); self.printfreq_spin.setRange(1, 1000000)
        self.printfreq_spin.setValue(1000)
        self.nruns_spin = QSpinBox(); self.nruns_spin.setRange(1, 8)
        self.nruns_spin.setValue(2)
        self.nchains_spin = QSpinBox(); self.nchains_spin.setRange(1, 20)
        self.nchains_spin.setValue(4)
        mform.addRow("Generations:", self.ngen_spin)
        mform.addRow("Sample frequency:", self.sample_spin)
        mform.addRow("Print frequency:", self.printfreq_spin)
        mform.addRow("Runs:", self.nruns_spin)
        mform.addRow("Chains:", self.nchains_spin)
        self.tabs.addTab(mcmc, "MCMC")

        self.preview = QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setStyleSheet("font-family: Consolas, monospace;")
        self.tabs.addTab(self.preview, "Preview")

        result_page = QWidget()
        result_layout = QVBoxLayout(result_page)
        self.output_preview = QPlainTextEdit()
        self.output_preview.setReadOnly(True)
        result_layout.addWidget(self.output_preview)
        self.tabs.addTab(result_page, "Result")

        actions = QHBoxLayout()
        refresh = QPushButton("Refresh preview")
        refresh.clicked.connect(self.refresh_preview)
        write = QPushButton("Write MRBAYES block…")
        write.clicked.connect(self._write_block)
        self.run_button = QPushButton("Run MrBayes")
        self.run_button.clicked.connect(self._request_run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._request_stop)
        self.import_button = QPushButton("Import result")
        self.import_button.setEnabled(False)
        self.import_button.clicked.connect(self.import_to_platform)
        self.mrbayes_button = QPushButton("Close && open MrBayes")
        self.mrbayes_button.setToolTip(
            "Write the current BEST NEXUS and open the standard MrBayes plugin."
        )
        self.mrbayes_button.clicked.connect(self._close_and_open_mrbayes)
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        actions.addWidget(refresh)
        actions.addWidget(write)
        actions.addWidget(self.run_button)
        actions.addWidget(self.stop_button)
        actions.addWidget(self.import_button)
        actions.addWidget(self.mrbayes_button)
        actions.addStretch(1)
        actions.addWidget(close)
        layout.addLayout(actions)

    def _browse_input(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select NEXUS file", "", "NEXUS (*.nex *.nexus);;All files (*)"
        )
        if path:
            self.input_edit.setText(path)
            self._load_partitions()

    def _load_partitions(self):
        path = self.input_edit.text().strip()
        if not path or not os.path.isfile(path):
            return
        try:
            partitions = parse_nexus_charsets(path)
        except OSError as exc:
            QMessageBox.warning(self, "Partitions", str(exc))
            return
        self.partition_table.setRowCount(len(partitions))
        for row, (name, range_text, group) in enumerate(partitions):
            self.partition_table.setItem(row, 0, QTableWidgetItem(name))
            self.partition_table.setItem(row, 1, QTableWidgetItem(range_text))
            self.partition_table.setItem(row, 2, QTableWidgetItem(str(group)))
        self.partition_table.resizeColumnsToContents()
        self.refresh_preview()

    def _auto_species_partition(self):
        """Fill one species per matrix row for single-individual datasets."""
        path = self.input_edit.text().strip()
        if not path or not os.path.isfile(path):
            QMessageBox.warning(
                self, "Species partition",
                "Select the NEXUS input file first.",
            )
            return
        taxa = parse_nexus_taxa(path)
        if not taxa:
            QMessageBox.warning(
                self, "Species partition",
                "No taxa were found in the NEXUS matrix.",
            )
            return
        self.species_edit.setPlainText(
            ", ".join("Sp%d: %d" % (index, index)
                      for index in range(1, len(taxa) + 1))
        )
        self.refresh_preview()

    def _write_launch_nexus(self) -> Path:
        """Write the data block plus the current MRBAYES block for MrBayes."""
        path = Path(self.input_edit.text().strip())
        if not path.is_file():
            raise ValueError("Select the NEXUS input file first.")
        if self.partition_table.rowCount() == 0:
            self._load_partitions()
        if self.partition_table.rowCount() == 0:
            raise ValueError(
                "No charsets were loaded. Load the input partitions first."
            )
        text = path.read_text(encoding="utf-8", errors="replace")
        text = re.sub(
            r"begin\s+mrbayes\s*;.*?end\s*;", "", text,
            flags=re.I | re.S,
        )
        # Charsets are emitted inside the MRBAYES block; a separate SETS
        # block is skipped by MrBayes and only causes confusion.
        text = re.sub(
            r"begin\s+sets\s*;.*?end\s*;", "", text,
            flags=re.I | re.S,
        )
        block = build_mrbayes_block(self._config())
        output = path.with_name(path.stem + "_best_mrbayes.nex")
        if output.resolve() == path.resolve():
            output = path.with_name(path.stem + "_best_launch.nex")
        output.write_text(text.rstrip() + "\n\n" + block, encoding="utf-8")
        return output

    def _close_and_open_mrbayes(self):
        if self.is_running:
            QMessageBox.warning(
                self, "MrBayes",
                "Stop the running BEST analysis before opening MrBayes.",
            )
            return
        try:
            output = self._write_launch_nexus()
        except Exception as exc:
            QMessageBox.warning(self, "MrBayes", str(exc))
            return
        parent = self.parent()
        opener = getattr(parent, "open_mrbayes_wrapper", None)
        if not callable(opener):
            QMessageBox.warning(
                self, "MrBayes",
                "The host cannot open the MrBayes plugin.",
            )
            return
        self.close()
        opener(import_path=str(output))

    def _config(self) -> BestConfig:
        partitions = []
        for row in range(self.partition_table.rowCount()):
            name_item = self.partition_table.item(row, 0)
            range_item = self.partition_table.item(row, 1)
            group_item = self.partition_table.item(row, 2)
            if not name_item:
                continue
            try:
                group = int(group_item.text()) if group_item else 1
            except ValueError:
                group = 1
            partitions.append((
                name_item.text(),
                range_item.text() if range_item else "",
                group,
            ))
        return BestConfig(
            partitions=partitions,
            species_partition=self.species_edit.toPlainText(),
            ploidy=self.ploidy_combo.currentText(),
            species_prior="yule" if self.prior_combo.currentIndex() == 0 else "birth-death",
            clock_model=self.clock_combo.currentText(),
            clock_rate_prior=self.clock_rate_edit.text().strip() or "fixed(1.0)",
            partition_rate_prior=self.partition_rate_combo.currentText(),
            gene_rate_prior=self.gene_rate_combo.currentText(),
            speciation_prior=self.speciation_edit.text().strip() or "exponential(10.0)",
            extinction_prior=self.extinction_edit.text().strip() or "beta(1,1)",
            sample_prob=self.sample_edit.text().strip() or "0.5",
            popsize_prior=self.popsize_edit.text().strip() or "gamma(1,100)",
            ngen=self.ngen_spin.value(),
            samplefreq=self.sample_spin.value(),
            printfreq=self.printfreq_spin.value(),
            nruns=self.nruns_spin.value(),
            nchains=self.nchains_spin.value(),
            extra_lines=list(self._extra_lines),
        )

    def refresh_preview(self):
        self.preview.setPlainText(build_mrbayes_block(self._config()))

    # ------------------------------------------------------------------
    # agent / execution surface
    # ------------------------------------------------------------------
    def set_task_queue(self, queue, window=None):
        self.task_queue = queue
        self.task_window = window
        if getattr(self, "run_button", None) is not None:
            self.run_button.setToolTip(
                "Run now, or queue while another analysis is active."
            )

    def _request_run(self):
        if self.is_running:
            return None
        queue = getattr(self, "task_queue", None)
        if queue is not None and callable(getattr(queue, "submit", None)):
            self.queue_job_id = queue.submit(
                self, window=self.task_window
            )
            return self.queue_job_id
        return self.run_analysis()

    def _request_stop(self):
        queue = getattr(self, "task_queue", None)
        job_id = getattr(self, "queue_job_id", "")
        if queue is not None and job_id and callable(getattr(queue, "get_job", None)):
            job = queue.get_job(job_id)
            if job is not None and getattr(job, "state", "") in (
                "queued", "running"
            ):
                queue.cancel(job_id)
                return
        self.stop_analysis()

    def run_analysis(self):
        if self.is_running:
            return
        input_path = self.input_edit.text().strip()
        if not input_path or not os.path.isfile(input_path):
            QMessageBox.warning(
                self, "BEST input", "Select a NEXUS data file first."
            )
            return
        if not self.tool_path or not os.path.isfile(self.tool_path):
            QMessageBox.critical(
                self, "BEST", "MrBayes executable was not found in config.json."
            )
            return
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.import_button.setEnabled(False)
        self.output_preview.setPlainText("Running BEST / MrBayes...")
        self._console = []
        self._last_error = ""
        self.current_output_files = []
        self.analysis_thread = BestThread(
            self.tool_path,
            input_path,
            build_mrbayes_block(self._config()),
            workdir=self.workdir,
        )
        self.analysis_thread.progress.connect(self._on_progress)
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.finished.connect(self._on_finished)
        self.analysis_thread.error.connect(self._on_error)
        self.analysis_thread.console_output.connect(self._on_console)
        self.analysis_thread.start()

    def stop_analysis(self):
        thread = getattr(self, "analysis_thread", None)
        if thread is not None and thread.isRunning():
            for name in ("progress", "result_ready", "finished", "error",
                         "console_output"):
                signal = getattr(thread, name, None)
                if signal is not None:
                    try:
                        signal.disconnect()
                    except (TypeError, RuntimeError):
                        pass
            thread.stop()
            thread.wait(3000)
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.output_preview.setPlainText("Stopped by user.")

    def _on_progress(self, text):
        self._console.append(str(text))

    def _on_console(self, message, kind="info"):
        self._console.append(str(message))

    def _on_result(self, payload):
        if isinstance(payload, dict):
            self.current_output_files = list(
                payload.get("output_files", []) or []
            )

    def _on_finished(self, output_files, html_files):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.current_output_files = list(output_files or [])
        tree_file = next(
            (path for path in self.current_output_files
             if path.endswith(".con.tre")), ""
        )
        if tree_file:
            try:
                from .mrbayes_plugin import read_nexus_trees
                trees = read_nexus_trees(tree_file)
                self.output_preview.setPlainText("\n".join(trees))
            except Exception:
                self.output_preview.setPlainText(
                    "Result files:\n" + "\n".join(self.current_output_files)
                )
        else:
            self.output_preview.setPlainText(
                "Result files:\n" + "\n".join(self.current_output_files)
            )
        self.import_button.setEnabled(bool(self.current_output_files))
        QMessageBox.information(
            self, "BEST completed",
            "BEST / MrBayes run completed. Review the Result tab or import "
            "the tree and chains into the workspace.",
        )

    def _on_error(self, message):
        self.is_running = False
        self._last_error = str(message)
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.output_preview.setPlainText("BEST failed:\n%s" % message)
        QMessageBox.critical(self, "BEST failed", str(message))

    def _build_result_spec(self):
        files = list(self.current_output_files or [])
        if not files:
            raise ValueError("No BEST output is available.")
        from .mrbayes_plugin import read_nexus_trees
        items = []
        artifacts = []
        for path in files:
            name = os.path.basename(path)
            if name.endswith(".con.tre"):
                trees = read_nexus_trees(path)
                if trees:
                    items.append(DatasetItemSpec(
                        loci_name=os.path.splitext(name)[0],
                        item_type="phylogeny",
                        file_path=str(path),
                        format="nexus",
                        data={"newick": trees[0], "trees": trees},
                        metadata={"source": "best", "tool": "MrBayes",
                                  "analysis": "BEST/MSC"},
                    ))
            else:
                artifacts.append(str(path))
        p_files = [path for path in files if path.endswith(".p")]
        by_run = {}
        for path in p_files:
            name = os.path.basename(path)
            try:
                run_number = int(name.split(".run")[1].split(".")[0])
            except (IndexError, ValueError):
                run_number = 1
            by_run.setdefault(run_number, []).append(str(path))
        for run_number in sorted(by_run):
            paths = sorted(by_run[run_number])
            items.append(DatasetItemSpec(
                loci_name="BEST run %d" % run_number,
                item_type="chain",
                file_path=paths[0],
                format="mrbayes",
                data={"file_paths": paths, "run_number": run_number,
                      "chain_count": len(paths), "tool": "mrbayes"},
                metadata={"source": "best", "tool": "MrBayes",
                          "analysis": "BEST/MSC"},
            ))
        if not items:
            raise ValueError("BEST produced no tree or chain item.")
        dataset = DatasetSpec(
            name="BEST / MrBayes MSC result",
            items=items,
            description="%d tree/chain item(s)" % len(items),
            source="best",
            summary={"items": len(items)},
            metadata={
                "input": self.input_edit.text().strip(),
                "config": self.agent_config(),
            },
        )
        return ResultSpec(source="best", datasets=[dataset], artifacts=artifacts)

    def import_to_platform(self):
        if not self.result_has_consumer():
            QMessageBox.warning(
                self, "No workspace",
                "This BEST window is not connected to a YR-MPE workspace."
            )
            return
        try:
            result = self._build_result_spec()
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        if self.export_result(result):
            self.import_button.setEnabled(False)
            QMessageBox.information(
                self, "Success", "BEST result sent to the current workspace."
            )

    def result_has_consumer(self):
        try:
            return self.receivers(self.result_ready) > 0
        except Exception:
            return False

    def export_result(self, result):
        problems = result.validate() if hasattr(result, "validate") else []
        if problems:
            return False
        self.result_ready.emit(result)
        return True

    def agent_config(self):
        config = self._config()
        return {
            "partitions": [
                {"name": name, "range": range_text, "group": group}
                for name, range_text, group in config.partitions
            ],
            "species_partition": config.species_partition,
            "ploidy": config.ploidy,
            "species_prior": config.species_prior,
            "clock_model": config.clock_model,
            "clock_rate_prior": config.clock_rate_prior,
            "partition_rate_prior": config.partition_rate_prior,
            "gene_rate_prior": config.gene_rate_prior,
            "speciation_prior": config.speciation_prior,
            "extinction_prior": config.extinction_prior,
            "sample_prob": config.sample_prob,
            "popsize_prior": config.popsize_prior,
            "ngen": config.ngen,
            "samplefreq": config.samplefreq,
            "printfreq": config.printfreq,
            "nruns": config.nruns,
            "nchains": config.nchains,
            "extra_lines": list(config.extra_lines),
        }

    def _apply_config(self, config):
        self.partition_table.setRowCount(len(config.partitions))
        for row, (name, range_text, group) in enumerate(config.partitions):
            self.partition_table.setItem(row, 0, QTableWidgetItem(name))
            self.partition_table.setItem(row, 1, QTableWidgetItem(range_text))
            self.partition_table.setItem(row, 2, QTableWidgetItem(str(group)))
        self.species_edit.setPlainText(config.species_partition)
        self.ploidy_combo.setCurrentText(config.ploidy)
        self.prior_combo.setCurrentIndex(
            0 if config.species_prior == "yule" else 1
        )
        self.clock_rate_edit.setText(config.clock_rate_prior)
        self.partition_rate_combo.setCurrentText(config.partition_rate_prior)
        self.gene_rate_combo.setCurrentText(config.gene_rate_prior)
        self.speciation_edit.setText(config.speciation_prior)
        self.extinction_edit.setText(config.extinction_prior)
        self.sample_edit.setText(config.sample_prob)
        self.popsize_edit.setText(config.popsize_prior)
        self.ngen_spin.setValue(config.ngen)
        self.sample_spin.setValue(config.samplefreq)
        self.printfreq_spin.setValue(config.printfreq)
        self.nruns_spin.setValue(config.nruns)
        self.nchains_spin.setValue(config.nchains)
        self._extra_lines = list(config.extra_lines)
        self.refresh_preview()

    def apply_options(self, options):
        config = self._config()
        if not isinstance(options, dict):
            raise ValueError("options must be an object")
        if "partitions" in options:
            partitions = []
            for row in options["partitions"] or []:
                if isinstance(row, (list, tuple)) and len(row) >= 1:
                    partitions.append((str(row[0]), str(row[1]) if len(row) > 1 else "",
                                       int(row[2]) if len(row) > 2 else 1))
                elif isinstance(row, dict):
                    partitions.append((str(row.get("name", "")),
                                       str(row.get("range", "")),
                                       int(row.get("group", 1))))
            config.partitions = partitions
        mapping = {
            "species_partition": "species_partition",
            "ploidy": "ploidy",
            "species_prior": "species_prior",
            "clock_rate_prior": "clock_rate_prior",
            "partition_rate_prior": "partition_rate_prior",
            "gene_rate_prior": "gene_rate_prior",
            "speciation_prior": "speciation_prior",
            "extinction_prior": "extinction_prior",
            "sample_prob": "sample_prob",
            "popsize_prior": "popsize_prior",
            "ngen": "ngen",
            "samplefreq": "samplefreq",
            "printfreq": "printfreq",
            "nruns": "nruns",
            "nchains": "nchains",
            "extra_lines": "extra_lines",
        }
        for key, attr in mapping.items():
            if key in options:
                value = options[key]
                if key == "extra_lines":
                    value = [str(line) for line in value]
                setattr(config, attr, value)
        self._apply_config(config)
        return self.agent_config()

    def apply_block(self, block):
        config = parse_mrbayes_block(block)
        self._apply_config(config)
        canonical = build_mrbayes_block(self._config())
        return {
            "config": self.agent_config(),
            "canonical_block": canonical,
            "unparsed_lines": list(config.extra_lines),
            "roundtrip_equal": (
                re.sub(r"\s+", " ", block.strip().lower())
                == re.sub(r"\s+", " ", canonical.strip().lower())
            ),
        }

    def control_bindings(self):
        return [
            {"id": "input", "widget": "input_edit", "label": "NEXUS input",
             "kind": "path", "target": "mr_bayes", "cli": [],
             "notes": "Data file; the generated MRBAYES block is appended."},
            {"id": "partitions", "widget": "partition_table",
             "label": "Partitions", "kind": "table", "target": "mr_bayes",
             "cli": [], "notes": "charset name, range and link group."},
            {"id": "species_partition", "widget": "species_edit",
             "label": "Species partition", "kind": "text", "target": "mr_bayes",
             "cli": []},
            {"id": "ploidy", "widget": "ploidy_combo", "label": "Ploidy",
             "kind": "choice", "target": "mr_bayes", "cli": []},
            {"id": "species_prior", "widget": "prior_combo",
             "label": "Species-tree prior", "kind": "choice",
             "target": "mr_bayes", "cli": []},
            {"id": "clock_rate_prior", "widget": "clock_rate_edit",
             "label": "Clock-rate prior", "kind": "text",
             "target": "mr_bayes", "cli": []},
            {"id": "partition_rate_prior", "widget": "partition_rate_combo",
             "label": "Partition rate prior", "kind": "choice",
             "target": "mr_bayes", "cli": []},
            {"id": "gene_rate_prior", "widget": "gene_rate_combo",
             "label": "Gene rate prior", "kind": "choice",
             "target": "mr_bayes", "cli": []},
            {"id": "speciation_prior", "widget": "speciation_edit",
             "label": "Speciation prior", "kind": "text",
             "target": "mr_bayes", "cli": []},
            {"id": "extinction_prior", "widget": "extinction_edit",
             "label": "Extinction prior", "kind": "text",
             "target": "mr_bayes", "cli": []},
            {"id": "sample_prob", "widget": "sample_edit",
             "label": "Sampling probability", "kind": "text",
             "target": "mr_bayes", "cli": []},
            {"id": "popsize_prior", "widget": "popsize_edit",
             "label": "Population-size prior", "kind": "text",
             "target": "mr_bayes", "cli": []},
            {"id": "ngen", "widget": "ngen_spin", "label": "Generations",
             "kind": "int", "target": "mr_bayes", "cli": []},
            {"id": "samplefreq", "widget": "sample_spin",
             "label": "Sample frequency", "kind": "int",
             "target": "mr_bayes", "cli": []},
            {"id": "printfreq", "widget": "printfreq_spin",
             "label": "Print frequency", "kind": "int",
             "target": "mr_bayes", "cli": []},
            {"id": "nruns", "widget": "nruns_spin", "label": "Runs",
             "kind": "int", "target": "mr_bayes", "cli": []},
            {"id": "nchains", "widget": "nchains_spin", "label": "Chains",
             "kind": "int", "target": "mr_bayes", "cli": []},
        ]

    def cli_preview(self):
        return ["MRBAYES block generated for MrBayes-MPI-BEAGLE"]

    def help_doc(self):
        return {
            "schema": "yrmpe.skill.v1",
            "id": "best-mrbayes",
            "title": "BEST / MrBayes Multi-Species Coalescent",
            "summary": "Full-likelihood multispecies-coalescent inference "
                       "through MrBayes, with block import/export.",
            "markdown": (
                "Use `best.apply_block` to send a complete MRBAYES block. "
                "The response contains the canonical block regenerated from "
                "the form and the lines that were not parsed. Use "
                "`best.modify_options` for individual controls, then "
                "`best.run` to execute MrBayes."
            ),
            "controls": self.control_bindings(),
            "cli": {"preview": self.cli_preview(), "target": "mr_bayes"},
        }

    def agent_status(self):
        job_id = str(getattr(self, "queue_job_id", "") or "")
        job_state = ""
        queue = getattr(self, "task_queue", None)
        if queue is not None and job_id and callable(getattr(queue, "get_job", None)):
            job = queue.get_job(job_id)
            if job is not None:
                job_state = str(getattr(job, "state", "") or "")
        return {
            "plugin": self.plugin_name,
            "tool": self.tool_name,
            "tool_path": self.tool_path,
            "input": self.input_edit.text().strip(),
            "is_running": bool(self.is_running),
            "job_id": job_id,
            "job_state": job_state,
            "outputs": list(self.current_output_files),
            "last_error": self._last_error,
        }

    def agent_contract(self):
        return {
            "schema": "yrmpe.plugin.v1",
            "id": "best-mrbayes",
            "status": self.agent_status(),
            "actions": {"actions": [
                {"id": "best.inspect_options"},
                {"id": "best.modify_options"},
                {"id": "best.apply_block"},
                {"id": "best.preview_block"},
                {"id": "best.write_block"},
                {"id": "best.set_input"},
                {"id": "best.auto_species"},
                {"id": "best.close_and_open_mrbayes"},
                {"id": "best.run"},
                {"id": "best.stop"},
                {"id": "best.export_to_workspace"},
            ]},
        }

    def agent_invoke(self, action_id, params=None):
        params = params or {}
        name = str(action_id or "").rsplit(".", 1)[-1]
        try:
            if name == "help_doc":
                return {"ok": True, "result": self.help_doc()}
            if name == "inspect_status":
                return {"ok": True, "result": self.agent_status()}
            if name in ("inspect_options", "inspect_controls"):
                return {"ok": True, "result": {
                    "options": self.agent_config(),
                    "block": build_mrbayes_block(self._config()),
                    "controls": self.control_bindings(),
                }}
            if name == "modify_options":
                options = params.get("options", params)
                result = self.apply_options(options)
                return {"ok": True, "result": {
                    "options": result,
                    "canonical_block": build_mrbayes_block(self._config()),
                }}
            if name == "apply_block":
                block = params.get("block", params.get("text", ""))
                if not str(block).strip():
                    return {"ok": False, "error": "provide a block string"}
                return {"ok": True, "result": self.apply_block(str(block))}
            if name == "auto_species":
                self._auto_species_partition()
                return {"ok": True, "result": {"species_partition": self.species_edit.toPlainText()}}
            if name == "close_and_open_mrbayes":
                self._close_and_open_mrbayes()
                return {"ok": True, "result": {"closed": True}}
            if name == "preview_block":
                return {"ok": True, "result": {
                    "block": build_mrbayes_block(self._config()),
                    "unparsed_lines": list(self._extra_lines),
                }}
            if name == "write_block":
                path = str(params.get("path", "") or "")
                if not path:
                    return {"ok": False, "error": "provide path"}
                with open(path, "w", encoding="utf-8") as handle:
                    handle.write(build_mrbayes_block(self._config()))
                return {"ok": True, "result": {"path": path}}
            if name == "set_input":
                path = str(params.get("path", "") or "")
                if not path or not os.path.isfile(path):
                    return {"ok": False, "error": "input file not found"}
                self.input_edit.setText(path)
                if params.get("load_partitions"):
                    self._load_partitions()
                return {"ok": True, "result": {"input": path}}
            if name == "run":
                parsed = None
                if "block" in params:
                    parsed = self.apply_block(str(params["block"]))
                if "options" in params:
                    self.apply_options(params["options"])
                if not self.input_edit.text().strip():
                    return {"ok": False, "error": "set an input NEXUS first"}
                job_id = self._request_run()
                result = {
                    "job_id": job_id or "",
                    "status": self.agent_status(),
                }
                if parsed is not None:
                    # Return the form-regenerated block even when the run is
                    # started, so the caller can diff what was understood.
                    result.update(parsed)
                return {"ok": True, "result": result}
            if name == "stop":
                self._request_stop()
                return {"ok": True, "result": self.agent_status()}
            if name == "inspect_console":
                limit = int(params.get("limit", 200) or 200)
                return {"ok": True, "result": {
                    "lines": self._console[-limit:],
                    "total_lines": len(self._console),
                }}
            if name == "export_to_workspace":
                if not self.result_has_consumer():
                    return {"ok": False, "error": "no workspace result channel"}
                result = self._build_result_spec()
                if not self.export_result(result):
                    return {"ok": False, "error": "ResultSpec validation failed"}
                return {"ok": True, "result": result.to_settings_dict()}
            if name == "close":
                if self.is_running:
                    return {"ok": False, "error": "BEST is running"}
                self.close()
                return {"ok": True, "result": {"closed": True}}
            return {"ok": False, "error": "unknown BEST action %r" % action_id}
        except Exception as exc:
            return {"ok": False, "error": str(exc)}

    def _write_block(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Write MRBAYES block", "best_mrbayes.nex",
            "NEXUS (*.nex *.nexus);;All files (*)",
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(build_mrbayes_block(self._config()))
        except OSError as exc:
            QMessageBox.warning(self, "Write failed", str(exc))
            return
        QMessageBox.information(
            self, "Written",
            "MRBAYES block written. Append it to the data file (or execute it "
            "as a separate NEXUS block) and run it with MrBayes.",
        )


class BestPluginEntry:
    def run(self, **kwargs):
        return BestDialog()
