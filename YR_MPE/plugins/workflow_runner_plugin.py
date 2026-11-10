"""Execute a compiled ``yrmpe.workflow.v1`` plan step by step.

The plugin is a thin host around the existing plugin windows:

* Input & Parameters: the JSONL spec, an output root and one path field per
  declared workflow input role.
* Procedures: the compiled DAG (``graph_from_plan``).
* Output Preview: one entry per produced artifact, with a text preview.
* Console: the command line and console output of every step.

Execution is sequential and stays on the Qt main thread. Each step creates
the ordinary plugin window with ``headless=True``, applies the step
parameters through ``control_bindings``, calls ``run_analysis()`` and waits
for ``is_running`` to fall. The typed result is collected with
``result_spec()`` and the produced files become the bindings for later
placeholders. A single merged ResultSpec is exported at the end.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from PyQt5.QtCore import QObject, QTimer
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from ..runtime.plugin_agent import apply_control
from ..runtime.results import DatasetItemSpec, DatasetSpec, ResultSpec
from ..runtime.workflow_compiler import (
    WorkflowCompileError,
    WorkflowPlan,
    _topological_order,
    compile_jsonl,
    compile_workflow,
)
from ..runtime.workflow_graph import graph_from_plan
from ..runtime.workflow_tools import create_plugin, entry_for_op
from ..templates.base_plugin_ui import BasePlugin


class WorkflowRunError(RuntimeError):
    pass


class WorkflowRunController(QObject):
    """Sequential executor for one ``WorkflowPlan``."""

    POLL_MS = 200

    def __init__(self, host: "WorkflowRunnerPlugin", plan: WorkflowPlan,
                 bindings: Dict[str, str], output_root: str):
        super().__init__(host)
        self.host = host
        self.plan = plan
        self.bindings = dict(bindings)
        self.output_root = os.path.abspath(output_root)
        # Keep the raw steps: placeholders must be re-resolved before each
        # step, because an earlier step replaces its declared output path
        # with the file the plugin actually produced.
        self.ordered = [dict(step) for step in _topological_order(plan.steps)]
        self.index = 0
        self.current_plugin = None
        self.current_step: Optional[Dict[str, Any]] = None
        self.started_at = 0.0
        self.step_log: List[Dict[str, Any]] = []
        self.items: List[DatasetItemSpec] = []
        self.outputs: List[Dict[str, str]] = []
        self._cancelled = False
        self._deadline = 0.0

    # ------------------------------------------------------------------
    # lifecycle
    # ------------------------------------------------------------------
    def start(self):
        self.host.add_console_message(
            "Workflow %s: %d step(s)" % (self.plan.name or self.plan.id,
                                          len(self.ordered)),
            "info",
        )
        QTimer.singleShot(0, self._run_next)

    def stop(self):
        self._cancelled = True
        plugin = self.current_plugin
        if plugin is not None:
            for name in ("stop_analysis", "stop"):
                method = getattr(plugin, name, None)
                if callable(method):
                    try:
                        method()
                    except Exception:
                        pass
                    break

    def _run_next(self):
        if self._cancelled:
            self.host._workflow_finished(False, "stopped")
            return
        if self.index >= len(self.ordered):
            self._finish()
            return

        step = self.ordered[self.index]
        step["params"] = self.plan.resolve(step.get("params") or {}, self.bindings)
        step["outputs"] = self.plan.resolve(step.get("outputs") or {}, self.bindings)
        self.current_step = step
        step_id = str(step.get("id") or "step_%d" % (self.index + 1))
        op = str(step.get("op") or step.get("method") or "")
        self.host.add_console_message(
            "[%s] %s" % (step_id, op), "info"
        )

        if entry_for_op(op) is None:
            self._fail(step_id, "no plugin entry for operation %r" % op)
            return

        step_dir = os.path.join(self.output_root, _safe_dir(step_id))
        os.makedirs(step_dir, exist_ok=True)
        try:
            plugin = create_plugin(op, workdir=step_dir, headless=True)
        except Exception as exc:
            self._fail(step_id, "could not create plugin for %s: %s" % (op, exc))
            return
        self.current_plugin = plugin
        plugin.setParent(self.host)
        plugin.hide()

        try:
            self._configure(plugin, step)
        except Exception as exc:
            self._fail(step_id, "could not configure %s: %s" % (op, exc))
            return

        command = _cli_preview(plugin)
        self.add_step_log(step_id, op, "running", command)
        if command:
            self.host.add_console_message(" ".join(command), "command")

        self.started_at = time.time()
        self._deadline = self.started_at + 6 * 3600
        try:
            plugin.run_analysis()
        except Exception as exc:
            self._fail(step_id, "run_analysis failed: %s" % exc)
            return
        QTimer.singleShot(self.POLL_MS, self._poll)

    def _configure(self, plugin, step):
        params = dict(step.get("params") or {})
        # A workflow input may be declared under a role that is not the
        # control id; map it to the plugin's path control when possible.
        for control in plugin.help_doc().get("controls", []):
            if control.get("kind") != "path":
                continue
            key = control.get("id", "")
            if key in params:
                continue
            for value in (step.get("params") or {}).values():
                if isinstance(value, str) and value in self.bindings.values():
                    params[key] = value
                    break
        for key, value in params.items():
            apply_control(plugin, key, value)

    # ------------------------------------------------------------------
    # polling / result collection
    # ------------------------------------------------------------------
    def _poll(self):
        if self._cancelled:
            self.host._workflow_finished(False, "stopped")
            return
        plugin = self.current_plugin
        if plugin is None:
            return
        if time.time() > self._deadline:
            self._fail(self.current_step.get("id", ""), "step timed out")
            return
        if getattr(plugin, "is_running", False):
            QTimer.singleShot(self.POLL_MS, self._poll)
            return
        QTimer.singleShot(100, self._collect)

    def _collect(self):
        plugin = self.current_plugin
        step = self.current_step or {}
        step_id = str(step.get("id") or "step_%d" % (self.index + 1))
        op = str(step.get("op") or step.get("method") or "")
        try:
            spec = plugin.result_spec()
        except Exception as exc:
            self._fail(step_id, "result construction failed: %s" % exc)
            return
        if spec is None:
            self._fail(step_id, "the plugin produced no typed result")
            return

        files = _spec_files(spec)
        for role, declared in (step.get("outputs") or {}).items():
            actual = _match_output(files, role, declared)
            if not actual:
                continue
            placeholder = _artifact_placeholder(self.plan, step_id, role)
            if placeholder:
                self.bindings[placeholder] = actual
            self.outputs.append({
                "step": step_id, "role": str(role),
                "path": actual,
            })
            self.host.add_output_entry(step_id, role, actual)

        for dataset in spec.datasets:
            for item in dataset.items:
                item.metadata = dict(item.metadata or {})
                item.metadata.setdefault("workflow_step", step_id)
                item.metadata.setdefault("workflow_op", op)
                if not item.loci_name.startswith(step_id + ":"):
                    item.loci_name = "%s:%s" % (step_id, item.loci_name)
                self.items.append(item)

        self.add_step_log(step_id, op, "completed", [], files)
        self.host.add_console_message(
            "[%s] completed (%d output file(s))" % (step_id, len(files)),
            "success",
        )
        self._close_current()
        self.index += 1
        QTimer.singleShot(0, self._run_next)

    def _fail(self, step_id: str, message: str):
        self.add_step_log(str(step_id), "", "failed", [], [], message)
        self.host.add_console_message(
            "[%s] failed: %s" % (step_id, message), "error"
        )
        self._close_current()
        self.host._workflow_finished(False, message)

    def _close_current(self):
        plugin = self.current_plugin
        self.current_plugin = None
        if plugin is not None:
            try:
                # Do not call close(): closing the last hidden plugin can
                # trigger QApplication.quitOnLastWindowClosed while the
                # workflow still has steps to run.
                plugin.hide()
                plugin.deleteLater()
            except Exception:
                pass

    def _finish(self):
        outputs = list(self.outputs)
        self.host._workflow_finished(True, "", self.items, outputs, self.step_log)

    def add_step_log(self, step_id, op, status, command=None,
                     outputs=None, error=""):
        self.step_log.append({
            "step": step_id,
            "op": op,
            "status": status,
            "command": list(command or []),
            "outputs": list(outputs or []),
            "error": error,
        })


def _safe_dir(value: str) -> str:
    text = "".join(ch if ch.isalnum() or ch in "._-" else "_"
                   for ch in str(value))
    return text or "step"


def _cli_preview(plugin) -> List[str]:
    try:
        doc = plugin.help_doc()
        preview = (doc.get("cli") or {}).get("preview") or []
        return [str(part) for part in preview]
    except Exception:
        return []


def _spec_files(spec: ResultSpec) -> List[str]:
    files = []
    for dataset in spec.datasets:
        for item in dataset.items:
            path = str(getattr(item, "file_path", "") or "")
            if path and path not in files:
                files.append(path)
    return files


def _match_output(files: List[str], role: str, declared: str) -> str:
    """Pick the file that best matches an output role."""
    if not files:
        return ""
    wanted = Path(str(declared or "")).suffix.lower()
    role_text = str(role).lower()
    for path in files:
        suffix = Path(path).suffix.lower()
        if wanted and suffix == wanted:
            return path
        if role_text and role_text in Path(path).stem.lower():
            return path
    return files[0]


def _artifact_placeholder(plan: WorkflowPlan, step_id: str,
                          role: str) -> str:
    for placeholder, artifact in plan.artifacts.items():
        if str(artifact.get("step")) == step_id and str(artifact.get("role")) == str(role):
            return placeholder
    return ""


class WorkflowRunnerPlugin(BasePlugin):
    """Plugin shell for a compiled workflow."""

    def __init__(self, import_from=None, import_data=None, **kwargs):
        self.plan: Optional[WorkflowPlan] = None
        self._controller: Optional[WorkflowRunController] = None
        self._outputs: List[Dict[str, str]] = []
        self._step_log: List[Dict[str, Any]] = []
        self._last_workflow_error = ""
        self._initial_workflow = kwargs.pop("workflow_path", "") or ""
        super().__init__(import_from=import_from, import_data=import_data,
                         **kwargs)
        self.setup_procedures_tab()
        if hasattr(self, "run_button"):
            self.run_button.setText("Run workflow")
        if self._initial_workflow:
            self.workflow_path_edit.setText(str(self._initial_workflow))
            self.load_workflow()

    # ------------------------------------------------------------------
    # plugin metadata
    # ------------------------------------------------------------------
    def init_plugin_info(self):
        super().init_plugin_info()
        self.plugin_name = "Workflow Runner"
        self.tool_name = "YR-MPE workflow engine"
        self.citation_keys = []
        self.input_types = {"JSONL workflow": ["jsonl", "json"]}
        self.output_types = {"Workspace items": ""}

    def config(self):
        """The runner has no external binary of its own."""
        return True

    def control_bindings(self):
        return [
            {"id": "workflow_file", "widget": "workflow_path_edit",
             "label": "Workflow JSONL", "kind": "path", "cli": [],
             "notes": "A yrmpe.workflow.v1 specification or an agent JSONL transcript."},
            {"id": "output_dir", "widget": "output_dir_edit",
             "label": "Output directory", "kind": "path", "cli": [],
             "notes": "Step artifacts are written below this directory."},
        ]

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def setup_input_tab(self):
        from PyQt5.QtWidgets import QVBoxLayout

        outer = QVBoxLayout(self.input_tab)
        outer.setContentsMargins(12, 12, 12, 12)
        title = QLabel("Workflow Runner")
        title.setStyleSheet("font-size: 17px; font-weight: 600;")
        outer.addWidget(title)
        subtitle = QLabel(
            "Compile a JSONL workflow, map its declared inputs, and execute "
            "the steps in order."
        )
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: palette(mid);")
        outer.addWidget(subtitle)

        group = QGroupBox("Workflow")
        form = QFormLayout(group)
        self.workflow_path_edit = QLineEdit()
        self.workflow_path_edit.setPlaceholderText(
            "workflow.jsonl or an agent transcript"
        )
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse_workflow)
        row = QHBoxLayout()
        row.addWidget(self.workflow_path_edit, 1)
        row.addWidget(browse)
        form.addRow("Workflow file:", row)

        self.output_dir_edit = QLineEdit()
        self.output_dir_edit.setPlaceholderText(
            "Directory for step outputs (default: workdir/workflow)"
        )
        out_browse = QPushButton("Browse...")
        out_browse.clicked.connect(self._browse_output)
        out_row = QHBoxLayout()
        out_row.addWidget(self.output_dir_edit, 1)
        out_row.addWidget(out_browse)
        form.addRow("Output directory:", out_row)

        self.validate_button = QPushButton("Load && Validate")
        self.validate_button.clicked.connect(self.load_workflow)
        form.addRow("", self.validate_button)
        outer.addWidget(group)

        self.inputs_group = QGroupBox("Workflow inputs")
        self.inputs_form = QFormLayout(self.inputs_group)
        outer.addWidget(self.inputs_group)

        self.validation_label = QLabel("No workflow loaded.")
        self.validation_label.setWordWrap(True)
        self.validation_label.setStyleSheet("color: palette(mid);")
        outer.addWidget(self.validation_label)
        outer.addStretch(1)

        if self.workdir:
            self.output_dir_edit.setText(
                str(Path(self.workdir) / "workflow")
            )

    def setup_output_tab(self):
        layout = QVBoxLayout(self.output_tab)
        layout.setContentsMargins(12, 12, 12, 12)
        self.output_combo = QComboBox()
        self.output_combo.currentIndexChanged.connect(self._show_output)
        layout.addWidget(self.output_combo)
        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)
        font = QFont("Cascadia Code")
        font.setStyleHint(QFont.Monospace)
        self.output_text.setFont(font)
        layout.addWidget(self.output_text, 1)
        self.export_button = QPushButton("Send results to workspace")
        self.export_button.clicked.connect(self._export_workflow_result)
        self.export_button.setEnabled(False)
        layout.addWidget(self.export_button)

    def setup_procedures_tab(self):
        from ..platforms.methods.workflow_preview import WorkflowPreviewWidget

        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 0, 0)
        self.procedures_note = QLabel("Load a workflow to see its steps.")
        self.procedures_note.setWordWrap(True)
        self.procedures_note.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.procedures_note)
        self.procedures_preview = WorkflowPreviewWidget()
        layout.addWidget(self.procedures_preview, 1)
        self.tab_widget.insertTab(1, page, "Procedures")

    # ------------------------------------------------------------------
    # loading / validation
    # ------------------------------------------------------------------
    def _browse_workflow(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select workflow", "",
            "Workflow files (*.jsonl *.json);;All files (*)",
        )
        if path:
            self.workflow_path_edit.setText(path)
            self.load_workflow()

    def _browse_output(self):
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", self.output_dir_edit.text()
        )
        if path:
            self.output_dir_edit.setText(path)

    def load_workflow(self):
        path = self.workflow_path_edit.text().strip()
        if not path:
            self._set_validation("Select a workflow file first.", error=True)
            return None
        try:
            if path.lower().endswith(".jsonl"):
                plan = compile_jsonl(path)
            else:
                payload = json.loads(Path(path).read_text(encoding="utf-8"))
                if (payload.get("schema") == "yrmpe.workflow.v1"
                        and "steps" in payload and "artifacts" in payload):
                    plan = WorkflowPlan.from_dict(payload)
                else:
                    plan = compile_workflow(payload)
        except (WorkflowCompileError, Exception) as exc:
            self.plan = None
            self.procedures_preview.set_graph(graph_from_plan(_empty_plan()))
            self._set_validation("Workflow compile failed: %s" % exc, error=True)
            self.add_console_message("Workflow compile failed: %s" % exc, "error")
            return None

        self.plan = plan
        self.procedures_preview.set_graph(graph_from_plan(plan))
        self.procedures_note.setText(
            "%s: %d step(s), %d declared input(s)"
            % (plan.name or plan.id, len(plan.steps), len(plan.inputs))
        )
        self._populate_inputs(plan)
        problems = self.validate_workflow(plan)
        if problems:
            self._set_validation("; ".join(problems), error=True)
            self.add_console_message("Workflow validation: " + "; ".join(problems), "warning")
        else:
            self._set_validation(
                "Workflow loaded: %d step(s), ready to run." % len(plan.steps)
            )
            self.add_console_message("Workflow validated.", "success")
        return plan

    def _populate_inputs(self, plan: WorkflowPlan):
        while self.inputs_form.rowCount():
            self.inputs_form.removeRow(0)
        self.input_edits = {}
        for name, meta in plan.inputs.items():
            edit = QLineEdit(str(meta.get("original", "")))
            browse = QPushButton("Browse...")
            browse.clicked.connect(
                lambda _=False, field=edit: self._browse_input(field)
            )
            row = QHBoxLayout()
            row.addWidget(edit, 1)
            row.addWidget(browse)
            self.inputs_form.addRow("%s:" % name, row)
            self.input_edits[name] = edit

    def _browse_input(self, field):
        path, _ = QFileDialog.getOpenFileName(self, "Select input file")
        if path:
            field.setText(path)

    def validate_workflow(self, plan: Optional[WorkflowPlan] = None):
        plan = plan or self.plan
        problems: List[str] = []
        if plan is None:
            return ["no workflow loaded"]
        if not plan.steps:
            problems.append("workflow has no steps")
        for step in plan.steps:
            op = step.get("op") or step.get("method") or ""
            if entry_for_op(op) is None:
                problems.append("no plugin entry for %r" % op)
        for name, edit in getattr(self, "input_edits", {}).items():
            path = edit.text().strip()
            if not path:
                problems.append("input %r has no file" % name)
            elif not os.path.exists(path):
                problems.append("input %r does not exist: %s" % (name, path))
        output_root = self._output_root()
        try:
            os.makedirs(output_root, exist_ok=True)
        except Exception as exc:
            problems.append("output directory is not writable: %s" % exc)
        return problems

    def _output_root(self) -> str:
        text = self.output_dir_edit.text().strip()
        if text:
            return os.path.abspath(text)
        base = self.workdir or os.getcwd()
        return os.path.abspath(os.path.join(base, "workflow"))

    def _set_validation(self, text, error=False):
        self.validation_label.setText(text)
        self.validation_label.setStyleSheet(
            "color: %s;" % ("#b91c1c" if error else "#166534")
        )

    # ------------------------------------------------------------------
    # execution
    # ------------------------------------------------------------------
    def run_analysis(self):
        if getattr(self, "is_running", False):
            return
        plan = self.plan or self.load_workflow()
        if plan is None:
            return
        problems = self.validate_workflow(plan)
        if problems:
            self._set_validation("; ".join(problems), error=True)
            self.add_console_message(
                "Workflow not started: " + "; ".join(problems), "error"
            )
            return

        bindings = plan.bindings(self._output_root())
        for name, edit in getattr(self, "input_edits", {}).items():
            meta = plan.inputs.get(name, {})
            placeholder = meta.get("placeholder", "")
            if placeholder:
                bindings[placeholder] = edit.text().strip()

        if getattr(self, "headless", False):
            from PyQt5.QtWidgets import QApplication
            app = QApplication.instance()
            if app is not None:
                app.setQuitOnLastWindowClosed(False)

        self._outputs = []
        self._step_log = []
        self.output_combo.clear()
        self.output_text.clear()
        self.export_button.setEnabled(False)
        self.is_running = True
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)

        self._controller = WorkflowRunController(
            self, plan, bindings, self._output_root()
        )
        self._controller.start()

    def stop_analysis(self):
        if self._controller is not None:
            self._controller.stop()

    def _workflow_finished(self, ok: bool, message: str = "",
                           items=None, outputs=None, step_log=None):
        self.is_running = False
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.progress_bar.setVisible(False)
        self._controller = None
        if outputs is not None:
            self._outputs = list(outputs)
        if step_log is not None:
            self._step_log = list(step_log)
        if not ok:
            self._last_workflow_error = message
            self._set_validation("Workflow stopped: %s" % message, error=True)
            self.add_console_message("Workflow stopped: %s" % message, "error")
            self._write_manifest(False)
            return

        self._last_workflow_error = ""
        self._pending_items = list(items or [])
        completed = sum(
            1 for entry in self._step_log if entry.get("status") == "completed"
        )
        self._set_validation(
            "Workflow completed: %d step(s), %d output file(s)."
            % (completed, len(self._outputs))
        )
        self.add_console_message("Workflow completed.", "success")
        self._write_manifest(True)
        if self._pending_items:
            self.export_button.setEnabled(True)
            self._export_workflow_result()

    def _write_manifest(self, ok: bool):
        if self.plan is None:
            return
        payload = {
            "schema": "yrmpe.workflow.run.v1",
            "plan": self.plan.to_dict(),
            "ok": bool(ok),
            "error": self._last_workflow_error,
            "outputs": list(self._outputs),
            "steps": list(self._step_log),
        }
        try:
            path = Path(self._output_root()) / "workflow_run.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                json.dumps(payload, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            self.add_console_message("Run manifest: %s" % path, "info")
        except Exception as exc:
            self.add_console_message(
                "Could not write workflow manifest: %s" % exc, "warning"
            )

    def _build_result_spec(self):
        items = list(getattr(self, "_pending_items", []) or [])
        if not items:
            raise ValueError("No workflow output is available.")
        return ResultSpec(
            source="workflow-runner",
            datasets=[DatasetSpec(
                name=(self.plan.name if self.plan else "Workflow run")
                     or "Workflow run",
                items=items,
                description="Outputs from a compiled workflow",
                source="workflow-runner",
                metadata={"workflow_id": getattr(self.plan, "id", "")},
            )],
            artifacts=[entry["path"] for entry in self._outputs],
            summary={"steps": len(self._step_log)},
        )

    def _export_workflow_result(self):
        try:
            spec = self._build_result_spec()
        except ValueError:
            return False
        if self.export_result(spec):
            self.add_console_message(
                "Workflow results sent to the workspace.", "success"
            )
            self.export_button.setEnabled(False)
            return True
        return False

    # ------------------------------------------------------------------
    # output preview
    # ------------------------------------------------------------------
    def add_output_entry(self, step_id: str, role: str, path: str):
        label = "%s.%s  —  %s" % (step_id, role, Path(path).name)
        self.output_combo.addItem(label, path)
        self._outputs.append({
            "step": step_id, "role": str(role), "path": path,
        })
        self.analysis_output_files = [entry["path"] for entry in self._outputs]

    def _show_output(self, index: int):
        if index < 0:
            self.output_text.clear()
            return
        path = self.output_combo.itemData(index)
        if not path:
            self.output_text.clear()
            return
        try:
            text = Path(str(path)).read_text(
                encoding="utf-8", errors="replace"
            )
        except Exception as exc:
            self.output_text.setPlainText("Cannot read %s: %s" % (path, exc))
            return
        if len(text) > 200000:
            text = text[:200000] + "\n... (truncated)"
        self.output_text.setPlainText(text)

    # ------------------------------------------------------------------
    # agent surface
    # ------------------------------------------------------------------
    def inspect_plan(self):
        if self.plan is None:
            return {"loaded": False}
        return {
            "loaded": True,
            "id": self.plan.id,
            "name": self.plan.name,
            "steps": [dict(step) for step in self.plan.steps],
            "inputs": {name: dict(meta)
                       for name, meta in self.plan.inputs.items()},
            "outputs": list(self._outputs),
            "status": self.inspect_status(),
        }

    def inspect_status(self):
        return {
            "loaded": self.plan is not None,
            "is_running": bool(getattr(self, "is_running", False)),
            "steps_total": len(self.plan.steps) if self.plan else 0,
            "steps_done": sum(
                1 for entry in self._step_log
                if entry.get("status") == "completed"
            ),
            "outputs": list(self._outputs),
            "error": self._last_workflow_error,
        }

    def agent_contract(self):
        return {
            "schema": "yrmpe.ui.workflow_runner.v1",
            "id": "workflow_runner",
            "plan": self.inspect_plan(),
            "status": self.inspect_status(),
            "actions": [
                {"id": "workflow.inspect_plan",
                 "description": "Return the compiled plan, inputs and outputs"},
                {"id": "workflow.validate",
                 "description": "Validate the current workflow"},
                {"id": "workflow.run",
                 "description": "Compile and execute the workflow"},
                {"id": "workflow.stop",
                 "description": "Stop the running workflow"},
                {"id": "workflow.inspect_status",
                 "description": "Return step progress"},
                {"id": "workflow.inspect_console",
                 "description": "Read the workflow console tail"},
            ],
        }

    def agent_invoke(self, action_id: str, params=None):
        params = params or {}
        name = str(action_id or "")
        if name.startswith("workflow."):
            name = name.split(".", 1)[1]
        if name == "inspect_plan":
            return {"ok": True, "result": self.inspect_plan()}
        if name == "inspect_status":
            return {"ok": True, "result": self.inspect_status()}
        if name == "validate":
            plan = self.plan or self.load_workflow()
            problems = self.validate_workflow(plan) if plan else ["no plan"]
            return {"ok": not problems, "result": {"problems": problems}}
        if name == "run":
            self.run_analysis()
            return {"ok": bool(getattr(self, "is_running", False)),
                    "status": self.inspect_status()}
        if name == "stop":
            self.stop_analysis()
            return {"ok": True, "status": self.inspect_status()}
        if name == "inspect_console":
            text = ""
            if hasattr(self, "console_text"):
                text = self.console_text.toPlainText()
            return {"ok": True, "result": {
                "lines": text.splitlines()[-200:],
                "total_lines": len(text.splitlines()),
            }}
        return super().agent_invoke(action_id, params)


def _empty_plan() -> WorkflowPlan:
    return WorkflowPlan(steps=[], inputs={}, artifacts={})


class WorkflowRunnerPluginEntry:
    def __init__(self, config=None, plugin_path=None):
        self.config = config
        self.plugin_path = plugin_path

    def run(self, import_from=None, import_data=None, workdir=None, **kwargs):
        return WorkflowRunnerPlugin(
            import_from=import_from,
            import_data=import_data,
            workdir=workdir,
            **kwargs,
        )
