"""About / help dialog and the copyable agent-bridge prompt."""

from __future__ import annotations

import json
import platform
import sys

from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QApplication,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from .. import __version__


def build_agent_prompt(bridge) -> str:
    """Return a self-contained instruction block for an external agent."""
    session = getattr(bridge, "session_dir", "")
    token = getattr(bridge, "token", "") or ""
    actions = getattr(getattr(bridge, "registry", None), "list_actions", lambda: [])()
    action_lines = []
    for action in actions:
        if not action.get("enabled", True):
            continue
        line = "- %s" % action["id"]
        if action.get("description"):
            line += ": %s" % action["description"]
        action_lines.append(line)
    if not action_lines:
        action_lines = ["- (read control.json for the current action list)"]

    command = {
        "schema": "yrmpe.bridge.command.v1",
        "id": "<unique-id>",
        "client": "<your-name>",
        "op": "<action-id>",
        "params": {},
        "token": token,
    }
    return (
        "You are controlling YPhylo Workbench through its local file bridge.\n"
        "Session directory: %s\n"
        "Agent token: %s\n\n"
        "Read/write these files:\n"
        "- Append one JSON object per line to %s/commands.jsonl.\n"
        "- Read matching responses from %s/responses.jsonl (match the \"id\").\n"
        "- Follow progress in %s/events.jsonl.\n"
        "- %s/READY exists while the bridge is alive; control.json lists actions "
        "and token_required.\n\n"
        "Command format (one line):\n%s\n\n"
        "Rules:\n"
        "- Include the exact token in every command; never write it to another file.\n"
        "- Generate a fresh id per command and wait for the response with that id.\n"
        "- Do not call arbitrary shell commands; use the semantic actions below.\n"
        "- Call workbench.list_windows once, then control a plugin through "
        "workbench.invoke_plugin with its window_id and native action id.\n"
        "- plugin.run returns a job_id; poll workbench.jobs or "
        "workbench.job_status instead of sleeping on a window. Use "
        "workbench.cancel_job to stop queued or running work.\n"
        "- A changed token means the Workbench restarted: copy the prompt again.\n\n"
        "Available actions:\n%s"
        % (
            session,
            token,
            session,
            session,
            session,
            session,
            json.dumps(command, ensure_ascii=False),
            "\n".join(action_lines),
        )
    )


class AboutDialog(QDialog):
    def __init__(
        self,
        parent=None,
        bridge=None,
        workspace_path="",
        repo_path="",
        start_session=None,
    ):
        super().__init__(parent)
        self.bridge = bridge
        self.start_session = start_session
        self.bridge_body = None
        self.bridge_error = ""
        self.setWindowTitle("About YPhylo Workbench")
        self.setMinimumWidth(760)
        outer = QVBoxLayout(self)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        content = QWidget()
        layout = QVBoxLayout(content)
        self.scroll.setWidget(content)
        outer.addWidget(self.scroll, 1)

        title = QLabel("<h2>YPhylo Workbench</h2>")
        layout.addWidget(title)
        subtitle = QLabel(
            "A workbench for YR-MPE plugins: alignment, model selection, "
            "distance, phylogeny and species-tree analysis. "
            "Plugin set version %s." % __version__
        )
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)

        help_group = QGroupBox("Quick start")
        help_layout = QVBoxLayout(help_group)
        help_text = QLabel(
            "1. WORKSPACE → import sequences, a dataset, or run One-Step OrthoFinder.\n"
            "2. ALIGN → align and trim; MODEL → choose a substitution model.\n"
            "3. DISTANCE / PHYLOGENY / SPECIES TREE → build and compare trees.\n"
            "4. Click a workspace item, or a dataset, to make it the active input; "
            "drag it onto Selected to queue it for a tool.\n"
            "5. Details shows metadata, workflow provenance and previews."
        )
        help_text.setWordWrap(True)
        help_layout.addWidget(help_text)
        layout.addWidget(help_group)

        environment = QGroupBox("Environment")
        form = QFormLayout(environment)
        form.addRow("Application", QLabel("YPhylo Workbench"))
        form.addRow("Plugin set", QLabel("YR-MPE %s" % __version__))
        form.addRow("Python", QLabel(sys.version.split()[0]))
        form.addRow("Platform", QLabel(platform.platform()))
        form.addRow("Workspace", QLabel(workspace_path or "Temporary"))
        if repo_path:
            form.addRow("Install path", QLabel(repo_path))
        layout.addWidget(environment)

        self.bridge_group = QGroupBox("Agent bridge")
        self.bridge_layout = QVBoxLayout(self.bridge_group)
        layout.addWidget(self.bridge_group)
        self._refresh_bridge_group()

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

        width, height = 820, 720
        screen = QApplication.primaryScreen()
        if screen is not None:
            available = screen.availableGeometry()
            width = min(width, max(760, available.width() - 80))
            height = min(height, max(480, available.height() - 80))
        self.resize(width, height)

    def _refresh_bridge_group(self):
        if self.bridge_body is not None:
            self.bridge_layout.removeWidget(self.bridge_body)
            self.bridge_body.hide()
            self.bridge_body.deleteLater()
            self.bridge_body = None

        body = QWidget()
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        self.bridge_body = body

        bridge = self.bridge
        token = getattr(bridge, "token", None) if bridge is not None else None
        if bridge is None:
            idle_note = QLabel(
                "The agent bridge is not running in this window. Start a local "
                "session to let an external agent drive the Workbench through "
                "semantic actions; no network port or shell access is opened."
            )
            idle_note.setWordWrap(True)
            body_layout.addWidget(idle_note)
            start_button = QPushButton("Start agent session")
            start_button.setEnabled(self.start_session is not None)
            start_button.clicked.connect(self._start_bridge_session)
            body_layout.addWidget(start_button, 0, Qt.AlignLeft)
            if self.start_session is None:
                hint = QLabel(
                    "This host did not provide a session starter."
                )
                hint.setWordWrap(True)
                body_layout.addWidget(hint)
        else:
            note = QLabel(
                "Copy this prompt to an external agent. It contains the session "
                "token and the exact JSONL files to use."
            )
            note.setWordWrap(True)
            body_layout.addWidget(note)
            prompt = build_agent_prompt(bridge)
            prompt_edit = QPlainTextEdit()
            prompt_edit.setReadOnly(True)
            prompt_edit.setPlainText(prompt)
            prompt_edit.setFont(QFont("Consolas", 9))
            prompt_edit.setMinimumHeight(220)
            body_layout.addWidget(prompt_edit)
            copy_row = QHBoxLayout()
            copy_button = QPushButton("Copy agent prompt")
            copy_button.clicked.connect(
                lambda: QApplication.clipboard().setText(prompt)
            )
            copy_row.addWidget(copy_button)
            copy_row.addStretch(1)
            body_layout.addLayout(copy_row)
            if not token:
                token_note = QLabel(
                    "This bridge was started without a token; enable token "
                    "authentication before exposing it to other scripts."
                )
                token_note.setWordWrap(True)
                body_layout.addWidget(token_note)

        if self.bridge_error:
            error = QLabel("Could not start the agent session: %s"
                           % self.bridge_error)
            error.setWordWrap(True)
            error.setStyleSheet("color: palette(text); font-weight: 600;")
            body_layout.addWidget(error)

        self.bridge_layout.addWidget(body)
        self.bridge_group.updateGeometry()

    def _scroll_to_bridge(self):
        try:
            self.scroll.ensureWidgetVisible(self.bridge_group)
        except RuntimeError:
            pass

    def _start_bridge_session(self):
        if self.start_session is None:
            return
        try:
            bridge = self.start_session()
        except Exception as exc:
            self.bridge_error = str(exc)
        else:
            self.bridge_error = ""
            if bridge is not None:
                self.bridge = bridge
        self._refresh_bridge_group()
        QTimer.singleShot(0, self._scroll_to_bridge)
