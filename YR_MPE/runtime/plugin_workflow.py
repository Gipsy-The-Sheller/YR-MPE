"""Workflow-step configuration mode for BasePlugin windows.

When a plugin is opened from a workflow preview, its path inputs are locked to
placeholder tokens and the Run button becomes ``Close & Apply``. The plugin
does not execute the tool in this mode; it emits a configuration dictionary
that the workflow executor can apply to the step.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List


def _placeholders_in(value: Any) -> Iterable[str]:
    import re
    pattern = re.compile(r"%[A-Za-z0-9_.-]+%")
    if isinstance(value, dict):
        for item in value.values():
            yield from _placeholders_in(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _placeholders_in(item)
    elif isinstance(value, str):
        for match in pattern.finditer(value):
            yield match.group(0)


def _first_placeholder(step: Dict[str, Any]) -> str:
    for value in (step.get("inputs") or {}).values():
        for placeholder in _placeholders_in(value):
            return placeholder
    for value in (step.get("params") or {}).values():
        for placeholder in _placeholders_in(value):
            return placeholder
    return ""


def _placeholder_for_binding(binding: Dict[str, Any], step: Dict[str, Any],
                             fallback: str) -> str:
    inputs = step.get("inputs") or {}
    for key in (binding.get("id"), binding.get("widget")):
        if key and key in inputs:
            values = list(_placeholders_in(inputs[key]))
            if values:
                return values[0]
    return fallback


def bind_workflow_step(plugin, step: Dict[str, Any]) -> Dict[str, Any]:
    """Lock path controls to placeholders and switch to Close & Apply."""
    context = dict(step or {})
    plugin.workflow_mode = True
    plugin.workflow_context = context
    fallback = _first_placeholder(context)
    controls = plugin.help_doc().get("controls", [])
    locked = []
    for row in controls:
        if row.get("kind") != "path":
            continue
        widget = getattr(plugin, row.get("widget", ""), None)
        if widget is None:
            continue
        placeholder = _placeholder_for_binding(row, context, fallback)
        if not placeholder:
            continue
        if callable(getattr(widget, "setText", None)):
            widget.setText(placeholder)
            if callable(getattr(widget, "setReadOnly", None)):
                widget.setReadOnly(True)
        if callable(getattr(widget, "setEnabled", None)):
            widget.setEnabled(False)
        locked.append({
            "id": row.get("id", ""),
            "widget": row.get("widget", ""),
            "placeholder": placeholder,
        })
    if not locked and fallback and getattr(plugin, "file_path_edit", None) is not None:
        widget = plugin.file_path_edit
        if callable(getattr(widget, "setText", None)):
            widget.setText(fallback)
        if callable(getattr(widget, "setReadOnly", None)):
            widget.setReadOnly(True)
        if callable(getattr(widget, "setEnabled", None)):
            widget.setEnabled(False)
        locked.append({
            "id": "input",
            "widget": "file_path_edit",
            "placeholder": fallback,
        })
    if fallback:
        plugin.imported_files = []
    plugin.workflow_locked_inputs = locked
    params = context.get("params") or {}
    apply_control = getattr(plugin, "_apply_agent_control", None)
    if callable(apply_control) and isinstance(params, dict):
        for row in controls:
            if row.get("kind") == "path":
                continue
            control_id = row.get("id", "")
            if control_id and control_id in params:
                try:
                    apply_control(control_id, params[control_id])
                except Exception:
                    pass
    for attribute in (
        "file_browse_btn", "sequence_text", "save_to_cwd", "save_to_tmp",
        "import_to_platform_btn",
    ):
        widget = getattr(plugin, attribute, None)
        if widget is not None and callable(getattr(widget, "setEnabled", None)):
            widget.setEnabled(False)
    if getattr(plugin, "stop_button", None) is not None:
        plugin.stop_button.setVisible(False)
    if getattr(plugin, "progress_bar", None) is not None:
        plugin.progress_bar.setVisible(False)
    if getattr(plugin, "run_button", None) is not None:
        plugin.run_button.setText("Close && Apply")
        plugin.run_button.setToolTip(
            "Apply these step parameters and close the configuration window."
        )
        try:
            plugin.run_button.clicked.disconnect()
        except (TypeError, RuntimeError):
            pass
        plugin.run_button.clicked.connect(plugin._apply_workflow_step)
    plugin.add_console_message(
        "Workflow step mode: inputs are locked to placeholders; Run is "
        "replaced by Close & Apply.",
        "info",
    )
    return context


def workflow_step_config(plugin) -> Dict[str, Any]:
    context = getattr(plugin, "workflow_context", {}) or {}
    parameters = {}
    if callable(getattr(plugin, "get_parameters", None)):
        try:
            parameters = plugin.get_parameters()
        except Exception:
            parameters = {}
    help_doc = {}
    if callable(getattr(plugin, "help_doc", None)):
        try:
            help_doc = plugin.help_doc()
        except Exception:
            help_doc = {}
    return {
        "id": context.get("step_id", context.get("id", "")),
        "op": context.get("op", ""),
        "inputs": context.get("inputs", {}),
        "outputs": context.get("outputs", {}),
        "parameters": parameters,
        "controls": help_doc.get("controls", []),
        "cli": help_doc.get("cli", {}).get("preview", []),
        "locked_inputs": list(getattr(plugin, "workflow_locked_inputs", [])),
    }


def apply_workflow_step(plugin) -> None:
    """Emit the step configuration and close the parent dialog."""
    config = workflow_step_config(plugin)
    plugin.workflow_configured.emit(config)
    window = plugin.window()
    if window is not None and window is not plugin:
        window.close()
