"""Generic agent surface shared by all ``BasePlugin`` subclasses.

Plugin-specific APIs (for example ``orthofinder.*``) override these methods.
Generic actions use the ``plugin.*`` namespace and are addressed together with
the window id through ``workbench.invoke_plugin``.
"""

from __future__ import annotations

from typing import Any, Dict, List


def agent_actions():
    return [
        {"id": "plugin.help_doc",
         "description": "Return the skill document and live GUI to CLI map"},
        {"id": "plugin.inspect_status",
         "description": "Return run state and capabilities"},
        {"id": "plugin.inspect_controls",
         "description": "Return current control values and CLI preview"},
        {"id": "plugin.inspect_console",
         "description": "Read the plugin console tail"},
        {"id": "plugin.modify_options",
         "description": "Set controls from an options object",
         "parameters": {"options": {"type": "object", "required": True}}},
        {"id": "plugin.run",
         "description": "Start the plugin using current controls"},
        {"id": "plugin.stop",
         "description": "Stop the running plugin thread"},
        {"id": "plugin.export_to_workspace",
         "description": "Build and export the plugin ResultSpec"},
        {"id": "plugin.close",
         "description": "Close the plugin when no run is active"},
    ]


def agent_status(plugin) -> Dict[str, Any]:
    job_id = str(getattr(plugin, "queue_job_id", "") or "")
    job_state = ""
    queue = getattr(plugin, "task_queue", None)
    if job_id and queue is not None and callable(getattr(queue, "get_job", None)):
        job = queue.get_job(job_id)
        if job is not None:
            job_state = str(getattr(job, "state", "") or "")
    return {
        "plugin": getattr(plugin, "plugin_name", ""),
        "tool": getattr(plugin, "tool_name", ""),
        "tool_path": getattr(plugin, "tool_path", "") or "",
        "is_running": bool(getattr(plugin, "is_running", False)),
        "job_id": job_id,
        "job_state": job_state,
        "can_run": bool(getattr(plugin, "run_button", None)
                        and plugin.run_button.isEnabled()),
        "can_stop": bool(getattr(plugin, "stop_button", None)
                         and plugin.stop_button.isEnabled()),
        "can_export": callable(getattr(plugin, "_build_result_spec", None))
                      or callable(getattr(plugin, "import_to_platform", None)),
        "import_from": getattr(plugin, "import_from", ""),
        "outputs": list(getattr(plugin, "analysis_output_files", []) or []),
    }


def agent_contract(plugin) -> Dict[str, Any]:
    return {
        "schema": "yrmpe.plugin.v1",
        "id": getattr(plugin, "plugin_name", plugin.__class__.__name__),
        "status": agent_status(plugin),
        "actions": {"actions": agent_actions()},
    }


def _set_widget_value(widget, value) -> bool:
    try:
        if callable(getattr(widget, "setCurrentText", None)) and callable(
            getattr(widget, "currentText", None)
        ):
            widget.setCurrentText(str(value))
        elif callable(getattr(widget, "setValue", None)):
            widget.setValue(value)
        elif callable(getattr(widget, "setChecked", None)):
            widget.setChecked(bool(value))
        elif callable(getattr(widget, "setText", None)) and callable(
            getattr(widget, "text", None)
        ):
            widget.setText(str(value))
        elif callable(getattr(widget, "setPlainText", None)):
            widget.setPlainText(str(value))
        else:
            return False
        return True
    except Exception:
        return False


def apply_control(plugin, control_id: str, value: Any) -> bool:
    doc = plugin.help_doc()
    binding = next(
        (row for row in doc.get("controls", [])
         if row.get("id") == control_id or row.get("widget") == control_id),
        None,
    )
    if binding is None:
        widget = getattr(plugin, control_id, None)
        if widget is None:
            return False
        binding = {"widget": control_id, "choices": {}}
    widget = getattr(plugin, binding.get("widget", ""), None)
    if widget is None:
        return False
    choices = binding.get("choices") or {}
    if isinstance(choices, dict):
        for gui_value, cli_value in choices.items():
            if value == cli_value or str(value) == str(cli_value):
                value = gui_value
                break
    ok = _set_widget_value(widget, value)
    if ok and binding.get("kind") == "path":
        paths = value if isinstance(value, (list, tuple)) else [value]
        paths = [str(path) for path in paths if str(path).strip()]
        if paths:
            plugin.imported_files = paths
            if hasattr(plugin, "file_path_edit"):
                plugin.file_path_edit.setText(
                    paths[0] if len(paths) == 1
                    else "%d files selected" % len(paths)
                )
    return ok


def agent_invoke(plugin, action_id: str, params=None) -> Dict[str, Any]:
    params = params or {}
    name = str(action_id or "").rsplit(".", 1)[-1]
    if name == "help_doc":
        return {"ok": True, "result": plugin.help_doc()}
    if name == "inspect_status":
        return {"ok": True, "result": agent_status(plugin)}
    if name == "inspect_actions":
        return {"ok": True, "result": {
            "schema": "yrmpe.plugin.actions.v1",
            "actions": agent_actions(),
            "status": agent_status(plugin),
        }}
    if name == "inspect_controls":
        doc = plugin.help_doc()
        return {"ok": True, "result": {
            "controls": doc.get("controls", []),
            "cli": doc.get("cli", {}),
            "warnings": doc.get("warnings", []),
        }}
    if name == "inspect_console":
        limit = int(params.get("limit", 200) or 200)
        text = ""
        if hasattr(plugin, "console_text") and hasattr(
            plugin.console_text, "toPlainText"
        ):
            text = plugin.console_text.toPlainText()
        lines = text.splitlines()
        return {"ok": True, "result": {
            "lines": lines[-limit:],
            "total_lines": len(lines),
        }}
    if name == "modify_options":
        options = params.get("options", params)
        if not isinstance(options, dict):
            return {"ok": False, "error": "options must be an object"}
        applied, ignored = [], []
        for key, value in options.items():
            key = str(key)
            if key in ("files", "imported_files") and isinstance(
                value, (list, tuple)
            ):
                paths = [str(path) for path in value if str(path).strip()]
                plugin.imported_files = paths
                if hasattr(plugin, "file_path_edit") and paths:
                    plugin.file_path_edit.setText(
                        paths[0] if len(paths) == 1
                        else "%d files selected" % len(paths)
                    )
                applied.append(key)
            elif apply_control(plugin, key, value):
                applied.append(key)
            else:
                ignored.append(key)
        return {"ok": True, "result": {
            "applied": applied,
            "ignored": ignored,
            "status": agent_status(plugin),
        }}
    if name == "run":
        if getattr(plugin, "is_running", False):
            return {"ok": False, "error": "plugin is already running",
                    "status": agent_status(plugin)}
        if not callable(getattr(plugin, "run_analysis", None)):
            return {"ok": False, "error": "plugin has no run_analysis"}
        starter = getattr(plugin, "_request_run", None)
        job_id = starter() if callable(starter) else plugin.run_analysis()
        result = {"ok": True, "status": agent_status(plugin)}
        if job_id:
            result["job_id"] = str(job_id)
        return result
    if name == "stop":
        queue = getattr(plugin, "task_queue", None)
        job_id = str(getattr(plugin, "queue_job_id", "") or "")
        job = None
        if queue is not None and job_id and callable(getattr(queue, "get_job", None)):
            job = queue.get_job(job_id)
        if job is not None and getattr(job, "state", "") in ("queued", "running"):
            queue.cancel(job_id)
            return {"ok": True, "job_id": job_id,
                    "status": agent_status(plugin)}
        if callable(getattr(plugin, "stop_analysis", None)):
            plugin.stop_analysis()
        return {"ok": True, "status": agent_status(plugin)}
    if name == "export_to_workspace":
        if not callable(getattr(plugin, "_build_result_spec", None)):
            return {"ok": False, "error": "plugin has no ResultSpec builder"}
        if not plugin.result_has_consumer():
            return {"ok": False, "error": "no workspace result channel"}
        try:
            result = plugin._build_result_spec()
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        if not plugin.export_result(result):
            return {"ok": False, "error": "ResultSpec validation failed"}
        return {"ok": True, "result": result.to_settings_dict()}
    if name == "preview_results":
        if not callable(getattr(plugin, "_build_result_spec", None)):
            return {"ok": False, "error": "plugin has no ResultSpec builder"}
        try:
            result = plugin._build_result_spec()
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        return {"ok": True, "result": result.to_settings_dict()}
    if name == "close":
        if getattr(plugin, "is_running", False):
            return {"ok": False, "error": "plugin is running"}
        queue = getattr(plugin, "task_queue", None)
        job_id = str(getattr(plugin, "queue_job_id", "") or "")
        job = None
        if queue is not None and job_id and callable(getattr(queue, "get_job", None)):
            job = queue.get_job(job_id)
        if job is not None and getattr(job, "state", "") in ("queued", "running"):
            return {"ok": False, "error": "plugin has an active queue job",
                    "job_id": job_id}
        plugin.close()
        return {"ok": True, "result": {"closed": True}}
    return {"ok": False, "error": "unknown plugin action %r" % action_id}
