"""Semantic host actions for an agent session started from the Workbench.

The file bridge never executes shell commands; every action calls a Python
method on the workbench or on a plugin window. The host exposes four discovery
surfaces before execution:

``workbench.menu``      current toolbar action tree with stable ids
``workbench.workspace`` datasets, items, selections and run history
``workbench.sidebar``   basket contents and supported send-to operations
``workbench.jobs``      queued/running/finished/error task records
``workbench.skills``    agent skill catalog; ``workbench.help_doc`` fetches one
"""

from __future__ import annotations

import os
import platform
import sys

from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication

from .. import __version__
from .menu_inspector import find_action, menu_tree
from .skills import SkillLibrary
from .ui_actions import UIActionRegistry


def _plugin_windows(host):
    """Return ``(window_id, window)`` pairs for agent-capable widgets."""
    windows = []
    seen = set()
    for widget in QApplication.allWidgets():
        if widget is host or id(widget) in seen:
            continue
        if not callable(getattr(widget, "agent_invoke", None)):
            continue
        seen.add(id(widget))
        window_id = getattr(widget, "_yrmpe_agent_window_id", None)
        if not window_id:
            counter = getattr(host, "_agent_window_counter", 0) + 1
            host._agent_window_counter = counter
            window_id = "plugin-%d" % counter
            widget._yrmpe_agent_window_id = window_id
        windows.append((window_id, widget))
    return windows


def _find_window(host, window_id):
    for candidate_id, widget in _plugin_windows(host):
        if candidate_id == window_id:
            return widget
    return None


def _help_widgets(host):
    """Widgets with a help_doc method, whether or not they expose actions."""
    rows = []
    seen = set()
    for widget in QApplication.allWidgets():
        if widget is host or id(widget) in seen:
            continue
        if not callable(getattr(widget, "help_doc", None)):
            continue
        seen.add(id(widget))
        rows.append(widget)
    return rows


def _find_help_widget(host, *, window_id="", plugin="", class_name=""):
    if window_id:
        widget = _find_window(host, window_id)
        if widget is not None and callable(getattr(widget, "help_doc", None)):
            return widget
    plugin = str(plugin or "").strip().lower()
    class_name = str(class_name or "").strip().lower()
    for widget in _help_widgets(host):
        names = {
            str(getattr(widget, "plugin_name", "")).lower(),
            str(getattr(widget, "tool_name", "")).lower(),
            widget.__class__.__name__.lower(),
        }
        if plugin and plugin in names:
            return widget
        if class_name and class_name == widget.__class__.__name__.lower():
            return widget
    return None


def _item_summary(item):
    data = item.data if isinstance(getattr(item, "data", None), dict) else {}
    data_keys = list(data.keys()) if isinstance(data, dict) else []
    best_model = (
        data.get("best_model")
        or data.get("model")
        or data.get("Model")
        or ""
    )
    if not best_model and data.get("type") == "model_table":
        rows = data.get("data") or []
        if rows:
            best_model = rows[0].get("Model", "")
    tree_count = 0
    if item.item_type == "treeset":
        tree_count = len(data.get("data") or [])
    return {
        "id": item.id,
        "dataset_id": getattr(item, "dataset_id", ""),
        "name": getattr(item, "loci_name", "") or getattr(item, "name", ""),
        "type": item.item_type,
        "format": getattr(item, "format", ""),
        "file": getattr(item, "file_path", "") or getattr(item, "file", ""),
        "length": getattr(item, "length", 0),
        "sequence_count": getattr(item, "sequence_count", 0),
        "is_aligned": bool(getattr(item, "is_aligned", False)),
        "model_sub_type": getattr(item, "model_sub_type", ""),
        "parent_id": getattr(item, "parent_id", ""),
        "import_id": getattr(item, "import_id", ""),
        "run_id": getattr(item, "run_id", ""),
        "selection_state": getattr(item, "selection_state", ""),
        "best_model": best_model,
        "tree_count": tree_count,
        "data_keys": data_keys[:24],
    }


def _dataset_rows(host):
    manager = getattr(host, "dataset_selection_manager", None)
    if manager is None:
        return []
    rows = []
    for dataset in manager.get_all_datasets():
        items = manager.get_items_by_dataset(dataset.id)
        settings = getattr(dataset, "settings", {}) or {}
        rows.append({
            "id": dataset.id,
            "name": dataset.name,
            "description": getattr(dataset, "description", ""),
            "selection_state": getattr(dataset, "selection_state", ""),
            "item_count": len(items),
            "settings": {
                key: settings.get(key)
                for key in ("topo_linked", "edge_linked")
                if key in settings
            },
            "items": [_item_summary(item) for item in items],
        })
    return rows


def _workspace_payload(host):
    manager = getattr(host, "dataset_selection_manager", None)
    datasets = _dataset_rows(host)
    items = []
    if manager is not None:
        items = [_item_summary(item) for item in manager.items.values()]
    green = [item for item in items
             if item["selection_state"] == "green"]
    basket_ids = []
    basket = []
    if manager is not None:
        basket_ids = list(getattr(manager, "basket_items", set()))
        basket = [_item_summary(item) for item in manager.get_basket_items()]
    by_type = {}
    for item in items:
        by_type.setdefault(item["type"], []).append(item["id"])

    runs = []
    streams = []
    imports = []
    store = getattr(manager, "store", None) if manager is not None else None
    if store is not None:
        for run in store.index.runs.values():
            runs.append({
                "id": run.id,
                "name": run.name,
                "method": run.method,
                "tool": run.tool,
                "status": run.status,
                "input_ids": list(run.input_ids),
                "output_ids": list(run.output_ids),
                "artifacts": list(run.artifacts),
                "stream_id": run.stream_id,
                "created_at": run.created_at,
                "completed_at": run.completed_at,
            })
        for stream in store.index.streams.values():
            streams.append({
                "id": stream.id,
                "name": stream.name,
                "dataset_id": stream.dataset_id,
                "run_ids": list(stream.run_ids),
            })
        for record in store.index.imports.values():
            imports.append({
                "id": record.id,
                "importer": record.importer,
                "source_paths": list(record.source_paths),
                "dataset_id": record.dataset_id,
                "item_ids": list(record.item_ids),
                "created_at": record.created_at,
            })

    return {
        "workspace_path": getattr(
            getattr(host, "workspace_manager", None), "workspace_path", ""
        ) or "",
        "current_dataset_id": getattr(host, "current_dataset_id", ""),
        "dataset_count": len(datasets),
        "item_count": len(items),
        "datasets": datasets,
        "items": items,
        "items_by_type": by_type,
        "selected_item_ids": [item["id"] for item in green],
        "selected_items": green,
        "basket_item_ids": basket_ids,
        "basket_items": basket,
        "runs": runs,
        "streams": streams,
        "imports": imports,
    }


def _sidebar_payload(host):
    manager = getattr(host, "dataset_selection_manager", None)
    panel = getattr(host, "selected_panel", None)
    basket = []
    green = []
    if manager is not None:
        basket = [_item_summary(item) for item in manager.get_basket_items()]
        green = [_item_summary(item) for item in manager.get_selected_items()]
    recommendations = []
    if panel is not None and callable(getattr(panel, "refresh", None)):
        try:
            panel.refresh()
        except Exception:
            pass
    if panel is not None and callable(getattr(panel, "recommended_operations", None)):
        try:
            recommendations = panel.recommended_operations()
        except Exception:
            recommendations = []
    inspector = {}
    inspector_panel = getattr(host, "inspector", None)
    if inspector_panel is not None and callable(getattr(inspector_panel, "state", None)):
        try:
            inspector = inspector_panel.state()
        except Exception:
            inspector = {}
    queue = getattr(host, "task_queue", None)
    jobs = {}
    if queue is not None and callable(getattr(queue, "state_payload", None)):
        try:
            jobs = queue.state_payload()
        except Exception:
            jobs = {}
    return {
        "selected_item_ids": [item["id"] for item in green],
        "selected_items": green,
        "basket_item_ids": [item["id"] for item in basket],
        "basket_items": basket,
        "recommended_operations": recommendations,
        "inspector": inspector,
        "jobs": jobs.get("jobs", []),
        "queue": {
            "running": jobs.get("running", 0),
            "queued": jobs.get("queued", 0),
            "max_concurrent": jobs.get("max_concurrent", 0),
        },
    }


def build_workbench_registry(host) -> UIActionRegistry:
    """Create the default action registry for a session started in the GUI."""
    registry = UIActionRegistry()

    def workbench_info(params):
        bridge = getattr(host, "agent_bridge", None)
        return {
            "ok": True,
            "result": {
                "application": "YPhylo Workbench",
                "plugin_set": "YR-MPE %s" % __version__,
                "workspace": getattr(
                    getattr(host, "workspace_manager", None),
                    "workspace_path", "",
                ) or "",
                "cwd": os.getcwd(),
                "python": sys.version.split()[0],
                "platform": platform.platform(),
                "plugin_windows": len(_plugin_windows(host)),
                "agent_session": str(getattr(bridge, "session_dir", "") or ""),
                "skills": len(SkillLibrary().list_skills()),
            },
        }

    registry.register(
        "workbench.info",
        workbench_info,
        description="Show application, workspace and runtime information",
    )

    def list_windows(params):
        rows = []
        for window_id, widget in _plugin_windows(host):
            row = {
                "id": window_id,
                "class": widget.__class__.__name__,
                "title": widget.windowTitle(),
                "visible": bool(widget.isVisible()),
            }
            if callable(getattr(widget, "agent_contract", None)):
                try:
                    contract = widget.agent_contract()
                    row["contract"] = {
                        key: value for key, value in contract.items()
                        if key != "options"
                    }
                except Exception as exc:
                    row["contract_error"] = str(exc)
            rows.append(row)
        return {"ok": True, "result": {"window_count": len(rows), "windows": rows}}

    registry.register(
        "workbench.list_windows",
        list_windows,
        description="List open plugin windows that expose agent actions",
    )

    def invoke_plugin(params):
        window_id = params.get("window_id")
        action = params.get("action") or params.get("action_id")
        if not window_id:
            return {"ok": False, "error": "missing window_id"}
        if not action:
            return {"ok": False, "error": "missing action"}
        widget = _find_window(host, window_id)
        if widget is None:
            return {"ok": False, "error": "unknown window_id %r" % window_id}
        outcome = widget.agent_invoke(action, params.get("params") or {})
        if isinstance(outcome, dict) and "ok" in outcome:
            return outcome
        return {"ok": True, "result": outcome}

    registry.register(
        "workbench.invoke_plugin",
        invoke_plugin,
        description=(
            "Call a plugin action by window_id and action id; use "
            "workbench.list_windows first"
        ),
        parameters={
            "window_id": {"type": "string", "required": True},
            "action": {"type": "string", "required": True},
            "params": {"type": "object", "required": False},
        },
    )

    def list_datasets(params):
        rows = _dataset_rows(host)
        return {
            "ok": True,
            "result": {
                "dataset_count": len(rows),
                "item_count": sum(row["item_count"] for row in rows),
                "datasets": rows,
            },
        }

    registry.register(
        "workbench.list_datasets",
        list_datasets,
        description="List workspace datasets, items and sequence counts",
    )

    def menu_state(params):
        payload = menu_tree(host)
        query = str(params.get("query", "")).strip().lower()
        if query:
            payload["actions"] = [
                row for row in payload["actions"]
                if query in row["id"].lower()
                or query in row["label"].lower()
                or query in str(row.get("description", "")).lower()
            ]
        return {"ok": True, "result": payload}

    registry.register(
        "workbench.menu",
        menu_state,
        description=(
            "Return the toolbar menu tree and a flat action list with stable "
            "ids and short descriptions; optional query filters the action list"
        ),
    )

    def invoke_menu(params):
        action_id = params.get("id") or params.get("action_id") or ""
        path = params.get("path")
        if action_id == "align_open_sequence_files":
            return {
                "ok": False,
                "error": (
                    "this action opens a human FileDialog; use "
                    "workbench.import_files with explicit paths instead"
                ),
            }
        if action_id == "align_align_by_mafft":
            opener = getattr(host, "open_mafft_wrapper", None)
            if callable(opener):
                QTimer.singleShot(0, lambda: opener(modal=False))
                return {
                    "ok": True,
                    "result": {
                        "id": action_id,
                        "label": "MAFFT",
                        "scheduled": True,
                        "non_modal": True,
                        "hint": "call workbench.list_windows for the new window_id",
                    },
                }
        action, node = find_action(host, action_id, path)
        if action is None:
            return {"ok": False, "error": "unknown menu action %r" % action_id}
        if not action.isEnabled():
            return {
                "ok": False,
                "error": "menu action is disabled",
                "action": node,
            }
        if action.menu() is not None:
            return {
                "ok": False,
                "error": "menu entry is a submenu; choose a child action",
                "action": node,
            }
        blocking = bool(params.get("blocking", False))
        if blocking:
            action.trigger()
        else:
            QTimer.singleShot(0, action.trigger)
        return {
            "ok": True,
            "result": {
                "id": node.get("id", action_id),
                "label": node.get("label", action.text()),
                "path": node.get("path", path or []),
                "scheduled": not blocking,
            },
        }

    registry.register(
        "workbench.invoke_menu",
        invoke_menu,
        description=(
            "Trigger a toolbar action by id or path. By default the trigger "
            "is scheduled after the response so modal dialogs do not block "
            "the bridge"
        ),
        parameters={
            "id": {"type": "string", "required": False},
            "path": {"type": "array", "required": False},
            "blocking": {"type": "boolean", "required": False},
        },
    )

    def workspace_state(params):
        return {"ok": True, "result": _workspace_payload(host)}

    registry.register(
        "workbench.workspace",
        workspace_state,
        description=(
            "Return datasets, items, selection states, runs, streams and "
            "imports in the current workspace"
        ),
    )

    def workflow_graph(params):
        from .workflow_graph import graph_from_plan, graph_from_workspace
        plan = params.get("plan")
        if plan:
            graph = graph_from_plan(plan)
        else:
            graph = graph_from_workspace(
                getattr(host, "dataset_selection_manager", None)
            )
        payload = graph.to_dict()
        payload["node_count"] = len(graph.nodes)
        payload["edge_count"] = len(graph.edges)
        return {"ok": True, "result": payload}

    registry.register(
        "workbench.workflow_graph",
        workflow_graph,
        description=(
            "Return the current workspace workflow as data/tool nodes and "
            "edges for previews or agent reasoning"
        ),
    )

    def show_workflow_preview(params):
        opener = getattr(host, "show_workflow_preview", None)
        if not callable(opener):
            return {"ok": False, "error": "workflow preview is unavailable"}
        QTimer.singleShot(0, opener)
        return {"ok": True, "result": {"scheduled": True}}

    registry.register(
        "workbench.show_workflow_preview",
        show_workflow_preview,
        description="Open the graphical workflow preview dialog",
    )

    def preview_workflow(params):
        opener = getattr(host, "show_workflow_plan", None)
        if not callable(opener):
            return {"ok": False, "error": "workflow plan preview is unavailable"}
        plan = params.get("plan")
        if not isinstance(plan, dict):
            return {"ok": False, "error": "params.plan must be an object"}
        QTimer.singleShot(0, lambda: opener(plan))
        return {"ok": True, "result": {"scheduled": True}}

    registry.register(
        "workbench.preview_workflow",
        preview_workflow,
        description=(
            "Open a compiled yrmpe.workflow.v1 plan in the graphical preview; "
            "tool nodes can be clicked to configure parameters"
        ),
        parameters={"plan": {"type": "object", "required": True}},
    )

    def sidebar_state(params):
        return {"ok": True, "result": _sidebar_payload(host)}

    registry.register(
        "workbench.sidebar",
        sidebar_state,
        description=(
            "Return Selected basket contents, recommended send-to operations "
            "and the current Inspector item/operations"
        ),
    )

    def skills_list(params):
        library = SkillLibrary()
        query = str(params.get("query", "")).strip()
        if query:
            entries = library.search(query)
        else:
            entries = library.list_skills()
        return {"ok": True, "result": {"skills": entries, "count": len(entries)}}

    registry.register(
        "workbench.skills",
        skills_list,
        description=(
            "List agent skills (one per plugin/tool) with summaries and menu "
            "paths; optional query searches ids, titles and summaries"
        ),
    )

    def help_doc(params):
        library = SkillLibrary()
        window_id = str(params.get("window_id", "")).strip()
        plugin = str(params.get("plugin", "")).strip()
        class_name = str(params.get("class", "")).strip()
        skill_id = str(params.get("id") or params.get("skill") or "").strip()

        widget = _find_help_widget(
            host, window_id=window_id, plugin=plugin, class_name=class_name
        )
        if widget is not None:
            doc = widget.help_doc()
        elif skill_id:
            doc = library.get_skill(skill_id)
            if doc is None:
                return {"ok": False, "error": "unknown skill %r" % skill_id}
        elif plugin:
            entry = library.find_for_plugin(plugin_name=plugin)
            if entry is None:
                return {"ok": False, "error": "unknown plugin %r" % plugin}
            doc = entry
        else:
            return {
                "ok": False,
                "error": "provide id, plugin, or window_id; use workbench.skills for ids",
            }
        if not callable(getattr(widget, "help_doc", None)) and not doc:
            return {"ok": False, "error": "help document unavailable"}
        return {"ok": True, "result": doc}

    registry.register(
        "workbench.help_doc",
        help_doc,
        description=(
            "Return a plugin skill document with Markdown and live GUI to CLI "
            "control bindings; use workbench.skills for ids"
        ),
        parameters={
            "id": {"type": "string", "required": False},
            "plugin": {"type": "string", "required": False},
            "window_id": {"type": "string", "required": False},
            "class": {"type": "string", "required": False},
        },
    )

    def import_files(params):
        importer = getattr(host, "import_sequence_files", None)
        if not callable(importer):
            return {"ok": False, "error": "file import API is unavailable"}
        paths = params.get("paths", params.get("path", []))
        if isinstance(paths, str):
            paths = [paths]
        try:
            imported = importer(
                paths, file_format=str(params.get("format", "") or "")
            )
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        return {
            "ok": True,
            "result": {
                "imported": imported,
                "file_count": len(imported),
                "sequence_count": sum(item["sequences"] for item in imported),
            },
        }

    registry.register(
        "workbench.import_files",
        import_files,
        description=(
            "Import local sequence files by path without opening a file "
            "dialog; supported formats: FASTA, PHYLIP, NEXUS"
        ),
        parameters={
            "paths": {"type": "array", "required": True},
            "format": {"type": "string", "required": False},
        },
    )

    def open_tool(params):
        tool = str(params.get("tool") or params.get("id") or "").lower()
        if tool in ("mafft", "align_mafft", "align_align_by_mafft"):
            opener = getattr(host, "open_mafft_wrapper", None)
            if callable(opener):
                QTimer.singleShot(0, lambda: opener(modal=False))
                return {
                    "ok": True,
                    "result": {
                        "tool": "mafft",
                        "scheduled": True,
                        "hint": "call workbench.list_windows for the new window_id",
                    },
                }
        if tool in ("best", "best_mrbayes", "mrbayes_best", "msc"):
            opener = getattr(host, "open_best_wrapper", None)
            if callable(opener):
                QTimer.singleShot(0, opener)
                return {
                    "ok": True,
                    "result": {
                        "tool": "best",
                        "scheduled": True,
                        "hint": "call workbench.list_windows, then use "
                                "plugin.inspect_controls / best.apply_block",
                    },
                }
        return {
            "ok": False,
            "error": "unsupported non-modal tool %r" % tool,
        }

    registry.register(
        "workbench.open_tool",
        open_tool,
        description=(
            "Open a tool in a non-modal window for agent use; supports "
            "MAFFT and BEST/MrBayes MSC. Use workbench.list_windows "
            "afterwards"
        ),
        parameters={"tool": {"type": "string", "required": True}},
    )

    def open_orthofinder(params):
        opener = getattr(host, "create_dataset_by_orthofinder", None)
        if not callable(opener):
            return {"ok": False, "error": "OrthoFinder is not available"}
        opener()
        return {
            "ok": True,
            "result": {
                "opened": True,
                "hint": "call workbench.list_windows to get its window_id",
            },
        }

    registry.register(
        "workbench.open_orthofinder",
        open_orthofinder,
        description="Open the One-Step OrthoFinder plugin window",
    )

    def jobs_state(params):
        queue = getattr(host, "task_queue", None)
        if queue is None or not callable(getattr(queue, "state_payload", None)):
            return {"ok": False, "error": "task queue is unavailable"}
        return {"ok": True, "result": queue.state_payload()}

    registry.register(
        "workbench.jobs",
        jobs_state,
        description=(
            "List queued/running/finished/error tasks with job_id, status, "
            "progress and artifacts"
        ),
    )

    def job_status(params):
        queue = getattr(host, "task_queue", None)
        job_id = str(params.get("job_id", "") or "")
        if queue is None or not job_id:
            return {"ok": False, "error": "provide job_id"}
        getter = getattr(queue, "job_payload", None)
        payload = getter(job_id) if callable(getter) else None
        if payload is None:
            return {"ok": False, "error": "unknown job %r" % job_id}
        return {"ok": True, "result": payload}

    registry.register(
        "workbench.job_status",
        job_status,
        description="Return one task's state, progress and artifacts",
        parameters={"job_id": {"type": "string", "required": True}},
    )

    def cancel_job(params):
        queue = getattr(host, "task_queue", None)
        job_id = str(params.get("job_id", "") or "")
        if queue is None or not job_id:
            return {"ok": False, "error": "provide job_id"}
        cancel = getattr(queue, "cancel", None)
        if not callable(cancel):
            return {"ok": False, "error": "task cancellation is unavailable"}
        if not cancel(job_id):
            return {"ok": False, "error": "job is unknown or already finished"}
        return {"ok": True, "result": {"job_id": job_id, "cancelled": True}}

    registry.register(
        "workbench.cancel_job",
        cancel_job,
        description="Cancel a queued or running task by job_id",
        parameters={"job_id": {"type": "string", "required": True}},
    )

    def inspect_job(params):
        queue = getattr(host, "task_queue", None)
        job_id = str(params.get("job_id", "") or "")
        if queue is None or not job_id:
            return {"ok": False, "error": "provide job_id"}
        getter = getattr(queue, "job_payload", None)
        payload = getter(job_id) if callable(getter) else None
        if payload is None:
            return {"ok": False, "error": "unknown job %r" % job_id}
        return {"ok": True, "result": payload}

    registry.register(
        "workbench.inspect_job",
        inspect_job,
        description=(
            "Inspect one task's status, command, console summary and result "
            "artifacts; this is the data view behind the queue preview"
        ),
        parameters={"job_id": {"type": "string", "required": True}},
    )

    return registry
