"""Agent-facing skill documents and GUI-to-CLI control maps.

A plugin inherits ``BasePlugin.help_doc()``. The default implementation reads
its catalog entry from ``YR_MPE/skills/manifest.json``, optionally loads a
Markdown companion, and adds live control bindings declared by the plugin's
``control_bindings()`` method.

The important contract is that every CLI-capable control must declare the
exact option it maps to. A control may also be GUI-only (``cli: []``) when it
only changes post-processing or presentation.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

SKILL_SCHEMA = "yrmpe.skill.v1"
MANIFEST_SCHEMA = "yrmpe.skills.v1"
SKILLS_DIR = Path(__file__).resolve().parents[1] / "skills"


def slug(text: str) -> str:
    value = re.sub(r"[^a-z0-9]+", "-", str(text or "").lower()).strip("-")
    return value or "skill"


class SkillLibrary:
    """Read-only catalog of skill metadata and Markdown files."""

    def __init__(self, root: Optional[str | Path] = None):
        self.root = Path(root) if root else SKILLS_DIR
        self.manifest_path = self.root / "manifest.json"
        self._manifest: Optional[Dict[str, Any]] = None

    def manifest(self) -> Dict[str, Any]:
        if self._manifest is None:
            if self.manifest_path.exists():
                with self.manifest_path.open("r", encoding="utf-8") as handle:
                    payload = json.load(handle)
                if payload.get("schema") != MANIFEST_SCHEMA:
                    raise ValueError(
                        "unsupported skills manifest schema %r"
                        % payload.get("schema")
                    )
                self._manifest = payload
            else:
                self._manifest = {"schema": MANIFEST_SCHEMA, "skills": []}
        return self._manifest

    def entries(self) -> List[Dict[str, Any]]:
        return list(self.manifest().get("skills", []))

    def list_skills(self) -> List[Dict[str, Any]]:
        rows = []
        for entry in self.entries():
            rows.append({
                "id": entry.get("id", ""),
                "title": entry.get("title", entry.get("id", "")),
                "summary": entry.get("summary", ""),
                "menu_path": entry.get("menu_path", []),
                "plugin_names": entry.get("plugin_names", []),
                "tool_names": entry.get("tool_names", []),
                "actions": entry.get("actions", []),
                "has_document": bool(entry.get("file")),
            })
        return rows

    def get_entry(self, skill_id: str) -> Optional[Dict[str, Any]]:
        skill_id = str(skill_id or "").strip()
        for entry in self.entries():
            if entry.get("id") == skill_id:
                return dict(entry)
        return None

    def get_skill(self, skill_id: str) -> Optional[Dict[str, Any]]:
        entry = self.get_entry(skill_id)
        return self.hydrate(entry) if entry is not None else None

    def hydrate(self, entry: Dict[str, Any]) -> Dict[str, Any]:
        result = dict(entry)
        result["schema"] = SKILL_SCHEMA
        file_name = entry.get("file", "")
        result["markdown"] = ""
        if file_name:
            path = self.root / file_name
            if path.exists():
                result["markdown"] = path.read_text(encoding="utf-8")
                result["document_path"] = str(path)
        if not result["markdown"]:
            lines = ["# %s" % result.get("title", result.get("id", "Skill")), ""]
            if result.get("summary"):
                lines.extend([result["summary"], ""])
            if result.get("menu_path"):
                lines.extend([
                    "Menu: %s" % " -> ".join(result["menu_path"]), ""
                ])
            lines.extend([
                "This catalog entry has no long-form document yet. Open the "
                "plugin and call `workbench.help_doc` with its window id to "
                "get live controls and CLI mapping.",
                "",
            ])
            result["markdown"] = "\n".join(lines)
        return result

    def find_for_plugin(
        self,
        plugin_name: str = "",
        tool_name: str = "",
        class_name: str = "",
    ) -> Optional[Dict[str, Any]]:
        candidates = {
            "plugin": str(plugin_name or "").strip().lower(),
            "tool": str(tool_name or "").strip().lower(),
            "class": str(class_name or "").strip().lower(),
        }
        # Exact plugin/class matches win over tool aliases. Several plugins
        # share tool_name "IQ-TREE 3"; checking tool names in the same pass
        # would make ModelFinder capture IQ-TREE/Model Parameter/RF requests.
        for entry in self.entries():
            plugin_names = {
                str(value).lower() for value in entry.get("plugin_names", [])
            }
            class_names = {
                str(value).lower() for value in entry.get("class_names", [])
            }
            if candidates["plugin"] and candidates["plugin"] in plugin_names:
                return self.hydrate(entry)
            if candidates["class"] and candidates["class"] in class_names:
                return self.hydrate(entry)

        for entry in self.entries():
            tool_names = {
                str(value).lower() for value in entry.get("tool_names", [])
            }
            if candidates["tool"] and candidates["tool"] in tool_names:
                return self.hydrate(entry)
        return None

    def search(self, query: str) -> List[Dict[str, Any]]:
        needle = str(query or "").strip().lower()
        if not needle:
            return self.list_skills()
        rows = []
        for entry in self.entries():
            haystack = " ".join([
                entry.get("id", ""),
                entry.get("title", ""),
                entry.get("summary", ""),
                " ".join(entry.get("alias", [])),
            ]).lower()
            if needle in haystack:
                rows.append(entry)
        return rows


def _widget_value(widget) -> Any:
    if widget is None:
        return None
    # QComboBox exposes currentText; QSpinBox must use value() rather than
    # text(), because its text may be the special value ("AUTO").
    for method in ("currentText", "value", "isChecked", "text", "toPlainText"):
        if callable(getattr(widget, method, None)):
            try:
                return getattr(widget, method)()
            except Exception:
                pass
    return None


def _resolved_value(binding: Dict[str, Any], current: Any) -> Any:
    choices = binding.get("choices") or {}
    if isinstance(choices, dict) and current in choices:
        return choices[current]
    if isinstance(choices, dict):
        for key, value in choices.items():
            if str(key).lower() == str(current).lower():
                return value
    return current


def _format_cli(binding: Dict[str, Any], value: Any) -> List[str]:
    cli = binding.get("cli") or []
    if not cli:
        return []
    if value in ("", None):
        return []
    if binding.get("kind") == "bool" and not value:
        return []
    if value in (binding.get("omit_values") or []):
        return []
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    result = []
    for token in cli:
        if token == "{value}":
            result.append(str(value))
        elif "{value}" in str(token):
            result.append(str(token).replace("{value}", str(value)))
        else:
            result.append(str(token))
    return result


def _resolve_widget(plugin, name: str):
    """Resolve a widget by attribute name or dotted path.

    Dotted paths let a plugin expose controls that live in a sub-panel or
    advanced dialog, for example ``advanced_dialog.round_spin`` or
    ``control_panel.methods_combo``.
    """
    name = str(name or "").strip()
    if not name:
        return None
    obj = plugin
    for part in name.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return None
    return obj


def _binding_rows(plugin) -> List[Dict[str, Any]]:
    rows = []
    bindings = []
    if callable(getattr(plugin, "control_bindings", None)):
        try:
            bindings = list(plugin.control_bindings() or [])
        except Exception:
            bindings = []
    for binding in bindings:
        row = dict(binding)
        widget_name = row.get("widget", "")
        widget = _resolve_widget(plugin, widget_name)
        current = _widget_value(widget)
        resolved = _resolved_value(row, current)
        active = True
        enabled_by = row.get("enabled_by", "")
        if enabled_by:
            controller = _resolve_widget(plugin, enabled_by)
            active = bool(_widget_value(controller))
        enabled_when = row.get("enabled_when") or {}
        if enabled_when:
            controller = _resolve_widget(plugin, enabled_when.get("widget", ""))
            controller_value = _widget_value(controller)
            expected = enabled_when.get("equals")
            if "in" in enabled_when:
                active = controller_value in (enabled_when.get("in") or [])
            else:
                active = controller_value == expected
        row["current"] = current
        row["resolved"] = resolved
        row["active"] = active
        row["cli_args"] = _format_cli(row, resolved) if active else []
        row["widget_found"] = widget is not None
        rows.append(row)
    return rows


def _generated_markdown(plugin, entry: Dict[str, Any],
                        controls: List[Dict[str, Any]]) -> str:
    title = entry.get("title") or getattr(plugin, "plugin_name", "Plugin")
    summary = entry.get("summary", "")
    lines = ["# %s" % title, ""]
    if summary:
        lines.extend([summary, ""])
    tool_path = getattr(plugin, "tool_path", "") or ""
    if tool_path:
        lines.extend(["Executable: `%s`" % tool_path, ""])
    menu_path = entry.get("menu_path") or []
    if menu_path:
        lines.extend(["Menu: %s" % " -> ".join(menu_path), ""])
    if controls:
        lines.extend([
            "## GUI to CLI",
            "",
            "| Control | CLI | Current |",
            "|---|---|---|",
        ])
        for row in controls:
            cli = " ".join(row.get("cli_args") or []) or "GUI-only"
            lines.append("| %s | `%s` | %s |" % (
                row.get("label", row.get("id", "")), cli, row.get("current", "")
            ))
        lines.append("")
    lines.extend([
        "## Discovery",
        "",
        "Call `workbench.menu` for the current action tree and "
        "`workbench.workspace` for available data.",
        "",
    ])
    return "\n".join(lines)


def build_plugin_help_doc(plugin) -> Dict[str, Any]:
    """Build the live help/skill payload for a plugin instance."""
    library = SkillLibrary()
    entry = library.find_for_plugin(
        plugin_name=getattr(plugin, "plugin_name", ""),
        tool_name=getattr(plugin, "tool_name", ""),
        class_name=plugin.__class__.__name__,
    )
    if entry is None:
        entry = {
            "id": slug(getattr(plugin, "plugin_name", "plugin")),
            "title": getattr(plugin, "plugin_name", "Plugin"),
            "summary": "",
            "menu_path": [],
            "actions": [],
        }
    controls = _binding_rows(plugin)
    actions = []
    contract = getattr(plugin, "agent_contract", None)
    if callable(contract):
        try:
            payload = contract()
            actions = payload.get("actions", {}).get("actions", [])
        except Exception:
            actions = []
    markdown = entry.get("markdown", "")
    if not markdown:
        markdown = _generated_markdown(plugin, entry, controls)
    cli = {
        "executable": getattr(plugin, "tool_path", "") or "",
        "help_option": "--help",
        "preview": [],
    }
    preview_method = getattr(plugin, "cli_preview", None)
    if callable(preview_method):
        try:
            preview = preview_method()
            cli["preview"] = preview if isinstance(preview, list) else [str(preview)]
        except Exception as exc:
            cli["preview_error"] = str(exc)
    warnings = []
    if getattr(plugin, "tool_path", "") and not controls:
        warnings.append(
            "No GUI to CLI control mapping declared; add control_bindings()."
        )
    return {
        "schema": SKILL_SCHEMA,
        "id": entry.get("id", slug(getattr(plugin, "plugin_name", "plugin"))),
        "title": entry.get("title", getattr(plugin, "plugin_name", "Plugin")),
        "summary": entry.get("summary", ""),
        "plugin": getattr(plugin, "plugin_name", ""),
        "tool": getattr(plugin, "tool_name", ""),
        "class": plugin.__class__.__name__,
        "menu_path": entry.get("menu_path", []),
        "controls": controls,
        "cli": cli,
        "actions": actions,
        "markdown": markdown,
        "warnings": warnings,
        "document_path": entry.get("document_path", ""),
    }
