"""Machine-readable view of the Workbench toolbar menus.

The toolbar is built from ``QToolButton`` instances rather than a QMenuBar.
This module walks those buttons and their submenus, assigns stable path-based
ids, attaches a short description when plugin metadata provides one, and can
resolve an id back to its ``QAction`` for invocation.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from PyQt5.QtWidgets import QAction, QToolBar, QToolButton, QWidget


_DESCRIPTION_CACHE = None


def _normalise_label(text: str) -> str:
    value = str(text or "").replace("&", "").strip().lower()
    return re.sub(r"[^a-z0-9]+", " ", value).strip()


def _description_maps():
    """Load short descriptions from ``settings.json`` and the skill manifest.

    ``settings.json`` is the canonical plugin manifest and its ``description``
    is preferred when it matches a menu label. The skill manifest's ``summary``
    fills entries that have no plugin description.
    """
    global _DESCRIPTION_CACHE
    if _DESCRIPTION_CACHE is not None:
        return _DESCRIPTION_CACHE

    by_path: Dict[Tuple[str, ...], str] = {}
    by_label: Dict[str, str] = {}
    settings_by_label: Dict[str, str] = {}
    root = Path(__file__).resolve().parents[1]
    manifest_path = root / "skills" / "manifest.json"
    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            for entry in payload.get("skills", []):
                summary = str(entry.get("summary", "") or "").strip()
                if not summary:
                    continue
                path = tuple(str(value) for value in entry.get("menu_path", []))
                if path:
                    by_path[path] = summary
                    by_label.setdefault(_normalise_label(path[-1]), summary)
                by_label.setdefault(_normalise_label(entry.get("title", "")), summary)
        except Exception:
            pass

    settings_path = Path(__file__).resolve().parents[2] / "settings.json"
    if settings_path.exists():
        try:
            payload = json.loads(settings_path.read_text(encoding="utf-8"))
            for entry in payload.get("plugins", []):
                name = str(entry.get("name", "") or "").strip()
                description = str(entry.get("description", "") or "").strip()
                if name and description:
                    settings_by_label[_normalise_label(name)] = description
        except Exception:
            pass

    _DESCRIPTION_CACHE = (by_path, by_label, settings_by_label)
    return _DESCRIPTION_CACHE


def _settings_description(label: str, settings_by_label: Dict[str, str]) -> str:
    key = _normalise_label(label)
    if not key:
        return ""
    exact = settings_by_label.get(key)
    if exact:
        return exact
    key_tokens = set(key.split())
    best = ""
    best_score = 0.0
    best_has_overlap = False
    for name, description in settings_by_label.items():
        name_tokens = set(name.split())
        if not name_tokens or not key_tokens:
            continue
        overlap = key_tokens & name_tokens
        has_overlap = bool(overlap)
        # Prefer a real token overlap over a generic prefix match. This keeps
        # "MrBayes3-MPI-BEAGLE" on the MPI entry rather than on "MrBayes".
        if has_overlap:
            score = len(overlap) / len(key_tokens | name_tokens)
        elif name.startswith(key) or key.startswith(name):
            score = 0.75
        else:
            continue
        # A single distinctive token such as "SeqDBG" or "SeqMatrix" can
        # identify a tool even when the menu label has extra words.
        for token in overlap:
            if len(token) < 6:
                continue
            appearances = sum(
                1 for other in settings_by_label if token in other.split()
            )
            if appearances == 1:
                score = max(score, 0.60)
                break
        if (len(name_tokens) == 1 and name_tokens.issubset(key_tokens)
                and len(next(iter(name_tokens))) >= 4):
            token = next(iter(name_tokens))
            appearances = sum(
                1 for other in settings_by_label if token in other.split()
            )
            if appearances == 1:
                score = max(score, 0.60)
        if (has_overlap, score) > (best_has_overlap, best_score):
            best_has_overlap = has_overlap
            best_score = score
            best = description
    return best if best_score >= 0.5 else ""


def _action_description(action: QAction, path: List[str], label: str) -> str:
    by_path, by_label, settings_by_label = _description_maps()
    # The plugin manifest is the canonical short description for a tool.
    settings_text = _settings_description(label, settings_by_label)
    if settings_text:
        return settings_text
    exact = by_path.get(tuple(path))
    if exact:
        return exact
    known = by_label.get(_normalise_label(label))
    if known:
        return known
    for value in (action.toolTip(), action.statusTip(), action.whatsThis()):
        text = str(value or "").strip()
        if text and text != label:
            return text
    return ""


def _clean_label(text: str) -> str:
    value = re.sub(r"(?<!&)&(?!&)", "", str(text or ""))
    value = value.replace("&&", "&").strip()
    return re.sub(r"\s+", " ", value)


def _slug(text: str) -> str:
    value = str(text or "").replace("&", "").strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_") or "item"


def _unique_id(base: str, used: set) -> str:
    candidate = base
    index = 2
    while candidate in used:
        candidate = "%s_%d" % (base, index)
        index += 1
    used.add(candidate)
    return candidate


def _walk_actions(
    actions,
    path: List[str],
    used_ids: set,
    nodes: List[Dict[str, Any]],
    flat: List[Dict[str, Any]],
) -> None:
    for action in actions:
        if action.isSeparator():
            nodes.append({"type": "separator"})
            continue
        label = _clean_label(action.text())
        if not label:
            continue
        current_path = path + [label]
        action_id = _unique_id(
            _slug(".".join([path[0]] + current_path[1:])), used_ids
        )
        menu = action.menu()
        description = "" if menu is not None else _action_description(
            action, current_path, label
        )
        node = {
            "id": action_id,
            "label": label,
            "path": current_path,
            "type": "submenu" if menu is not None else "action",
            "enabled": bool(action.isEnabled()),
            "checkable": bool(action.isCheckable()),
            "checked": bool(action.isChecked()),
            "shortcut": action.shortcut().toString(),
            "tooltip": action.toolTip(),
            "description": description,
        }
        if menu is not None:
            node["items"] = []
            _walk_actions(menu.actions(), current_path, used_ids,
                          node["items"], flat)
        nodes.append(node)
        flat.append({
            key: node[key]
            for key in ("id", "label", "path", "type", "enabled",
                        "checkable", "checked", "shortcut", "description")
        })


def menu_tree(host: QWidget) -> Dict[str, Any]:
    """Return ``{"menus": [...], "actions": [...], "action_count": n}``."""
    menus: List[Dict[str, Any]] = []
    flat: List[Dict[str, Any]] = []
    used_ids: set = set()
    seen_buttons = set()
    for toolbar in host.findChildren(QToolBar):
        for button in toolbar.findChildren(QToolButton):
            if id(button) in seen_buttons:
                continue
            seen_buttons.add(id(button))
            label = button.text().strip()
            if not label:
                continue
            button_id = _unique_id(_slug(label), used_ids)
            items: List[Dict[str, Any]] = []
            actions = button.menu().actions() if button.menu() is not None else button.actions()
            _walk_actions(actions, [label], used_ids, items, flat)
            menus.append({
                "id": button_id,
                "label": label,
                "type": "menu",
                "tooltip": button.toolTip(),
                "items": items,
            })
    return {
        "menus": menus,
        "actions": flat,
        "action_count": len(flat),
    }


def find_action(host: QWidget, action_id: str = "",
                path: Optional[List[str]] = None) -> Tuple[Optional[QAction], Dict[str, Any]]:
    """Resolve an action by its generated id or by its label path."""
    tree = menu_tree(host)
    wanted_path = [str(item) for item in (path or [])]

    def walk(actions, current_path):
        for action in actions:
            if action.isSeparator():
                continue
            label = _clean_label(action.text())
            if not label:
                continue
            current = current_path + [label]
            node_id = _slug(".".join([current_path[0]] + current[1:]))
            if action_id and (node_id == action_id or
                              node_id.startswith(action_id + "_")):
                return action, current
            if wanted_path and current == wanted_path:
                return action, current
            if action.menu() is not None:
                found = walk(action.menu().actions(), current)
                if found[0] is not None:
                    return found
        return None, []

    seen = set()
    for toolbar in host.findChildren(QToolBar):
        for button in toolbar.findChildren(QToolButton):
            if id(button) in seen:
                continue
            seen.add(id(button))
            label = button.text().strip()
            if not label:
                continue
            actions = button.menu().actions() if button.menu() is not None else button.actions()
            action, found_path = walk(actions, [label])
            if action is not None:
                node = next(
                    (entry for entry in tree["actions"]
                     if entry["path"] == found_path),
                    None,
                )
                return action, node or {"path": found_path}
    return None, {}
