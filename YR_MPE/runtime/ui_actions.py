"""Registry for semantic UI actions.

Panels and plugins register named actions with a parameter schema. GUI buttons
and external agents call the same handler; the registry only deals in
JSON-serializable dicts so it can be driven by the file bridge or a future
YRTools agent.

Example::

    registry = UIActionRegistry()
    registry.register("inspector.show_item", panel.show_item,
                      parameters={"item_id": {"type": "string", "required": True}})
    registry.invoke("inspector.show_item", {"item_id": "item_1"})
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional


@dataclass
class UIAction:
    id: str
    handler: Callable[[Dict[str, Any]], Any]
    description: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)
    enabled: bool = True

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "description": self.description,
            "parameters": self.parameters,
            "enabled": self.enabled,
        }


class UIActionRegistry:
    """A small, transport-agnostic action registry."""

    def __init__(self):
        self._actions: Dict[str, UIAction] = {}

    def register(
        self,
        action_id: str,
        handler: Callable[[Dict[str, Any]], Any],
        description: str = "",
        parameters: Optional[Dict[str, Any]] = None,
        enabled: bool = True,
    ) -> UIAction:
        if action_id in self._actions:
            raise ValueError("duplicate UI action %r" % action_id)
        action = UIAction(
            id=action_id,
            handler=handler,
            description=description,
            parameters=parameters or {},
            enabled=enabled,
        )
        self._actions[action_id] = action
        return action

    def unregister(self, action_id: str) -> None:
        self._actions.pop(action_id, None)

    def get(self, action_id: str) -> Optional[UIAction]:
        return self._actions.get(action_id)

    def list_actions(self) -> List[Dict[str, Any]]:
        return [action.to_dict() for action in self._actions.values()]

    def invoke(self, action_id: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        action = self._actions.get(action_id)
        if action is None:
            return {"ok": False, "error": "unknown action %r" % action_id}
        if not action.enabled:
            return {"ok": False, "error": "action %r is disabled" % action_id}
        try:
            result = action.handler(params or {})
        except Exception as exc:
            return {"ok": False, "error": str(exc)}
        if isinstance(result, dict) and "ok" in result:
            return result
        return {"ok": True, "result": result}
