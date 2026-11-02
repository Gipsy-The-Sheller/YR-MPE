# Semantic UI actions (v1)

> Goal: every user-visible operation should also be invokable by an agent,
> with named parameters instead of synthetic mouse/keyboard events.

## Problem with widget-only automation

Clicking a QPushButton by coordinates or object name is brittle: it depends on
layout, focus, enabled state, and locale. An agent should call the same intent
the button calls, with a typed payload:

```python
panel.agent_invoke("inspector.show_item", {"item_id": "item_..."})
panel.agent_invoke("inspector.preview", {"item_id": "item_...", "kind": "tree"})
panel.agent_invoke("inspector.rename_run", {"run_id": "run_...", "name": "IQ-TREE final"})
```

## Inspector contract

`InspectorPanel.agent_contract()` returns:

```json
{
  "schema": "yrmpe.ui.inspector.v1",
  "id": "inspector",
  "state": {"item_id": "...", "run_id": "...", "dataset_id": "..."},
  "actions": [
    {"id": "inspector.show_item", "parameters": {"item_id": {"type": "string", "required": true}}},
    {"id": "inspector.show_run", "parameters": {"run_id": {"type": "string", "required": true}}},
    {"id": "inspector.show_dataset", "parameters": {"dataset_id": {"type": "string", "required": true}}},
    {"id": "inspector.preview", "parameters": {"item_id": {"type": "string"}, "kind": {"type": "string"}}},
    {"id": "inspector.rename_run", "parameters": {"run_id": {"type": "string"}, "name": {"type": "string"}}},
    {"id": "inspector.invoke_operation", "parameters": {"operation": {"type": "string"}, "target_id": {"type": "string"}}}
  ]
}
```

`agent_invoke()` returns `{"ok": true}` or `{"ok": false, "error": "..."}`.
GUI buttons call the same handlers; the agent path is not a separate code path.

## Run parameters and rerun

`RunRecord` now has both `parameters` (values) and `parameter_schema` (semantic
descriptions: type, required, default, choices, label, help). The inspector can
render a form from the schema, and an agent can call:

```python
panel.agent_invoke("inspector.rerun_run", {
    "run_id": "run_...",
    "parameters": {"model": "GTR+F+I+G4", "boot": 1000},
})
```

The host maps `(method, tool)` to a plugin entry and passes the parameter dict;
it does not replay the old command line. This is why semantic parameters make
"Modify & rerun" simpler: the same dict drives the form, validation, agent call,
and the rerun.

## Future host registry

A YRTools-level `UIActionRegistry` should:

- register panels by id (`inspector`, `workspace`, `dataset_manager`, ...);
- expose `list_actions()` and `invoke(action_id, params)`;
- map each action to an object name and a Qt widget for focus/disabled state;
- log every invocation with `run_id` so the message bus can replay it;
- reject actions whose widget is disabled, unless the agent explicitly requests
  a headless equivalent.

Do not synthesize Qt mouse events for normal agent use. Keep a separate
screenshot/geometry API for visual testing only.
