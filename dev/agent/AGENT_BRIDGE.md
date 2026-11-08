# Local file bridge for external agents

> Goal: any CLI (or an LLM with shell access) can drive YR-MPE without an
> embedded model interface. YR-MPE watches a command file; the agent tails the
> event/log files. No network port, no shell execution.

## Session layout

```
<session>/
  commands.jsonl    # agent appends one JSON command per line
  events.jsonl      # YR-MPE appends structured events and responses
  responses.jsonl   # request/response pairs only
  stdout.log        # raw stdout/app lines
  stderr.log        # raw stderr/error lines
  control.json      # session id, pid, action list (never the token)
  READY             # present while the bridge is alive
```

`control.json` has `session_dir`, `schema`, `pid`, `token_required` and the
registered action list. It never contains the token itself; the token is
copied from About → Agent bridge. `<temp>/yrmpe-agent-*/latest.txt` points at
the most recent session.

## Command

One JSON object per line. `id` correlates the response; `op` is a registered
semantic action; `params` is an object.

```json
{"schema":"yrmpe.bridge.command.v1","id":"c1","client":"agent","op":"orthofinder.inspect_status","params":{}}
{"id":"c2","op":"orthofinder.modify_options","params":{"options":{"run_mode":"parse_only","threads":4}}}
{"id":"c3","op":"orthofinder.run","params":{"show_dialogs":false}}
{"id":"c4","op":"bridge.ping"}
```

Built-ins: `bridge.ping`, `bridge.actions`, `bridge.stop`.

## Response

```json
{"schema":"yrmpe.bridge.response.v1","id":"c1","op":"orthofinder.inspect_status",
 "client":"agent","ok":true,"ts":"2026-11-02T22:30:00+0800",
 "duration_ms":0.42,"result":{...}}
```

`responses.jsonl` is append-only; `events.jsonl` also carries a
`kind=response` event so a single tail can see everything.

## Event

```json
{"schema":"yrmpe.bridge.event.v1","seq":42,"ts":"...","stream":"stdout",
 "level":"info","kind":"log","source":"orthofinder",
 "message":"Parsed 4 species / 335 CDS","code":null,"run_id":null,"payload":{}}
```

`stream` is `stdout|stderr|app`; `level` is `debug|info|notice|warning|error|
fatal`. They are independent: a tool writing progress to stderr is not an
error. Raw `stdout.log` / `stderr.log` mirror the corresponding streams for
tools that cannot parse JSON.

## Tailing from any CLI

Bash:

```bash
tail -f "$SESSION/events.jsonl"
```

PowerShell:

```powershell
Get-Content "$SESSION\events.jsonl" -Wait -Tail 0
```

Python:

```python
with open(session / "events.jsonl", "r", encoding="utf-8") as handle:
    handle.seek(0, 2)  # start at EOF
    while True:
        line = handle.readline()
        if not line:
            time.sleep(0.1); continue
        event = json.loads(line)
        if event["kind"] == "response" and event["payload"]["id"] == command_id:
            ...
```

Appending a command:

```bash
printf '%s\n' '{"id":"c1","op":"bridge.ping"}' >> "$SESSION/commands.jsonl"
```

## Concurrency and durability

- **One writer** for `commands.jsonl` in v1. Multiple agents should use
  distinct `client` fields and serialize, or wait for a directory-queue mode.
- The bridge polls by byte offset and buffers partial lines; appending half a
  line is safe as long as the rest arrives later.
- If the file is truncated/rotated, the offset resets to 0.
- A command is processed at most once per bridge lifetime; offsets are not
  persisted across restarts.
- `bridge.stop` stops the bridge and removes `READY`. The host process may
  still be alive.

## Security

- Local files only; no TCP port, no shell execution.
- Only registered semantic actions can run. There is no `eval`, `exec`, or
  arbitrary command action.
- Sessions started from About require a random token; commands must echo it.
  `control.json` only reports `token_required`, never the token value.
- The session directory inherits the user's OS permissions; do not place it in
  a world-writable directory.
- Paths supplied to actions are validated by the action handler (for example,
  OrthoFinder validates input/runtime before running).

## QTimer vs a tail thread

`poll_once()` is not only "read agent input": it reads new command bytes,
dispatches registered actions, writes responses/events and emits heartbeats.
In the GUI it should run on the Qt main thread via a `QTimer` (100–200 ms).
Qt widgets must not be touched from a worker thread; dispatching actions in
the timer keeps `modify_options`, `run`, `show_item`, etc. on the GUI thread.

A `tail -f`-style thread is possible, but split it correctly:

- the thread only does file I/O and parses JSON;
- it emits a Qt signal with the command dict (queued connection);
- a main-thread slot calls `registry.invoke()`.

Do **not** call widget handlers directly from the tail thread. Regular files
do not support blocking `readline()` like a pipe, so a thread still needs a
poll/sleep loop or an OS watcher; `QFileSystemWatcher` can miss append/truncate
events on Windows, so byte-offset polling is more reliable. For true blocking
reads, add a named pipe / Unix FIFO / localhost socket transport later and
reuse the same command/response schema.

Start with `QTimer`; add a thread only if the GUI stalls or sub-50 ms command
latency is required. Tool runs, not file I/O, dominate the workflow.

## Starting a session from the Workbench

When YPhylo Workbench runs inside a host such as YRTools, no bridge exists
unless the user starts one. Use **WORKSPACE → About YPhylo Workbench →
Agent bridge → Start agent session**. The host creates a token-protected
`AgentBridge` in a temporary session directory and polls it from a `QTimer`
on the GUI thread; the About dialog then shows the ready-to-copy prompt.

The default registry exposes:

- `workbench.info` — application, workspace and runtime paths.
- `workbench.list_datasets` — workspace datasets, items and sequence counts.
- `workbench.list_windows` — open plugin windows that expose `agent_invoke`.
- `workbench.invoke_plugin` — call a plugin action by `window_id` and action id.
- `workbench.open_orthofinder` — open the One-Step OrthoFinder window.

Call `workbench.list_windows` first, then pass its `window_id` plus the
plugin's native action id (for example `orthofinder.run`) to
`workbench.invoke_plugin`. This keeps one stable host API while plugins keep
their own action namespaces.

A launcher that already created a bridge (for example
`E:/YPhylo-Workbench-AgentTest/start_workbench.py`) sets
`window.agent_bridge`; About reuses that session and does not start a second
one.

## Host wiring

1. Create `UIActionRegistry` at app startup.
2. Register the active plugins/panels:
   - OrthoFinder: `orthofinder.*` methods from `dev/agent/ORTHOFINDER_API.md`.
   - Inspector: `inspector.*` from `UI_SEMANTICS.md`.
3. Create `AgentBridge(session_dir=..., registry=...)`.
4. Use a `QTimer` (100–200 ms) in the Qt main thread to call `poll_once()`;
   never dispatch Qt widget actions from a worker thread.
5. Route plugin console output to `bridge.emit_stdout/emit_stderr` (or the
   message bus once it lands).
6. On app close, call `bridge.close()`.

A headless CLI can use `run_forever()` without Qt. The GUI must stay alive for
actions that manipulate widgets.

## Relation to other protocols

- This is the first transport for the message protocol: `events.jsonl` uses
  the same `stream` / `level` / `kind` separation.
- The `op` values are exactly `UIActionRegistry` ids, not invented commands.
- A future named-pipe / Unix-socket / localhost JSON-RPC transport can reuse
  the same command/response/event schema.

## Discovery and skills

An agent should not guess plugin names. The host exposes:

- `workbench.menu` - toolbar menu tree plus a flat list of generated action
  ids, labels, paths, enabled/checked state and a short `description` taken
  from the plugin manifest / skill catalog; `query` matches id, label or
  description. `workbench.invoke_menu` triggers
  a leaf action by id or label path; triggers are deferred so a modal dialog
  does not block the bridge response.
- `workbench.workspace` - datasets, items, selection states, basket ids and
  run/stream/import history.
- `workbench.sidebar` - Selected basket contents, recommended send-to
  operations (available and unavailable with reasons), the current Inspector
  operations and the task-queue summary.
- `workbench.jobs` - list queued/running/finished/error tasks. Each row has a
  stable `job_id`, `state`, `status`, progress text, command preview and
  artifact paths.
- `workbench.job_status` / `workbench.inspect_job` - inspect one job by
  `job_id`.
- `workbench.cancel_job` - cancel a queued or running job. A plugin's own
  Stop button uses the same queue cancellation path when the job is queued.
- `workbench.skills` - the skill catalog under `YR_MPE/skills/manifest.json`.
- `workbench.help_doc` - fetch a skill by id, plugin name, or open window id.
  A live plugin document includes `controls`: one row per GUI control with
  `current`, `resolved`, `cli_args`, and any GUI-only note.

Every `BasePlugin` also exposes a generic window API through
`workbench.invoke_plugin`:

- `plugin.help_doc`
- `plugin.inspect_status`
- `plugin.inspect_controls`
- `plugin.inspect_console`
- `plugin.modify_options` (`options` object keyed by control id/widget)
- `plugin.run` - submits through the host task queue and returns `job_id`
  when the workbench started a queue; poll `workbench.job_status`
- `plugin.stop` - cancels the queue job when one is active
- `plugin.export_to_workspace`
- `plugin.close`

Plugin-specific APIs such as `orthofinder.*` override the generic surface.

The skill contract requires `control_bindings()` to map each CLI-capable GUI
control to the exact option printed by the tool's `--help`. GUI-only controls
use `"cli": []` and a note. Tests for ModelTest-NG and RAxML-NG assert every
declared flag occurs in the executable's `--help` output.

## Workflow preview

- `workbench.workflow_graph` returns the current workspace (or a supplied
  compiled `plan`) as data/tool nodes and edges.
- `workbench.show_workflow_preview` opens the workspace graph.
- `workbench.preview_workflow` opens a compiled `yrmpe.workflow.v1` plan.
- Clicking a tool node opens its plugin in placeholder-locked workflow mode;
  the plugin exposes `workflow_configured` instead of running immediately.

## BEST / MrBayes MSC

- `workbench.open_tool` with `tool="best"` opens the BEST window.
- `best.apply_block` accepts a complete MRBAYES block; the form is updated
  from the parsed values and the response returns `canonical_block`,
  `unparsed_lines` and `roundtrip_equal`.
- `best.modify_options` addresses individual controls; `best.run` starts the
  single-process MrBayes runner; `best.stop` cancels it.
