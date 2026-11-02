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
  control.json      # session id, pid, token, action list
  READY             # present while the bridge is alive
```

`control.json` has `session_dir`, `schema`, `pid`, `token_required`,
`token`, and the registered action list. `<temp>/yrmpe-agent-*/latest.txt`
points at the most recent session.

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
- Optional `token` in `control.json`; commands must echo it.
- The session directory inherits the user's OS permissions; do not place it in
  a world-writable directory.
- Paths supplied to actions are validated by the action handler (for example,
  OrthoFinder validates input/runtime before running).

## Host wiring (next step)

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
