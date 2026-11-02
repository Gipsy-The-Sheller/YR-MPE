# YR-MPE / YRTools message protocol (v1 proposal)

> Pilot: the One-Step OrthoFinder plugin. If it works here, promote the same
> protocol to `BasePlugin` / `BaseProcessThread` and the YRTools Console.

## Why

Today plugins emit two different things:

- `console_output(str, str)` — a free-form `(message, type)` pair.
- ad-hoc `QMessageBox` dialogs and `print()` calls.

That is hard to filter, replay, or aggregate across plugins. A single
structured message queue gives the Console one source of truth, makes
`stderr` a first-class semantic channel, and lets errors be handled by code
instead of by dialog strings.

## Core rule: stream ≠ level

| Field | Meaning | Examples |
|---|---|---|
| `stream` | Where the bytes came from | `stdout`, `stderr`, `app`, `host` |
| `level` | How severe the message is | `debug`, `info`, `notice`, `warning`, `error`, `fatal` |
| `kind` | What the message represents | `log`, `progress`, `result`, `state`, `metric` |

A line written by a tool to stderr is **`stream=stderr`**, but its default
level is `warning`, not `error`. It becomes `error` only when the process
exits nonzero, an exception is raised, or a known error pattern matches.
Progress is **not** an error even if the tool prints it to stderr.

## Message schema

One JSON object per line (JSONL), UTF-8, no embedded newlines.

```json
{
  "schema": "yrmpe.message.v1",
  "id": "2f5c...",
  "run_id": "20261102_211500_orthofinder",
  "seq": 42,
  "ts": "2026-11-02T21:15:03.123+08:00",
  "source": "orthofinder",
  "stream": "stderr",
  "level": "warning",
  "kind": "log",
  "message": "WARNING: species name contains a dot",
  "code": "ORTHOFINDER_SPECIES_NAME",
  "context": {
    "plugin": "one_step_orthofinder",
    "tool": "OrthoFinder 2.5.5",
    "command": ["orthofinder", "-f", "..."],
    "input": "testdata/gb",
    "output": "testdata/work"
  },
  "payload": {}
}
```

Field rules:

- `schema`: fixed `yrmpe.message.v1`; bump to `v2` on breaking changes.
- `id`: UUID4 per message.
- `run_id`: one per plugin run; used to group, replay, and link Console tabs.
- `seq`: monotonic integer per run; ordering even if timestamps collide.
- `ts`: ISO-8601 with timezone.
- `source`: `orthofinder`, `blastp`, `makeblastdb`, `plugin`, `host`, ...
- `stream`: `stdout|stderr|app|host`.
- `level`: `debug|info|notice|warning|error|fatal`.
- `kind`: `log|progress|result|state|metric`.
- `message`: human-readable one-line text.
- `code`: stable machine-readable code, e.g. `ORTHOFINDER_EXIT_1`,
  `PLUGIN_EXCEPTION`, `MISSING_BLASTP`; optional for plain logs.
- `context`: small dict (never large payloads); command, input, output, exit
  code, tool version.
- `payload`: structured data for `progress`/`result`/`metric`.

## Message taxonomy

| Event | stream | level | kind | payload |
|---|---|---|---|---|
| Tool stdout line | `stdout` | `info` | `log` | — |
| Tool stderr line | `stderr` | `warning` by default; `error` on nonzero exit | `log` | — |
| Plugin progress | `app` | `info` | `progress` | `{"current": 3, "total": 10}` |
| Result summary | `app` | `info` | `result` | `{"species": 4, "cds": 335, ...}` |
| State transition | `app` | `info` | `state` | `{"state": "running"}` |
| Recoverable problem | `app`/`stderr` | `warning` | `log` | — |
| Fatal error | `stderr`/`app` | `error`/`fatal` | `state` | `{"state": "failed"}` |
| Metric | `app` | `info` | `metric` | `{"elapsed_s": 19.2}` |

## API (pilot)

```python
from YR_MPE.runtime.messages import MessageBus

bus = MessageBus(
    run_id="20261102_211500_orthofinder",
    sinks=[JsonlSink("workspace/logs/20261102_211500_orthofinder.jsonl"),
           QtSignalSink(self.message_signal)],
)

bus.stdout("orthofinder", line)
bus.stderr("orthofinder", line)
bus.progress("plugin", current=3, total=10, label="Running BLAST")
bus.result("plugin", {"species": 4, "cds": 335, "single_copy_ogs": 62})
bus.state("plugin", "running")
bus.error("plugin", "OrthoFinder exited 1", code="ORTHOFINDER_EXIT_1",
          context={"exit_code": 1})
```

Implementation:

- Core bus is **stdlib only** (`queue.Queue`, `json`, `uuid`, `datetime`).
  It must work headless/CLI without PyQt.
- `QtSignalSink` is a thin adapter (`pyqtSignal(object)`); the GUI consumes
  dicts, not formatted strings.
- Sinks are pluggable: JSONL file, Console ring buffer, pytest capture,
  future YRTools host bus.
- Bounded queue (`maxsize=10000`): if full, drop `debug`/`info` first and emit
  a `memory.dropped` metric; never block a running tool.
- The Console/UI is a **consumer**, not a writer. It renders:
  `[stream badge] [level] source: message`, with filters for run_id/source/
  stream/level/kind, copy/export, and a "stderr" badge that is not red unless
  level is error.

## Mapping existing code

`BaseProcessThread.execute_command` currently emits
`console_output(line, "output"|"error")`. Keep that for compatibility, but
add an optional `message_bus`:

```python
if self.message_bus:
    self.message_bus.stdout(self.get_tool_name(), line)
else:
    self.console_output.emit(line, "output")
```

`BasePlugin.add_console_message(message, type)` remains as a legacy adapter:

```python
if self.message_bus:
    self.message_bus.log(level=..., stream=..., source="plugin", message=message)
else:
    self.console_output.emit(message, type)
```

That gives an incremental migration: new plugins use the bus; old plugins keep
working; the Console can read JSONL from both.

## Console / YRTools integration

- **In-process**: YRTools Console subscribes to the global bus (same process)
  and shows tabs/filters per `run_id`/`source`.
- **Cross-process / replay**: tail the JSONL file under
  `workspace/logs/<run_id>.jsonl`; the Console can open it directly.
- **No raw `print()`** in new plugins; every user-visible line is a message.
- **Dialogs**: a fatal error is a message with `level=error`/`fatal` and a
  stable `code`; the plugin may show a dialog as a *view* of that message,
  but the message is the source of truth.
- **Progress**: `kind=progress` messages are coalesced by the UI (e.g. keep the
  latest per `source`/`label`) to avoid flooding the queue.

## Acceptance for the pilot

- OrthoFinder stdout/stderr are both visible; stderr is not all red.
- `run_id` groups one pipeline; JSONL replays the Console exactly.
- A nonzero exit produces one `error`/`fatal` message with a code and exit
  code, plus a final `state=failed`.
- Legacy `console_output` plugins still work unchanged.
- Unit tests: schema validation, seq ordering, JSONL round-trip, queue
  backpressure, QtSignalSink (headless).

## Open questions

1. Should the bus live in `YR_MPE.runtime` (plugin package) or in YRTools
   host? Proposal: start in `YR_MPE.runtime`, promote later.
2. JSONL retention policy: per-run file, keep N days or N runs?
3. Should `stderr` default to `warning` or `notice`? Proposal: `warning`, with
   `error` escalation on nonzero exit/exception.
4. Do we want a global `run_id` across chained YR-MPEA steps, or one per
   plugin invocation? Proposal: one per plugin run, plus a `session_id` later.
