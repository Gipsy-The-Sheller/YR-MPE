"""Local file bridge for external agents.

YR-MPE does not embed an LLM. A GUI (or headless) process can expose a session
directory containing append-only JSONL files and raw stdout/stderr logs:

    <session>/
      commands.jsonl    # external agent appends one JSON command per line
      events.jsonl      # YR-MPE appends structured events/responses
      responses.jsonl   # request/response pairs only
      stdout.log        # raw stdout/app lines
      stderr.log        # raw stderr/error lines
      control.json      # session, pid, action list (never the token)
      READY             # present while the bridge is alive

The host polls ``poll_once()`` from a Qt timer (GUI thread) or
``run_forever()`` in a headless process. Commands are dispatched through a
``UIActionRegistry``; no shell is involved.

Protocol: see ``dev/agent/AGENT_BRIDGE.md``.
"""

from __future__ import annotations

import json
import os
import secrets
import tempfile
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

from .ui_actions import UIActionRegistry

COMMAND_SCHEMA = "yrmpe.bridge.command.v1"
RESPONSE_SCHEMA = "yrmpe.bridge.response.v1"
EVENT_SCHEMA = "yrmpe.bridge.event.v1"
CONTROL_SCHEMA = "yrmpe.bridge.control.v1"


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def _atomic_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temporary, path)


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


class AgentBridge:
    """Append-only file command/event bridge."""

    def __init__(
        self,
        session_dir: Optional[str] = None,
        registry: Optional[UIActionRegistry] = None,
        token: Optional[str] = None,
        allow_anonymous: bool = False,
        heartbeat_interval: float = 5.0,
    ):
        if session_dir:
            self.session_dir = Path(session_dir).expanduser().resolve()
        else:
            self.session_dir = Path(
                tempfile.mkdtemp(prefix="yrmpe-agent-", dir=tempfile.gettempdir())
            )
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.registry = registry or UIActionRegistry()
        if token is None and not allow_anonymous:
            token = secrets.token_urlsafe(24)
        self.token = token
        self.heartbeat_interval = float(heartbeat_interval)

        self.commands_path = self.session_dir / "commands.jsonl"
        self.events_path = self.session_dir / "events.jsonl"
        self.responses_path = self.session_dir / "responses.jsonl"
        self.stdout_path = self.session_dir / "stdout.log"
        self.stderr_path = self.session_dir / "stderr.log"
        self.control_path = self.session_dir / "control.json"
        self.ready_path = self.session_dir / "READY"

        for path in (self.commands_path, self.events_path, self.responses_path,
                     self.stdout_path, self.stderr_path):
            path.touch(exist_ok=True)

        self.session_id = uuid.uuid4().hex
        self._offset = 0
        self._buffer = ""
        self._seq = 0
        self._io_lock = threading.RLock()
        self._closed = False
        self._started_at = _now()
        self._last_heartbeat = 0.0

        self._write_control(state="running")
        self.ready_path.write_text(self.session_id, encoding="utf-8")
        self._write_latest_pointer()
        self.emit_event(
            stream="app", level="info", kind="state", source="bridge",
            message="Bridge ready", payload={"state": "ready"},
        )

    # ------------------------------------------------------------------
    # discovery
    # ------------------------------------------------------------------
    def _write_latest_pointer(self) -> None:
        try:
            (self.session_dir.parent / "latest.txt").write_text(
                str(self.session_dir), encoding="utf-8"
            )
        except Exception:
            pass

    def _write_control(self, state: str) -> None:
        _atomic_json(self.control_path, {
            "schema": CONTROL_SCHEMA,
            "session_id": self.session_id,
            "pid": os.getpid(),
            "state": state,
            "started_at": self._started_at,
            "last_heartbeat": _now(),
            "session_dir": str(self.session_dir),
            "commands": self.commands_path.name,
            "events": self.events_path.name,
            "responses": self.responses_path.name,
            "stdout": self.stdout_path.name,
            "stderr": self.stderr_path.name,
            "token_required": bool(self.token),
            "actions": self.registry.list_actions(),
        })

    # ------------------------------------------------------------------
    # output
    # ------------------------------------------------------------------
    def _append_jsonl(self, path: Path, payload: Dict[str, Any]) -> None:
        line = json.dumps(_jsonable(payload), ensure_ascii=False)
        with self._io_lock:
            with path.open("a", encoding="utf-8") as handle:
                handle.write(line + "\n")

    def emit_event(
        self,
        stream: str = "app",
        level: str = "info",
        kind: str = "log",
        source: str = "app",
        message: str = "",
        payload: Optional[Dict[str, Any]] = None,
        run_id: Optional[str] = None,
        code: Optional[str] = None,
    ) -> Dict[str, Any]:
        with self._io_lock:
            self._seq += 1
            event = {
                "schema": EVENT_SCHEMA,
                "seq": self._seq,
                "ts": _now(),
                "stream": stream,
                "level": level,
                "kind": kind,
                "source": source,
                "message": message,
                "code": code,
                "run_id": run_id,
                "payload": _jsonable(payload or {}),
            }
            self._append_jsonl(self.events_path, event)
            raw = message if message else json.dumps(event["payload"], ensure_ascii=False)
            if raw:
                target = self.stderr_path if stream == "stderr" else self.stdout_path
                with target.open("a", encoding="utf-8") as handle:
                    handle.write(raw + "\n")
        return event

    def emit_stdout(self, message: str, source: str = "app", level: str = "info") -> Dict[str, Any]:
        return self.emit_event(stream="stdout", level=level, kind="log",
                               source=source, message=message)

    def emit_stderr(self, message: str, source: str = "app", level: str = "warning") -> Dict[str, Any]:
        return self.emit_event(stream="stderr", level=level, kind="log",
                               source=source, message=message)

    def _respond(self, response: Dict[str, Any]) -> Dict[str, Any]:
        self._append_jsonl(self.responses_path, response)
        self.emit_event(
            stream="app", level="info" if response.get("ok") else "error",
            kind="response", source="bridge",
            message="response to %s" % response.get("id"),
            payload=response,
        )
        return response

    # ------------------------------------------------------------------
    # command polling
    # ------------------------------------------------------------------
    def _read_new_lines(self) -> list[str]:
        try:
            size = self.commands_path.stat().st_size
        except OSError:
            return []
        if size < self._offset:
            # Truncated/rotated: start from the beginning.
            self._offset = 0
            self._buffer = ""
        if size == self._offset:
            return []
        with self.commands_path.open("rb") as handle:
            handle.seek(self._offset)
            data = handle.read()
            self._offset = handle.tell()
        text = self._buffer + data.decode("utf-8", errors="replace")
        self._buffer = ""
        lines = text.splitlines(keepends=True)
        if lines and not lines[-1].endswith(("\n", "\r")):
            self._buffer = lines.pop()
        return [line.strip() for line in lines if line.strip()]

    def _dispatch(self, command: Dict[str, Any]) -> Dict[str, Any]:
        command_id = str(command.get("id") or uuid.uuid4().hex)
        op = command.get("op") or command.get("action")
        params = command.get("params") or {}
        client = command.get("client", "")
        started = time.time()

        if not op:
            return self._respond({
                "schema": RESPONSE_SCHEMA, "id": command_id, "op": "",
                "client": client, "ok": False, "error": "missing op",
                "ts": _now(),
            })
        if self.token and command.get("token") != self.token:
            return self._respond({
                "schema": RESPONSE_SCHEMA, "id": command_id, "op": op,
                "client": client, "ok": False, "error": "invalid token",
                "ts": _now(),
            })

        stop_after = False
        if op == "bridge.ping":
            outcome = {"ok": True, "result": {"pong": True, "session_id": self.session_id}}
        elif op == "bridge.actions":
            outcome = {"ok": True, "result": {"actions": self.registry.list_actions()}}
        elif op == "bridge.stop":
            outcome = {"ok": True, "result": {"stopping": True}}
            stop_after = True
        else:
            outcome = self.registry.invoke(op, params)

        response = {
            "schema": RESPONSE_SCHEMA,
            "id": command_id,
            "op": op,
            "client": client,
            "ok": bool(outcome.get("ok")),
            "ts": _now(),
            "duration_ms": round((time.time() - started) * 1000, 2),
        }
        if response["ok"]:
            response["result"] = _jsonable(outcome.get("result"))
        else:
            response["error"] = outcome.get("error", "unknown error")
        self._respond(response)
        if stop_after:
            self.close()
        return response

    def poll_once(self, max_commands: int = 200) -> int:
        """Process buffered commands. Call from the GUI thread."""
        if self._closed:
            return 0
        now = time.time()
        if self.heartbeat_interval > 0 and now - self._last_heartbeat >= self.heartbeat_interval:
            self._last_heartbeat = now
            self._write_control(state="running")
            self.emit_event(stream="app", level="info", kind="heartbeat",
                            source="bridge", message="heartbeat",
                            payload={"pid": os.getpid(), "state": "running"})

        processed = 0
        for line in self._read_new_lines():
            if processed >= max_commands:
                break
            try:
                command = json.loads(line)
                if not isinstance(command, dict):
                    raise ValueError("command must be a JSON object")
            except Exception as exc:
                self.emit_event(stream="stderr", level="error", kind="parse_error",
                                source="bridge", message="bad command: %s" % exc,
                                code="BRIDGE_BAD_COMMAND")
                continue
            self._dispatch(command)
            processed += 1
            if self._closed:
                break
        return processed

    def run_forever(self, poll_interval: float = 0.1) -> None:
        """Headless loop; not used by the GUI (which uses a QTimer)."""
        while not self._closed:
            self.poll_once()
            time.sleep(poll_interval)
        self.close()

    # ------------------------------------------------------------------
    # lifecycle
    # ------------------------------------------------------------------
    def close(self) -> None:
        if self._closed and not self.ready_path.exists():
            return
        self._closed = True
        self.emit_event(stream="app", level="info", kind="state", source="bridge",
                        message="Bridge closed", payload={"state": "closed"})
        self._write_control(state="closed")
        try:
            self.ready_path.unlink()
        except OSError:
            pass
