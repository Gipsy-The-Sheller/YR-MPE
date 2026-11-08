"""Tests for the file-based agent bridge (stdlib only)."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# Load the two runtime modules without importing YR_MPE/__init__ (PyQt).
ui_actions = _load("yrmpe_ui_actions", REPO / "YR_MPE" / "runtime" / "ui_actions.py")
sys.modules["yrmpe_runtime_ui_actions"] = ui_actions

# agent_bridge imports .ui_actions relatively; patch its package context by
# loading it under a tiny package-like namespace.
import types  # noqa: E402

package = types.ModuleType("yrmpe_runtime")
package.__path__ = [str(REPO / "YR_MPE" / "runtime")]
sys.modules["yrmpe_runtime"] = package
sys.modules["yrmpe_runtime.ui_actions"] = ui_actions
bridge = _load("yrmpe_runtime.agent_bridge", REPO / "YR_MPE" / "runtime" / "agent_bridge.py")


class AgentBridgeTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.registry = ui_actions.UIActionRegistry()
        self.registry.register(
            "math.add", lambda p: p["a"] + p["b"],
            description="add two numbers",
            parameters={"a": {"type": "number"}, "b": {"type": "number"}},
        )
        self.registry.register("demo.echo", lambda p: {"echo": p.get("value")})
        self.session = self.root / "session"
        self.bridge = bridge.AgentBridge(
            str(self.session), registry=self.registry,
            allow_anonymous=True, heartbeat_interval=0,
        )

    def tearDown(self):
        self.bridge.close()

    def _append(self, command):
        with (self.session / "commands.jsonl").open("a", encoding="utf-8") as handle:
            if isinstance(command, str):
                handle.write(command + "\n")
            else:
                handle.write(json.dumps(command) + "\n")

    def _responses(self):
        path = self.session / "responses.jsonl"
        return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]

    def test_command_response_and_events(self):
        self._append({"id": "c1", "op": "math.add", "params": {"a": 2, "b": 3}})
        self.assertEqual(self.bridge.poll_once(), 1)
        response = self._responses()[0]
        self.assertTrue(response["ok"])
        self.assertEqual(response["id"], "c1")
        self.assertEqual(response["result"], 5)
        self.assertIn("duration_ms", response)

        events = [
            json.loads(line) for line in
            (self.session / "events.jsonl").read_text(encoding="utf-8").splitlines()
            if line
        ]
        self.assertTrue(any(event["kind"] == "response" for event in events))

    def test_unknown_action_and_bad_json(self):
        self._append({"id": "c1", "op": "nope", "params": {}})
        self._append("{not json")
        self.bridge.poll_once()
        responses = self._responses()
        self.assertFalse(responses[0]["ok"])
        self.assertIn("unknown action", responses[0]["error"])
        self.assertIn("bad command", (self.session / "stderr.log").read_text(encoding="utf-8"))

    def test_partial_line_and_token(self):
        command = json.dumps({"id": "c1", "op": "math.add", "params": {"a": 1, "b": 1}})
        with (self.session / "commands.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(command[:10])
        self.assertEqual(self.bridge.poll_once(), 0)
        with (self.session / "commands.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(command[10:] + "\n")
        self.assertEqual(self.bridge.poll_once(), 1)
        self.assertEqual(self._responses()[0]["result"], 2)

        secured = bridge.AgentBridge(
            str(self.root / "secured"), registry=self.registry, token="secret",
            heartbeat_interval=0,
        )
        try:
            with (self.root / "secured" / "commands.jsonl").open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({"id": "c2", "op": "math.add", "params": {"a": 1, "b": 2}}) + "\n")
            secured.poll_once()
            response = [json.loads(line) for line in
                        (self.root / "secured" / "responses.jsonl").read_text(encoding="utf-8").splitlines()
                        if line][0]
            self.assertFalse(response["ok"])
            self.assertIn("token", response["error"])
        finally:
            secured.close()

    def test_token_is_generated_and_required_by_default(self):
        secured = bridge.AgentBridge(
            str(self.root / "default_token"), registry=self.registry,
            heartbeat_interval=0,
        )
        try:
            self.assertTrue(secured.token)
            control = json.loads(
                (self.root / "default_token" / "control.json").read_text(encoding="utf-8")
            )
            self.assertTrue(control["token_required"])
            self.assertNotIn("token", control)
            self._append_to(secured, {"id": "t1", "op": "math.add", "params": {"a": 1, "b": 2}})
            secured.poll_once()
            response = self._last_response(secured)
            self.assertFalse(response["ok"])
            self.assertIn("token", response["error"])
            self._append_to(secured, {
                "id": "t2", "op": "math.add", "params": {"a": 1, "b": 2},
                "token": secured.token,
            })
            secured.poll_once()
            self.assertTrue(self._last_response(secured)["ok"])
        finally:
            secured.close()

    def _append_to(self, bridge_obj, command):
        path = Path(bridge_obj.session_dir) / "commands.jsonl"
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(command) + "\n")

    @staticmethod
    def _last_response(bridge_obj):
        path = Path(bridge_obj.session_dir) / "responses.jsonl"
        lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
        return json.loads(lines[-1])

    def test_builtin_actions_and_control(self):
        self._append({"id": "p", "op": "bridge.ping"})
        self._append({"id": "a", "op": "bridge.actions"})
        self.bridge.poll_once()
        responses = self._responses()
        self.assertTrue(responses[0]["result"]["pong"])
        ids = [action["id"] for action in responses[1]["result"]["actions"]]
        self.assertIn("math.add", ids)

        control = json.loads((self.session / "control.json").read_text(encoding="utf-8"))
        self.assertTrue((self.session / "READY").exists())
        self.assertEqual(control["schema"], "yrmpe.bridge.control.v1")
        self.assertEqual(control["pid"], __import__("os").getpid())

    def test_stop_and_close(self):
        self._append({"id": "s", "op": "bridge.stop"})
        self.bridge.poll_once()
        self.assertFalse((self.session / "READY").exists())
        control = json.loads((self.session / "control.json").read_text(encoding="utf-8"))
        self.assertEqual(control["state"], "closed")

    def test_stdout_stderr_logs(self):
        self.bridge.emit_stdout("hello")
        self.bridge.emit_stderr("problem", level="error")
        self.assertIn("hello", (self.session / "stdout.log").read_text(encoding="utf-8"))
        self.assertIn("problem", (self.session / "stderr.log").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
