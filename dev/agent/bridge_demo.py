"""Demo of the local file bridge with a few semantic actions.

Run:

    python dev/agent/bridge_demo.py

Then append commands to <session>/commands.jsonl and tail events.jsonl.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from YR_MPE.runtime.agent_bridge import AgentBridge  # noqa: E402
from YR_MPE.runtime.ui_actions import UIActionRegistry  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--session-dir", default=None)
    parser.add_argument("--token", default=None)
    args = parser.parse_args()

    registry = UIActionRegistry()
    registry.register(
        "demo.echo",
        lambda params: {"echo": params},
        description="Echo the parameters",
        parameters={"value": {"type": "string"}},
    )
    registry.register(
        "demo.add",
        lambda params: params["a"] + params["b"],
        description="Add two numbers",
        parameters={"a": {"type": "number"}, "b": {"type": "number"}},
    )

    bridge = AgentBridge(
        session_dir=args.session_dir,
        registry=registry,
        token=args.token,
        heartbeat_interval=2.0,
    )

    def status(params):
        return {
            "session_id": bridge.session_id,
            "commands": str(bridge.commands_path),
            "events": str(bridge.events_path),
            "seq": bridge._seq,
        }

    registry.register("demo.status", status, description="Bridge status")
    bridge._write_control(state="running")

    print("session:", bridge.session_dir)
    print("commands:", bridge.commands_path)
    print("events:", bridge.events_path)
    print("stdout:", bridge.stdout_path)
    print("stderr:", bridge.stderr_path)
    print()
    print('append a command:')
    print('  printf \'%%s\\n\' \'{"id":"c1","op":"demo.add","params":{"a":1,"b":2}}\' >> "%s"'
          % bridge.commands_path)
    print("tail events:")
    print('  tail -f "%s"' % bridge.events_path)
    print()
    bridge.run_forever()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
